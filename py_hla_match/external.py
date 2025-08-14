# external.py
import requests
import logging
from enum import Enum
from dataclasses import dataclass
from typing import Dict, Optional

# Set up logging
logger = logging.getLogger(__name__)


@dataclass
class DPB1TCEConfig:
    """Configuration for DPB1 TCE API endpoints and response keys"""
    endpoints: Dict[str, str] = None
    response_keys: Dict[str, str] = None

    def __post_init__(self):
        if self.endpoints is None:
            self.endpoints = {
                "2.0": "https://www.ebi.ac.uk/cgi-bin/ipd/matching/"
                "dpb1_tce_v2",
                "2.1": "https://www.ebi.ac.uk/cgi-bin/ipd/matching/"
                "dpb1_tce_v21",
                "3.0": "https://www.ebi.ac.uk/cgi-bin/ipd/matching/"
                "dpb1_tce_v3"
            }

        if self.response_keys is None:
            self.response_keys = {
                "2.0": "HLA-DPB1_TCE_report_V2.0",
                "2.1": "HLA-DPB1_TCE_report_V2.1",
                "3.0": "HLA-DPB1_TCE_report_V3.0"
            }


class DPB1TCEStatus(Enum):
    """TCE prediction results from the DPB1 TCE API"""
    API_ERROR = "Unknown API Error"
    INVALID_ALLELES = "API Invalid Alleles"
    PERMISSIVE = "Permissive"
    NON_PERMISSIVE_GVH = "Non-Permissive GvH"
    NON_PERMISSIVE_HVG = "Non-Permissive HvG"
    TIMEOUT_ERROR = "Timeout Error"
    REQUEST_ERROR = "Request Error"
    VALUE_ERROR = "Value Error"
    CONFIGURATION_ERROR = "Configuration Error"
    UNEXPECTED_ERROR = "Unexpected Error"


def query_dpb1_tce(
        patient_dpb1: str,
        patient_dpb2: str,
        donor_dpb1: str,
        donor_dpb2: str,
        version: str = "2.1",
        timeout: int = 10,
        config: Optional[DPB1TCEConfig] = None
) -> DPB1TCEStatus:
    """
    Query the EBI DPB1 TCE API for T-Cell Epitope matching.

    Args:
        patient_dpb1: Patient's first DPB1 allele (e.g. "01:01")
        patient_dpb2: Patient's second DPB1 allele
        donor_dpb1: Donor's first DPB1 allele
        donor_dpb2: Donor's second DPB1 allele
        version: API version - "2.0", "2.1", or "3.0"
        timeout: API request timeout in seconds
        config: Configuration object (default if None)

    Returns:
        DPB1TCEStatus enum
    """
    if config is None:
        config = DPB1TCEConfig()

    # Keys vary by version
    DONORS_KEY = 'donors'
    RESULT_KEY = 'result'

    # Validate API version
    if version not in config.endpoints or version not in config.response_keys:
        logger.error(
            f"Unsupported or misconfigured API version: {version}"
            f" must be endpoints key: {list(config.endpoints.keys())}"
            f" and response key: {list(config.response_keys.keys())}"
        )
        return DPB1TCEStatus.CONFIGURATION_ERROR

    # Prepare request
    url = config.endpoints[version]
    query_params = {
        'pid': 1,
        'patdpb1': patient_dpb1,
        'patdpb2': patient_dpb2,
        'did': 2,
        'dondpb1': donor_dpb1,
        'dondpb2': donor_dpb2
    }

    try:
        logger.debug(
            f"Querying EBI DPB1 TCE API v{version} with parameters: "
            f"{query_params}"
        )

        # Make request
        response = requests.get(url, params=query_params, timeout=timeout)
        response.raise_for_status()

        # Parse response
        response_json = response.json()

        # Navigate to the TCE prediction
        report = response_json.get(config.response_keys[version], {})
        if not report:
            logger.error(
                f"Key '{config.response_keys[version]}' not found in API "
                "response_json."
            )
            return DPB1TCEStatus.API_ERROR

        donors_list = report.get(DONORS_KEY, [])

        if not donors_list:
            logger.warning(
                f"No '{DONORS_KEY}' data (list) in EBI API response report."
            )
            return DPB1TCEStatus.API_ERROR

        first_donor_object = donors_list[0]

        result_details = first_donor_object.get(RESULT_KEY, {})
        tce_prediction = result_details.get('tce_prediction')

        # Check for invalid alleles
        if tce_prediction \
                and "Not possible as typing contains non-existent allele" \
                    in tce_prediction:
            logger.warning(f"EBI API reports invalid allele: {tce_prediction}")
            return DPB1TCEStatus.INVALID_ALLELES

        if not tce_prediction:
            logger.warning("No TCE prediction in EBI API response")
            return DPB1TCEStatus.API_ERROR

        # normalize
        normalized_prediction = tce_prediction.strip().lower()

        EXACT_TCE_MAPPINGS = {
            "permissive": DPB1TCEStatus.PERMISSIVE,
            "permissive (core)": DPB1TCEStatus.PERMISSIVE,
            "ard matched": DPB1TCEStatus.PERMISSIVE,

            "non-permissive gvh": DPB1TCEStatus.NON_PERMISSIVE_GVH,
            "non permissive gvh": DPB1TCEStatus.NON_PERMISSIVE_GVH,

            "non-permissive hvg": DPB1TCEStatus.NON_PERMISSIVE_HVG,
            "non permissive hvg": DPB1TCEStatus.NON_PERMISSIVE_HVG,
        }

        if "not possible as typing contains non-existent allele"\
                in normalized_prediction:
            logger.warning(f"EBI API reports invalid allele: {tce_prediction}")
            return DPB1TCEStatus.INVALID_ALLELES

        # exact matching for safety
        if normalized_prediction in EXACT_TCE_MAPPINGS:
            status = EXACT_TCE_MAPPINGS[normalized_prediction]
            logger.info(
                f"Mapping '{tce_prediction}' (normalized: "
                f"'{normalized_prediction}') "
                f"to {status.name}"
            )
            return status
        else:
            # unknown response
            logger.error(
                "CRITICAL: Unknown TCE prediction from EBI API: "
                f"'{tce_prediction}' (normalized: '{normalized_prediction}'). "
                "API may changed. Known values: "
                f"{list(EXACT_TCE_MAPPINGS.keys())}"
            )
            return DPB1TCEStatus.API_ERROR

    except requests.Timeout:
        logger.error(f"Timeout connecting to EBI DPB1 TCE API (>{timeout}s)")
        return DPB1TCEStatus.TIMEOUT_ERROR
    except requests.RequestException as e:
        logger.error(f"Error connecting to EBI DPB1 TCE API: {str(e)}")
        return DPB1TCEStatus.REQUEST_ERROR
    except ValueError as e:
        logger.error(f"Value error in EBI DPB1 TCE API response: {str(e)}")
        return DPB1TCEStatus.VALUE_ERROR
    except Exception as e:
        logger.error(f"Unexpected error in DPB1 TCE prediction: {str(e)}")
        return DPB1TCEStatus.UNEXPECTED_ERROR
