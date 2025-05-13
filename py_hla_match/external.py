import requests


def is_permissive_dpb1_match(patient_dpb1: str, patient_dpb2: str, donor_dpb1: str, donor_dpb2: str) -> str:

    # Construct the query URL
    base_url = "https://www.ebi.ac.uk/cgi-bin/ipd/pl/hla/dpb_v2.cgi"
    query_params = {
        'pid': 1,
        'patdpb1': patient_dpb1,
        'patdpb2': patient_dpb2,
        'did': 2,
        'dondpb1': donor_dpb1,
        'dondpb2': donor_dpb2
    }

    response = requests.get(base_url, params=query_params)

    # Check if the request was successful
    if response.status_code != 200:
        raise Exception(f"API request failed with status code {response.status_code}")

    response_json = response.json()

    return response_json['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['results']['tce_prediction']

