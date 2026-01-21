import logging
from pyard import ard
import threading
from py_hla_match.config import get_config

logger = logging.getLogger(__name__)

# singleton instance of ARD
_ard_instance = None
_config_read = False
_ard_lock = threading.RLock()


def _reset_ard_instance():
    """
    Reset ARD singleton so subsequent get_ard_instance calls re-read
    config and re-initialize ARD.
    """
    global _ard_instance, _config_read
    with _ard_lock:
        _ard_instance = None
        _config_read = False


def get_ard_instance(imgt_version=None, data_dir=None, **config):
    global _ard_instance, _config_read
    with _ard_lock:
        if _ard_instance is None:
            # Read config only once, only if no explicit params provided
            if imgt_version is None and data_dir is None and not _config_read:
                cfg = get_config()
                imgt_version = cfg.ard_imgt_version
                data_dir = cfg.ard_data_dir
                _config_read = True

            # Fall back to original defaults if still None
            if imgt_version is None:
                imgt_version = "Latest"
            _ard_instance = ard.ARD(
                imgt_version=imgt_version,
                data_dir=data_dir,
                load_mac=False,
                **config
            )

    return _ard_instance
