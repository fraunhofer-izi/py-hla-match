from pyard import ard

# singleton instance of ARD
_ard_instance = None


def get_ard_instance(imgt_version="Latest", data_dir=None, **config):
    global _ard_instance
    if _ard_instance is None:
        _ard_instance = ard.ARD(
            imgt_version=imgt_version, data_dir=data_dir, **config
        )
    return _ard_instance
