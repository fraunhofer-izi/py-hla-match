import time
import logging
from typing import Dict, Union
import psutil

from py_hla_match.parser import HLADataSource
from py_hla_match.export import PairwiseMatch


def benchmark_hla_match(
    patient_path: str,
    donor_path: str,
    output_path: str,
    n_patients: int,
    n_donors: int,
    *,
    chunk_size: int = 10000,
    col_idx_start: int = 1,
    col_idx_stop: int = 13,
    row_idx_start: int = 1,
    verbose: bool = False
) -> Dict[str, Union[int, float, bool]]:
    """
    Lightweight benchmark for py_hla_match.
    """

    # logging
    if verbose:
        logging.basicConfig(format="%(message)s", level=logging.INFO)
    log = logging.getLogger("hla-bench")

    proc = psutil.Process()
    t0 = time.perf_counter()

    # core library work
    src = HLADataSource(patient_path,
                        col_idx_start=col_idx_start,
                        col_idx_stop=col_idx_stop,
                        row_idx_start=row_idx_start)

    tgt = HLADataSource(donor_path,
                        col_idx_start=col_idx_start,
                        col_idx_stop=col_idx_stop,
                        row_idx_start=row_idx_start)

    PairwiseMatch(
        source=src,
        target=tgt,
        storage_filename=output_path,
        stream=True,
        chunk_size=chunk_size
    ).run()

    # stats
    runtime = time.perf_counter() - t0
    mem_peak_mb = round(proc.memory_info().rss / 1024 / 1024, 1)

    comparisons = min(n_patients, n_donors)
    throughput = int(comparisons / runtime) if runtime else 0

    stats: Dict[str, int | float | bool] = {
        "patients": n_patients,
        "donors": n_donors,
        "total_comparisons": comparisons,
        "time_s": round(runtime, 2),
        "comparisons_per_s": throughput,
        "mem_peak_mb": mem_peak_mb,
        "success": True,
        "output": output_path,
    }

    log.info("Finished in %.1fs – %s cmp/s – mem %.1f MB",
             runtime, f"{throughput:,}", mem_peak_mb)

    return stats
