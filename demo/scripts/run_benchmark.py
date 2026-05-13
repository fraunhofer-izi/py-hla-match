import os
import pandas as pd

from demo.scripts.performance_benchmark import benchmark_hla_match

dataset_sizes = [1000, 10000, 100000, 200000, 500000, 1000000]

# run synthetic data generation (-> "synthetic_data_generator.py")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "synthetic_data")

results = []
for n in dataset_sizes:
    patient_path = os.path.join(DATA_DIR, f"recipients_{n}.csv")
    donor_path = os.path.join(DATA_DIR, f"donors_{n}.csv")
    output_path = os.path.join(DATA_DIR, f"match_results_{n}.csv")

    stats = benchmark_hla_match(
        donor_path,
        patient_path,
        output_path,
        n_patients=n,
        n_donors=n,
        chunk_size=100000,  # the bigger the better, because results are cached
        verbose=False
    )
    results.append(stats)

df = pd.DataFrame(results)[
    ["patients", "donors", "time_s", "comparisons_per_s", "mem_peak_mb"]
]

df.rename(columns={
    "patients": "patients",
    "donors": "donors",
    "time_s": "time [s]",
    "comparisons_per_s": "throughput [/s]",
    "mem_peak_mb": "mem peak [MB]"
}, inplace=True)

print(df)