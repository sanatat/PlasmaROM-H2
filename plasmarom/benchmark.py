"""Small serial CPU walltime check, NOT HPC scaling."""
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from .flow import FlowConfig,run_flow


def run_benchmark(out,grid_sizes=(16,24,32),repeats=3):
    rows=[]
    for n in grid_sizes:
        timings=[]
        for _ in range(repeats):
            begin=perf_counter()
            run_flow(FlowConfig(n=n,t_end=.12,dt=.012))
            timings.append(perf_counter()-begin)
        rows.append({'n':n,'cpu_run_seconds_median':float(np.median(timings)),
                     'repeats':repeats})
    result={'status':'single-process serial Python benchmark; NOT parallel or HPC scaling', 'results':rows}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'BENCHMARK.json').write_text(json.dumps(result,indent=2))
    return result
