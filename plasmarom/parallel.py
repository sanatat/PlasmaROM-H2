"""Local multicore parameter sweep and serial/parallel equality check.

Benchmarks only independent 0D tasks. It is NOT an HPC solver benchmark,
and no speedup is assumed. Operating-system scheduling affects runtimes.
"""
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from time import perf_counter
from pathlib import Path
import json
import numpy as np
from .electron_energy import run_coupled_electron_energy,ElectronEnergyConfig


def _compute_power(power):
    d=run_coupled_electron_energy(config=ElectronEnergyConfig(input_power_W_m3=float(power)))
    return (float(d['te_eV'][-1]),float(d['fractions'][1,-1]),float(d['relative_energy_error']))


def compare_parallel(out, powers=(100.,500.,1000.,3000.), workers=2):
    if workers<1 or not powers or any(not np.isfinite(x) or x<0 for x in powers):
        raise ValueError('Invalid worker count or sweep parameters')
    start=perf_counter();serial=list(map(_compute_power,powers));tserial=perf_counter()-start
    start=perf_counter()
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
        parallel=list(pool.map(_compute_power,powers))
    tparallel=perf_counter()-start
    maxdiff=float(np.max(np.abs(np.asarray(serial)-np.asarray(parallel))))
    result={'status':'local CPU process-pool task-parallel demonstration, NOT HPC scaling',
            'powers_W_m3':list(powers),'worker_processes':workers,
            'serial_wall_s':tserial,'parallel_wall_s':tparallel,
            'measured_speedup':tserial/max(tparallel,1e-12),
            'max_absolute_serial_parallel_difference':maxdiff,
            'repetitions':1,'warning':'Do not extrapolate one-run walltime to HPC performance.'}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'parallel_comparison.json').write_text(json.dumps(result,indent=2))
    return result


def fixed_workload_scaling(out, powers=(100.,500.,1000.,3000.,5000.,7000.),
                           workers=(1,2,4), repetitions=3):
    """Repeated local fixed-workload task-parallel scaling, NOT DNS/LES HPC.

    Every worker count solves the SAME independent 0D cases. Uses median
    walltime and compares numerical outputs with a serial reference.
    """
    if (not powers or repetitions<1 or any(type(k) is not int or k<1 for k in workers)
        or not workers or 1 not in workers or len(set(workers)) != len(workers)
        or any(not np.isfinite(p) or p<0 for p in powers)):
        raise ValueError('Invalid fixed-workload scaling parameters')
    reference=list(map(_compute_power,powers))
    timings={}
    for n in workers:
        elapsed=[]
        for _ in range(repetitions):
            start=perf_counter()
            if n==1:
                values=list(map(_compute_power,powers))
            else:
                with ProcessPoolExecutor(max_workers=n,mp_context=get_context('spawn')) as pool:
                    values=list(pool.map(_compute_power,powers))
            elapsed.append(perf_counter()-start)
            if not np.allclose(values,reference,rtol=0,atol=0):
                raise RuntimeError('Serial/parallel mismatch')
        timings[n]={'median_wall_s':float(np.median(elapsed)),
                    'wall_times_s':[float(x) for x in elapsed]}
    baseline=timings[1]['median_wall_s']
    for n in workers:
        timings[n]['speedup_vs_serial']=baseline/timings[n]['median_wall_s']
        timings[n]['parallel_efficiency']=timings[n]['speedup_vs_serial']/n
    result={'kind':'local fixed-workload independent-task scaling, NOT solver HPC scaling',
            'workload_powers_W_m3':list(powers),'repetitions':repetitions,
            'by_workers':{str(n):timings[n] for n in workers},
            'warning':'Timing not portable. No MPI, GPU, CFD scaling or weak-scaling assessment.'}
    output=Path(out);output.mkdir(parents=True,exist_ok=True)
    (output/'fixed_workload_scaling.json').write_text(json.dumps(result,indent=2))
    return result
