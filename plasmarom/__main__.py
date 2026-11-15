"""Command-line interface: python -m plasmarom [run|kinetics|flow] [--quick]."""
import argparse
from .diagnostics import zero_d,spatial

def main():
    ap=argparse.ArgumentParser(description='Educational plasma-inspired chemical reduction and 2D reactive-flow prototype')
    ap.add_argument('command',choices=('run','kinetics','flow','cantera','skeletal','verify','plasma-demo','plasma','thermal','benchmark','electron-energy','eedf-rate','channel-reduce','flame-1d','parallel-check','report-v5','assess-skeletal','eedf-field','scaling','flame-grid','ammonia-h2','boltzmann','transport-1d'),nargs='?',default='run')
    ap.add_argument('--quick',action='store_true',help='Use 32x32 rather than 64x64 spatial grid')
    ap.add_argument('--out',default='results',help='Output directory')
    ap.add_argument('--mechanism', default='h2o2.yaml', help='Cantera mechanism (skeletal)')
    ap.add_argument('--protocol', default=None, help='Pinned v10 scientific evaluation JSON/id; required for verified v10 acceptance')
    ap.add_argument('--cases-file',default=None,
                    help='JSON with explicit training/holdout [T_K, P_atm, composition] cases')
    ap.add_argument('--reactor-t-end',type=float,default=None,
                    help='Optional zero-D reactor horizon in seconds for skeletal')
    ap.add_argument('--snapshots',type=int,default=None,
                    help='Optional number of reactor snapshots for skeletal')
    ap.add_argument('--sigma01',default=None,help='Cross section CSV for H2 v0->v1')
    ap.add_argument('--sigma12',default=None,help='Cross section CSV for H2 v1->v2')
    ap.add_argument('--source01',default=None,help='Cross section citation/source for v0->v1')
    ap.add_argument('--source12',default=None,help='Cross section citation/source for v1->v2')
    ap.add_argument('--te',type=float,default=2.0,help='Electron Te [eV]')
    ap.add_argument('--ne',type=float,default=1e17,help='Electron density [m^-3]')
    ap.add_argument('--power',type=float,default=3000.,help='Electron input power [W/m3] for 0D coupled model')
    ap.add_argument('--eedf',default=None,help='External EEDF CSV; probability-per-eV format')
    ap.add_argument('--source-eedf',default=None,help='EEDF provenance/citation')
    ap.add_argument('--field-Td',type=float,default=None,help='Target reduced field E/N [Townsend]')
    ap.add_argument('--eedf-manifest',default=None,help='CSV of externally computed EEDFs at E/N values')
    ap.add_argument('--summary',default=None,help='Input skeletal_summary.json for assess-skeletal')
    ap.add_argument('--require-holdout-pass',action='store_true',help='Fail command if no validated skeletal reduction')
    ap.add_argument('--parallel-workers',type=int,default=2,help='Local processes for parameter-sweep check')
    a=ap.parse_args()
    if a.command == 'transport-1d':
        from .reactive_transport_1d import run_transport_demo
        print('Verified scalar conservation only:',run_transport_demo(a.out))
        return
    if a.command == 'boltzmann':
        from .cantera_boltzmann import run_phelps_eedf
        en=(a.field_Td,) if a.field_Td is not None else (50.,100.,200.)
        print('True two-term Boltzmann EEDF in Cantera (N2/O2):',run_phelps_eedf(a.out,fields_Td=en))
        return
    if a.command == 'ammonia-h2':
        from .ammonia_h2 import run_ammonia_h2
        print('Real NH3/H2 Cantera (NO plasma):',run_ammonia_h2(a.out))
        return
    if a.command == 'assess-skeletal':
        from pathlib import Path
        from .acceptance import assess_skeletal_file
        source=a.summary or str(Path(a.out)/'skeletal_summary.json')
        assessment=assess_skeletal_file(source,Path(a.out)/'acceptance_assessment.json')
        print('Verified numerical acceptance evidence:',assessment)
        if not assessment['accepted']:
            raise SystemExit(1)
        return
    if a.command == 'eedf-field':
        from .eedf_field import run_field_rate
        from .plasma_states import levels_eV
        if a.eedf_manifest is None or a.field_Td is None or a.sigma01 is None or not a.source01:
            ap.error('eedf-field requires --eedf-manifest --field-Td --sigma01 --source01')
        result=run_field_rate(a.out,a.eedf_manifest,a.sigma01,a.source01,float(levels_eV()[1]),a.field_Td)
        print('Imported field-dependent EEDF convolution (NOT Boltzmann solver):',result)
        return
    if a.command == 'skeletal':
        from .cantera_demo import run_skeletal
        from pathlib import Path
        import json
        kwargs={}
        if a.cases_file:
            spec=json.loads(Path(a.cases_file).read_text(encoding='utf-8'))
            if not isinstance(spec,dict) or set(spec)-{'training','holdout'}:
                ap.error('--cases-file must contain only training and holdout arrays')
            kwargs.update(training=spec.get('training'),holdout=spec.get('holdout'))
        if a.protocol is not None:
            kwargs['protocol_path']=a.protocol
        if a.reactor_t_end is not None:
            kwargs['t_end']=a.reactor_t_end
        if a.snapshots is not None:
            kwargs['n_snapshots']=a.snapshots
        result=run_skeletal(a.out,mechanism=a.mechanism,**kwargs)
        print('Cantera skeleton:',result['status'])
        if a.require_holdout_pass and not result['acceptance_evidence']['accepted']:
            raise SystemExit('Strict skeletal validation failed: NO confirmed HOLDOUT_PASS')
        return
    if a.command=='report-v5':
        from .v5_report import generate_v5_report
        print('Computed report:',generate_v5_report(a.out))
        return
    if a.command=='electron-energy':
        from .electron_energy import run_energy_demo,ElectronEnergyConfig
        d=run_energy_demo(a.out,ElectronEnergyConfig(input_power_W_m3=a.power,ne_m3=a.ne,te_initial_eV=a.te))
        print('Power-coupled electron energy (synthetic cross sections):',d)
        return
    if a.command=='channel-reduce':
        from .channel_reduction import evaluate_channel_pruning
        d=evaluate_channel_pruning(a.out)
        print('Synthetic channel reduction:',d['status'],d['worst_holdout_error'])
        return
    if a.command=='eedf-rate':
        from .eedf import run_eedf_rate
        from .plasma_states import levels_eV
        if not all((a.eedf,a.sigma01,a.source_eedf,a.source01)):
            ap.error('eedf-rate requires --eedf --sigma01 --source-eedf --source01')
        d=run_eedf_rate(a.out,a.eedf,a.sigma01,a.source_eedf,a.source01,float(levels_eV()[1]))
        print('External EEDF collision rate [m3/s]:',d['rate_coefficient_m3_s'])
        return
    if a.command=='flame-grid':
        from .cantera_flame import run_flame_refinement_study
        d=run_flame_refinement_study(a.out,mechanism=a.mechanism)
        print('Cantera 1D flame mesh sensitivity, NOT experimental validation:',d['status'])
        if d['status'] != 'GRID_SENSITIVITY_PASS':
            raise SystemExit(1)
        return
    if a.command=='flame-1d':
        from .cantera_flame import run_flame
        d=run_flame(a.out,mechanism=a.mechanism)
        print('Cantera 1D chemical flame (NO plasma):',d)
        return
    if a.command=='scaling':
        from .parallel import fixed_workload_scaling
        d=fixed_workload_scaling(a.out,workers=tuple(sorted({1,2,a.parallel_workers})))
        print('Repeated local independent-task scaling (NOT HPC):',d['by_workers'])
        return
    if a.command=='parallel-check':
        from .parallel import compare_parallel
        d=compare_parallel(a.out,workers=a.parallel_workers)
        print('Local CPU parallel sweep:',d)
        return
    if a.command in ('plasma','plasma-demo'):
        from .plasma_states import PlasmaConfig,run_plasma
        kw={} if a.command=='plasma-demo' else dict(path01=a.sigma01,path12=a.sigma12,source01=a.source01,source12=a.source12)
        if a.command=='plasma' and not all(kw.values()):
            ap.error('plasma needs --sigma01 --sigma12 --source01 --source12')
        result=run_plasma(a.out,config=PlasmaConfig(te_eV=a.te,ne_m3=a.ne),**kw)
        print('Three-level H2 demonstration:',result['data_status'])
        return
    if a.command == 'verify':
        from .verification import convergence_study
        print('Observed periodic FV convergence:',{k:round(v['observed_order'],3) for k,v in convergence_study(out=a.out).items()})
        return
    if a.command == 'thermal':
        from .thermal_flow import run_thermal
        print('Dimensionless heat-release proxy:',run_thermal(a.out))
        return
    if a.command == 'benchmark':
        from .benchmark import run_benchmark
        print('Serial-only benchmark:',run_benchmark(a.out))
        return
    if a.command == 'cantera':
        from .cantera_demo import run_cantera
        d=run_cantera(out=a.out+'/cantera')
        print('Cantera chemistry (NO plasma); species candidate ranking:',d['candidate_species'])
        return
    if a.command in ('run','kinetics'):
        d=zero_d(a.out)
        print('0D max absolute errors [F0,B=F1+F2,R,P]:',d['max_abs_error_F0_B_R_P'])
        print('0D median runtime (s):',d['runtimes'])
    if a.command in ('run','flow'):
        d=spatial(a.out,a.quick)
        print('2D simulation diagnostics:',d)
    if a.command == 'run':
        from .report import write_report
        print('Previous-series report:',write_report(a.out))
        from .v4_report import generate_v4_report
        generate_v4_report(a.out)
        print('Version 4 computed report:',a.out+'/REPORT_V4.md')
    print('Figures and machine-readable data:',a.out)

if __name__=='__main__':main()
