"""CLI: compare model predictions against independently sourced measurements."""
import argparse
from plasmarom.experiment_compare import evaluate_experiment
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--measurements',required=True)
    ap.add_argument('--predictions',required=True)
    ap.add_argument('--out',default='results_v11/experiment_comparison.json')
    a=ap.parse_args()
    print(evaluate_experiment(a.measurements,a.predictions,a.out))
