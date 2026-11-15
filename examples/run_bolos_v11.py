"""Provide YOUR sourced BOLSIG+ cross-section set and exact target mix."""
import argparse
from plasmarom.bolos_bridge import run_bolos
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cross-sections',required=True)
    p.add_argument('--source',required=True,help='Citation/DOI/URL to collision dataset')
    p.add_argument('--gas',required=True,help='One BOLOS target species, e.g., H2')
    p.add_argument('--field-Td',type=float,default=100.)
    args=p.parse_args()
    print(run_bolos(args.cross_sections,'results_v11/bolos',source_url=args.source,mixture={args.gas:1.},field_Td=args.field_Td))
