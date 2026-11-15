from pathlib import Path
import csv
import numpy as np
import pytest
from scipy.integrate import trapezoid
from plasmarom.eedf import maxwellian_pdf
from plasmarom.eedf_field import read_field_manifest,interpolate_field_eedf,run_field_rate
from plasmarom.plasma_states import synthetic_table


def make_example(tmp_path):
    energy=np.linspace(0,80,1601)
    for field,te in [(50,1.),(200,3.)]:
        pdf=maxwellian_pdf(energy,te)
        pdf/=trapezoid(pdf,energy)
        with (tmp_path/f'e_{field}.csv').open('w',newline='') as out:
            w=csv.writer(out);w.writerow(['energy_eV','probability_per_eV']);w.writerows(zip(energy,pdf))
    path=tmp_path/'manifest.csv'
    path.write_text('EoverN_Td,eedf_csv,source\n50,e_50.csv,SYNTHETIC_50_Td\n200,e_200.csv,SYNTHETIC_200_Td\n')
    sig=synthetic_table(.5)
    with (tmp_path/'sigma.csv').open('w',newline='') as out:
        w=csv.writer(out);w.writerow(['energy_eV','sigma_m2']);w.writerows(sig)
    return path,tmp_path/'sigma.csv'


def test_interpolation_from_externally_supplied_field_grid(tmp_path):
    manifest,sigma=make_example(tmp_path)
    entries=read_field_manifest(manifest)
    mid,source=interpolate_field_eedf(entries,100.)
    assert len(source)==2
    assert trapezoid(mid[:,1],mid[:,0])==pytest.approx(1,abs=0.005)
    result=run_field_rate(tmp_path/'out',manifest,sigma,'SYNTHETIC_SIGMA',.5,100)
    assert result['rate_coefficient_m3_s']>0
    assert result['EoverN_Td']==100
    assert (tmp_path/'out'/'field_dependent_rate.json').exists()


def test_field_extrapolation_and_duplicate_rejected(tmp_path):
    manifest,sigma=make_example(tmp_path)
    entries=read_field_manifest(manifest)
    for field in (10,300,0,np.nan):
        with pytest.raises(ValueError,match='outside'):
            interpolate_field_eedf(entries,field)
    manifest.write_text('EoverN_Td,eedf_csv,source\n50,e_50.csv,SYNTHETIC\n50,e_200.csv,SYNTHETIC\n')
    with pytest.raises(ValueError,match='unique'):
        read_field_manifest(manifest)
