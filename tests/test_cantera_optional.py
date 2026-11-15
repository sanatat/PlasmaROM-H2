import json
from pathlib import Path
import pytest


def test_cantera_full_chemistry_optional(tmp_path):
    pytest.importorskip('cantera',minversion='3.2')
    from plasmarom.cantera_demo import run_cantera
    output=run_cantera(tmp_path,conditions=((1100.,1.),),t_end=0.00025,n_snapshots=8)
    assert output['n_species_full']>=8
    assert output['n_reactions_full']>10
    assert output['max_relative_H_atom_drift']<1e-5
    assert output['max_relative_O_atom_drift']<1e-5
    assert output['max_relative_enthalpy_drift']<1e-5
    assert Path(tmp_path/'states.csv').exists()
    assert len(output['candidate_species'])>0


def test_cantera_parameter_validation_without_optional_install(tmp_path):
    from plasmarom.cantera_demo import run_cantera
    for kw in ({'n_snapshots':1}, {'t_end':0}, {'threshold':1.5},
               {'conditions':()}, {'conditions':((0.,1.),)}):
        with pytest.raises(ValueError):
            run_cantera(tmp_path,**kw)


def test_cantera_executable_skeletal_optional(tmp_path):
    pytest.importorskip('cantera',minversion='3.2')
    from plasmarom.cantera_demo import run_skeletal
    result=run_skeletal(tmp_path,mechanism='h2o2.yaml',
                        training=((1100.,1.),),holdout=((1150.,1.),),
                        t_end=0.001,n_snapshots=15,thresholds=(0.4,0.1))
    assert result['status'] in ('NO REDUCTION','HOLDOUT_PASS','FAILED_OR_INCONCLUSIVE_HOLDOUT')
    assert (tmp_path/'skeletal_summary.json').exists()
    if result['status']!='NO REDUCTION':
        assert (tmp_path/'candidate_skeletal.yaml').exists()
    if result['status']=='HOLDOUT_PASS':
        assert result['acceptance_evidence']['accepted'] is True
        assert (tmp_path/'holdout_pass_skeletal.yaml').exists()
    else:
        assert result['acceptance_evidence']['accepted'] is False


def test_cantera_rejects_hash_matching_fake_yaml(tmp_path):
    """Integration: only an actual Cantera-loadable phase can pass the artifact gate."""
    pytest.importorskip('cantera',minversion='3.2')
    import hashlib
    from plasmarom.acceptance import assess_skeletal_summary
    artifact=tmp_path/'holdout_pass_skeletal.yaml'
    artifact.write_text('NOT CANTERA YAML',encoding='utf-8')
    fake={'source_mechanism':'h2o2.yaml', 'artifact':'holdout_pass_skeletal.yaml',
          'artifact_sha256':hashlib.sha256(artifact.read_bytes()).hexdigest()}
    # Its hash matches but its contents are NOT valid Cantera chemistry.
    assert hashlib.sha256((tmp_path/'holdout_pass_skeletal.yaml').read_bytes()).hexdigest()==fake['artifact_sha256']
    decision=assess_skeletal_summary(fake,tmp_path)
    assert decision['accepted'] is False
    assert 'invalid_or_unreadable_cantera_mechanism' in decision['reasons']


def test_v9_cantera_species_thermo_and_transport_fingerprint_roundtrip(tmp_path):
    ct=pytest.importorskip('cantera',minversion='3.2')
    from plasmarom.acceptance import _species_fingerprint
    parent=ct.Solution('h2o2.yaml')
    saved=tmp_path/'source_roundtrip.yaml'
    parent.write_yaml(str(saved))
    reloaded=ct.Solution(str(saved))
    assert _species_fingerprint(parent.species('H2O')) == _species_fingerprint(reloaded.species('H2O'))


def test_v9_cantera_reactor_accepts_distinct_mixtures(tmp_path):
    ct=pytest.importorskip('cantera',minversion='3.2')
    from plasmarom.cantera_demo import _reactor_curve
    gas=ct.Solution('h2o2.yaml')
    cases=[]
    for mixture in ('H2:2,O2:1,N2:4','H2:1,O2:1,N2:4'):
        case=_reactor_curve(ct,gas,1200.,1.,.0001,8,composition=mixture)
        assert case['T_K'].shape==(8,)
        cases.append(case)
    assert cases[0]['T_K'][0]==pytest.approx(1200.)
    assert cases[1]['T_K'][0]==pytest.approx(1200.)
