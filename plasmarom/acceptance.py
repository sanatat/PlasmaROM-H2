"""Fail-closed, independently recomputed evidence gate for optional Cantera runs.

This is NOT experimental validation. Independent JSON checks prevent accidental
acceptance of missing, contradictory, malformed or duplicated evidence. Actual
Cantera loading is mandatory for chemical-structure validation; hashes alone
cannot attest that a file is a chemically meaningful mechanism.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Mapping
from pathlib import Path


_TRAJECTORY_METRICS = (
    'T_K_normalized_max_abs_error',
    'OH_Y_normalized_max_abs_error',
    'H2O_Y_normalized_max_abs_error',
)


def _finite(value, *, strictly_positive=False):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and (value > 0 if strictly_positive else value >= 0))


def _integer(value, *, minimum=0):
    return type(value) is int and value >= minimum


DEFAULT_COMPOSITION = 'H2:2,O2:1,N2:4'


def _composition_signature(raw):
    """Canonical mole-fraction key: distinguish mixtures, not text ordering/scale.

    Reject malformed / negative / nonfinite / repeated species without depending
    on optional Cantera; its species-availability check happens during replay.
    """
    if not isinstance(raw,str) or not raw.strip():
        return None
    parts={}
    for token in raw.split(','):
        pair=token.strip().split(':')
        if len(pair)!=2 or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_+()\-]*',pair[0].strip()):
            return None
        sp=pair[0].strip()
        try:
            amount=float(pair[1].strip())
        except ValueError:
            return None
        if sp in parts or not math.isfinite(amount) or amount<=0:
            return None
        parts[sp]=amount
    total=sum(parts.values())
    if not math.isfinite(total) or total<=0:
        return None
    return tuple((sp,round(v/total,12)) for sp,v in sorted(parts.items()))


def _conditions(raw):
    """Parse cases [T_K, P_atm] or [T_K, P_atm, mole-fraction string].

    The third element permits explicit off-stoichiometric mixture evaluation.
    Legacy two-entry cases retain the documented H2/O2/N2 composition.
    """
    if not isinstance(raw,list) or not raw:
        return None
    parsed=[]
    for item in raw:
        if not isinstance(item,(list,tuple)) or len(item) not in (2,3):
            return None
        if not all(_finite(v,strictly_positive=True) for v in item[:2]):
            return None
        composition=DEFAULT_COMPOSITION if len(item)==2 else item[2]
        if _composition_signature(composition) is None:
            return None
        point=(float(item[0]),float(item[1]),str(composition).strip())
        if _overlap([point],parsed):
            return None
        parsed.append(point)
    return parsed


def _overlap(a,b):
    return any(math.isclose(x[0],y[0],rel_tol=1e-10,abs_tol=1e-12) and
               math.isclose(x[1],y[1],rel_tol=1e-10,abs_tol=1e-12) and
               _composition_signature(x[2] if len(x)>2 else DEFAULT_COMPOSITION)==
               _composition_signature(y[2] if len(y)>2 else DEFAULT_COMPOSITION)
               for x in a for y in b)


def _source_path(source, ct):
    """Locate the exact built-in or file-backed Cantera mechanism in use."""
    candidates=[Path(source)]
    try:
        candidates.extend(Path(d)/source for d in ct.get_data_directories())
    except Exception:
        pass
    for path in candidates:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError(f'Cannot establish source mechanism bytes: {source}')


def mechanism_sha256(source, ct):
    return hashlib.sha256(_source_path(source,ct).read_bytes()).hexdigest()


def _json_kinetic_data(x):
    """Canonicalize Cantera rate/YAML input maps without dropping kinetic data."""
    if isinstance(x,Mapping):
        return {str(k):_json_kinetic_data(v) for k,v in sorted(x.items(),key=lambda p:str(p[0]))}
    if isinstance(x,(list,tuple)):
        return [_json_kinetic_data(v) for v in x]
    if hasattr(x,'tolist'):
        return _json_kinetic_data(x.tolist())
    if isinstance(x,(int,float,str,bool)) or x is None:
        if isinstance(x,float) and not math.isfinite(x):
            raise ValueError('Non-finite reaction parameter')
        return x
    raise TypeError(f'Unknown Cantera parameter representation {type(x).__name__}')


def _reaction_fingerprint(reaction):
    """Compare all rate, stoichiometry, reversibility and collider fields.

    input_data contains Cantera's normalized YAML-equivalent rate description.
    Repeated identical reactions are counted using a multiset, not a set.
    """
    third=getattr(reaction,'third_body',None)
    fields={
        'equation':str(reaction.equation),
        'reactants':dict(reaction.reactants),
        'products':dict(reaction.products),
        'reversible':bool(reaction.reversible),
        'duplicate':bool(reaction.duplicate),
        'orders':dict(getattr(reaction,'orders',{})),
        'reaction_type':str(getattr(reaction,'reaction_type','')),
        'rate_input':reaction.input_data,
        'third_body':None if third is None else {
            'collider':str(getattr(third,'collider','M')),
            'efficiencies':dict(getattr(third,'efficiencies',{})),
            'default_efficiency':float(getattr(third,'default_efficiency',1.)),
        },
    }
    return json.dumps(_json_kinetic_data(fields),sort_keys=True,separators=(',',':'),allow_nan=False)


def _species_fingerprint(species):
    """Includes Cantera's full species thermo and transport declarations.

    Comparing only elemental compositions misses changed NASA polynomials,
    transport parameters, charge and temperature bounds. These identifiers are
    calculated independently for the source and written skeletal YAML.
    """
    metadata={
        'name':species.name,
        'composition':dict(species.composition),
        'charge':getattr(species,'charge',0),
        'size':getattr(species,'size',1),
        'input_data':species.input_data,
    }
    thermo=getattr(species,'thermo',None)
    if thermo is not None:
        metadata['thermo_coeffs']=getattr(thermo,'coeffs',None)
        metadata['thermo_temperature_range']=[getattr(thermo,'min_temp',None),
                                              getattr(thermo,'max_temp',None)]
    return json.dumps(_json_kinetic_data(metadata),sort_keys=True,
                      separators=(',',':'),allow_nan=False)


def _load_mechanism_metadata(artifact: Path, source: str):
    """Independent Cantera load and full species/reaction identity comparison."""
    import cantera as ct
    reduced=ct.Solution(str(artifact))
    full=ct.Solution(source)
    return {
        'full_species':full.n_species,
        'full_reactions':full.n_reactions,
        'reduced_species':reduced.n_species,
        'reduced_reactions':reduced.n_reactions,
        'source_sha256':mechanism_sha256(source,ct),
        'source_species_names':tuple(full.species_names),
        'reduced_species_names':tuple(reduced.species_names),
        'source_species_composition':{sp.name:dict(sp.composition) for sp in full.species()},
        'reduced_species_composition':{sp.name:dict(sp.composition) for sp in reduced.species()},
        'source_species_fingerprints':{sp.name:_species_fingerprint(sp) for sp in full.species()},
        'reduced_species_fingerprints':{sp.name:_species_fingerprint(sp) for sp in reduced.species()},
        'source_reaction_signatures':tuple(_reaction_fingerprint(r) for r in full.reactions()),
        'reduced_reaction_signatures':tuple(_reaction_fingerprint(r) for r in reduced.reactions()),
    }


def _replay_training_observations(source,conditions,t_end,snapshots):
    """Recompute each training reference independently, rather than trust JSON."""
    import cantera as ct
    from .cantera_demo import _reactor_curve
    gas=ct.Solution(source)
    observations=[]
    for T,P,composition in conditions:
        curve=_reactor_curve(ct,gas,T,P,t_end,snapshots,composition=composition)
        observations.append({
            'reference_ignition_observed':bool(curve['ignition_observed']),
            'threshold_crossing_time_s':curve['threshold_crossing_time_s'],
            'peak_dTdt_time_s':curve['peak_dTdt_time_s'],
            'max_time_step_s':curve['max_time_step_s'],
        })
    return observations


def _replay_holdout_metrics(artifact: Path,source: str,conditions,t_end,snapshots):
    """Re-simulate held-out reference and candidate curves independently."""
    import cantera as ct
    from .cantera_demo import _reactor_curve,_curve_error
    full=ct.Solution(source)
    reduced=ct.Solution(str(artifact))
    return [_curve_error(_reactor_curve(ct,full,T,P,t_end,snapshots,composition=comp),
                         _reactor_curve(ct,reduced,T,P,t_end,snapshots,composition=comp))
            for T,P,comp in conditions]


def _same_metric(a,b):
    if type(a) is bool or type(b) is bool:
        return type(a) is type(b) and a==b
    if isinstance(a,(int,float)) or isinstance(b,(int,float)):
        return _finite(a) and _finite(b) and math.isclose(a,b,rel_tol=1e-5,abs_tol=1e-10)
    return a==b


def assess_skeletal_summary(summary, base_dir=None, *, require_artifact=True, require_pinned=False):
    """Strictly recompute acceptance from per-case evidence and a Cantera phase.

    Strict v10 records require fixed protocol and replay of training and holdout
    candidate and reference trajectories; this is NOT experimental validation.
    Legacy direct-call fixtures remain supported for historical unit tests only;
    CLI/file acceptance requires a pinned v10 report.
    A missing Cantera installation necessarily produces `accepted=False`.
    """
    reasons = []
    if not isinstance(summary, dict):
        return {'accepted':False,'reasons':['not_a_json_object']}

    def fail(code):
        if code not in reasons:
            reasons.append(code)

    counts = {}
    for name in ('full_species','reduced_species','full_reactions','reduced_reactions'):
        value = summary.get(name)
        if not _integer(value, minimum=1):
            fail('invalid_'+name)
        else:
            counts[name] = value
    if len(counts) == 4 and not (counts['reduced_species'] < counts['full_species'] and
                                 counts['reduced_reactions'] < counts['full_reactions']):
        fail('no_nontrivial_mechanism_reduction')

    if summary.get('status') != 'HOLDOUT_PASS':
        fail('status_not_holdout_pass')
    v10=summary.get('evidence_schema')=='v10-pinned'
    if require_pinned and not v10:
        fail('pinned_v10_protocol_required_for_scientific_claim')
    if summary.get('evidence_schema','').startswith('v10'):
        if not v10:
            fail('unlocked_exploratory_protocol_cannot_be_certified')
        else:
            try:
                from .protocol import load_pinned_protocol,verify_summary_protocol
                declared_protocol=load_pinned_protocol(
                    summary.get('protocol_id',''))
                verify_summary_protocol(declared_protocol,summary)
            except Exception:
                fail('protocol_not_pinned_or_mismatched')
    tol = summary.get('holdout_tolerance')
    if not _finite(tol, strictly_positive=True) or tol >= 1:
        fail('invalid_holdout_tolerance')
        tol = None

    t_end=summary.get('reactor_t_end_s')
    snapshots=summary.get('reactor_n_snapshots')
    if not _finite(t_end,strictly_positive=True) or not _integer(snapshots,minimum=6):
        fail('invalid_replay_protocol')
    # Mixed-composition scenarios are allowed, but their declaration must agree
    # with every per-case mixture. Legacy two-field cases have fixed composition.
    declared_composition=summary.get('initial_composition')
    if declared_composition not in (DEFAULT_COMPOSITION,'scenario_specific'):
        fail('unknown_initial_composition')
    if summary.get('ignition_definition')!='peak_dTdt_if_threshold_resolved_with_20K_rise':
        fail('unknown_ignition_definition')

    training = _conditions(summary.get('training_conditions'))
    holdout = _conditions(summary.get('holdout_conditions'))
    if training is None:
        fail('invalid_training_conditions')
    if holdout is None:
        fail('invalid_holdout_conditions')
    if training is not None and holdout is not None and _overlap(training,holdout):
        fail('training_holdout_overlap')
    if training is not None and holdout is not None:
        varies=any(_composition_signature(x[2])!=_composition_signature(DEFAULT_COMPOSITION)
                   for x in training+holdout)
        expected='scenario_specific' if varies else DEFAULT_COMPOSITION
        if declared_composition != expected:
            fail('composition_protocol_mismatch')

    observed_training = 0
    records = summary.get('training_ignition_observations')
    if not isinstance(records, list) or training is None or len(records) != len(training):
        fail('missing_training_ignition_evidence')
    else:
        for condition, record in zip(training, records):
            if not isinstance(record, dict) or record.get('reference_ignition_observed') not in (True,False) or type(record.get('reference_ignition_observed')) is not bool:
                fail('invalid_training_ignition_record')
                continue
            case = _conditions([record.get('condition_T_K_P_atm')])
            if case is None or not _overlap(case, [condition]):
                fail('training_condition_record_mismatch')
            observed_training += int(record['reference_ignition_observed'])
    if not _integer(summary.get('ignition_training_cases_observed')) or summary.get('ignition_training_cases_observed') != observed_training or observed_training < 1:
        fail('inconsistent_training_ignition_count')

    metrics = summary.get('holdout_metrics')
    observed_holdout=0
    all_match=True
    all_errors=[]
    valid_case_set=(isinstance(metrics,list) and holdout is not None and
                    len(metrics)==len(holdout))
    if not valid_case_set:
        fail('missing_per_case_holdout_evidence')
    else:
        for condition, row in zip(holdout, metrics):
            if not isinstance(row, dict):
                fail('invalid_holdout_case');all_match=False;continue
            case=_conditions([row.get('condition_T_K_P_atm')])
            if case is None or not _overlap(case,[condition]):
                fail('holdout_condition_record_mismatch')
            ref=row.get('reference_ignition_observed')
            cand=row.get('candidate_ignition_observed')
            matched=row.get('ignition_classification_matches')
            evaluable=row.get('ignition_error_evaluable')
            if not all(type(v) is bool for v in (ref,cand,matched,evaluable)):
                fail('invalid_ignition_case_status');all_match=False
            else:
                actual_match = ref == cand
                actual_evaluable = ref and cand
                if matched is not actual_match or evaluable is not actual_evaluable:
                    fail('contradictory_ignition_classification')
                if not actual_match:
                    all_match=False
                    fail('ignition_classification_mismatch')
                observed_holdout += int(ref)
            ign_err = row.get('ignition_proxy_rel_error')
            if evaluable is True:
                if not _finite(ign_err):
                    fail('invalid_ignition_error')
                else:
                    all_errors.append(float(ign_err))
            elif ign_err is not None:
                # A censored or mismatched observation must not sneak in a
                # fabricated ignition proxy error.
                if not (ref is not cand and ign_err == 1):
                    fail('non_evaluable_ignition_error_present')
                else:
                    all_errors.append(1.0)
            diagnostic_metrics=(_TRAJECTORY_METRICS+
                                (('OH_Y_trace_relative_max_error',) if v10 else ()))
            for key in diagnostic_metrics:
                val=row.get(key)
                if not _finite(val):
                    fail('invalid_'+key)
                else:
                    all_errors.append(float(val))

    training_metrics=summary.get('training_metrics')
    if v10:
        if (not isinstance(training_metrics,list) or training is None or
                len(training_metrics)!=len(training)):
            fail('missing_training_candidate_evidence')
        else:
            for condition,record in zip(training,training_metrics):
                if not isinstance(record,dict):
                    fail('invalid_training_candidate_record')
                    continue
                case=_conditions([record.get('condition_T_K_P_atm')])
                if case is None or not _overlap(case,[condition]):
                    fail('training_candidate_condition_mismatch')
                for metric in (*_TRAJECTORY_METRICS,'OH_Y_trace_relative_max_error'):
                    if not _finite(record.get(metric)):
                        fail('invalid_training_'+metric)
            vals=[row.get(k) for row in training_metrics if isinstance(row,dict)
                  for k in (*_TRAJECTORY_METRICS,'OH_Y_trace_relative_max_error')]
            claimed=summary.get('training_max_normalized_error')
            if (not vals or not all(_finite(x) for x in vals) or
                    not _finite(claimed) or
                    not math.isclose(float(claimed),max(vals),rel_tol=1e-9,abs_tol=1e-12)):
                fail('training_aggregate_error_inconsistent')
            elif (not _finite(summary.get('training_tolerance'),strictly_positive=True)
                  or claimed>summary['training_tolerance']):
                fail('training_error_exceeds_tolerance')

    declared_count=summary.get('ignition_holdout_cases_observed')
    if not _integer(declared_count) or declared_count != observed_holdout or observed_holdout < 1:
        fail('inconsistent_holdout_ignition_count')
    if type(summary.get('holdout_ignition_classifications_match')) is not bool or summary.get('holdout_ignition_classifications_match') != all_match:
        fail('inconsistent_aggregate_ignition_status')

    claimed_err=summary.get('holdout_max_normalized_error')
    actual_err=max(all_errors) if all_errors else None
    if (not _finite(claimed_err) or actual_err is None or
        not math.isclose(claimed_err,actual_err,rel_tol=1e-9,abs_tol=1e-12)):
        fail('aggregate_error_inconsistent_with_cases')
    if tol is not None and (actual_err is None or actual_err > tol):
        fail('holdout_error_exceeds_tolerance')

    # The flag may not supply its own evidence; it must match case-level truth.
    evidence=(observed_training > 0 and observed_holdout > 0 and all_match and
              actual_err is not None and tol is not None and actual_err <= tol and
              training is not None and holdout is not None and not _overlap(training,holdout))
    if summary.get('accepted_with_ignition_evidence') is not evidence:
        fail('inconsistent_claimed_acceptance')
    if not evidence:
        fail('insufficient_ignition_and_holdout_evidence')

    if not require_artifact:
        fail('artifact_validation_required')
    else:
        artifact=summary.get('artifact')
        if artifact != 'holdout_pass_skeletal.yaml':
            fail('missing_validated_artifact_name')
        elif base_dir is None or not (Path(base_dir)/artifact).is_file():
            fail('missing_holdout_pass_skeletal_yaml')
        else:
            artifact_path=Path(base_dir)/artifact
            sha=summary.get('artifact_sha256')
            if not isinstance(sha,str) or not re.fullmatch('[a-f0-9]{64}',sha):
                fail('invalid_skeletal_yaml_hash')
            elif hashlib.sha256(artifact_path.read_bytes()).hexdigest()!=sha:
                fail('skeletal_yaml_hash_mismatch')
            source=summary.get('source_mechanism')
            if not isinstance(source,str) or not source.strip():
                fail('missing_source_mechanism')
            else:
                try:
                    metadata=_load_mechanism_metadata(artifact_path,source)
                except (ImportError,ModuleNotFoundError):
                    fail('cantera_unavailable_cannot_validate_mechanism')
                except Exception:
                    fail('invalid_or_unreadable_cantera_mechanism')
                else:
                    for field in counts:
                        if metadata.get(field) != counts[field]:
                            fail('cantera_'+field+'_count_mismatch')
                    source_names=set(metadata.get('source_species_names',()))
                    reduced_names=set(metadata.get('reduced_species_names',()))
                    if not reduced_names or not reduced_names <= source_names:
                        fail('reduced_species_not_subset_of_source')
                    srccomp=metadata.get('source_species_composition',{})
                    redcomp=metadata.get('reduced_species_composition',{})
                    if any(srccomp.get(sp) != redcomp.get(sp) for sp in reduced_names):
                        fail('modified_species_composition')
                    srcfp=metadata.get('source_species_fingerprints')
                    redfp=metadata.get('reduced_species_fingerprints')
                    if not isinstance(srcfp,dict) or not isinstance(redfp,dict):
                        fail('missing_full_species_identity_evidence')
                    elif any(srcfp.get(sp) != redfp.get(sp) or sp not in srcfp
                             for sp in reduced_names):
                        fail('modified_species_thermo_or_transport')
                    if summary.get('source_mechanism_sha256') != metadata.get('source_sha256'):
                        fail('source_mechanism_hash_mismatch')
                    src_rxns=metadata.get('source_reaction_signatures')
                    red_rxns=metadata.get('reduced_reaction_signatures')
                    if src_rxns is None or red_rxns is None:
                        fail('missing_reaction_identity_evidence')
                    elif not Counter(red_rxns) <= Counter(src_rxns):
                        fail('modified_or_foreign_reaction_parameters')
                    if (training is not None and isinstance(records,list)
                            and len(records)==len(training) and
                            _finite(t_end,strictly_positive=True) and
                            _integer(snapshots,minimum=6)):
                        try:
                            training_replay=_replay_training_observations(
                                source,training,t_end,snapshots)
                        except Exception:
                            fail('independent_training_replay_failed')
                        else:
                            if len(training_replay)!=len(records):
                                fail('independent_training_replay_case_count_mismatch')
                            else:
                                required=('reference_ignition_observed',
                                          'threshold_crossing_time_s',
                                          'peak_dTdt_time_s','max_time_step_s')
                                for record, recomputed in zip(records,training_replay):
                                    if any(key not in record or key not in recomputed or
                                           not _same_metric(record[key],recomputed[key])
                                           for key in required):
                                        fail('independent_training_replay_mismatch')
                    if (v10 and training is not None and
                            isinstance(training_metrics,list) and
                            len(training_metrics)==len(training) and
                            _finite(t_end,strictly_positive=True) and
                            _integer(snapshots,minimum=6)):
                        try:
                            train_replayed=_replay_holdout_metrics(
                                artifact_path,source,training,t_end,snapshots)
                        except Exception:
                            fail('independent_training_candidate_replay_failed')
                        else:
                            if len(train_replayed)!=len(training_metrics):
                                fail('independent_training_candidate_replay_case_count_mismatch')
                            else:
                                for record,computed in zip(training_metrics,train_replayed):
                                    if any(k not in record or not _same_metric(record[k],v)
                                           for k,v in computed.items()):
                                        fail('independent_training_candidate_replay_mismatch')
                    if (holdout is not None and _finite(t_end,strictly_positive=True)
                            and _integer(snapshots,minimum=6) and valid_case_set):
                        try:
                            replay=_replay_holdout_metrics(artifact_path,source,holdout,t_end,snapshots)
                        except Exception:
                            fail('independent_reactor_replay_failed')
                        else:
                            if len(replay)!=len(metrics):
                                fail('independent_replay_case_count_mismatch')
                            else:
                                for row, recomputed in zip(metrics,replay):
                                    if any(k not in row or not _same_metric(row[k],v)
                                           for k,v in recomputed.items()):
                                        fail('independent_reactor_replay_mismatch')
                                    # Every holdout is required to expose both ignition clocks
                                    required=('reference_threshold_crossing_time_s',
                                              'candidate_threshold_crossing_time_s',
                                              'reference_peak_dTdt_time_s',
                                              'candidate_peak_dTdt_time_s',
                                              'reference_max_time_step_s',
                                              'candidate_max_time_step_s')
                                    if any(k not in row for k in required):
                                        fail('missing_explicit_ignition_clock_evidence')
    return {'accepted':not reasons,'reasons':reasons}


def assess_skeletal_file(path, output=None):
    source=Path(path)
    data=json.loads(source.read_text(encoding='utf-8'))
    report=assess_skeletal_summary(data,source.parent,require_pinned=True)
    if output is not None:
        Path(output).write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report
