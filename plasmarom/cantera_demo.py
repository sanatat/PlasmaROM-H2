"""OPTIONAL true H2/O2 chemistry demo using Cantera's supplied h2o2.yaml.

Run only if installed with: pip install 'cantera>=3.1,<4'
This code does not model plasma; it exercises thermochemistry, ignition, and
real chemical DRGEP ranking. No publication-level validation is claimed.
"""
import csv
import json
from pathlib import Path
import numpy as np
from .reduction import (multi_state_importance, choose_species,
                        gross_flux_interaction_coefficients,propagate_importance)


def _ct():
    try:
        import cantera as ct
    except ImportError as exc:
        raise RuntimeError('Cantera unavailable. Install: pip install "cantera>=3.1,<4"') from exc
    return ct


def run_cantera(out='results/cantera', conditions=((1000.,1.),(1100.,1.),(1200.,1.)),
                t_end=0.005, n_snapshots=81, threshold=0.08):
    """Adiabatic constant-pressure reactor + chemical DRGEP reference-state ranking.

    P supplied in atmospheres; T in kelvin; t in seconds. Output selection is
    only a candidate roster, NOT an executable validated skeletal mechanism.
    """
    # Validate configuration BEFORE importing the optional dependency.
    if n_snapshots < 2 or t_end <= 0 or not 0 <= threshold <= 1:
        raise ValueError('Need >=2 snapshots, t_end>0, and 0<=threshold<=1')
    if not conditions or any(not (np.isfinite(t) and np.isfinite(p) and t>0 and p>0)
                             for t,p in conditions):
        raise ValueError('Conditions must contain positive finite T [K] and P [atm]')
    ct=_ct()
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    gas=ct.Solution('h2o2.yaml')
    names=gas.species_names
    nu=gas.product_stoich_coeffs - gas.reactant_stoich_coeffs
    participation=_reaction_participation_matrix(gas)
    targets=[gas.species_index(s) for s in ('H2','O2','OH','H2O')]
    recorded_rates=[]
    gross_importance=[]
    rows=[]
    max_H_drift=0.
    max_O_drift=0.
    max_energy_drift=0.
    for T,press_atm in conditions:
        gas.TPX=float(T),press_atm*ct.one_atm,'H2:2,O2:1,N2:4'
        reactor=ct.IdealGasConstPressureReactor(gas,clone=True)
        net=ct.ReactorNet([reactor]);net.rtol=1e-9;net.atol=1e-17
        initial=reactor.thermo
        y_initial=initial.Y.copy()
        atomic_h=np.array([gas.n_atoms(k,'H') for k in range(gas.n_species)])
        atomic_o=np.array([gas.n_atoms(k,'O') for k in range(gas.n_species)])
        molecular_weights=gas.molecular_weights
        # kmol atom/kg gas: mass fractions -> kmol/kg -> atom abundance
        atoms_H0=float(np.dot(y_initial/molecular_weights, atomic_h))
        atoms_O0=float(np.dot(y_initial/molecular_weights, atomic_o))
        h0=float(initial.enthalpy_mass)
        for t in np.linspace(0,t_end,n_snapshots):
            net.advance(float(t))
            state=reactor.thermo
            rates=np.asarray(state.net_rates_of_progress, dtype=float)
            recorded_rates.append(rates)
            gross_importance.append(propagate_importance(
                gross_flux_interaction_coefficients(nu,state.forward_rates_of_progress,
                                                    state.reverse_rates_of_progress,participation),targets))
            y=state.Y.copy()
            a_H=float(np.dot(y/molecular_weights, atomic_h))
            a_O=float(np.dot(y/molecular_weights, atomic_o))
            max_H_drift=max(max_H_drift,abs(a_H-atoms_H0)/max(atoms_H0,1e-30))
            max_O_drift=max(max_O_drift,abs(a_O-atoms_O0)/max(atoms_O0,1e-30))
            max_energy_drift=max(max_energy_drift,abs(state.enthalpy_mass-h0)/max(abs(h0),1.))
            rows.append({'initial_T_K':T,'P_atm':press_atm,'time_s':t,
                         'T_K':state.T,'H2_mass_fraction':y[gas.species_index('H2')],
                         'OH_mass_fraction':y[gas.species_index('OH')],
                         'H2O_mass_fraction':y[gas.species_index('H2O')]})
    importance=multi_state_importance(nu,np.asarray(recorded_rates),participation,targets)
    mandatory=('H2','O2','OH','H2O','N2')
    selected=choose_species(importance,names,threshold=threshold,required=mandatory)
    candidate_reactions=sum(all(names[k] in selected for k in np.flatnonzero(participation[:,i]))
                            for i in range(gas.n_reactions))
    with (out/'states.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={'status':'true chemical kinetics, NO plasma','mechanism':'Cantera h2o2.yaml',
             'n_species_full':gas.n_species,'n_reactions_full':gas.n_reactions,
             'reference_conditions_T_K_P_atm':[list(q) for q in conditions],
             'drgep_threshold':threshold,'drgep_targets':['H2','O2','OH','H2O'],
             'species_importance':dict(zip(names,importance.tolist())),
             'gross_flux_diagnostic_importance':dict(zip(names,np.max(gross_importance,axis=0).tolist())),
             'gross_flux_note':'Conservative collider-aware diagnostic; selection is still a species roster',
             'candidate_species':selected,'candidate_reactions_by_simple_membership':candidate_reactions,
             'max_relative_H_atom_drift':max_H_drift,'max_relative_O_atom_drift':max_O_drift,
             'max_relative_enthalpy_drift':max_energy_drift,
             'warning':'NO plasma-specific energy targets, no validated skeletal reactor and no DNS/LES. Reaction membership count conservatively includes colliders.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    return summary


def _reaction_participation_matrix(phase):
    """Stoichiometric and collider participation shared by net and gross screens.

    Explicit collision efficiencies, collider names and generic M collisions
    participate in the kinetics even though stoichiometry alone omits them.
    For generic M with nonzero default efficiency, mark every available
    species: conservative screening can prevent reduction, but cannot
    misrepresent a third-body-dependent rate as collider-independent.
    """
    participation=np.asarray((phase.product_stoich_coeffs +
                               phase.reactant_stoich_coeffs)>0, dtype=bool).copy()
    species_names=tuple(phase.species_names)
    for index,reaction in enumerate(phase.reactions()):
        mentioned=set(getattr(reaction,'orders',{}))
        third=getattr(reaction,'third_body',None)
        if third is not None:
            collider=str(getattr(third,'collider','M')).strip('()+')
            efficiencies=getattr(third,'efficiencies',{}) or {}
            mentioned.update(efficiencies)
            if collider not in ('M',''):
                mentioned.add(collider)
            elif float(getattr(third,'default_efficiency',1.))!=0.:
                # Every species is an implicit collider for M (conservative).
                mentioned.update(species_names)
        for name in mentioned:
            if name in species_names:
                participation[phase.species_index(name),index]=True
    return participation


def _safe_reaction_subset(phase, species_names):
    """Conservatively keep reactions only if all declared participants exist.

    Protect third-body efficiencies/colliders, not only stoichiometric species.
    This may reject a valid pruned reaction and reduce less aggressively; it
    must never silently change its collision chemistry.
    """
    selected=set(species_names)
    result=[]
    for reaction in phase.reactions():
        mentioned=set(reaction.reactants)|set(reaction.products)
        mentioned.update(getattr(reaction,'orders',{}))
        third=getattr(reaction,'third_body',None)
        if third is not None:
            collider=getattr(third,'collider','M')
            if collider not in ('M','(+M)'):
                mentioned.add(collider.strip('()+'))
            mentioned.update(getattr(third,'efficiencies',{}))
        if mentioned <= selected:
            result.append(reaction)
    return result


def _ignition_assessment(time_s, temperature_K, *, rise_K=20., endpoint_margin=2):
    """Disclose BOTH 20 K threshold time and peak dT/dt ignition proxy.

    The two definitions need not agree: peak dT/dt can occur earlier than
    threshold crossing. Threshold crossing is interpolated linearly between
    samples. No physical ignition is inferred from a short/censored window.
    """
    t = np.asarray(time_s,dtype=float)
    T = np.asarray(temperature_K,dtype=float)
    if (t.ndim!=1 or T.shape!=t.shape or t.size<6 or
        not np.isfinite(t).all() or not np.isfinite(T).all() or
        not np.all(np.diff(t)>0) or not np.isfinite(rise_K) or rise_K<=0 or
        type(endpoint_margin) is not int or endpoint_margin<1):
        raise ValueError('Invalid ignition time/temperature trajectory')
    delta=T-T[0]
    cross=np.flatnonzero(delta>=rise_K)
    grad=np.gradient(T,t)
    peak=int(np.argmax(grad))
    peak_time=float(t[peak]) if 0<peak<len(t)-1 and grad[peak]>0 else None
    crossing_time=None
    first_index=None
    if len(cross):
        first_index=int(cross[0])
        if first_index==0:
            crossing_time=float(t[0])
        else:
            j=first_index
            y0,y1=float(delta[j-1]),float(delta[j])
            frac=(rise_K-y0)/(y1-y0) if y1!=y0 else 1.
            crossing_time=float(t[j-1]+np.clip(frac,0,1)*(t[j]-t[j-1]))
    result={
        'threshold_crossing_time_s':crossing_time,
        'peak_dTdt_time_s':peak_time,
        'peak_dTdt_K_s':float(grad[peak]),
        'max_time_step_s':float(np.max(np.diff(t))),
        'ignition_time_definition':'peak_dTdt_if_threshold_resolved',
        'temperature_threshold_rise_K':float(rise_K),
        'temperature_rise_K':float(np.max(delta)),
        'ignition_observed':False,'ignition_time_s':None,
        'ignition_status':'right_censored_no_threshold_crossing',
    }
    if first_index is None:
        return result
    if first_index>=len(t)-endpoint_margin:
        result['ignition_status']='right_censored_crossing_near_end'
    elif peak_time is None:
        result['ignition_status']='unresolved_temperature_rise'
    else:
        result['ignition_observed']=True
        result['ignition_time_s']=peak_time  # backward-compatible explicit alias
        result['ignition_status']='observed_numerical_temperature_rise'
    return result


def _reactor_curve(ct, gas, temperature, pressure_atm, t_end, snapshots, *, composition='H2:2,O2:1,N2:4'):
    """Independent constant-pressure, adiabatic homogeneous reactor case."""
    gas.TPX=temperature,pressure_atm*ct.one_atm,composition
    reactor=ct.IdealGasConstPressureReactor(gas,clone=True)
    net=ct.ReactorNet([reactor]);net.rtol=1e-9;net.atol=1e-17
    tt=np.linspace(0,t_end,snapshots)
    fields={'T_K':[],'OH_Y':[],'H2O_Y':[]}
    for t in tt:
        net.advance(float(t))
        g=reactor.thermo
        fields['T_K'].append(g.T)
        for sp,col in [('OH','OH_Y'),('H2O','H2O_Y')]:
            fields[col].append(float(g.Y[g.species_index(sp)]))
    result={'time_s':tt}
    result.update({name:np.array(values) for name,values in fields.items()})
    result.update(_ignition_assessment(tt,result['T_K']))
    return result


def _curve_error(ref, candidate):
    """Compare curves AND disclose the two distinct numerical ignition clocks."""
    metrics={}
    for name,scale in [('T_K',max(20.,float(max(ref['T_K'])-ref['T_K'][0]))),
                       ('OH_Y',max(1e-7,float(max(ref['OH_Y'])))),
                       ('H2O_Y',max(1e-7,float(max(ref['H2O_Y']))))]:
        metrics[name+'_normalized_max_abs_error']=float(max(abs(ref[name]-candidate[name]))/scale)
    # In addition to the legacy 1e-7 normalization, expose OH errors at
    # trace-radical levels. A candidate cannot hide a 5000-fold spurious OH
    # increase merely because the reference concentration is below 1e-7.
    r_oh=np.asarray(ref['OH_Y'],dtype=float)
    c_oh=np.asarray(candidate['OH_Y'],dtype=float)
    floor=1.e-12  # MASS FRACTION; fixed by the versioned protocol
    difference=np.abs(r_oh-c_oh)
    active=(r_oh>=floor)|(c_oh>=floor)
    metrics['OH_Y_max_abs_difference']=float(np.max(difference))
    metrics['OH_Y_activity_threshold']=floor
    metrics['OH_Y_active_points']=int(np.count_nonzero(active))
    metrics['OH_Y_trace_relative_max_error']=(float(np.max(
        difference[active]/np.maximum(np.abs(r_oh[active]),floor)))
        if np.any(active) else 0.)
    obs_ref=bool(ref['ignition_observed'])
    obs_red=bool(candidate['ignition_observed'])
    metrics['reference_ignition_observed']=obs_ref
    metrics['candidate_ignition_observed']=obs_red
    metrics['ignition_classification_matches']=obs_ref==obs_red
    metrics['ignition_error_evaluable']=obs_ref and obs_red
    metrics['ignition_time_definition']='peak_dTdt_if_threshold_resolved'
    metrics['temperature_threshold_rise_K']=float(ref['temperature_threshold_rise_K'])
    for prefix,curve in (('reference',ref),('candidate',candidate)):
        for label in ('threshold_crossing_time_s','peak_dTdt_time_s','max_time_step_s'):
            metrics[f'{prefix}_{label}']=curve[label]
    if obs_ref and obs_red:
        metrics['ignition_proxy_rel_error']=float(abs(ref['ignition_time_s']-candidate['ignition_time_s'])/
                                                  max(ref['ignition_time_s'],1e-12))
    elif obs_ref != obs_red:
        metrics['ignition_proxy_rel_error']=float(1.)
    else:
        metrics['ignition_proxy_rel_error']=None
    return metrics


def _max_curve_error(metrics):
    """True errors only; absent ignition comparisons cannot imply success."""
    return max(float(v) for k,v in metrics.items()
               if k.endswith('_error') and v is not None)


def _finalize_scientific_acceptance(result, out):
    """Downgrade any provisional curve-level pass rejected by independent gate.

    Safe to test without Cantera by mocking the evidence gate; production uses
    `assess_skeletal_summary` to reload the saved mechanism with Cantera.
    """
    from .acceptance import assess_skeletal_summary
    result['acceptance_evidence'] = assess_skeletal_summary(result, out, require_pinned=True)
    if result['status'] == 'HOLDOUT_PASS' and not result['acceptance_evidence']['accepted']:
        result['scientific_acceptance_failures'] = result['acceptance_evidence']['reasons']
        result['status'] = 'FAILED_OR_INCONCLUSIVE_HOLDOUT'
        result['accepted_with_ignition_evidence'] = False
        result['artifact'] = 'candidate_skeletal.yaml'
        result.pop('artifact_sha256', None)
        (Path(out)/'holdout_pass_skeletal.yaml').unlink(missing_ok=True)
    return result


def run_skeletal(out='results/cantera_reduction',mechanism='h2o2.yaml',
                 training=((1050.,1.),(1150.,1.),(1250.,1.)),
                 holdout=((1100.,1.),(1300.,1.)),
                 t_end=.006,n_snapshots=80,thresholds=(.4,.2,.1,.05),
                 train_tolerance=.10,holdout_tolerance=.12, protocol_path=None):
    """Select, construct, simulate, and TEST a genuine chemical submechanism.

    Optional Cantera. DRGEP importance is estimated using positive and signed
    rates from reference trajectories. Only pre-registered training cases are
    used for selection; holdout cases NEVER influence candidate choice.
    Results can be NO REDUCTION or FAILED HOLDOUT, reported as such.
    No plasma modelling or experimental validation is implied.
    """
    if (not training or not holdout or n_snapshots<6 or not np.isfinite(t_end) or t_end<=0 or not thresholds
        or not 0<train_tolerance<=1 or not 0<holdout_tolerance<=1 or
        any(not np.isfinite(x) or not 0<=x<=1 for x in thresholds)):
        raise ValueError('Invalid train/holdout configuration')
    # Reject duplicate cases BEFORE importing optional Cantera or writing files.
    # A case is (T_K, P_atm [, gas composition string]).
    from .acceptance import _conditions, _overlap
    train_pairs=_conditions([list(x) for x in training])
    holdout_pairs=_conditions([list(x) for x in holdout])
    if train_pairs is None or holdout_pairs is None or _overlap(train_pairs,holdout_pairs):
        raise ValueError('Training and holdout conditions must be disjoint and unique')
    protocol_info=None
    if protocol_path is not None:
        from .protocol import load_pinned_protocol, verify_protocol_parameters
        protocol_info=load_pinned_protocol(protocol_path)
        verify_protocol_parameters(protocol_info, mechanism, train_pairs, holdout_pairs,
                                   t_end, n_snapshots, thresholds, train_tolerance,
                                   holdout_tolerance)
    ct=_ct()
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    # Discard stale YAMLs from earlier passes: results must match THIS run.
    for old in ('candidate_skeletal.yaml','holdout_pass_skeletal.yaml'):
        (out/old).unlink(missing_ok=True)
    gas=ct.Solution(mechanism)
    from .acceptance import mechanism_sha256
    source_digest=mechanism_sha256(mechanism,ct)
    # Preserve *every* input species, including a fuel or diluent added in a
    # composition sweep. Omitting an initial species would invalidate the test.
    from .acceptance import _composition_signature
    required_base=('H2','O2','H','O','OH','H2O','N2')
    initial_species={sp for case in train_pairs+holdout_pairs
                     for sp,_ in _composition_signature(case[2])}
    mandatory=tuple(dict.fromkeys([*required_base,*sorted(initial_species)]))
    if not set(mandatory).issubset(gas.species_names):
        missing=sorted(set(mandatory)-set(gas.species_names))
        raise ValueError(f'Mechanism is missing mandatory or supplied mixture species: {missing}')
    train_refs={};net_rates=[];gross_importance=[]
    nu_dummy=gas.product_stoich_coeffs-gas.reactant_stoich_coeffs
    participation_dummy=_reaction_participation_matrix(gas)
    targets_dummy=[gas.species_index(sp) for sp in ('H2','O2','OH','H2O')]
    for case in train_pairs:
        T,P,composition=case
        gas.TPX=T,P*ct.one_atm,composition
        reactor=ct.IdealGasConstPressureReactor(gas,clone=True)
        net=ct.ReactorNet([reactor]);net.rtol=1e-9;net.atol=1e-17
        times=np.linspace(0,t_end,n_snapshots)
        history={'time_s':times,'T_K':[],'OH_Y':[],'H2O_Y':[]}
        for t in times:
            net.advance(float(t))
            local=reactor.thermo
            history['T_K'].append(local.T)
            history['OH_Y'].append(float(local.Y[local.species_index('OH')]))
            history['H2O_Y'].append(float(local.Y[local.species_index('H2O')]))
            net_rates.append(np.array(local.net_rates_of_progress))
            gross_importance.append(propagate_importance(
                gross_flux_interaction_coefficients(nu_dummy,local.forward_rates_of_progress,
                                                    local.reverse_rates_of_progress,participation_dummy),targets_dummy))
        history={k:np.asarray(v,float) for k,v in history.items()}
        history.update(_ignition_assessment(times, history['T_K']))
        train_refs[case]=history
    nu=gas.product_stoich_coeffs-gas.reactant_stoich_coeffs
    participation=_reaction_participation_matrix(gas)
    targets=[gas.species_index(sp) for sp in ('H2','O2','OH','H2O')]
    net_importance=multi_state_importance(nu,np.array(net_rates),participation,targets)
    gross_importance_arr=np.max(gross_importance,axis=0)
    # Fail-conservative: high gross bidirectional flux prevents omission even
    # near equilibrium when net rates accidentally cancel.
    importance=np.maximum(net_importance,gross_importance_arr)
    n_igniting_train=sum(bool(v['ignition_observed']) for v in train_refs.values())
    tests=[];chosen=None
    for th in thresholds:
        roster=choose_species(importance,gas.species_names,th,required=mandatory)
        reactions=_safe_reaction_subset(gas,roster)
        attempt={'threshold':th,'species':len(roster),'reactions':len(reactions)}
        if len(roster)>=gas.n_species or not reactions:
            attempt['status']='not reduced';tests.append(attempt);continue
        try:
            skel=ct.Solution(thermo='ideal-gas',kinetics='gas',
                             species=[gas.species(s) for s in roster],reactions=reactions,
                             name='h2_skeletal_candidate')
            errors=[]
            for T,P,composition in train_pairs:
                curve=_reactor_curve(ct,skel,T,P,t_end,n_snapshots,composition=composition)
                errors.append(_curve_error(train_refs[T,P,composition],curve))
            worst=max(_max_curve_error(row) for row in errors)
            attempt['train_max_normalized_error']=worst
            attempt['ignition_train_cases_observed']=n_igniting_train
            attempt['status']='train_pass' if worst<=train_tolerance and n_igniting_train>0 else 'train_inconclusive_or_fail'
            if worst<=train_tolerance and n_igniting_train>0 and (chosen is None or len(roster)<chosen[0]):
                chosen=(len(roster),skel,attempt,errors)
        except Exception as ex:
            attempt['status']='construction_or_integration_failed'
            attempt['error']=str(ex)[:300]
        tests.append(attempt)
    from .acceptance import _composition_signature, DEFAULT_COMPOSITION
    mixture_varies=any(_composition_signature(case[2]) != _composition_signature(DEFAULT_COMPOSITION)
                       for case in train_pairs+holdout_pairs)
    result={'model':'homogeneous adiabatic H2/O2 chemical kinetics; NO plasma',
            'cantera_version':getattr(ct,'__version__','unknown'),
            'source_mechanism':mechanism,'source_mechanism_sha256':source_digest,
            'reactor_t_end_s':float(t_end),
            'reactor_n_snapshots':int(n_snapshots),'initial_composition':'scenario_specific' if mixture_varies else DEFAULT_COMPOSITION,
            'ignition_definition':'peak_dTdt_if_threshold_resolved_with_20K_rise',
            'full_species':gas.n_species,
            'full_reactions':gas.n_reactions,'training_conditions':[list(x) for x in train_pairs],
            'holdout_conditions':[list(x) for x in holdout_pairs],
            'mechanism_trials':tests,'status':'NO REDUCTION',
            'ignition_training_cases_observed':n_igniting_train,
            'training_ignition_observations':[
                {'condition_T_K_P_atm':[float(T),float(P),composition],
                 'reference_ignition_observed':bool(train_refs[T,P,composition]['ignition_observed']),
                 'threshold_crossing_time_s':train_refs[T,P,composition]['threshold_crossing_time_s'],
                 'peak_dTdt_time_s':train_refs[T,P,composition]['peak_dTdt_time_s'],
                 'max_time_step_s':train_refs[T,P,composition]['max_time_step_s']}
                for T,P,composition in train_pairs],
            'ignition_required_for_acceptance':True,
            'candidate_selection':'max(net DRGEP, gross bidirectional-flux screen)',
            'net_species_importance':dict(zip(gas.species_names,net_importance.tolist())),
            'experimental_validation':False,
            'holdout_tolerance':holdout_tolerance,
            'training_tolerance':train_tolerance,
            'gross_flux_species_importance':dict(zip(gas.species_names,gross_importance_arr.tolist())),
            'warning':'This is conservative CHEMICAL DRGEP screening; NOT plasma P-DRGEP or experimental validation.'}
    if protocol_info is not None:
        result['evidence_schema']='v10-pinned'
        result['protocol_id']=protocol_info['id']
        result['protocol_sha256']=protocol_info['sha256']
    else:
        result['evidence_schema']='v10-unpinned-exploratory'
    if chosen is not None:
        _,skel,attempt,training_errors=chosen
        training_per_case=[
            {'condition_T_K_P_atm':[T,P,composition],**error}
            for (T,P,composition),error in zip(train_pairs,training_errors)]
        result['training_metrics']=training_per_case
        result['training_max_normalized_error']=max(_max_curve_error(e) for e in training_errors)
        skel.write_yaml(str(out/'candidate_skeletal.yaml'))
        # Re-parse the saved artifact with Cantera. A generated YAML is not
        # evidence until the solver itself can load it as a chemical phase.
        reread=ct.Solution(str(out/'candidate_skeletal.yaml'))
        if (reread.n_species!=skel.n_species or reread.n_reactions!=skel.n_reactions):
            raise RuntimeError('Round-trip skeletal YAML species/reaction count changed')
        errors=[]
        for T,P,composition in holdout_pairs:
            ref=_reactor_curve(ct,gas,T,P,t_end,n_snapshots,composition=composition)
            red=_reactor_curve(ct,skel,T,P,t_end,n_snapshots,composition=composition)
            errors.append({'condition_T_K_P_atm':[T,P,composition],**_curve_error(ref,red)})
        worst=max(_max_curve_error(case) for case in errors)
        n_igniting_holdout=sum(bool(x['reference_ignition_observed']) for x in errors)
        matched=all(x['ignition_classification_matches'] for x in errors)
        passed=worst<=holdout_tolerance and n_igniting_holdout>0 and matched
        result.update({'status':'HOLDOUT_PASS' if passed else 'FAILED_OR_INCONCLUSIVE_HOLDOUT',
                       'reduced_species':skel.n_species,'reduced_reactions':skel.n_reactions,
                       'species_reduction_fraction':1-skel.n_species/gas.n_species,
                       'candidate_threshold':attempt['threshold'],
                       'holdout_metrics':errors,'holdout_max_normalized_error':worst,
                       'ignition_holdout_cases_observed':n_igniting_holdout,
                       'holdout_ignition_classifications_match':matched,
                       'accepted_with_ignition_evidence':passed,
                       'artifact':'candidate_skeletal.yaml'})
        if passed:
            skel.write_yaml(str(out/'holdout_pass_skeletal.yaml'))
            import hashlib
            candidate_file=out/'holdout_pass_skeletal.yaml'
            result['artifact']='holdout_pass_skeletal.yaml'
            result['artifact_sha256']=hashlib.sha256(candidate_file.read_bytes()).hexdigest()
    _finalize_scientific_acceptance(result, out)
    (out/'skeletal_summary.json').write_text(json.dumps(result,indent=2))
    return result
