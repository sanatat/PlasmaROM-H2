"""Versioned, pinned evaluation plans for *numerical* Cantera reductions.

These repository-committed hashes are reproducibility guards, NOT cryptographic
proof of pre-registration with an independent third party. A different plan
requires a reviewed source change rather than a permissive CLI tolerance.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

PINNED_PROTOCOLS = {
    'v10_h2o2_default': '9aa8551b7d20308564723fce939a77f8d2d56efdfb5f11d78d4980be0d856981',
    'v10_h2o2_composition': '3af58dace503f30346471ae90489610afb49fa15757d9d99ebe27858571d6b69',
}


def load_pinned_protocol(path_or_id):
    if not isinstance(path_or_id,str):
        raise ValueError('Protocol must be a pinned id or repository protocol JSON path')
    pid=Path(path_or_id).stem
    if pid not in PINNED_PROTOCOLS:
        raise ValueError('Unknown / unlocked evaluation protocol')
    reference=Path(__file__).resolve().parent/'protocols'/f'{pid}.json'
    contents=reference.read_bytes()
    sha=hashlib.sha256(contents).hexdigest()
    if sha!=PINNED_PROTOCOLS[pid]:
        raise ValueError('Pinned evaluation protocol changed without review')
    if Path(path_or_id).is_file() and Path(path_or_id).resolve()!=reference.resolve():
        if hashlib.sha256(Path(path_or_id).read_bytes()).hexdigest()!=sha:
            raise ValueError('Supplied protocol differs from pinned repository plan')
    plan=json.loads(contents)
    if plan.get('id')!=pid:
        raise ValueError('Wrong protocol identifier')
    return {**plan,'sha256':sha}


def verify_protocol_parameters(plan,mechanism,training,holdout,t_end,snapshots,
                               thresholds,train_tolerance,holdout_tolerance):
    from .acceptance import _conditions,_overlap
    ref_train=_conditions(plan['training'])
    ref_holdout=_conditions(plan['holdout'])
    if ref_train is None or ref_holdout is None or _overlap(ref_train,ref_holdout):
        raise ValueError('Pinned protocol has invalid or leaked conditions')
    if (mechanism != plan['source_mechanism'] or
        _conditions([list(row) for row in training]) != ref_train or
        _conditions([list(row) for row in holdout]) != ref_holdout or
        not math.isclose(float(t_end),float(plan['reactor_t_end_s']),rel_tol=1e-12) or
        snapshots != plan['reactor_n_snapshots'] or
        list(thresholds)!=plan['thresholds'] or
        not math.isclose(float(train_tolerance),plan['train_tolerance'],rel_tol=1e-12) or
        not math.isclose(float(holdout_tolerance),plan['holdout_tolerance'],rel_tol=1e-12)):
        raise ValueError('Run parameters differ from locked scientific evaluation protocol')
    if plan.get('oh_mass_fraction_relative_floor')!=1e-12:
        raise ValueError('Unexpected OH evaluation floor')


def verify_summary_protocol(plan,summary):
    from .acceptance import _conditions
    if summary.get('protocol_sha256')!=plan['sha256']:
        raise ValueError('Invalid protocol checksum')
    if summary.get('protocol_id')!=plan['id']:
        raise ValueError('Protocol identifier changed')
    if summary.get('source_mechanism')!=plan['source_mechanism']:
        raise ValueError('Source mechanism changed')
    # Recorded thresholds can be present among the trials, but must be from
    # the fixed list; generated trials should exactly match plan length.
    if [t.get('threshold') for t in summary.get('mechanism_trials',[]) if isinstance(t,dict)]!=plan['thresholds']:
        raise ValueError('Threshold search differed from protocol')
    verify_protocol_parameters(
        plan, summary.get('source_mechanism'),
        _conditions(summary.get('training_conditions')) or [],
        _conditions(summary.get('holdout_conditions')) or [],
        summary.get('reactor_t_end_s',float('nan')),
        summary.get('reactor_n_snapshots'), plan['thresholds'],
        summary.get('training_tolerance',float('nan')),
        summary.get('holdout_tolerance',float('nan')),
    )
