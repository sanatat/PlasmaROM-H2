"""Simple integral electron-energy channel screening; NOT plasma P-DRGEP.

Given measured/derived rate history, assign gross source power to each
collisional channel, rank and select until a fraction of gross power is
retained. No inter-channel kinetic importance or post-reduction guarantee.
"""
import numpy as np
from scipy.integrate import trapezoid


def screen_power_channels(times, rates, gaps_eV, names, fraction=0.95):
    times=np.asarray(times,float)
    rates=np.asarray(rates,float)
    gaps=np.asarray(gaps_eV,float)
    if (times.ndim!=1 or times.size<2 or rates.shape!=(times.size,len(names))
        or gaps.shape!=(len(names),) or not np.all(np.isfinite(times))
        or not np.all(np.diff(times)>0) or not np.all(np.isfinite(rates))
        or np.any(rates<0) or not np.all(np.isfinite(gaps)) or np.any(gaps<0)
        or not 0<fraction<=1):
        raise ValueError('Invalid power-screening data')
    power=np.array([trapezoid(rates[:,i]*gaps[i],times) for i in range(len(names))])
    total=float(power.sum())
    if total==0:
        return {'selected':[],'retained_fraction':0.,'integrated_power_eV_flux':dict(zip(names,power.tolist()))}
    rank=np.argsort(-power,kind='stable')
    picked=[];accum=0.
    for idx in rank:
        if power[idx]==0:continue
        picked.append(names[idx]);accum+=power[idx]
        if accum/total>=fraction-1e-12:break
    return {'selected':picked,'retained_fraction':float(accum/total),
            'integrated_power_eV_flux':dict(zip(names,power.tolist()))}
