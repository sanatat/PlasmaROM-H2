"""2D periodic incompressible Navier-Stokes and advected surrogate species.

Fourier pseudo-spectral RK2 projection for velocity, first-order monotone
upwind + centered diffusion + local reaction for scalar chemical equivalents.
Optional Smagorinsky-*style* strain stress divergence; not validated DNS/LES.
"""
from dataclasses import dataclass
import numpy as np
from .kinetics import rhs_full, rhs_reduced, Parameters

@dataclass
class FlowConfig:
    n: int = 40
    dt: float = 0.02
    t_end: float = 4.2
    nu: float = 0.01
    diffusion: float = 0.006
    cs: float = 0.0
    reduced: bool = False

class SpectralFlow:
    def __init__(self, n, nu=0.01, cs=0.0):
        self.n = n
        self.nu = nu
        self.cs = cs
        self.dx = 2*np.pi/n
        k = np.fft.fftfreq(n, d=1/n)
        self.kx = k[:,None]
        self.ky = k[None,:]
        self.k2 = self.kx**2 + self.ky**2
        self.k2safe = np.where(self.k2 == 0, 1, self.k2)
        self.filter = ((abs(self.kx) < n/3) & (abs(self.ky) < n/3)).astype(float)
        x = 2*np.pi*np.arange(n)/n
        self.x, self.y = np.meshgrid(x,x,indexing='ij')
        self.u = np.sin(self.x)*np.cos(self.y)
        self.v = -np.cos(self.x)*np.sin(self.y)

    def deriv(self, a, axis):
        ik = 1j*(self.kx if axis == 0 else self.ky)
        return np.fft.ifft2(ik*np.fft.fft2(a)).real

    def lap(self, a):
        return np.fft.ifft2(-self.k2*np.fft.fft2(a)).real

    def project(self,u,v):
        uh,vh = np.fft.fft2(u),np.fft.fft2(v)
        dp = (self.kx*uh + self.ky*vh)/self.k2safe
        uh = (uh-self.kx*dp)*self.filter
        vh = (vh-self.ky*dp)*self.filter
        return np.fft.ifft2(uh).real,np.fft.ifft2(vh).real

    def divergence(self):
        return self.deriv(self.u,0)+self.deriv(self.v,1)

    def tendencies(self,u,v):
        ux,uy = self.deriv(u,0),self.deriv(u,1)
        vx,vy = self.deriv(v,0),self.deriv(v,1)
        tu = -(u*ux+v*uy)
        tv = -(u*vx+v*vy)
        if self.cs:
            strain = np.sqrt(2*(ux**2+vy**2+2*(0.5*(uy+vx))**2))
            nut = (self.cs*self.dx)**2 * strain
            viscosity = self.nu+nut
            # Full symmetric strain tensor tau_ij = 2(nu+nu_t) S_ij.
            # This implements a variable-eddy-viscosity demonstration,
            # NOT a verified LES closure for turbulent reacting flows.
            sxy = 0.5*(uy+vx)
            tu += self.deriv(2*viscosity*ux,0)+self.deriv(2*viscosity*sxy,1)
            tv += self.deriv(2*viscosity*sxy,0)+self.deriv(2*viscosity*vy,1)
        else:
            tu += self.nu*self.lap(u)
            tv += self.nu*self.lap(v)
        return self.project(tu,tv)

    def step(self,dt):
        u0,v0 = self.u,self.v
        du0,dv0 = self.tendencies(u0,v0)
        ue,ve = self.project(u0+dt*du0,v0+dt*dv0)
        du1,dv1 = self.tendencies(ue,ve)
        self.u,self.v = self.project(u0+dt*(du0+du1)/2,v0+dt*(dv0+dv1)/2)

    def taylor_green_error(self,t):
        amp = np.exp(-2*self.nu*t)
        utrue = amp*np.sin(self.x)*np.cos(self.y)
        vtrue = -amp*np.cos(self.x)*np.sin(self.y)
        return np.sqrt(np.mean((self.u-utrue)**2+(self.v-vtrue)**2))


def divergence_free_face_velocities(u, v, dx):
    """Convert periodic Fourier-projected cell-center velocities to face fluxes.

    Reconstruct a streamfunction at cell corners; its discrete curl gives
    face-normal velocities with zero *finite-volume* divergence (up to
    floating-point rounding). Constant mean flow is retained separately.
    uface[i,j] is x-flux through face (i+1/2,j); vface[i,j] is y-flux
    through (i,j+1/2). No density variation or compressibility model.
    """
    u=np.asarray(u,dtype=float);v=np.asarray(v,dtype=float)
    if u.shape!=v.shape or u.ndim!=2 or u.shape[0]!=u.shape[1]:
        raise ValueError('Expected square collocated 2D velocity arrays')
    n=u.shape[0]
    k=np.fft.fftfreq(n,d=1/n)
    kx,ky=k[:,None],k[None,:]
    k2=kx**2+ky**2
    uh=np.fft.fft2(u-u.mean());vh=np.fft.fft2(v-v.mean())
    psi_hat=np.divide(1j*(kx*vh-ky*uh),k2,
                      out=np.zeros_like(uh),where=k2>0)
    # Fourier evaluate streamfunction at (x+dx/2,y+dx/2) corners.
    psi_corner=np.fft.ifft2(psi_hat*np.exp(0.5j*dx*(kx+ky))).real
    uface=u.mean()+(psi_corner-np.roll(psi_corner,1,axis=1))/dx
    vface=v.mean()-(psi_corner-np.roll(psi_corner,1,axis=0))/dx
    return uface,vface


def scalar_transport(y,u,v,dx,diffusivity,faces=None):
    """Conservative first-order upwind finite-volume flux divergence.

    Face-centered flow is divergence-free under the SAME discrete operator.
    The flux divergence telescopes on a periodic domain, preserving the
    global amount of every individual scalar to floating-point accuracy.
    Upwind + centered diffusion preserves positivity under the CFL check.
    """
    if faces is None:
        faces=divergence_free_face_velocities(u,v,dx)
    uf,vf=faces
    xnext=np.roll(y,-1,axis=-2)
    ynext=np.roll(y,-1,axis=-1)
    flux_x=np.maximum(uf,0)*y+np.minimum(uf,0)*xnext - diffusivity*(xnext-y)/dx
    flux_y=np.maximum(vf,0)*y+np.minimum(vf,0)*ynext - diffusivity*(ynext-y)/dx
    return -((flux_x-np.roll(flux_x,1,axis=-2))+
             (flux_y-np.roll(flux_y,1,axis=-1)))/dx


def finite_volume_cfl(faces,dx,diffusivity):
    """Maximum positive outgoing finite-volume rate, including diffusion."""
    uf,vf=faces
    outgoing=(np.maximum(uf,0)+np.maximum(-np.roll(uf,1,axis=0),0)+
              np.maximum(vf,0)+np.maximum(-np.roll(vf,1,axis=1),0))/dx
    return float(np.max(outgoing)+4*diffusivity/dx**2)


def reaction_outflow_bound(p=Parameters(), reduced=False):
    """Guaranteed bound for the prescribed nonnegative Gaussian drive (mask <= 1)."""
    if not np.isfinite(p.pulse_strength) or not np.isfinite(p.pulse_width) or p.pulse_strength < 0 or p.pulse_width <= 0:
        raise ValueError('Pulse strength must be >= 0 and width > 0')
    raw_rates=(p.excitation,p.second_excitation,p.relaxation10,
               p.relaxation21,p.activated_oxidation,p.background_oxidation,
               p.radical_oxidation,p.radical_quench)
    if any((not np.isfinite(v)) or v<0 for v in raw_rates):
        raise ValueError('All chemical transfer coefficients must be finite and >= 0')
    drive_max=1.75*p.pulse_strength  # sum of two Gaussians <= 1.75
    # Linear transition network with positive, conservative generator.
    outgoing = [p.excitation*drive_max+p.background_oxidation,
                p.radical_oxidation+p.radical_quench]
    if reduced:
        # B->F0 and B->R are weighted by fractions (1-f2) and f2.
        # Sum is <= max(relaxation10, activated_oxidation).
        if p.relaxation21+p.activated_oxidation<=0:
            raise ValueError('QSSA fast-state denominator must be positive')
        outgoing.append(max(p.relaxation10,p.activated_oxidation))
    else:
        outgoing.extend((p.relaxation10+p.second_excitation*drive_max,
                         p.relaxation21+p.activated_oxidation))
    if any(v<0 for v in outgoing):
        raise ValueError('All transition rates must be nonnegative')
    return max(outgoing)


def reaction_ssprk2(y,t,dt,p,mask,reduced=False):
    """Reaction operator: SSPRK2 with conservative positive subcycling.

    Each forward Euler stage is nonnegative for a linear transfer network
    if h*outflow_max <= 1; SSPRK2 is a convex combination of Euler stages.
    """
    if dt<=0:
        return y,0
    bound=reaction_outflow_bound(p,reduced)
    nsteps=max(1,int(np.ceil(dt*bound/0.7)))
    h=dt/nsteps
    f=rhs_reduced if reduced else rhs_full
    out=y.copy()
    for j in range(nsteps):
        tj=t+j*h
        stage=out+h*f(tj,out,p,mask)
        out=0.5*out+0.5*(stage+h*f(tj+h,stage,p,mask))
    return out,nsteps


def run_flow(cfg=FlowConfig(), p=Parameters()):
    """2D idealized advection-diffusion-reaction demonstration.

    Strang split: local positive conservative chemical SSPRK2 half steps,
    upwind-forward-Euler transport, chemical half step. Every species uses
    the same transport, so sum-of-pools remains approximately one locally.
    Reaction auto-subcycles; spatial advection/diffusion CFL is checked.
    Thermal energy, density change, real plasma kinetics and flames absent.
    """
    if cfg.n < 12 or cfg.n % 2 != 0 or cfg.dt<=0 or cfg.t_end<0:
        raise ValueError('n must be even >=12, dt >0, t_end >=0')
    if cfg.nu<0 or cfg.diffusion<0 or cfg.cs<0:
        raise ValueError('Viscosities and diffusivity must be nonnegative')
    flow = SpectralFlow(cfg.n,cfg.nu,cfg.cs)
    n = cfg.n
    coords = flow.x,flow.y
    mask = np.exp(-((flow.x-np.pi)**2+(flow.y-np.pi)**2)/(2*0.9**2))
    m = 4 if cfg.reduced else 5
    y = np.zeros((m,n,n)); y[0] = 1.
    t = 0.
    count = int(np.ceil(cfg.t_end/cfg.dt))
    peak_div=0.
    min_pop=1.
    max_substeps=0
    for i in range(count):
        dt=min(cfg.dt,cfg.t_end-t)
        if dt <= 0:break
        # SSP-compatible positivity bound for first-order upwind transport
        # plus explicit 5-point diffusion; unlike two independent CFL checks
        # this enforces a combined stability condition.
        faces=divergence_free_face_velocities(flow.u,flow.v,flow.dx)
        cfl=dt*finite_volume_cfl(faces,flow.dx,cfg.diffusion)
        if cfl>0.9:
            raise ValueError('Transport CFL exceeded: decrease dt (combined advection/diffusion)')
        # Conservative bound on the largest retained Fourier diffusion mode.
        # The eddy viscosity is an instantaneous estimate, not predictive.
        if cfg.cs:
            ux=flow.deriv(flow.u,0);uy=flow.deriv(flow.u,1)
            vx=flow.deriv(flow.v,0);vy=flow.deriv(flow.v,1)
            max_strain=float(np.max(np.sqrt(2*(ux**2+vy**2+2*(0.5*(uy+vx))**2))))
            nu_max=cfg.nu+(cfg.cs*flow.dx)**2*max_strain
        else:
            nu_max=cfg.nu
        if dt*nu_max*2*(cfg.n/2)**2>1.5:
            raise ValueError('Spectral velocity viscous limit exceeded: reduce dt')
        y,ns1=reaction_ssprk2(y,t,dt/2,p,mask,cfg.reduced)
        y=y+dt*scalar_transport(y,flow.u,flow.v,flow.dx,cfg.diffusion,faces=faces)
        y,ns2=reaction_ssprk2(y,t+dt/2,dt/2,p,mask,cfg.reduced)
        max_substeps=max(max_substeps,ns1,ns2)
        flow.step(dt)
        t += dt
        if not np.all(np.isfinite(y)) or not np.all(np.isfinite(flow.u)) or not np.all(np.isfinite(flow.v)):
            raise FloatingPointError('Nonfinite field; decrease dt or inspect parameters')
        if y.min() < -1e-10:
            raise FloatingPointError('Negative chemical pool; integration is not positivity preserving')
        peak_div=max(peak_div,float(np.max(np.abs(flow.divergence()))))
        min_pop=min(min_pop,float(y.min()))
    return {'x':coords[0],'y':coords[1],'u':flow.u,'v':flow.v,
            'species':y,'t':t,'peak_divergence':peak_div,'min_pool':min_pop,
            'sum_deviation':float(np.max(np.abs(y.sum(axis=0)-1))),
            'global_species_integrals':(y.sum(axis=(1,2))*flow.dx**2).tolist(),
            'taylor_green_l2':flow.taylor_green_error(t), 'cs':cfg.cs,
            'n':cfg.n,'reduced':cfg.reduced,'max_chemical_substeps_per_half_step':max_substeps}
