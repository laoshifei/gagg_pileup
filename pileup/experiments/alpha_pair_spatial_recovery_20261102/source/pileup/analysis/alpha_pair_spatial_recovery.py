"""Spatially partitioned recovery with a separately learned illumination prior.

The input profile is fixed from development mothers. No evaluation event's
realized positions enter E/t prediction. Regionwise mean-field prompt XT omits
explicit inter-region neighbour transport; this remains a model approximation.
"""
from __future__ import annotations

import numpy as np

from .alpha_pair_causal_recovery import CausalRecovery, DifferentiableAge


def symmetric_profile(cell_probabilities, divisions):
    """Square D4 average on exact integer blocks of a 120x120-cell sensor."""
    p = np.asarray(cell_probabilities, float)
    if p.shape != (120,120) or 120 % divisions or divisions < 1:
        raise ValueError("Expected120x120 probabilities and an exact partition")
    if np.any(p < 0) or not np.isfinite(p).all() or not np.isclose(p.sum(),1.):
        raise ValueError("Normalized nonnegative cell probabilities required")
    width = 120 // divisions
    coarse = p.reshape(divisions,width,divisions,width).sum((1,3))
    images = [np.rot90(coarse,k) for k in range(4)]
    images += [np.fliplr(im) for im in images]
    return np.mean(images,axis=0)


def symmetry_groups(profile):
    """Merge equal-illumination D4 tiles without changing independent states."""
    p = np.asarray(profile,float)
    n = len(p)
    if p.shape != (n,n) or 120 % n or not np.isclose(p.sum(),1.) or np.any(p<0):
        raise ValueError("Invalid square spatial profile")
    groups = {}
    for y in range(n):
        for x in range(n):
            key = tuple(sorted((min(x,n-1-x),min(y,n-1-y))))
            groups.setdefault(key,[]).append((y,x))
    cells, probabilities = [], []
    for indices in groups.values():
        values = np.array([p[y,x] for y,x in indices])
        if np.ptp(values)>1e-14:
            raise ValueError("Profile is not D4 symmetric")
        cells.append(len(indices)*(120//n)**2)
        probabilities.append(float(values.sum()))
    return np.asarray(cells,float),np.asarray(probabilities,float)


class SpatialAge:
    """Independent recovering-age distributions for frozen illumination groups."""
    def __init__(self,sensor,library,profile,dt_ns=1.):
        if (sensor['cells_x'],sensor['cells_y']) != (120,120):
            raise ValueError("This spatial partition is defined for120x120cells")
        self.profile=np.asarray(profile,float).copy()
        self.cell_counts,self.probabilities=symmetry_groups(self.profile)
        self.cells=float(self.cell_counts.sum())
        self.dt=float(dt_ns)
        self.engines=[]
        for cells in self.cell_counts:
            engine=DifferentiableAge(sensor,library,self.dt)
            engine.cells=float(cells)
            self.engines.append(engine)

    def forward(self,light):
        x=np.asarray(light,float)
        if x.ndim!=1 or np.any(x<0) or not np.isfinite(x).all():
            raise ValueError("Expected nonnegative finite light counts")
        output=np.zeros_like(x)
        for engine,probability in zip(self.engines,self.probabilities):
            output+=engine.forward(x*probability)
        return output

    def with_gradient(self,light,derivatives):
        x=np.asarray(light,float)
        d=np.asarray(derivatives,float)
        output=np.zeros_like(x)
        gradient=np.zeros_like(d)
        for engine,probability in zip(self.engines,self.probabilities):
            y,dy=engine.with_gradient(x*probability,d*probability)
            output+=y
            gradient+=dy
        return output,gradient


class SpatialRecovery(CausalRecovery):
    """Calibrated isolated anchor plus a spatial shared-minus-isolated correction.

Charge-to-light inversion uses the SAME spatial operator and frozen profile as
the pair response. Single-event voltage is inherited exactly by construction;
an additional truth-time raw response check is required to test spatial physics.
"""
    def __init__(self,anchor,configuration,library,charge_map,profile):
        super().__init__(anchor,configuration,library,charge_map,strength=1.)
        self.age=SpatialAge(configuration['sensor'],library,profile,self.dt)
