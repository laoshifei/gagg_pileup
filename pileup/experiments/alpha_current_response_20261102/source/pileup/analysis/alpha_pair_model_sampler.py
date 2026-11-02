"""Exact Gaussian draws for the frozen joint likelihood, not native photons.

The model specifies the marginal distribution of FIR-transformed samples.
The first FIR-order voltage samples must therefore be specified separately.
No parameter clipping, ADC saturation, or amplitude rescaling is applied.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from scipy.signal import lfilter


def lower_band_product(band, values):
    """Multiply a lower-triangular band matrix along the final array axis."""
    x = np.asarray(values, float)
    n = band.shape[1]
    if x.shape[-1] != n:
        raise ValueError('Incompatible band and sample dimensions')
    out = np.zeros_like(x)
    for k in range(len(band)):
        out[..., k:] += band[k, :n-k] * x[..., :n-k]
    return out


@dataclass
class GaussianFactor:
    """A square root and inverse for C = B - V V.T in O(n b + n k)."""
    band_cholesky: np.ndarray
    U: np.ndarray
    subtract: np.ndarray
    add: np.ndarray
    eigenvalues: np.ndarray

    @classmethod
    def from_covariance(cls, band, V):
        factor = cholesky_banded(band, lower=True, check_finite=True)
        U = solve_banded((len(factor)-1, 0), factor, V, check_finite=True) if V.shape[1] else V.copy()
        if U.shape[1]:
            eigenvalues, directions = np.linalg.eigh(U.T @ U)
            if eigenvalues.min() < -1e-10 or eigenvalues.max() >= 1.:
                raise np.linalg.LinAlgError('Covariance downdate is not positive definite')
            # The forms below remain finite when an eigenvalue tends to zero.
            root = np.sqrt(1. - np.maximum(eigenvalues, 0.))
            subtract = (directions / (1. + root)) @ directions.T
            add = (directions / (root * (1. + root))) @ directions.T
        else:
            eigenvalues = np.empty(0)
            subtract = add = np.empty((0, 0))
        return cls(factor, U, subtract, add, eigenvalues)

    @classmethod
    def from_noise(cls, noise, curves, charges):
        return cls.from_covariance(*noise.covariance(curves, charges))

    @property
    def size(self):
        return self.band_cholesky.shape[1]

    @property
    def logdet(self):
        return float(2. * np.log(self.band_cholesky[0]).sum() + np.log1p(-self.eigenvalues).sum())

    def color(self, standard):
        z = np.asarray(standard, float)
        w = z - ((z @ self.U) @ self.subtract) @ self.U.T
        return lower_band_product(self.band_cholesky, w)

    def draw(self, rng, count=None):
        shape = (self.size,) if count is None else (count, self.size)
        return self.color(rng.standard_normal(shape))

    def whiten(self, residual):
        r = np.asarray(residual, float)
        shape = r.shape
        x = solve_banded((len(self.band_cholesky)-1, 0), self.band_cholesky,
                         r.reshape(-1, self.size).T, check_finite=True).T.reshape(shape)
        return x + ((x @ self.U) @ self.add) @ self.U.T


def lift_innovations(innovations, coefficients, prefix=None):
    """Invert the fixed FIR conditional on explicitly supplied first samples.

    This supplies a right inverse of the rectangular likelihood transform;
    it does not invent an unmodelled stationary prior for the first samples.
    The experiment uses zero residual prefix and checks its transient decay.
    """
    z = np.asarray(innovations, float)
    a = np.asarray(coefficients, float)
    order = len(a)-1
    if z.ndim != 1 or a.ndim != 1 or not len(a) or a[0] == 0.:
        raise ValueError('Expected a causal nonzero-leading FIR and a vector')
    first = np.zeros(order) if prefix is None else np.asarray(prefix, float)
    if first.shape != (order,):
        raise ValueError('The prefix must contain exactly FIR-order samples')
    forcing = np.r_[lfilter(a, [1.], first) if order else [], z]
    return lfilter([1.], a, forcing)


def conditional_mean(response, noise, times, charges):
    return response.components(times, charges).sum(0) + response.baseline + noise.background_mV
