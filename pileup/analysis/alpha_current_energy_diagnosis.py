"""Truth-assisted residual bookkeeping and local log-energy diagnostics.

These calculations diagnose saved fits. They do not reconstruct an event or
estimate calibrated uncertainties, and must never be used as blind inputs.
"""
import numpy as np

TERM_NAMES = ('strong_isolated', 'weak_isolated', 'shared_difference', 'background_readout')


def residual_terms(voltage, shared, isolated, predicted, predicted_isolated, background, weak):
    y, s, i = map(lambda a: np.asarray(a, float), (voltage, shared, isolated))
    p, q = np.asarray(predicted, float), np.asarray(predicted_isolated, float)
    if s.shape != i.shape or s.shape != (len(y), 2, y.shape[-1]) or p.shape != q.shape or p.shape != s.shape[1:]:
        raise ValueError('Expected repeat/component/sample arrays and fixed predictions')
    errors = i - q
    terms = np.stack((errors[:, 1-weak], errors[:, weak],
                      (s-i).sum(1) - (p-q).sum(0),
                      y-s.sum(1)-background), axis=1)
    total = y-p.sum(0)-background
    if not np.allclose(terms.sum(1), total, rtol=1e-12, atol=1e-10):
        raise ValueError('Residual decomposition failed to close')
    return terms, total


def local_operator(jacobian):
    """Return a two-column SVD inverse; columns are strong/weak dV/dlnE."""
    d = np.asarray(jacobian, float)
    if d.ndim != 2 or d.shape[1] != 2 or not np.isfinite(d).all():
        raise ValueError('Expected finite sample-by-two derivatives')
    u, singular, vt = np.linalg.svd(d, full_matrices=False)
    if singular[-1] <= singular[0]*1e-12:
        raise ValueError('Local energy directions are numerically rank deficient')
    inverse = (vt.T/singular)@u.T
    norms = np.linalg.norm(d, axis=0)
    rho = float(d[:, 0]@d[:, 1]/np.prod(norms))
    perpendicular = d[:, 1]-d[:, 0]*(d[:, 0]@d[:, 1]/(norms[0]**2))
    weak_operator = perpendicular/(perpendicular@perpendicular)
    if not np.allclose(inverse[1], weak_operator, rtol=1e-8, atol=1e-10):
        raise ValueError('SVD and strong-profiled weak projection disagree')
    return inverse, dict(cosine=rho, normalized_condition=float(np.linalg.cond(d/norms)),
        condition=float(singular[0]/singular[-1]), strong_norm=float(norms[0]), weak_norm=float(norms[1]),
        weak_retained_fraction=float(np.linalg.norm(perpendicular)/norms[1]),
        weak_projection_norm=float(np.linalg.norm(inverse[1])))


def repeat_covariance(projected_terms):
    """Keep all cross terms: isolated and shared differences are correlated."""
    x = np.asarray(projected_terms, float)
    if x.ndim != 2 or x.shape[1] != 4 or len(x) < 2:
        raise ValueError('Expected repeats by four weak projections')
    covariance = np.cov(x, rowvar=False, ddof=1)
    total_variance = float(np.var(x.sum(1), ddof=1))
    if not np.isclose(covariance.sum(), total_variance, rtol=1e-10, atol=1e-10):
        raise ValueError('Full covariance does not reproduce total variance')
    return dict(covariance=covariance.tolist(), total_variance=total_variance,
        diagonal_variance=float(np.trace(covariance)), cross_variance=float(covariance.sum()-np.trace(covariance)),
        mean_mc_standard_error=float(np.sqrt(max(0., total_variance)/len(x))))
