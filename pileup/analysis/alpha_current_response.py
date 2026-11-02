"""Statistics for paired transport controls; repeats never count as mothers."""
import numpy as np


def paired_mother_interval(old, current, *, seed, draws=4000):
    """Equal mother weighting, paired cluster bootstrap of the mean change.

    Arrays have (mother, repeat) shape. The interval is descriptive for the
    selected engineering controls; resampling does not make them independent
    of earlier model development or establish population equivalence.
    """
    old, current = np.asarray(old, float), np.asarray(current, float)
    if old.shape != current.shape or old.ndim != 2 or old.shape[0] < 2:
        raise ValueError('Matched mother-by-repeat arrays with at least two mothers required')
    if not np.isfinite(old).all() or not np.isfinite(current).all():
        raise ValueError('Finite observations required')
    difference = (current - old).mean(axis=1)
    index = np.random.default_rng(seed).integers(len(difference), size=(draws, len(difference)))
    return dict(mothers=len(difference), repeats=old.shape[1], mean_change=float(difference.mean()),
                ci95=np.percentile(difference[index].mean(axis=1), [2.5, 97.5]).tolist())


def energy_metrics(errors, boundaries, converged):
    values = np.asarray(errors, float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('At least two finite relative energy errors required')
    if len(boundaries) != len(values) or len(converged) != len(values):
        raise ValueError('Fit status and energy errors must align')
    return dict(rows=len(values), bias_percent=float(values.mean()),
                rmse_percent=float(np.sqrt(np.mean(values**2))),
                error_sd_percent=float(values.std(ddof=1)),
                strict10_count=int((abs(values) <= 10).sum()),
                boundary=int(np.sum(boundaries)), nonconverged=int(np.sum(~np.asarray(converged, bool))))
