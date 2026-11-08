"""scipy.special.logsumexp for real NumPy arrays, computed step by step as scipy 1.16 computes it.

The DP sampler calls logsumexp several times per mutation and cluster on every iteration. scipy's array-API version
spends most of each call on dispatch and copies, which made it the largest cost of clustering. This follows the same
operations in the same order (scipy/special/_logsumexp.py: logsumexp and _logsumexp, b=None, return_sign=False), so
results are bitwise identical; test_exact_logsumexp.py checks this. scipy 1.15 differs from 1.16 only where the result
is not finite (-inf, +inf or NaN), and gives the same values there. Other scipy versions use other algorithms, and
there `logsumexp` below is scipy's own. Arguments other than `axis`, and arrays that are not native-order floats, go
to scipy.
"""
import numpy as np
import scipy
from scipy.special import logsumexp as _scipy_logsumexp

_SCIPY_MINOR = tuple(int(x) for x in scipy.__version__.split(".")[:2])
REPLICATES_SCIPY = _SCIPY_MINOR in ((1, 15), (1, 16))


def exact_logsumexp(a, axis=None, **kwargs):
    """scipy's result for a real, non-empty float array. Steps of scipy that cannot change a real result are left out:
    with b=None the sum s is non-negative or NaN and the tie count m at least 1, so scipy's sign tests never fire and
    its s == 0 branch equals s / m. Its fallback, log(sum(exp(a))), is computed only where the result is not finite."""
    a = np.asarray(a)
    if kwargs or a.dtype.kind != "f" or not a.dtype.isnative or a.size == 0:
        return _scipy_logsumexp(a, axis=axis, **kwargs)
    if a.ndim == 0:
        a = a.reshape(1)
    axis = tuple(range(a.ndim)) if axis is None else axis
    # scipy ignores overflow only in its fallback, so an overflowing a - a_max still warns as in scipy.
    with np.errstate(divide='ignore', invalid='ignore'):
        a_max = np.max(a, axis=axis, keepdims=True)
        i_max = a == a_max
        m = np.sum(i_max, axis=axis, keepdims=True, dtype=a.dtype)
        exp = np.exp(np.where(i_max, -np.inf, a) - a_max)
        s = np.sum(exp, axis=axis, keepdims=True, dtype=exp.dtype) / m
        out = np.log1p(s) + np.log(m) + a_max
        finite = np.isfinite(out)
    if not finite.all():
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            out = np.where(finite, out, np.log(np.sum(np.exp(a), axis=axis, keepdims=True)))
    out = np.squeeze(out, axis=axis)
    return out[()] if out.ndim == 0 else out


logsumexp = exact_logsumexp if REPLICATES_SCIPY else _scipy_logsumexp
