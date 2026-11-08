"""Check that exact_logsumexp matches scipy.special.logsumexp bit for bit on the inputs the DP sampler produces.

Run: python Cluster/test_exact_logsumexp.py
"""
import numpy as np
from scipy.special import logsumexp as scipy_logsumexp

from exact_logsumexp import REPLICATES_SCIPY, exact_logsumexp


def same(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return x.dtype == y.dtype and x.shape == y.shape and x.tobytes() == y.tobytes()


def cases(rng):
    for dtype in (np.float32, np.float64):
        for _ in range(4000):
            row = rng.normal(-50, 40, size=101).astype(dtype)
            kind = rng.integers(6)
            if kind == 1:
                row[rng.random(101) < 0.3] = -np.inf
            elif kind == 2:
                row[:] = -np.inf
            elif kind == 3:
                row[rng.integers(101, size=3)] = row.max()
            elif kind == 4:
                row = np.log(np.clip(rng.dirichlet(np.ones(101) * 0.05), 1e-40, 1)).astype(dtype)
            elif kind == 5:
                row = np.full(101, -1e4, dtype=dtype)
            yield row
        for _ in range(4000):
            # loglik vectors over clusters: a few entries, some -inf for skipped clusters
            vec = rng.normal(-2000, 800, size=rng.integers(1, 12)).astype(dtype)
            vec[rng.random(vec.size) < 0.2] = -np.inf
            yield vec


def main():
    if not REPLICATES_SCIPY:
        print("scipy version not replicated: exact_logsumexp is not used, nothing to test")
        return
    rng = np.random.default_rng(0)
    n = 0
    for row in cases(rng):
        assert same(exact_logsumexp(row), scipy_logsumexp(row)), row
        two_d = row.reshape(1, -1)
        assert same(exact_logsumexp(two_d[0]), scipy_logsumexp(two_d[0]))
        n += 1
    print(f"{n} inputs: exact_logsumexp matches scipy bit for bit")


if __name__ == "__main__":
    main()
