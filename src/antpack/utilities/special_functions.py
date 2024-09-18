"""Contains special functions implemented here to avoid having to
introduce additional dependencies (logsumexp from scipy)."""
import numpy as np


def logsumexp(a, axis):
    """Compute the log of the sum of exponentials of input elements
    in a numerically stable way. This is a simplified version of
    scipy's logsumexp function.

    Args:
        a (np.ndarray): The input array.
        axis (int): The axis over which to apply this operation.

    Returns:
        out (np.ndarray): The result of the logsumexp operation.
    """
    initial_value = -np.inf if np.size(a) == 0 else None
    a_max = np.amax(a, axis=axis, keepdims=True, initial=initial_value)

    if a_max.ndim > 0:
        a_max[~np.isfinite(a_max)] = 0
    elif not np.isfinite(a_max):
        a_max = 0

    tmp = np.exp(a - a_max)

    # suppress warnings about log of zero
    with np.errstate(divide='ignore'):
        s = np.sum(tmp, axis=axis, keepdims=False)
        out = np.log(s)

    a_max = np.squeeze(a_max, axis=axis)
    out += a_max

    return out
