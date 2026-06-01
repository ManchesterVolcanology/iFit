"""Contains functions to characterise and generate an instrument line shape."""
import numpy as np
from scipy.special import gamma
from scipy.interpolate import griddata


# =============================================================================
# Super Gaussian
# =============================================================================

def super_gaussian(grid, w, k, a_w, a_k, shift=0, amp=1, offset=0):
    """Return a super-Gaussian line shape."""
    # Compute A
    A = k / (2 * w * gamma(1/k))

    # Form empty array
    ils = np.zeros(len(grid))

    # Iterate over x grid. If negative do one thing, if positive do the other
    ils = np.array([left_func(x, w, k, a_w, a_k) if x <= 0
                    else right_func(x, w, k, a_w, a_k)
                    for x in grid])

    # Shift the lineshape
    if shift != 0:
        mod_grid = grid + shift

        # Interpolate onto the measurement grid
        ils = griddata(mod_grid, ils, grid, method='cubic', fill_value=0.0)

    return ils * A * amp + offset


def left_func(x, w, k, a_w, a_k):
    """Left function for asymetric Gaussian."""
    return np.exp(-np.power(np.abs((x) / (w - a_w)), k - a_k))


def right_func(x, w, k, a_w, a_k):
    """Right function for asymetric Gaussian."""
    return np.exp(-np.power(np.abs((x) / (w + a_w)), k + a_k))


# =============================================================================
# make_ils
# =============================================================================

def make_ils(interval, FWEM, k=2, a_w=0, a_k=0):
    """Generate a synthetic instrument line shape.

    Generates a lineshape based on the super-Gaussian function:

    .                { exp(-| x / (w-a_w) | ^ (k-a_k)) for x <= 0
    G(x) = A(w, k) * {
    .                { exp(-| x / (w+a_w) | ^ (k+a_k)) for x > 0

    where A(w, k) = k / (2 * w * Gamma(1/k)).

    See Beirle et al (2017) for more details: doi:10.5194/amt-10-581-2017

    Parameters
    ----------
    interval : int
        The spacing of the wavelength grid on which the ILS is built
    FWEM : float
        The Full Width eth Maximum of the lineshape, defined as FWEM = 2*w
    k : float, optional
        Controls the shape of the lineshape (default = 2):
            - k < 2 -> sharp point and wide tails
            - k = 2 -> normal Gaussian
            - k > 2 -> flat top, approaches boxcar at k -> inf
    a_w and a_k : float, optional
        Controls the asymetry of the lineshape. Defaults are 0

    Returns
    -------
    ils : numpy array
        The calculated ILS function on a wavelength grid of the given spacing
        and 5 times the width of the supplied FWEM
    """
    # Create a 4 nm grid
    grid = np.arange(-2, 2, interval)

    # Calculate w as half of the FWEM
    w = 0.5 * FWEM

    # Make the line shape
    ils = super_gaussian(grid, w, k, a_w, a_k)

    # Normalise
    ils = np.divide(ils, sum(ils))

    return ils
