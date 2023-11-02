"""Synthetic data generation functions.

Written by Alexandre Fournier, modified by Leonard Seydoux.
"""

import numpy as np


def time_range(
    time_start: float, time_end: float, time_step: float
) -> np.ndarray:
    """Generate a time range.

    Parameters:
    -----------
    time_start : float
        The start time of the time window.
    time_end : float
        The end time of the time window.
    time_step : float
        The time interval between two consecutive points in the time window.

    Returns:
    --------
    time_window : numpy.ndarray
        A 1D numpy array containing the time window.
    """
    num_points = int((time_end - time_start) / time_step)
    time_window = np.linspace(
        time_start, time_end, num=num_points, endpoint=True, dtype=float
    )
    return time_window


def sine(
    times: np.ndarray, period: float, amplitude: float, phase: float
) -> np.ndarray:
    """Sine function.

    Generate a sine function with the given parameters. The formula is:

    .. math::

            y(t) = A \sin(2 \pi t / T + \phi)

    where :math:`A` is the amplitude, :math:`T` is the period and :math:`\phi`
    is the phase.

    Parameters:
    -----------
    times : numpy.ndarray
        A 1D numpy array containing the time values.
    period : float
        The period of the cosine function.
    amplitude : float
        The amplitude of the cosine function.
    phase : float
        The phase of the cosine function.

    Returns:
    --------
    numpy.ndarray
        A 1D numpy array containing the sine time series.
    """
    return amplitude * np.sin(2 * np.pi * times / period + phase)
