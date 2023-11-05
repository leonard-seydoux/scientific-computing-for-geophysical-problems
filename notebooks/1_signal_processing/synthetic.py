"""Synthetic data generation functions.

Written by Alexandre Fournier, modified by Leonard Seydoux.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.signal.windows import tukey


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


def make_illustration():
    # Define parameters for the chirp signal
    f0 = 5  # starting frequency (Hz)
    f1 = 40  # ending frequency (Hz)
    T = 3  # duration (seconds)
    Fs = 200  # Sampling frequency (samples/second)
    t = np.arange(0, T, 1 / Fs)  # time array

    # Calculate rate of change of frequency k
    k = (f1 - f0) / T

    # Generate the chirp signal
    x = np.sin(2 * np.pi * (f0 * t + 0.5 * k * t**2))

    # Design a bandpass filter for the 20-30 Hz band
    b, a = butter(2, [12, 20], btype="band", fs=Fs)

    # Apply the bandpass filter to the chirp signal
    y_bandpassed = filtfilt(b, a, x, padlen=100, padtype="odd")

    # Compute the Fourier Transform of the bandpassed signal
    window = tukey(len(x), alpha=0.4)
    X = np.fft.rfft(x * window)
    Y_bandpassed = np.fft.rfft(y_bandpassed * window)
    frequencies = np.fft.rfftfreq(len(x), 1 / Fs)

    # Adjusting the color to 'C3' and regenerating the plots
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))

    # Making figure background transparent
    fig.patch.set_alpha(0.0)

    # Chirp signal and bandpassed signal
    ax[0].plot(t, x, color="k", label="Original", lw=0.5)
    ax[0].plot(t, y_bandpassed, color="C3", label="Bandpassed")
    ax[0].set_title("Chirp Signal")
    ax[0].set_xlabel("Time $t$ in seconds")
    ax[0].set_ylabel(r"Amplitude $x(t) = A_0 \sin t\omega(t)$")
    ax[0].legend(loc="upper right")
    ax[0].set_xticks(np.linspace(0, T, T + 1))
    ax[0].set_yticks(np.linspace(-1, 1, 3))
    ax[0].set_yticklabels(["$-A_0$", "0", "$A_0$"])
    ax[0].grid(True)

    # Fourier Transforms
    ax[1].plot(frequencies, np.abs(X), color="k", lw=0.5)
    ax[1].plot(frequencies, np.abs(Y_bandpassed), color="C3")
    ax[1].set_title("Spectrum")
    ax[1].set_xlabel(r"Angular Frequency $\omega$ in radians/second")
    ax[1].set_ylabel(r"Spectral amplitude $|\hat x(\omega)|$")
    ax[1].set_xticks(np.linspace(0, Fs / 2, 3))
    ax[1].set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])
    ax[1].set_yticks(np.linspace(0, np.abs(X).max(), 2))
    ax[1].set_yticklabels(["0", r"$\hat A_0$"])
    ax[1].grid(True)

    # Making axes background semi-transparent
    for axis in ax:
        axis.set_facecolor((1.0, 1.0, 1.0, 0.5))

    # Adjust layout for better spacing and save to SVG with transparent background
    svg_c3_path = "illustration_signal_processing.svg"
    fig.tight_layout()
    fig.savefig(svg_c3_path, facecolor="none")


if __name__ == "__main__":
    make_illustration()
