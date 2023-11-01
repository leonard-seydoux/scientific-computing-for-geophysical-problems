"""Ambient noise processing functions.

This module contains functions for processing ambient noise data ang generating
noise cross-correlations, and for generating synthetic noise.

Made in 2023 by Léonard Seydoux (seydoux@ipgp.fr) for the class "Scientific
Computing for Geophysical Problems" at the Institut de Physique du Globe de
Paris. Some functions were written by Celine Hadziioannou and Ashim Rijal.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt, correlate


def source_time_function(n_sources, n_samples, amplitudes=None, kind="noise"):
    """Generate a source time function.

    Creates a source time function for each source. By default, the source time
    function is a random noise. The amplitudes of the sources can be specified
    by the user.

    Parameters
    ----------
    n_sources : int
        Number of sources.
    n_samples : int
        Number of samples.
    kind : str, optional
        Kind of source time function. Default is "noise".

    Returns
    -------
    np.ndarray
        Source time function for each source.
    """
    # Set amplitudes to 1 if not specified
    amplitudes = np.ones(n_sources) if amplitudes is None else amplitudes

    # Generate source time function
    if kind == "noise":
        records = np.random.randn(n_sources, n_samples)

    # Scale source time function
    records *= amplitudes[:, np.newaxis]

    return records


def bandpass(records, freqmin, freqmax, sampling_rate=1.0, order=4, **kwargs):
    """Bandpass filter.

    Filter the records between the frequencies `freqmin` and `freqmax`
    using a Butterworth filter and a zero-phase filter.

    Parameters
    ----------
    records : np.ndarray
        Records to filter.
    freqmin : float
        Minimum frequency.
    freqmax : float
        Maximum frequency.

    Returns
    -------
    np.ndarray
        Filtered records.
    """
    # Compute Nyquist frequency
    nyquist = 0.5 * sampling_rate

    # Compute normalized frequencies
    freqmin /= nyquist
    freqmax /= nyquist

    # Compute filter coefficients
    butter_filter = butter(order, [freqmin, freqmax], btype="band", **kwargs)

    # Filter records
    return filtfilt(*butter_filter, records)


def whiten(trace, freqmin, freqmax):
    """Whiten a trace.

    Whiten a trace between the frequencies `freqmin` and `freqmax` using a
    cosine taper. The taper is applied in the frequency domain.

    Parameters
    ----------
    trace : obspy.Trace
        Trace to whiten.
    freqmin : float
        Minimum frequency.
    freqmax : float
        Maximum frequency.

    Returns
    -------
    obspy.Trace
        Whiten trace.

    Notes
    -----
    This function was written by Celine Hadziioannou and Ashim Rijal, and was
    modified by Léonard Seydoux in 2023.
    """
    # Get basic information
    sampling_rate = trace.stats.sampling_rate
    n_samples = len(trace.data)

    # Whiten trace
    bandwidth = freqmax - freqmin
    smoothing_frequency = min(0.01, 0.5 * bandwidth)
    smoothing_len = int(np.fix(smoothing_frequency * n_samples / sampling_rate))

    # Frequencies
    frequencies = np.arange(n_samples) * sampling_rate / (n_samples - 1)
    indices = ((frequencies > freqmin) & (frequencies < freqmax)).nonzero()[0]

    # Signal spectra
    spectra = np.fft.fft(trace.data)
    spectra_white = np.zeros(n_samples, dtype="complex")

    # Apodization to the left with cos^2 (to smooth the discontinuities)
    smooth_left = np.cos(np.linspace(np.pi / 2, np.pi, smoothing_len + 1)) ** 2
    left = slice(indices[0], indices[0] + smoothing_len + 1)
    spectra_white[left] = smooth_left * np.exp(1j * np.angle(spectra[left]))

    # Boxcar
    center = slice(indices[0] + smoothing_len + 1, indices[-1] - smoothing_len)
    spectra_white[center] = np.exp(1j * np.angle(spectra[center]))

    # Apodization to the right with cos^2 (to smooth the discontinuities)
    smooth_right = np.cos(np.linspace(0, np.pi / 2, smoothing_len + 1)) ** 2
    right = slice(indices[-1] - smoothing_len, indices[-1] + 1)
    spectra_white[right] = smooth_right * np.exp(1j * np.angle(spectra[right]))

    # Inverse Fourier transform
    trace.data = 2 * np.fft.ifft(spectra_white).real

    return trace


def normalize(trace, clip_factor=6, clip_weight=10, method="onebit"):
    """Normalize a trace.

    Normalize a trace using one of the following methods:

    - "clipping": clip the trace at `clip_factor` times the standard deviation
        of the trace.
    - "clipping_iter": clip the trace at `clip_factor` times the standard
        deviation of the absolute value of the trace. Repeat until no values
        are left above the waterlevel.
    - "onebit": normalize the trace using the 1-bit normalization method.

    Parameters
    ----------
    trace : obspy.Trace
        Trace to normalize.
    clip_factor : float, optional
        Factor to multiply the standard deviation of the trace with to get the
        clipping level. Default is 6.
    clip_weight : float, optional
        Factor to divide the clipped values with. Default is 10.
    method : str, optional
        Normalization method. Default is "onebit".

    Returns
    -------
    obspy.Trace
        Normalized trace.

    Notes
    -----
    This function was written by Celine Hadziioannou and Ashim Rijal, and was
    modified by Léonard Seydoux in 2023.
    """
    if method == "clipping":
        threshold = clip_factor * np.std(trace.data)
        trace.data[trace.data > threshold] = threshold
        trace.data[trace.data < -threshold] = -threshold

    elif method == "clipping_iter":
        threshold = clip_factor * np.std(np.abs(trace.data))
        while trace.data[np.abs(trace.data) > threshold] != []:
            trace.data[trace.data > threshold] /= clip_weight
            trace.data[trace.data < -threshold] /= clip_weight

    elif method == "onebit":
        trace.data = np.sign(trace.data)

    return trace


def plot_trace_and_spectrum(trace, title=None):
    """Plot trace and spectrum.

    Parameters
    ----------
    trace : obspy.Trace
        Trace to plot.
    title : str, optional
        Title of the figure.
    """
    # Prepare figure
    _, ax = plt.subplots(2)

    # Plot trace
    ax[0].plot(trace.times(), trace.data)

    # Plot spectrum
    spectrum = np.abs(np.fft.rfft(trace.data))
    frequency = np.fft.rfftfreq(trace.data.size, d=trace.stats.delta)
    ax[1].plot(frequency[1:], spectrum[1:])

    # Labels
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Amplitude")
    ax[0].grid()
    ax[0].set_title(title)
    ax[1].set_xlabel("Frequency (Hz)")
    ax[1].set_ylabel("Spectrum")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].grid()

    # Show
    plt.tight_layout()
    plt.show()


def crosscorrelate(stream, stations, window_duration):
    """Compute the correlation between the records of two stations.

    Parameters
    ----------
    stream : obspy.Stream
        Stream containing the records.
    stations : list of str
        Stations to correlate. Must be in the stream.
    window_duration : float
        Duration of the correlation window in seconds.

    Returns
    -------
    np.ndarray
        Cross-correlation.
    np.ndarray
        Lags.

    Notes
    -----
    This function was written by Celine Hadziioannou and Ashim Rijal, and was
    modified by Léonard Seydoux in 2023.
    """
    # Extract traces
    trace_1 = stream.select(station=stations[0])[0]
    trace_2 = stream.select(station=stations[1])[0]

    # Starttime and endtime
    starttime = max(trace_1.stats.starttime, trace_2.stats.starttime)
    endtime = min(trace_1.stats.endtime, trace_2.stats.endtime)

    # Slice streams
    cross_correlation = list()
    while starttime + window_duration < endtime:
        # Get records
        record_1 = trace_1.slice(starttime, starttime + window_duration)
        record_2 = trace_2.slice(starttime, starttime + window_duration)

        # Cross-correlate
        xcorr = correlate(record_1.data, record_2.data, "same")
        cross_correlation.append(xcorr)

        # Shift timewindow by one correlation window length
        starttime += window_duration

    # Get lags
    lag_max = window_duration / 2
    dt = trace_1.stats.delta
    lags = np.arange(-lag_max, lag_max + dt, dt)

    # Stack cross-correlations
    cross_correlation = np.array(cross_correlation)
    cross_correlation = np.mean(cross_correlation, axis=0)

    return cross_correlation, lags
