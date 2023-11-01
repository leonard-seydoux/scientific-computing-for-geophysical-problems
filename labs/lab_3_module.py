# Functions
# collection of functions used in noise correlation processing

import numpy as np


def normalize(
    trace, clip_factor=6, clip_weight=10, norm_win=None, norm_method="1bit"
):

    if norm_method == "clipping":
        lim = clip_factor * np.std(trace.data)
        trace.data[trace.data > lim] = lim
        trace.data[trace.data < -lim] = -lim

    elif norm_method == "clipping_iter":
        lim = clip_factor * np.std(np.abs(trace.data))

        # as long as still values left above the waterlevel, clip_weight
        while trace.data[np.abs(trace.data) > lim] != []:
            trace.data[trace.data > lim] /= clip_weight
            trace.data[trace.data < -lim] /= clip_weight

    elif norm_method == "ramn":
        lwin = trace.stats.sampling_rate * norm_win
        st = 0  # starting point
        N = lwin  # ending point

        while N < trace.stats.npts:
            win = trace.data[st:N]

            w = np.mean(np.abs(win)) / (2.0 * lwin + 1)

            # weight center of window
            trace.data[st + lwin / 2] /= w

            # shift window
            st += 1
            N += 1

        # taper edges
        taper = get_window(trace.stats.npts)
        trace.data *= taper

    elif norm_method == "1bit":
        trace.data = np.sign(trace.data)
        trace.data = np.float32(trace.data)

    return trace


def get_window(N, alpha=0.2):

    window = np.ones(N)
    x = np.linspace(-1.0, 1.0, N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window


def whiten(tr, freqmin, freqmax):

    nsamp = tr.stats.sampling_rate

    n = len(tr.data)
    if n == 1:
        return tr
    else:
        frange = float(freqmax) - float(freqmin)
        nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
        f = np.arange(n) * nsamp / (n - 1.0)
        JJ = ((f > float(freqmin)) & (f < float(freqmax))).nonzero()[0]

        # signal FFT
        FFTs = np.fft.fft(tr.data)
        FFTsW = np.zeros(n) + 1j * np.zeros(n)

        # Apodization to the left with cos^2 (to smooth the discontinuities)
        smo1 = np.cos(np.linspace(np.pi / 2, np.pi, nsmo + 1)) ** 2
        FFTsW[JJ[0] : JJ[0] + nsmo + 1] = smo1 * np.exp(
            1j * np.angle(FFTs[JJ[0] : JJ[0] + nsmo + 1])
        )

        # boxcar
        FFTsW[JJ[0] + nsmo + 1 : JJ[-1] - nsmo] = np.ones(
            len(JJ) - 2 * (nsmo + 1)
        ) * np.exp(1j * np.angle(FFTs[JJ[0] + nsmo + 1 : JJ[-1] - nsmo]))

        # Apodization to the right with cos^2 (to smooth the discontinuities)
        smo2 = np.cos(np.linspace(0.0, np.pi / 2.0, nsmo + 1)) ** 2.0
        espo = np.exp(1j * np.angle(FFTs[JJ[-1] - nsmo : JJ[-1] + 1]))
        FFTsW[JJ[-1] - nsmo : JJ[-1] + 1] = smo2 * espo

        whitedata = 2.0 * np.fft.ifft(FFTsW).real

        tr.data = np.require(whitedata, dtype="float32")

        return tr


def correlateNoise(st, stations, corrwin):

    # initialize sliding timewindow (length = corrwin) for correlation
    # start 1 corrwin after the start to account for different stream lengths
    timewin = st.select(station=stations[1])[0].stats.starttime + corrwin

    # loop over timewindows
    # stop 1 corrwin before the end to account for different stream lengths
    while (
        timewin < st.select(station=stations[0])[-1].stats.endtime - 2 * corrwin
    ):
        sig1 = st.select(station=stations[0]).slice(timewin, timewin + corrwin)
        sig1.merge(method=0, fill_value=0)
        sig2 = st.select(station=stations[1]).slice(timewin, timewin + corrwin)
        sig2.merge(method=0, fill_value=0)
        xcorr = np.correlate(sig1[0].data, sig2[0].data, "same")

        try:
            # build array with all correlations
            corr = np.vstack((corr, xcorr))
        except:
            # if corr doesn't exist yet
            corr = xcorr

        # shift timewindow by one correlation window length
        timewin += corrwin
        lags = np.arange(
            -corrwin / 2, corrwin / 2 + st[0].stats.delta, st[0].stats.delta
        )

        # stack the correlations; normalize
        # stack = np.mean(corr, 0)
        # stack = stack / float((np.abs(stack).max()))

    return corr, lags
