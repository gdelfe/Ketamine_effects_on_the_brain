import matplotlib.pyplot as plt


def spike_phase_hist(hist, ax=None, ticks=True, tickstride=5):
    if ax is not None:
        plt.sca(ax)
    plt.imshow(hist.values, interpolation='none')
    plt.axis('auto')
    if not ticks:
        return
    freqtick = range(len(hist.columns))[::tickstride]
    freqlabel = np.round(hist.columns, 1)[::tickstride]
    phasetick = np.arange(5) / 4 * (hist.shape[0] - 1)
    phaselabel = ['$-\pi$', '$-\pi/2$', 0, '$\pi/2$', '$-\pi$']
    plt.xticks(freqtick, freqlabel)
    plt.yticks(phasetick, phaselabel)
    return plt.gcf()
