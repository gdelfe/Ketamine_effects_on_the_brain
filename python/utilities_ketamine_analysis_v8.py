# This example imports functions from the DemoReadSGLXData module to read
# digital data. The metadata file must be present in the same directory as the binary file.
# Works with both imec and nidq digital channels.

import pickle
from scipy.stats.stats import pearsonr
from scipy import stats
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy import stats
from scipy.interpolate import interp1d
from scipy.signal import gaussian, savgol_filter
from scipy.stats import linregress, spearmanr, kendalltau, pearsonr, circstd
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
from DemoReadSGLXData.readSGLX import readMeta, SampRate, makeMemMapRaw, ExtractDigital, Int2Volts
from datetime import datetime
from scipy.io import savemat, loadmat
from sklearn.manifold import Isomap
from scipy.optimize import curve_fit
from scipy.signal import chirp, peak_widths
from scipy.signal import find_peaks, savgol_filter, medfilt
# load pickle module
import time, os, math
from nptdms import TdmsFile

from scipy import signal

from joblib import Parallel, delayed

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Bitstream Vera Sans', 'DejaVu Sans', 'Lucida Grande',
                                   'Verdana', 'Geneva', 'Lucid', 'Avant Garde', 'sans-serif']


# for publication quality plots
def set_pub_plots(pal=sns.blend_palette(["gray", "crimson", 'cyan', 'magenta', 'purple'], 5)):
    sns.set_style("white")
    sns.set_palette('pastel')
    sns.set_context("poster", font_scale=1.5,
                    rc={"lines.linewidth": 2.5, "axes.linewidth": 2.5, 'figure.facecolor': 'white'})
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
    # optional, makes markers bigger, too, axes.linewidth doesn't seem to work
    plt.rcParams['axes.linewidth'] = 2.5


rc_pub = {'font.size': 25, 'axes.labelsize': 25, 'legend.fontsize': 25.0,
          'axes.titlesize': 25, 'xtick.labelsize': 25, 'ytick.labelsize': 25,
          # 'axes.color_cycle':pal, # image.cmap - rewritesd the default colormap
          'axes.linewidth': 2.5, 'lines.linewidth': 2.5,
          'xtick.color': 'black', 'ytick.color': 'black', 'axes.edgecolor': 'black', 'axes.labelcolor': 'black',
          'text.color': 'black'}
# to restore the defaults, call plt.rcdefaults()

# set_pub_bargraphs()
set_pub_plots()
sns.set_context('poster')
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False


# smooth hist with boxcar # max of 1 second
def bin_smooth_spk(df, binSize):
    trace_length = round(df.spk.max())
    boxcar = np.ones(int(1 / binSize))
    edges = np.arange(0, trace_length, binSize)
    ls_df = []
    for unit in df.unit_id.unique():
        spkt = df[df.unit_id == unit].spk
        hist, bins = np.histogram(spkt, bins=edges)
        spkd = np.convolve(hist, boxcar, mode='same')
        df_tmp = pd.DataFrame({'cids': unit, 'hz': spkd, 'times': bins[:-1]})
        ls_df.append(df_tmp)
    #     print(spkt)
    df_out = pd.concat(ls_df)
    return df_out


def NormalizeData(data):
    return (data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data))


# from scipy.spatial.distance import pdist
from ripser import ripser as tda
# from persim import plot_diagrams
# from sklearn import neighbors
# from sklearn.manifold import Isomap
# from mpl_toolkits.mplot3d import axes3d

# decoding from manifold, velocity ratio, HD anchoring

from scipy.stats import linregress, spearmanr, kendalltau, pearsonr, circstd


## Kendall Tau correlations parallized version ## you can also use pearson or any kind of pairwise computation
### features [cells x time]
def kendall_tau_parallel(features):
    import itertools
    import numpy as np
    import scipy.stats as ss
    from tqdm.notebook import tqdm
    import dask
    from dask.diagnostics import ProgressBar

    ## all cells
    tau_corrs_list = []
    tauVec = []
    corr_mat = []
    units_common = []

    units_common = np.arange(0, features.shape[0], 1)
    #     print(units_common)

    # keep pairs id
    cell_pairs = []
    cell_pairs = [(x1, x2) for x1, x2 in itertools.combinations(units_common, 2)]

    print(len(cell_pairs))

    # delayed_results = [dask.delayed(ss.pearsonr)(features[icell1,:],features[icell2,:],) for cell1, cell2 in tqdm(cell_pairs)]
    delayed_results = [dask.delayed(ss.kendalltau)(features[cell1, :], features[cell2, :], ) for cell1, cell2 in
                       tqdm(cell_pairs)]

    with ProgressBar():
        tau_corrs_list = dask.compute(delayed_results, scheduler='processes')
    tauVec = np.array(tau_corrs_list).copy()
    tau_corrs_list = [corr[0] for corr in tau_corrs_list[0]]  # remove p-values from list
    # put corrs list into matrix
    corr_mat = np.zeros((features.shape[0], features.shape[0]))

    for i, (cell1, cell2) in enumerate(cell_pairs):
        corr_mat[cell1, cell2] = tau_corrs_list[i]
        corr_mat[cell2, cell1] = tau_corrs_list[i]
    np.fill_diagonal(corr_mat, 1)

    return corr_mat, cell_pairs, tauVec


## Kendall Tau correlations parallized version ## you can also use pearson or any kind of pairwise computation
### features [cells x time]
def pearson_parallel(features):
    import itertools
    import numpy as np
    import scipy.stats as ss
    from tqdm.notebook import tqdm
    import dask
    from dask.diagnostics import ProgressBar

    ## all cells
    tau_corrs_list = []
    tauVec = []
    corr_mat = []
    units_common = []

    units_common = np.arange(0, features.shape[0], 1)
    #     print(units_common)

    # keep pairs id
    cell_pairs = []
    cell_pairs = [(x1, x2) for x1, x2 in itertools.combinations(units_common, 2)]

    print(len(cell_pairs))

    delayed_results = [dask.delayed(ss.pearsonr)(features[cell1, :], features[cell2, :], ) for cell1, cell2 in
                       tqdm(cell_pairs)]

    with ProgressBar():
        tau_corrs_list = dask.compute(delayed_results, scheduler='processes')
    tauVec = np.array(tau_corrs_list).copy()
    tau_corrs_list = [corr[0] for corr in tau_corrs_list[0]]  # remove p-values from list
    # put corrs list into matrix
    corr_mat = np.zeros((features.shape[0], features.shape[0]))

    for i, (cell1, cell2) in enumerate(cell_pairs):
        corr_mat[cell1, cell2] = tau_corrs_list[i]
        corr_mat[cell2, cell1] = tau_corrs_list[i]
    np.fill_diagonal(corr_mat, 1)

    return corr_mat, cell_pairs, tauVec


## IsoMap projection
## units in [cells x time], note: isomap algorithm takes the input as time x cells!!
## ndim = dimensions to where you want to project (10 for betti numbers computation) (3 for visualization in 3d plot)
## nneigh = number of neigbors, you will have to systematically increase this parameter and look at the 3d projections and evalute
## wheter they are changing as you increase the number of neigbors, it will depend on the time resolution of your neural activity
## as well as the time scale of the cognitive processes you are interested in observing
def isomap_proj(behavior, units, ndim, nneigh):
    from sklearn.manifold import Isomap

    tmp_sim = []
    _inp_l = []
    time = []

    _inp_l = units

    Xl = np.sqrt(_inp_l.T)
    print(Xl.shape)

    embedding = []
    embedding = Isomap(n_neighbors=nneigh, n_components=ndim)
    tmp_sim = embedding.fit_transform(Xl)

    time = np.linspace(0, units.shape[1], units.shape[1]) / 5

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    p3 = ax.scatter(tmp_sim[:, 0], tmp_sim[:, 1], tmp_sim[:, 2], c=behavior, s=10, cmap='viridis')
    ax.set_xlabel('Dim 1')
    ax.set_ylabel('Dim 2')
    ax.set_zlabel("Dim 3")
    # ax.axis('equal')
    cb_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(p3, cax=cb_ax, label = 'Dose')
    plt.show()

    behavior = []
    units = []

    return tmp_sim


## compute persistent homology using betti numbers
## input is an ISOMAP projection of your spiking activity to a 10-Dimensional ISOmap space
def compute_betti(dataIso):
    from scipy.spatial.distance import pdist
    from sklearn import neighbors
    from ripser import ripser as tda

    dataBetti = []
    H1_rates = []
    H2_rates = []
    barcodes = []
    dist = []
    rad = []
    neigh = []
    num_nbrs = []
    threshold = []
    thrsh_rates = []
    results = []
    idx = []

    Bettithrsh = False  # treshold nb of data point to h2MaxDt
    doBetti = True
    h2MaxDt = 700

    if Bettithrsh:
        #         print('neighbor thresholding 2')
        # a) find number of neighbors of each point within radius of 1st percentile of all
        # pairwise dist.
        dist = pdist(dataIso, 'euclidean')
        rad = np.percentile(dist, 1)
        neigh = neighbors.NearestNeighbors()
        neigh.fit(dataIso)
        num_nbrs = [*map(len, neigh.radius_neighbors(X=dataIso, radius=rad, return_distance=False))]

        # b) threshold out points with low density
        thrsh_prcnt = 20
        threshold = np.percentile(num_nbrs, thrsh_prcnt)
        thrsh_rates = dataIso[num_nbrs > threshold]
        dataBetti = thrsh_rates
    else:
        dataBetti = dataIso

    #     print('computing Betti numbers...')
    results = {'h0': [], 'h1': [], 'h2': []}
    # Betti
    # H0 & H1
    H1_rates = dataBetti
    #     print('h0-1')
    barcodes = tda(H1_rates, maxdim=1, coeff=2)['dgms']
    results['h0'] = barcodes[0]
    results['h1'] = barcodes[1]
    if len(dataBetti) > h2MaxDt:
        #         print('shortening data for Betti nb 2')
        idx = np.random.choice(np.arange(len(dataBetti)), h2MaxDt, replace=False)
        H2_rates = dataBetti[idx]
    else:
        H2_rates = dataBetti
    #     print('h2')
    barcodes = tda(H2_rates, maxdim=2, coeff=2)['dgms']
    results['h2'] = barcodes[2]
    #     print('done')

    plot_barcode = True

    h0, h1, h2 = results['h0'], results['h1'], results['h2']
    # replace the infinity bar (-1) in H0 by a really large number (this will depend on the spread of your data in isomap space)
    h0[~np.isfinite(h0)] = 50
    h1[~np.isfinite(h1)] = 50
    h2[~np.isfinite(h2)] = 50
    # Plot the longest barcodes only
    plot_prcnt = [99, 98, 90]  # order is h0, h1, h2

    if plot_barcode:
        print('plotting...')
        col_list = ['r', 'g', 'm', 'c']
        to_plot = []
        for curr_h, cutoff in zip([h0, h1, h2], plot_prcnt):
            bar_lens = curr_h[:, 1] - curr_h[:, 0]
            if len(curr_h) > 0:
                plot_h = curr_h[bar_lens > np.percentile(bar_lens, cutoff)]
                to_plot.append(plot_h)

        for curr_betti, curr_bar in enumerate(to_plot):
            #             ax = fig.add_subplot()
            figBetti = plt.figure(figsize=(15, 15))
            ax = figBetti.add_subplot(3, 9, curr_betti * 9 + 1)
            for i, interval in enumerate(reversed(curr_bar)):
                ax.plot([interval[0], interval[1]], [i, i], '.-', color=col_list[curr_betti],
                        lw=1.5)
            # ax.set_xlim([0, xlim])
            # ax.set_xticks([0, xlim])
            ax.set_ylim([-1, 55])  # len(curr_bar)])
            ax.set_xlim([-1, 55])
            # ax.set_yticks([])
        plt.show()
        print('plotted')

    data_betti = h0, h1, h2

    betti_data = []
    h1_ratio = []
    h1_max = []
    h1 = []

    return data_betti


## Linear Discriminant projections
## units in [cells x time]
## behavior, in this case is position [X,Y]
def lda_projections(units, behavior):
    ## here I load many packages that I do not use, but helpful in case you want to switch dimensionality reduction method
    from sklearn.model_selection import train_test_split
    from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.neighbors import KNeighborsClassifier, NeighborhoodComponentsAnalysis
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import hamming_loss
    from sklearn.metrics import multilabel_confusion_matrix
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.metrics import r2_score
    from sklearn.metrics import median_absolute_error
    from sklearn.manifold import TSNE
    from scipy.ndimage import gaussian_filter1d

    lda_coeffs = []
    lda_cov = []

    LDA_data = []
    X_LDA = []
    X_LDA_tmp = []
    X_to_LDA = []

    X_LDA_tmp = gaussian_filter1d(units, 1, 1)  # gaussian filter to calcium imaging data
    X_to_LDA = np.transpose(X_LDA_tmp) / np.std(X_LDA_tmp, axis=1)
    X_LDA = np.nan_to_num(np.sqrt(X_to_LDA.copy()))

    X = []
    Y = []

    X = behavior[:, 0]
    Y = behavior[:, 1]

    occ = []
    occ = stats.binned_statistic_2d(X, Y, None, 'count', bins=8)

    Pos_2d = []
    Pos_2d = occ.binnumber

    print(X_LDA.shape)
    lda = []
    lda = LinearDiscriminantAnalysis(n_components=3, store_covariance=True)
    LDA_data = lda.fit_transform(X_LDA, Pos_2d)

    lda_coeffs = lda.coef_
    lda_cov = lda.covariance_

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    tmp_emb = []
    tmp_emb = LDA_data

    # Creating plot
    p2 = ax.scatter3D(tmp_emb[:, 0], tmp_emb[:, 1], tmp_emb[:, 2], c=Pos_2d, s=20, cmap='viridis')
    cb_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(p2, cax=cb_ax, label='Linearized bin position')
    ax.set_xlabel('Dim 1')
    ax.set_ylabel('Dim 2')
    ax.set_zlabel('Dim 3')
    plt.show()

    return LDA_data, lda_coeffs, lda_cov


def kendall_visualization_PCo(data1, data2):
    PCo_out = []

    try:
        a, b = np.polyfit(np.nan_to_num(data1[:, 0]), np.nan_to_num(data2[:, 0]), 1)

    except:
        a = 0
        b = 0

    PCo_out = pearsonr(np.nan_to_num(data1[:, 0]), np.nan_to_num(data2[:, 0]))

    fig = plt.figure(figsize=(5, 5))
    # Creating plot
    plt.scatter(data1[:, 0], data2[:, 0], s=50, color='orange', alpha=0.3)
    x = np.linspace(-1, 1, 10)
    plt.plot(x, x, '--k')
    plt.plot(x, a * x + b, linewidth=6, color='orange', alpha=0.5)
    plt.xlabel('Pairwise correlations - Condition 1')
    plt.ylabel('Pairwise correlations - Condition 2')
    # plt.xlim(0,1)
    plt.show()

    print('PCo', PCo_out.statistic)

    return PCo_out


def kendall_visualization_PCo_doses(data1, data2, data3, data4):
    PCo_out = []

    try:
        a, b = np.polyfit(np.nan_to_num(data1[:, 0]), np.nan_to_num(data2[:, 0]), 1)

    except:
        a = 0
        b = 0

    try:
        a2, b2 = np.polyfit(np.nan_to_num(data1[:, 0]), np.nan_to_num(data3[:, 0]), 1)

    except:
        a2 = 0
        b2 = 0

    try:
        a3, b3 = np.polyfit(np.nan_to_num(data1[:, 0]), np.nan_to_num(data4[:, 0]), 1)

    except:
        a3 = 0
        b3 = 0

    PCo_out = [pearsonr(np.nan_to_num(data1[:, 0]), np.nan_to_num(data2[:, 0])),
               pearsonr(np.nan_to_num(data1[:, 0]), np.nan_to_num(data3[:, 0])),
               pearsonr(np.nan_to_num(data1[:, 0]), np.nan_to_num(data4[:, 0]))]

    fig = plt.figure(figsize=(5, 5))
    # Creating plot
    plt.scatter(data1[:, 0], data2[:, 0], s=20, color='green', alpha=0.3)
    plt.scatter(data1[:, 0], data3[:, 0], s=20, color='orange', alpha=0.3)
    plt.scatter(data1[:, 0], data4[:, 0], s=20, color='blue', alpha=0.3)

    x = np.linspace(-1, 1, 10)
    plt.plot(x, x, '--k')
    plt.plot(x, a * x + b, linewidth=6, color='green', alpha=0.5, label='Dose 1')
    plt.plot(x, a2 * x + b2, linewidth=6, color='orange', alpha=0.5, label='Dose 2')
    plt.plot(x, a3 * x + b3, linewidth=6, color='blue', alpha=0.5, label='Dose 3')

    plt.xlabel('Pairwise correlations - Baseline')
    plt.ylabel('Pairwise correlations')
    plt.legend()
    # plt.xlim(0,1)
    plt.show()

    print('PCo', PCo_out[0].statistic, PCo_out[1].statistic, PCo_out[2].statistic)

    return PCo_out

def kendall_visualization_PCo_doses_3D(data1, data2, data3, data4, data5, data6, data7, data8, data9):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1, projection = '3d')
    # Creating plot
    ax.scatter(data1[:, 0], data2[:, 0],data3[:, 0], s=20, color='green', alpha=0.3)
    ax.scatter(data4[:, 0], data5[:, 0],data6[:, 0], s=20, color='orange', alpha=0.3)
    ax.scatter(data7[:, 0], data8[:, 0],data9[:, 0], s=20, color='red', alpha=0.3)
    ax.set_xlabel('Baseline')
    ax.set_ylabel('Dose 1')
    ax.set_zlabel('Dose 3')
    ax.legend()
    plt.show()

    return



def getMeanBinnedVector_2d(matrix, bins, return_format="array"):
    """
    Breaks up a matrix into bins along the time axis (axis=1). Takes the mean of each bin along the same axis and returns the new matrix
    """
    binned_matrices = np.array_split(matrix, bins, axis=0)
    means = [np.mean(bin, axis=0) for bin in binned_matrices]
    if return_format == "list":
        return means
    elif return_format == "array":
        return np.array(means).T  # returns as Cells x Time
    else:
        raise Exception(
            "Invalid 'return_format' argument. Must be one of 'list' or 'array'."
        )


## functions for place cell properties

## useful functions from CircleSquare manuscript
def neighbor_sum(rate_map):
    import numpy as np
    import scipy.ndimage

    mask = np.isnan(rate_map)
    rate_map[mask] = 0

    conv_rate_map = scipy.ndimage.convolve(rate_map, np.ones([3, 3]), mode='constant')

    conv_rate_map[mask] = np.nan
    rate_map[mask] = np.nan
    out_tmp = np.corrcoef(rate_map[~mask].flat, conv_rate_map[~mask].flat)[0, 1]

    # rate_map = (rate_map) / np.std(rate_map[~np.isnan(rate_map)])
    # out_tmp = []
    # for index, val in np.ndenumerate(rate_map):
    #     if not np.isnan(val):
    #         map_tmp = rate_map.copy()
    #         map_tmp[index] = 0
    #         map_tmp[np.isnan(rate_map)] = 0
    #         sum_neigh = scipy.ndimage.convolve(map_tmp, np.ones([3, 3]), mode='constant')[index]
    #
    #         map_norm = np.ones(rate_map.shape)
    #         map_norm[index] = 0
    #         map_norm[np.isnan(rate_map)] = 0
    #         norm = scipy.ndimage.convolve(map_norm, np.ones([3, 3]), mode='constant')[index]
    #         if norm > 0:
    #             out_tmp.append(abs(val - sum_neigh/norm))

    # print(out_tmp)
    return np.mean(out_tmp)


def place_info_content(occ, rate_map):
    import numpy as np
    T = np.sum(occ)  # total time spent
    R = (np.sum(occ[~np.isnan(rate_map)] * rate_map[~np.isnan(rate_map)])) / T  # mean rate
    info_map = (occ / T) * (rate_map / R) * np.log2(rate_map / R)

    return np.sum(info_map[~np.isnan(info_map)])


def nan_gaussian_filter(map, sig):
    import scipy.ndimage
    import numpy as np

    V = map.copy()
    V[np.isnan(map)] = 0
    V_filtered = scipy.ndimage.gaussian_filter(V, sig)

    W = 0 * map.copy() + 1
    W[np.isnan(map)] = 0
    W_filtered = scipy.ndimage.gaussian_filter(W, sig)

    mask = W_filtered < 0.4

    ratio = V_filtered / W_filtered
    out = scipy.ndimage.gaussian_filter(ratio, sig)
    out[mask] = np.nan
    return out


def read_timestamps(filename):
    file = open(filename, 'r')
    columns = (line.split(sep='\t') for line in file)
    next(columns)  # skip first line
    z = list(zip(*columns))
    frameNum = list(map(int, z[1]))
    sysClock = list(map(int, z[2]))
    file.close()
    return frameNum, sysClock


def Pos_process(filename):
    file = open(filename, 'r')
    uu = (line for line in file)
    arg = (str.split(sep='\t') for str in uu if not (str[0] == '%' or str[0] == '\t'))
    z = list(zip(*arg))

    NumFrame = list(map(int, z[0]))
    TimeFrame = list(map(int, z[1]))
    Xpos = list(map(int, z[2]))
    Ypos = list(map(int, z[3]))
    Sector = list(map(int, z[4]))
    Shock = list(map(int, z[5]))
    file.close()

    return NumFrame, TimeFrame, Xpos, Ypos, Sector, Shock


def center_gravity(rate_map):
    import numpy as np
    y = 0
    x = 0
    w = 0

    for index, val in np.ndenumerate(rate_map):
        if not (np.isnan(val)):
            x += index[0] * val
            y += index[1] * val
            w += val

    return x / w, y / w


def center_high_rate(occ, rate_map):
    import numpy as np
    y = 0
    x = 0
    w = 0

    T = np.sum(occ)  # total time spent
    R = (np.sum(occ[~np.isnan(rate_map)] * rate_map[~np.isnan(rate_map)])) / T  # mean rate

    thresh = R + 4 * np.nanstd(rate_map.flat)
    for index, val in np.ndenumerate(rate_map):
        if not (np.isnan(val)) and val > thresh:
            x += index[0] * val
            y += index[1] * val
            w += val

    if w == 0:
        return np.nan, np.nan
    else:
        return x / w, y / w


def center_high_rate_contiguous(rm):
    import numpy as np
    import scipy.ndimage

    rate_map = nan_gaussian_filter(rm.copy(), 0.5)

    thresh = 4 * np.nanstd(rate_map.flat)
    rate_map[rate_map > thresh] = 1
    rate_map[rate_map <= thresh] = 0

    mask = np.isnan(rate_map)
    rate_map[mask] = 0
    conv_rate_map = scipy.ndimage.convolve(rate_map, np.ones([3, 3]), mode='constant')
    conv_rate_map[mask] = np.nan
    rate_map[mask] = np.nan

    map_norm = np.ones(rate_map.shape)
    map_norm[mask] = 0
    norm = scipy.ndimage.convolve(map_norm, np.ones([3, 3]), mode='constant')
    norm[mask] = np.inf

    conv_rate_map == norm - 1

    y = 0
    x = 0
    w = 0

    for index, val in np.ndenumerate(zip(conv_rate_map, rate_map)):
        if val:
            x += index[0] * val
            y += index[1] * val
            w += val

    return x / w, y / w


def conv_nan(x, w):
    import numpy as np
    half_w = int(np.floor(w / 2))

    a = []
    for idx in range(len(x)):
        min_b = max(0, idx - half_w)
        max_b = min(len(x), idx + half_w + 1)
        a.append(np.nanmean(x[min_b:max_b]))

    return a


def fisherZ(r):
    import numpy as np
    z = 0.5 * (np.log(1 + r) - np.log(1 - r))

    return z


def loadPhysiology(day, animal, xlAll, pathAll, expt, frame):
    import numpy as np
    import h5py
    # set up path
    xl = xlAll[xlAll.day == day]
    date = np.unique(xl.date)[0]
    print(date)
    TrackerPath = pathAll + 'tracker_files/' + date + '/'
    MnsTSPath = pathAll + 'timestamps/' + date + '/'
    MatPath = pathAll + 'mat/' + animal + '_' + date + '_out_1b_ss3' + '.mat'

    # extract recording info
    VideoName = list(xl.videoname)
    Sessions = list(xl.session)
    DayCount = list(xl.day_idx)[0]
    TrackerName = list(xl.trackername)
    doAnalysis = list(xl.analyze)
    Training = list(xl.training)[0]
    if expt == 'PlaceAvoidance' and frame == 'Arena':
        TrackerName = list(xl.arena)

    if expt == 'ATN':
        AtnSplit = [list(xl.start), list(xl.soundStart),
                    list(xl.foodStart), list(xl.foodEnd),
                    list(xl.soundEnd), list(xl.end)]
    else:
        AtnSplit = []

    MnsTSFileName = [MnsTSPath + s + '_timestamp.dat' for s in VideoName]
    TrackerFileName = [TrackerPath + str(s) for s in TrackerName]

    # check ordering ?!

    # evaluate length of each recording to later extract from longer file
    rec_idx = []
    rec_len = 0
    for fn in MnsTSFileName:
        VideoFrameNum, sysClock = read_timestamps(fn)
        nExtraTiff = (len(VideoFrameNum) - 1) // 10000
        # because of segmentation : frame 2 to 10001 (first and last frame always lost) then remaining modulo 3
        rec_idx.append([rec_len, rec_len + ((len(VideoFrameNum) - 1 - 10000 * nExtraTiff) // 3) + 3333 * nExtraTiff])
        rec_len += ((len(VideoFrameNum) - 1 - 10000 * nExtraTiff) // 3) + 3333 * nExtraTiff

    # load physiology
    f = h5py.File(MatPath, 'r')
    S_all = np.array(f['output']['S'])
    # C_all = np.array(f['output']['C'])
    A_tmp = np.array(f['output']['A'])

    xlOut = [VideoName, DayCount, Sessions, doAnalysis, Training, AtnSplit]
    fileNames = [MnsTSFileName, TrackerFileName]
    physio = [S_all, A_tmp, rec_idx]

    return [xlOut, fileNames, physio]


def runsTest(series, s0):
    import numpy as np
    runs = 0
    for i in range(len(series)):
        if (series[i] > s0 > series[i - 1]) or (series[i] < s0 < series[i - 1]):
            runs += 1

    n1 = np.sum([val > s0 for val in series])
    n2 = np.sum([val < s0 for val in series])

    runs_exp = ((2 * n1 * n2) / (n1 + n2)) + 1
    stan_dev = np.sqrt((2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / (((n1 + n2) ** 2) * (n1 + n2 - 1)))

    z = (runs - runs_exp) / stan_dev

    return z


def place_cell_properties(spk_data, beh_data, binsxy, fsrate):
    occ, rate_maps_temp, p_info, coh = [], [], [], []
    nRepeat = 100  # number of shuffles

    # COMPUTING THE RATE MAPS
    occ = np.histogram2d(beh_data[0], beh_data[1], bins=binsxy)[0] / fsrate
    rate_maps_temp = [np.histogram2d(beh_data[0], beh_data[1], bins=binsxy, weights=s)[0]
                      / occ for s in spk_data]

    # analyze rate maps
    p_info = [place_info_content(occ, m) for m in rate_maps_temp]
    coh = [neighbor_sum(m) for m in rate_maps_temp]

    ## pinfo and coherence
    #     print('Place info', np.hstack(p_info))
    #     print('Coherence', np.hstack(coh))

    # Randomizing
    sTrackedCopy = spk_data.copy()
    sTrackedTimePairs = list(zip(*sTrackedCopy))

    pInfoShuffled = []
    cohShuffled = []
    for iRand in range(nRepeat):
        np.random.shuffle(sTrackedTimePairs)
        sShuffled = list(zip(*sTrackedTimePairs))

        RateMapsShuffled = [np.histogram2d(beh_data[0], beh_data[1], bins=binsxy, weights=s)[0]
                            / occ for s in sShuffled]
        pInfoShuffled.append([place_info_content(occ, m) for m in RateMapsShuffled])

        cohShuffled.append([neighbor_sum(m) for m in RateMapsShuffled])

    pInfoTest = []
    for p, z in zip(p_info, zip(*pInfoShuffled)):
        Out = [p, (p - np.mean(z)) / np.std(z)]
        pInfoTest.append([Out, Out[0] > 1 and Out[1] > 1.96])

    cohTest = []
    for p, z in zip(coh, zip(*cohShuffled)):
        Out = [p, (p - np.mean(z)) / np.std(z)]
        testTmp = Out[0] > 0.5 and Out[1] > 1.96
        cohTest.append([Out, testTmp])

    # determine whether cells are place cells
    PC_test = [c[0][1] > 1.96 and p[0][1] > 1.96 for c, p in zip(cohTest, pInfoTest)]
    print('number of PCs', np.sum(PC_test))

    if np.sum(PC_test) > 5:

        ## plot firing rate maps
        for icell in np.random.choice(np.hstack(np.argwhere(np.hstack(PC_test) == True)), 5):
            print('place information', p_info[icell])
            print('coherence', coh[icell])

            fig = plt.figure(figsize=(2, 2))
            ax = fig.add_subplot(1, 1, 1)
            sns.heatmap(rate_maps_temp[icell], cmap='viridis', ax=ax, cbar_kws={'label': 'Normalized rate'})
            ax.set_xlabel('X position')
            ax.set_ylabel('Y position')
            fname = 'Cell' + str(icell)
            ax.set_title(fname)
            plt.show()

    return occ, rate_maps_temp, p_info, coh, PC_test


#####

def bin_spk(df, binSize):
    trace_length = round(df.spk.max())
    edges = np.arange(0, trace_length, binSize)
    ls_df = []
    for unit in df.unit_id.unique():
        spkt = df[df.unit_id == unit].spk
        hist, bins = np.histogram(spkt, bins=edges)
        spkd = hist.copy()
        df_tmp = pd.DataFrame({'cids': unit, 'hz': spkd, 'times': bins[:-1]})
        ls_df.append(df_tmp)
    #     print(spkt)
    df_out = pd.concat(ls_df)
    return df_out

####

def get_theta_ch_v2(data_lfp,lim):

    from scipy.signal import butter, lfilter, freqz, iirnotch, welch 

    theta_power_all = []
    for ich in range(0,data_lfp.shape[1]):
        freqs, psd = welch(data_lfp[:,ich], fs=2500, nperseg=1024, scaling='density')

        theta_band_freqs = np.argwhere(((freqs < 12) & (freqs > 4)) == True)
        theta_power = np.sum(psd[theta_band_freqs])
    #     print('theta_freqs', freqs[theta_band_freqs])
    #     print(theta_power)
        theta_power_all.append(theta_power)
        
    thres = []
    thres = np.percentile(np.hstack(theta_power_all),lim)
#     print('psd threshold', thres)

    theta_ch = []
    theta_ch = np.hstack(np.argwhere(np.hstack(theta_power_all) > thres))
#     print('theta channels', np.arange(theta_ch[0],theta_ch[-1]))
      
    return np.arange(theta_ch[0],theta_ch[-1]), np.hstack(theta_power_all)[np.arange(theta_ch[0],theta_ch[-1])]


####

def get_swr_ch_v2(data_lfp_tmp,theta_channels):
    
    data_lfp = []    
    data_lfp  = data_lfp_tmp[:,theta_channels]

    from scipy.signal import butter, lfilter, freqz, iirnotch, welch 

    ripple_power_all = []
    for ich in range(0,data_lfp.shape[1]):
        freqs, psd = welch(data_lfp[:,ich], fs=2500, nperseg=1024, scaling='density')

        ripple_band_freqs = np.argwhere(((freqs > 250) & (freqs < 400)) == True)
        ripple_power = np.sum(psd[ripple_band_freqs])
    #     print('ripple_freqs', freqs[ripple_band_freqs])
#         print(ripple_power)
        ripple_power_all.append(ripple_power)

    thres = []
    thres = np.percentile(np.hstack(ripple_power_all),75)
#     print('swr - psd threshold', thres)

    ripple_ch = []
    ripple_ch = np.hstack(np.argwhere(np.hstack(ripple_power_all) > thres))
    non_ripple_ch = np.hstack(np.argwhere(np.hstack(ripple_power_all) <= thres))

    # print('ripple channels', np.arange(ripple_ch[0],ripple_ch[-1]))
#     print('ripple channels',theta_channels[ripple_ch])
    
    return np.hstack(ripple_power_all), [ripple_ch, theta_channels[ripple_ch]]

###

def theta_fir(ripple_data):
    
    n_channels = ripple_data.shape[1]
    print(n_channels)
    outputSignal = []
    fs_val = 2500
    
    theta_all = []
    raw_all = []
    ripple_all = []
    
    
    for ich in range(0,n_channels):
            
        signal_tmp = []
        signal_tmp = ripple_data[:,ich] - np.mean(ripple_data[:,ich])

        # reject any signal above 350 Hz
        sos_ini2 = signal.butter(4, 12, 'lp', fs=fs_val, output='sos')
        signal_zscore2 = signal.sosfiltfilt(sos_ini2, signal_tmp)

        # reject any signal below 200 Hz
        sos3 = signal.butter(4, 4, 'hp', fs=fs_val, output='sos')

        ripple_only = signal.sosfiltfilt(sos3, signal_zscore2)   
    
        ripple_all.append(ripple_only)
        
    
    ripple_all = np.vstack(ripple_all)    
        
    return ripple_all   

###

def csd_fir(ripple_data):
    
    n_channels = ripple_data.shape[1]
    print(n_channels)
    outputSignal = []
    fs_val = 2500
    
    theta_all = []
    raw_all = []
    ripple_all = []
    
    filter_60Hz = False
    outputSignal = []
    
    if filter_60Hz == True:
        
        for ich in range(0,n_channels):

            signal_tmp = []
            signal_tmp = ripple_data[:,ich] - np.mean(ripple_data[:,ich])

            # reject any signal above 350 Hz
            sos_ini2 = signal.butter(4, 100, 'lp', fs=fs_val, output='sos')
            signal_zscore2 = signal.sosfiltfilt(sos_ini2, signal_tmp)

            # reject any signal below 200 Hz
            sos3 = signal.butter(4, 5, 'hp', fs=fs_val, output='sos')
            ripple_only = signal.sosfiltfilt(sos3, signal_zscore2)   


            # Create/view notch filter
            notch_freq = 60.0  # Frequency to be removed from signal (Hz)
            quality_factor = 20.0  # Quality factor

            # Design a notch filter using signal.iirnotch
            b_notch, a_notch = signal.iirnotch(notch_freq, quality_factor, fs_val)

            # Apply notch filter to the noisy signal using signal.filtfilt
            outputSignal_tmp = []
            outputSignal_tmp = signal.filtfilt(b_notch, a_notch, ripple_only)
            
            ripple_all.append(outputSignal_tmp)
    
    if filter_60Hz == False:
        for ich in range(0,n_channels):

            signal_tmp = []
            signal_tmp = ripple_data[:,ich] - np.mean(ripple_data[:,ich])

            # reject any signal above 350 Hz
            sos_ini2 = signal.butter(4, 100, 'lp', fs=fs_val, output='sos')
            signal_zscore2 = signal.sosfiltfilt(sos_ini2, signal_tmp)

            # reject any signal below 200 Hz
            sos3 = signal.butter(4, 5, 'hp', fs=fs_val, output='sos')
            ripple_only = signal.sosfiltfilt(sos3, signal_zscore2)   
            ripple_all.append(ripple_only)
        
    
    ripple_all = np.vstack(ripple_all)    
        
    return ripple_all    
    
    
    ####
    
from scipy.interpolate import interp1d
import pdb

## speed aligned
## lfp aligned
## speed rate is Hz of behavior

def upsample_speed(speed, LFP, LFP_rate = 2500, speed_rate = 100):
    
    # time length for speed and LFP variables
    speed_T = np.linspace(0, len(speed)/ speed_rate, len(speed))
    LFP_T = np.linspace(0, len(LFP[:,0])/ LFP_rate, len(LFP[:,0])) # all channels have same length

    # interpolate speed variable based on LFP
    interpolator = interp1d(speed_T, speed, kind = 'linear', fill_value="extrapolate")
    speed_upsampled = interpolator(LFP_T) 
    
    return speed_upsampled

#####

# Create a speed mask. 
# True if speed < threshold at any point in a given time window of length win
# False if there is at least one value above threshold in the given time window

## input: upsampled speed
## period: 
# threshold in mm/s
# win is computation window for thresholding in samples (2500 is seconds)

def create_speed_mask(speed_up,win,thresh):
    mask = []
    speed_th = speed_up > thresh # speed mask above threshold
    cnt = 0 # count for high speed trials 
    for i in range(0,len(speed_th),win):
        data_win = speed_th[i:i+win]
        if np.mean(data_win): # if all the speed value are above threshold
            mask.extend([True]) # create True mask for that window
            cnt +=1
        else:
            mask.extend([False]) # create False mask for that window
            
    return np.array(mask)


####

def create_sta_mask(speed_up,win,thresh):
    mask = []
    speed_th = speed_up < thresh # speed mask above threshold
    cnt = 0 # count for high speed trials 
    for i in range(0,len(speed_th),win):
        data_win = speed_th[i:i+win]
        if np.all(data_win): # if all the speed value are above threshold
            mask.extend([True]) # create True mask for that window
            cnt +=1
        else:
            mask.extend([False]) # create False mask for that window
            
    return np.array(mask)

####

def HPC_lfp_CA1_DG_fisure_detection_v3(data_lfp_tmp,data_lfp_sta,samples_to_plot,scale):

    import matplotlib.pyplot as plt
    import numpy as np
    import icsd
    from scipy import io
    import neo
    import quantities as pq

    #patch quantities with the SI unit Siemens if it does not exist
    for symbol, prefix, definition, u_symbol in zip(
        ['siemens', 'S', 'mS', 'uS', 'nS', 'pS'],
        ['', '', 'milli', 'micro', 'nano', 'pico'],
        [pq.A/pq.V, pq.A/pq.V, 'S', 'mS', 'uS', 'nS'],
        [None, None, None, None, u'µS', None]):
        if type(definition) is str:
            definition = lastdefinition / 1000
        if not hasattr(pq, symbol):
            setattr(pq, symbol, pq.UnitQuantity(
                prefix + 'siemens',
                definition,
                symbol=symbol,
                u_symbol=u_symbol))
        lastdefinition = definition
    
    ## get theta channels
    theta_channels = []
    theta_channels = get_theta_ch_v2(data_lfp_tmp,50)

    data_lfp = []
    data_lfp = data_lfp_tmp[:,theta_channels]
    
    fs_val = 2500
    n_channels = data_lfp.shape[1]    
    #loading test data
    ## input is channel x time
    # test_data = io.loadmat(r'G:\My Drive\Fenton_Lab\Ketamine_analysis\Ketamine_data\Neurotar_HPC_ketamine\iCSD-master\test_data.mat')

    #prepare lfp data for use, by changing the units to SI and append quantities,
    #along with electrode geometry, conductivities and assumed source geometry
    lfp_data = csd_fir(data_lfp) * 1E-6 * pq.V 
    z_data = np.arange(20E-6,3860E-6,20E-6) * pq.m # [m]
    diam = 500E-6 * pq.m                              # [m]
    h = 20E-6 * pq.m                                 # [m]
    sigma = 0.3 * pq.S / pq.m                         # [S/m] or [1/(ohm*m)]
    sigma_top = 0.3 * pq.S / pq.m    

    # [S/m] or [1/(ohm*m)]

    z_data = z_data[theta_channels]

    # Input dictionaries for each method
    spline_input = {
        'lfp' : lfp_data,
        'coord_electrode' : z_data,
        'diam' : diam,
        'sigma' : sigma,
        'sigma_top' : sigma,
        'num_steps' : int(z_data.size),      # Spatial CSD upsampling to N steps
        'tol' : 1E-12,
        'f_type' : 'gaussian',
        'f_order' : (20, 5),
    }

    #Create the different CSD-method class instances. We use the class methods
    #get_csd() and filter_csd() below to get the raw and spatially filtered
    #versions of the current-source density estimates.
    csd_dict = dict(spline_icsd = icsd.SplineiCSD(**spline_input),)

    ## for plotting purposes 
    data_lfp2 = theta_fir(data_lfp)
    theta_filtered = np.hstack(data_lfp2[:,np.argmax(data_lfp2,axis=1)[0]])

    # fit spline to derivate of theta amplitude across channels, to find fisure
    from scipy.interpolate import UnivariateSpline

    data = theta_filtered.copy()
    spl = UnivariateSpline(np.arange(len(data)), np.gradient(data), k=5)
    spl.set_smoothing_factor(1000)

    # fig = plt.figure(figsize=(6,6))
    # ax = fig.add_subplot(1,1,1)
    # plt.plot(spl(np.arange(len(data))), label='Smooth Fct 1e3')
    # plt.legend(loc='lower left')
    # plt.show()

    max_idx = np.argmax(spl(np.arange(len(data))))    # local
    min_idx = np.argmin(spl(np.arange(len(data))))    # local
    max_abs = np.argmax(np.abs(spl(np.arange(len(data))))) # max max absolute, the one we use

    print('absolute max min method',np.hstack([max_idx, min_idx]),theta_channels[np.hstack([max_idx, min_idx])])    
    print('absolute abs max method',np.hstack([max_abs]),theta_channels[np.hstack([max_abs])])    

    # plot csd and raw lfp on top
    for method, csd_obj in list(csd_dict.items()):

        time_array = np.arange(0,data_lfp.shape[0])

        #plot
        fig = plt.figure(figsize=(24,12))
        ax = fig.add_subplot(1,1,1)

        for ich in range(0,n_channels):

            ax.plot(time_array,(data_lfp[:,ich]/scale + ich),color='black', linewidth = 0.5)

        #plot raw csd estimate
        csd = []
        csd = csd_obj.get_csd()

        #plot spatially filtered csd estimate
        csd = csd_obj.filter_csd(csd)

        im = ax.imshow(np.array(csd), origin='lower', vmin=-abs(csd).max()/3, \
              vmax=abs(csd).max()/3, cmap='jet', interpolation='nearest')
#         ax.plot
        ax.axis(ax.axis('tight'))
        ax.set_title(csd_obj.name + ', filtered')
        cb = plt.colorbar(im, ax=ax)
        cb.set_label('CSD (%s)' % csd.dimensionality.string)
        ax.set_ylabel('ch #')
        ax.set_xlabel('Samples @ 2500 Hz')        

        plt.show()

        max_idx = np.unravel_index(np.argmax(csd, axis=None), csd.shape)
        idx_plot = max_idx[1]

        if max_idx[1] < samples_to_plot:
            idx_plot = samples_to_plot

        if max_idx[1] > (data_lfp.shape[0]-samples_to_plot):
            idx_plot = data_lfp.shape[0] - samples_to_plot

        ## plot csd at a given time window, chosen at the beginning
        fig2 = plt.figure(figsize=(4,8))
        ax = fig2.add_subplot(1,1,1)

        for ich in range(0,n_channels):

            ax.plot(time_array[:samples_to_plot*2],(data_lfp2[ich,(idx_plot-samples_to_plot):(idx_plot+samples_to_plot)]/scale + ich),
                    color='black', linewidth = 0.25)


        im = ax.imshow(np.array(csd[:,(idx_plot-samples_to_plot):(idx_plot+samples_to_plot)]), origin='lower', vmin=-abs(csd).max()/3, \
              vmax=abs(csd).max()/3, cmap='jet', interpolation='nearest')
        ax.axis(ax.axis('tight'))
        ax.set_title(csd_obj.name + ', filtered')
        cb = plt.colorbar(im, ax=ax)
        cb.set_label('CSD (%s)' % csd.dimensionality.string)
        ax.set_ylabel('ch #')
        ax.set_xlabel('Samples @ 2500 Hz')      

        plt.show()

    # get fisure id
    fisure_id = []
    fisure_id = [np.hstack([max_abs]),theta_channels[np.hstack([max_abs])]]

    ###############

    ## get ripple band power

    swr_power, swr_band_ch = [], []
    swr_power, swr_band_ch = get_swr_ch_v2(data_lfp_sta,theta_channels)

    ## Print summary        

    print('Theta channels',theta_channels)
    print('Ripple channels',swr_band_ch[1])
    print('Brain region identification ...')
    print('CA1 = theta + high swr power, channels', theta_channels[:fisure_id[0][0]])
    print('CA1- swr power', 'whole layer',np.mean(swr_power[:np.hstack([max_abs])[0]]), 'pyramidale',
          np.mean(swr_power[swr_band_ch[0]]))

    print('Fisure = flip in theta sign, channel', fisure_id[1][0])

    print('DG, channels', theta_channels[fisure_id[0][0]:])
    print('DG- swr power', np.mean(swr_power[np.hstack([max_abs])[0]:]))

    return csd, fisure_id, data_lfp, theta_channels, [theta_channels[:fisure_id[0][0]], theta_channels[fisure_id[0][0]:]]

########


def plot_lfps_theta(ripple_data, time, window_idx,scale):
    
    from scipy import signal
    
    n_channels = ripple_data.shape[1]
    print(n_channels)
    outputSignal = []
    fs_val = 2500
    
    theta_all = []
    raw_all = []
    ripple_all = []
    
    
    for ich in range(0,n_channels):

        filter_60hz = False

        if filter_60hz:
            # Create/view notch filter
            notch_freq = 60.0  # Frequency to be removed from signal (Hz)
            quality_factor = 20.0  # Quality factor

            # Design a notch filter using signal.iirnotch
            b_notch, a_notch = signal.iirnotch(notch_freq, quality_factor, fs_val)

            # Apply notch filter to the noisy signal using signal.filtfilt
            outputSignal_tmp = signal.filtfilt(b_notch, a_notch, ripple_data[:,ich])

            # Create/view notch filter
            notch_freq = 120.0  # Frequency to be removed from signal (Hz)
            quality_factor = 20.0  # Quality factor

            # Design a notch filter using signal.iirnotch
            b_notch, a_notch = signal.iirnotch(notch_freq, quality_factor, fs_val)

            # Apply notch filter to the noisy signal using signal.filtfilt
            outputSignal.append(signal.filtfilt(b_notch, a_notch, outputSignal_tmp))

            # low pass butterworth filter at 280Hz
            #sos = signal.butter(10, 100, 'lp', fs=fs, output='sos', analog=False)        
            #outputSignal = signal.sosfilt(sos, outputSignal_tmp2)

        else:
            
            signal_tmp = []
            signal_tmp = ripple_data[:,ich] - np.mean(ripple_data[:,ich])

            # reject any signal above 12 Hz
            sos_ini = signal.butter(4, 12, 'lp', fs=fs_val, output='sos')
            signal_zscore = signal.sosfiltfilt(sos_ini, signal_tmp)

            # reject any signal above 12 Hz
            sos2 = signal.butter(4, 4, 'hp', fs=fs_val, output='sos')

            theta_only = signal.sosfiltfilt(sos2, signal_zscore)   
            
            
            # reject any signal above 350 Hz
            sos_ini2 = signal.butter(4, 350, 'lp', fs=fs_val, output='sos')
            signal_zscore2 = signal.sosfiltfilt(sos_ini2, signal_tmp)

            # reject any signal below 200 Hz
            sos3 = signal.butter(4, 200, 'hp', fs=fs_val, output='sos')

            ripple_only = signal.sosfiltfilt(sos3, signal_zscore2)   
    
        theta_all.append(theta_only) 
        ripple_all.append(ripple_only) 

        raw_all.append(signal_tmp)
        
    
    raw_all = np.vstack(raw_all)    
    print(raw_all.shape)
    
    raw_all -= np.mean(raw_all,axis=0)
    
    window_start = window_idx - time* fs_val #np.random.randint(0,outputSignal.size - fs*time,1)[0]
    # ripple_data[window_start[0]:window_start[0]+fs*60]
    fig = plt.figure(figsize=(48,40))
    
    time_array = np.linspace(0,time,fs_val*time)
    windows_end = window_start+fs_val*time
    
    for ich in range(0,n_channels):
    
        plt.plot(time_array,(raw_all[ich, window_start:windows_end]/scale + ich),color='blue', linewidth = 2)
#         plt.plot(time_array,(theta_all[ich][window_start:windows_end]/scale - ich),color='red', linewidth = 2)

#         plt.xlim(0,time)
#         label_y = 'Ch' + str(ich*4)
        
#         if ich == n_channels - 1:
#             plt.xlabel('Time [sec]')
#             #get current axes
#             ax = plt.gca()
#             #hide x-axis
#             ax.get_yaxis().set_visible(False)   
                        
#         else:
        plt.axis('off')    
    
    plt.show()
    
    fig = plt.figure(figsize=(48,40))
    ax = fig.add_subplot(1,1,1)
    
    for ich in range(0,n_channels):
    
#         plt.plot(time_array,(raw_all[ich, window_start:windows_end]/scale - ich),color='blue', linewidth = 2)
        ax.plot(time_array,theta_all[ich][window_start:windows_end]/scale + ich,color='red', linewidth = 2)
#         plt.plot(time_array,(ripple_all[ich][window_start:windows_end]/scale - ich),color='green', linewidth = 2)



#         plt.xlim(0,time)
#         label_y = 'Ch' + str(ich*4)
        
#         if ich == n_channels - 1:
#             plt.xlabel('Time [sec]')
#             #get current axes
#             ax = plt.gca()
#             #hide x-axis
#             ax.get_yaxis().set_visible(False)   
                        
#         else:
        ax.axis('off')
    
    
    ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)    
    ax_histy.plot(np.vstack(theta_all)[:,2000], np.arange(0,385,1), color = 'black', linewidth = 10)
    ax_histy.axis('off')
    
    plt.show()       
    
    fig = plt.figure(figsize=(48,40))
    ax = fig.add_subplot(1,1,1)

    for ich in range(0,n_channels):
    
#         plt.plot(time_array,(raw_all[ich, window_start:windows_end]/scale - ich),color='blue', linewidth = 2)
        ax.plot(time_array,ripple_all[ich][window_start:windows_end]/scale + ich,color='green', linewidth = 2)

#         plt.xlim(0,time)
#         label_y = 'Ch' + str(ich*4)
        
#         if ich == n_channels - 1:
#             plt.xlabel('Time [sec]')
#             #get current axes
#             ax = plt.gca()
#             #hide x-axis
#             ax.get_yaxis().set_visible(False)   
                        
#         else:
        ax.axis('off')

    ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)    
    ax_histy.plot(np.vstack(ripple_all)[:,1000], np.arange(0,385,1), color = 'black', linewidth = 10)
    ax_histy.axis('off')    
    
    plt.show()          
        
        
    return raw_all, ripple_all, theta_all    

####
    
    # prelimiary spectral analysis
from scipy.signal import butter, lfilter, freqz, iirnotch, welch 

def bandpass_filter(data, lowcut, highcut, fs=2500, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')

    filtered_data = np.empty_like(data)

    # Loop through each channel and apply the filter
    for i in range(data.shape[1]):
        filtered_data[:, i] = lfilter(b, a, data[:, i])

    return filtered_data

def notch_filter(data, notch_freq, fs=2500, Q=30):
    b, a = iirnotch(notch_freq, Q, fs)
    filtered_data = np.empty_like(data)

    # Loop through each channel and apply the filter
    for i in range(data.shape[1]):
        filtered_data[:, i] = lfilter(b, a, data[:, i])

    return filtered_data


    
    