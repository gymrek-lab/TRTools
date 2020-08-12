from typing import List

import matplotlib.pyplot as plt
import numpy as np
import sklearn.model_selection
import sklearn.neighbors

NOFFSETS = 4

def _get_bins(min_val: float,
              max_val: float,
              nbins: int,
              offset: int) -> List[float]:
    """
    Get an array of bin endpoints with the corresponding number and offset.

    The bins will extend a bit to either side of min_val and max_val
    to make sure they are covered. The offset determines how much to one
    side or the other.

    Parameters
    ----------
    min_val:
        The min value the histogram should cover.
    max_val:
        The max value the histogram should cover.
    nbins:
        The number of bins in the histogram
    offset:
        An int between 1 and NOFFSETS. 1 means most over the
        overhang (and subsequent offsetting of bins) is to the left of min_val,
        NOFFSETS means most of the overhang (and prior offsetting of bins) is to
        the right of max_val.

    Returns
    -------
    bins: List[float]
        A List of bin endpoints. length = len(nbins) + 1
    """

    assert 1 <= offset and offset <= NOFFSETS
    eps = (max_val - min_val)/10e3
    binrange = (max_val - min_val + 2*eps)*(nbins + 1)/nbins
    binsize = binrange/nbins
    start = min_val - eps - binsize*(1 - offset/NOFFSETS)
    return np.arange(start, start + binsize*nbins + eps, binsize)


def PlotHistogram(data: np.ndarray,
                  xlabel: str,
                  title: str,
                  fname: str,
                  random_state: int = 13):
    """
    Plot a histogram with learned bin sizes and alignment.

    Learn using 5 fold cross validation

    Parameters
    ----------
    data:
        Either a 1D array of statistics to create a histogram of,
        or a 2D array where each column represents a different stratification
        of the data that will get a different line in the plot
    xlabel:
        the x label for the graph
    title:
        the title for the graph
    fname:
        the file name to save the graph. Must include the extension,
        and one that matplotlib will recognize so that it produces
        a file of that type.
    randome_state:
        used to control the splitting in the cross validation.
    """
    # Handle case when data is all the same
    if np.all(data == data[0]):
        fig, ax = plt.subplots()
        ax.hist(data)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Counts")
        plt.savefig(fname)
        return

    # Get the best nbins and offset
    MAX_BIN = 30
    kfold = sklearn.model_selection.KFold(
        n_splits=5,
        random_state=random_state,
        shuffle=True
    )
    best = (1, 1)
    best_prob = 0
    for nbins in range(1, MAX_BIN + 1):
        for offset in range(1, NOFFSETS + 1):
            print("Testing fit {} of {}".format(
                (nbins-1)*NOFFSETS+ offset, MAX_BIN * NOFFSETS),
                end="\r", flush=True)
            prob = 0
            for train_idxs, test_idxs in kfold.split(data):
                train = data[train_idxs]
                test = data[test_idxs]
                min_val = np.min(train)
                max_val = np.max(train)
                bins = _get_bins(min_val, max_val, nbins, offset)
                train_hist, _ = np.histogram(train, bins=bins, density=True)
                test_hist, _ = np.histogram(test, bins=bins)
                prob += np.sum(np.multiply(train_hist, test_hist))
            if prob > best_prob:
                best = (nbins, offset)
                best_prob = prob

    print("Done fitting. Now plotting         ", end="\r", flush=True)
    min_val = np.min(data)
    max_val = np.max(data)
    best_bins = _get_bins(min_val, max_val, best[0], best[1])

    # plot using those parameters
    fig, ax = plt.subplots()
    ax.hist(data, bins=best_bins)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    plt.savefig(fname)
    print("Done plotting                  ", flush=True)


def PlotKDE(data: np.ndarray,
            xlabel: str,
            title: str,
            fname: str,
            random_state: int = 13):
    """
    Plots a kernel density estimation of the distribution.

    This is a smoother representation of the distribution
    than a historgram. Kernel bandwidth (which determins plot
    smoothness) is determined by cross validation

    Parameters
    ----------
    data:
        Either a 1D array of statistics to create a histogram of,
        or a 2D array where each column represents a different stratification
        of the data that will get a different line in the plot
    xlabel:
        the x label for the graph
    title:
        the title for the graph
    fname:
        the file name to save the graph. Must include the extension,
        and one that matplotlib will recognize so that it produces
        a file of that type.
    randome_state:
        used to control the splitting in the cross validation.
    """
    # Handle case when data is all the same
    if np.all(data == data[0]):
        fig, ax = plt.subplots()
        ax.hist(
            np.array(data[0]),
            bins=np.arange(data[0] - 1.5, data[0] + 2.5, 1)
        )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Percentage of loci")
        plt.savefig(fname)
        return

    # code from
    # https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
    data = data.reshape(-1, 1)

    # Use gridsearch to choose the bandwidth
    kfold = sklearn.model_selection.KFold(
        n_splits=5,
        random_state=random_state,
        shuffle=True
    )
    bandwidths = 10 ** np.linspace(-2, 2, 20)
    grid = sklearn.model_selection.GridSearchCV(
        sklearn.neighbors.KernelDensity(kernel='gaussian'),
        {'bandwidth': bandwidths},
        cv=kfold,
        verbose=1
    )
    grid.fit(data)
    bandwidth = grid.best_params_['bandwidth']

    #compute the kde with the best bandwidth
    kde = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                          bandwidth=bandwidth)
    kde.fit(data)
    min_val = np.min(data)
    max_val = np.max(data)
    eps = (max_val - min_val)/10e3
    xs = np.arange(min_val - eps, max_val + eps, eps)
    curve = np.exp(kde.score_samples(xs.reshape(-1, 1)))

    # plot
    fig, ax = plt.subplots()
    ax.fill_between(xs, curve)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")
    plt.savefig(fname)


