from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import numpy.random
import sklearn.model_selection
import sklearn.neighbors

import trtools.utils.common as common

NOFFSETS = 4

def _get_bins(min_val: float,
              max_val: float,
              nbins: int,
              offset: int) -> np.ndarray:
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
                  strata_labels: Optional[List[str]] = None,
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
    strata_labels:
        if data is 2D, then a different label for each column in data
    randome_state:
        used to control the splitting in the cross validation.
    """
    # Handle case when data is all the same
    if np.all(data[~np.isnan(data)] == data.reshape(-1)[0]):
        fig, ax = plt.subplots()
        if len(data.shape) > 1:
            ax.hist(data, label=strata_labels)
            ax.legend()
        else:
            ax.hist(data)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Counts")
        plt.savefig(fname)
        return

    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    elif strata_labels is not None:
        assert len(strata_labels) == data.shape[1]

    n_strata = data.shape[1]
    flat_data = data.reshape(-1)
    flat_data = flat_data[~np.isnan(flat_data)]

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
            prob = 0
            for train_idxs, test_idxs in kfold.split(flat_data):
                train = flat_data[train_idxs]
                test = flat_data[test_idxs]
                min_val = np.min(train)
                max_val = np.max(train)
                bins = _get_bins(min_val, max_val, nbins, offset)
                train_hist, _ = np.histogram(train, bins=bins, density=True)
                test_hist, _ = np.histogram(test, bins=bins)
                prob += np.sum(np.multiply(train_hist, test_hist))
            if prob > best_prob:
                best = (nbins, offset)
                best_prob = prob

    min_val = np.min(flat_data)
    max_val = np.max(flat_data)
    best_bins = _get_bins(min_val, max_val, best[0], best[1])

    # plot using those parameters
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    if n_strata == 1:
        ax.hist(data[~np.isnan(data)], bins=best_bins)
    else:
        for stratum in range(n_strata):
            y, _ = np.histogram(data[:, stratum][~np.isnan(data[:, stratum])], bins=best_bins)
            bincenters = 0.5*(best_bins[1:] + best_bins[:-1])
            ax.plot(bincenters, y, '-', label=strata_labels[stratum])
    if n_strata > 1:
        ax.legend()
    plt.savefig(fname)


def PlotKDE(data: np.ndarray,
            xlabel: str,
            title: str,
            fname: str,
            strata_labels: Optional[List[str]] = None,
            random_state: int = 13):
    """
    Plots a kernel density estimation of the distribution.

    This is a smoother representation of the distribution
    than a historgram. Kernel bandwidth (which determins plot
    smoothness) is determined by cross validation (if more
    than 1000 loci, training is done on a subset of 1000
    chosen at random).

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
    strata_labels:
        if data is 2D, then a different label for each column in data
    randome_state:
        used to control the splitting in the cross validation.
    """
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    elif strata_labels is not None:
        assert len(strata_labels) == data.shape[1]

    if np.all(data == data[0, 0]):
        common.WARNING("Omitting graph {} because all the data points equal"
                       " {}".format(title, data[0, 0]))
        return

    n_strata = data.shape[1]

    # only train on up to 1k loci for speed
    if data.shape[0] > 1e3:
        rng = numpy.random.default_rng(random_state)
        train_data = data[rng.choice(1e3), :]
    else:
        train_data = data

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")

    # Fit and plot each stratum individually
    for col in range(n_strata):

        stratum = train_data[:, col]
        stratum = stratum[~np.isnan(stratum)]
        # Handle case when data is all the same
        if np.all(stratum == stratum[0]):
            common.WARNING("Omitting strata {} from graph {} because all the "
                           "data points equal {}".format(strata_labels[col],
                                                         title, data[0, col]))
            continue

        # code from
        # https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
        stratum = stratum.reshape(-1, 1)

        # Use gridsearch to choose the bandwidth
        kfold = sklearn.model_selection.KFold(
            n_splits=5,
            random_state=random_state,
            shuffle=True
        )
        bandwidths = 10 ** np.linspace(-1.8, 1.8, 20)
        grid = sklearn.model_selection.GridSearchCV(
            sklearn.neighbors.KernelDensity(kernel='gaussian'),
            {'bandwidth': bandwidths},
            cv=kfold
        )
        grid.fit(stratum)
        bandwidth = grid.best_params_['bandwidth']

        #compute the kde with the best bandwidth
        stratum = data[:, col] # now use all the data
        stratum = stratum[~np.isnan(stratum)]
        stratum = stratum.reshape(-1, 1)
        kde = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                              bandwidth=bandwidth)
        kde.fit(stratum)
        min_val = np.min(stratum)
        max_val = np.max(stratum)
        eps = (max_val - min_val)/10e3
        xs = np.arange(min_val - eps, max_val + eps, eps)
        curve = np.exp(kde.score_samples(xs.reshape(-1, 1)))

        # plot
        if n_strata == 1:
            ax.fill_between(xs, curve)
        else:
            ax.plot(xs, curve, label=strata_labels[col])
    if n_strata > 1:
        ax.legend()
    plt.savefig(fname)

