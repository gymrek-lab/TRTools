import matplotlib.pyplot as plt
import numpy as np
import sklearn.model_selection
import sklearn.neighbors

NOFFSETS = 4

def _get_bins(min_val, max_val, nbins, offset):
    assert 1 <= offset and offset <= NOFFSETS
    eps = (max_val - min_val)/10e3
    binrange = (max_val - min_val + 2*eps)*(nbins + 1)/nbins
    binsize = binrange/nbins
    start = min_val - eps - binsize*(1 - offset/NOFFSETS)
    return np.arange(start, start + binsize*nbins + eps, binsize)


def PlotHistogram(data, xlabel, title, fname, random_state=13):
    """
    Plot a histogram with learned bin sizes and alignment.

    Learn parameters with 5 fold cross validation
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
            prob = 0
            for train, test in kfold.split(data):
                min_val = np.min(train)
                max_val = np.max(train)
                bins = _get_bins(min_val, max_val, nbins, offset)
                train_hist, _ = np.histogram(train, bins=bins, density=True)
                test_hist, _ = np.histogram(test, bins=bins)
                prob += np.sum(np.multiply(train_hist, test_hist))
            if prob > best_prob:
                best = (nbins, offset)
                best_prob = prob

    min_val = np.min(data)
    max_val = np.max(data)
    best_bins = _get_bins(min_val, max_val, best[0], best[1])

    fig, ax = plt.subplots()
    ax.hist(data, bins=best_bins)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    plt.savefig(fname)


def PlotKDE(data, xlabel, title, fname, random_state=13):
    """
    Plots a kernel density estimation of the distribution.
    
    This is a smoother representation of the distribution
    than a historgram.
    Kernel bandwidth (which determins plot smoothness)
    is determined by cross validation
    """
    if np.all(data == data[0]):
        fig, ax = plt.subplots()
        ax.hist(
            np.array(data[0]), 
            bins = np.arange(data[0] - 1.5, data[0] + 2.5, 1)
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
    bandwidths = 10 ** np.linspace(-1, 1, 20)
    grid = sklearn.model_selection.GridSearchCV(
        sklearn.neighbors.KernelDensity(kernel='gaussian'),
        {'bandwidth': bandwidths},
        cv=kfold
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


