import numpy as np
from fast_forward.interaction_distribution import BINS_DICT

def hellinger(p,q):
    return np.round(np.sqrt(np.sum(np.power((np.sqrt(p) - np.sqrt(q)),2))) / np.sqrt(2),2)

def std_dev_hist(hist, bins):
    bin_centers = (bins[:-1] + bins[1:]) / 2
    return np.sqrt(np.cov(bin_centers, aweights=hist, bias=True))

def calc_score(ref, test, bins=BINS_DICT['distances']):
    '''
    Compute the score between two distributions.
    The score is a weighted sum of the Hellinger distance and the difference in means of the distributions
    normalized by the standard deviation of the reference distribution. If mean is outside of the standard deviation
    of the reference distribution, it is capped at 1.

    ----------
    ref : numpy 1D-array
        Reference distribution.
    test : numpy 1D-array
        Test distribution.  
    bins : numpy 1D-array
        Bins for the distributions, used to calculate the mean difference.

    Returns
    -------
    float
        Score between the two distributions, between [0, 1].
    '''
    ref = ref / np.sum(ref) if np.sum(ref) > 0 else ref # normalize distributions
    test = test / np.sum(test) if np.sum(test) > 0 else test

    bin_centers = (bins[:-1] + bins[1:]) / 2

    mean_diff = np.average(bin_centers, weights=ref) - np.average(bin_centers, weights=test)
    mean_diff_norm = np.min([np.abs(mean_diff), 1]) # normalize mean difference by standard deviation of reference distribution, if it is zero, use 1 as maximum penalty

    score = hellinger(ref, test) * 0.7 + mean_diff_norm * 0.3 # score is a weighted sum of Hellinger distance and mean difference normalized by standard deviation
    return np.round(score, 2)