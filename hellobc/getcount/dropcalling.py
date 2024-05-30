import numpy as np
import numpy.matlib as npm
import pandas as pd
import scipy
import scipy.special
import scipy.stats
from typing import NamedTuple
import sys

class countMatrix(object):

    def __init__(self, bcs, features, mtx) -> None:

        self.bcs = bcs
        self.features = features
        self.mtx = mtx

        self.bcs_dim, = self.bcs.shape

        self.counts_per_bc = self.getCountsPerBc()

    def into_csc(self):

        if type(self.mtx) is not scipy.sparse.csc_matrix:
            self.mtx = scipy.sparse.csc_matrix(self.mtx)
    
    def getCountsPerBc(self):

        return np.squeeze(np.asarray(self.mtx.sum(axis=0)))
        
        # column_sum = self.mtx.sum(axis=0)
        # column_sum_array = np.array([item[0] for item in column_sum][0])[0]
        # return column_sum_array

    def getCountsBeforeMapping(self, unmapdict:dict):
        tmp_unmap = list()
        for bc in self.bcs:
            # if bc[:-2] in unmapdict:
            tmp_unmap.append(unmapdict[bc[:-2]])
        self.unmapped_counts = np.array(tmp_unmap)
        return np.array(tmp_unmap)
        
    @staticmethod
    def readLegMatrix(f_mtx:str):

        barcodes = pd.read_csv(f"{f_mtx}/barcodes.tsv", delimiter='\t', header=None, usecols=[0]).values.squeeze()
        features = pd.read_csv(f"{f_mtx}/features.tsv", delimiter='\t', header=None, usecols=[0, 1], names=['gene_id', 'gene_name'])
        matrix_mtx = scipy.io.mmread(f"{f_mtx}/matrix.mtx")

        matrix = countMatrix(bcs=barcodes, features=features, mtx=matrix_mtx)
        matrix.into_csc()

        return matrix


class Metrics:
    def update(self, other):
        for k, v in other.__dict__.items():
            if v is not None:
                setattr(self, k, v)


class BarcodeFilterResults(Metrics):
    def __init__(self, default_value: int = 0):
        self.filtered_bcs = default_value
        self.filtered_bcs_lb = default_value
        self.filtered_bcs_ub = default_value
        self.filtered_bcs_var = default_value
        self.filtered_bcs_cv = float(default_value)

    @staticmethod
    def init_with_constant_call(n_bcs: int):
        res = BarcodeFilterResults()
        res.filtered_bcs = n_bcs
        res.filtered_bcs_lb = n_bcs
        res.filtered_bcs_ub = n_bcs
        res.filtered_bcs_var = 0
        res.filtered_bcs_cv = 0
        return res

    def to_dict_with_prefix(
        self, i: int, sample, method: str
    ) -> dict:
        sample_prefix = "_" + sample if sample else ""
        return {
            "gem_group_%d%s_%s_%s" % (i, sample_prefix, key, method): value
            for key, value in self.__dict__.items()
        }        


def ordmag_algrithm(raw_bcs_nonzero:np.ndarray, upper_expect_limit:int=1<<18, appro_step_num:int=2000, boot_num:int=100):

    rs = np.random.RandomState(0)
    loss = list()
    estimate_cnums = list()

    for _ in range(boot_num):
    
        raw_bcs_nonzero_random = rs.choice(raw_bcs_nonzero, len(raw_bcs_nonzero))
        estimate_cnum = np.linspace(1, np.log2(upper_expect_limit), appro_step_num)    
        estimate_cnum = np.unique(np.round(np.power(2, estimate_cnum)).astype(int))   
        baseline_bc_idx = np.round(estimate_cnum * (1 - 0.99))                       
        baseline_bc_idx = np.minimum(baseline_bc_idx.astype(int), len(raw_bcs_nonzero_random) - 1)
        raw_bcs_nonzero_dec = np.sort(raw_bcs_nonzero_random)[::-1]
        baseline = raw_bcs_nonzero_dec[baseline_bc_idx]
        cutoff = np.maximum(1, np.round(0.1 * baseline)).astype(int)                  

        observed_cnum = len(raw_bcs_nonzero_random) - np.searchsorted(raw_bcs_nonzero_dec[::-1], cutoff)
        los = np.power(observed_cnum - estimate_cnum, 2) / estimate_cnum
        idx = np.argmin(los)

        estimate_cnums.append(estimate_cnum[idx])
        loss.append(los[idx])

    return np.array(estimate_cnums), np.array(loss)


def bootstrap_odmag(raw_bcs_nonzero:np.ndarray, baselind_bc_ind:np.ndarray, boot_num:int=100):
    rs = np.random.RandomState(0)
    top_n_boot = list()
    for _ in range(boot_num):
        raw_bcs_nonzero_random = rs.choice(raw_bcs_nonzero, len(raw_bcs_nonzero))
        raw_bcs_nonzero_dec = np.sort(raw_bcs_nonzero_random)[::-1]
        # Add +1 as we're getting from the other side
        baseline = raw_bcs_nonzero_dec[baselind_bc_ind]
        cutoff = np.maximum(1, np.round(0.1 * baseline)).astype(int)
        # Return the index corresponding to the cutoff in descending order
        top_n_boot.append(len(raw_bcs_nonzero_random) - np.searchsorted(raw_bcs_nonzero_dec[::-1], cutoff))
        top_n_boot = np.array(top_n_boot)
        return top_n_boot


def robust_divide(a, b) -> float:
    """Handles 0 division and conversion to floats automatically."""
    a = float(a)
    b = float(b)
    if b == 0:
        return float("NaN")
    else:
        return a / b
    

def summarize_bootstrapped_top_n(top_n_boot, nonzero_counts):
    top_n_bcs_mean = np.mean(top_n_boot)
    top_n_bcs_var = np.var(top_n_boot)
    top_n_bcs_sd = np.sqrt(top_n_bcs_var)
    result = BarcodeFilterResults()
    result.filtered_bcs_var = top_n_bcs_var
    result.filtered_bcs_cv = robust_divide(top_n_bcs_sd, top_n_bcs_mean)
    result.filtered_bcs_lb = np.round(scipy.stats.norm.ppf(0.025, top_n_bcs_mean, top_n_bcs_sd), 0)
    result.filtered_bcs_ub = np.round(scipy.stats.norm.ppf(0.975, top_n_bcs_mean, top_n_bcs_sd), 0)

    nbcs = int(np.round(top_n_bcs_mean))
    result.filtered_bcs = nbcs

    # make sure that if a barcode with count x is selected, we select all barcodes with count >= x
    # this is true for each bootstrap sample, but is not true when we take the mean

    if nbcs > 0:
        order = np.argsort(nonzero_counts, kind="stable")[::-1]
        sorted_counts = nonzero_counts[order]

        cutoff = sorted_counts[nbcs - 1]
        index = nbcs - 1
        if cutoff > 0:
            while (index + 1) < len(sorted_counts) and sorted_counts[index] == cutoff:
                index += 1
                # if we end up grabbing too many barcodes, revert to initial estimate
                if (index + 1 - nbcs) > 0.20 * nbcs:
                    return result
        result.filtered_bcs = index + 1

    return result


# Add References (emptyDrops)
def emptyDrops_init_calling(raw_mtx:countMatrix):
    
    metrics = BarcodeFilterResults(0)
    
    raw_bcs_nonzero = raw_mtx.counts_per_bc[raw_mtx.counts_per_bc > 0]
    
    if len(raw_bcs_nonzero) == 0:
        raise Warning("All barcodes have 0 counts.")
        return []
    
    estimate_cnums, loss = ordmag_algrithm(raw_bcs_nonzero=raw_bcs_nonzero)
    mean_estimate_cnum = np.mean(estimate_cnums)
    mean_los = np.mean(loss)
    mean_estimate_cnum = max(int(np.round(mean_estimate_cnum)), 50)
    print(f"Found recovered_cells = {mean_estimate_cnum} with loss = {mean_los}")

    baseline_bc_idx = int(np.round(float(mean_estimate_cnum) * (1 - 0.99)))
    baseline_bc_idx = min(baseline_bc_idx, len(raw_bcs_nonzero) - 1)

    # Bootstrap sampling; run algo with many random samples of the data
    top_n_boot = bootstrap_odmag(raw_bcs_nonzero=raw_bcs_nonzero, baselind_bc_ind=baseline_bc_idx, boot_num=100)

    metrics.update(summarize_bootstrapped_top_n(top_n_boot, raw_bcs_nonzero))

    # Get the filtered barcodes
    top_n = metrics.filtered_bcs
    top_bc_idx = np.sort(np.argsort(raw_mtx.counts_per_bc, kind="stable")[::-1][0:top_n])             
    assert top_n <= len(raw_bcs_nonzero), "Invalid selection of 0-count barcodes!"
    
    return top_bc_idx, metrics
    
# ---------------------------------------------------------- init knee method --------------------------------------------------- #
def getKneeDistance(values):
    '''
    This function is based on
    https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

    and https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1

    The idea is to draw a line from the first to last point on the
    cumulative counts curve and then find the point on the curve
    which is the maximum distance away from this line
    '''

    # get coordinates of all the points
    nPoints = len(values)
    allCoord = np.vstack((range(nPoints), values)).T

    # get the first point
    firstPoint = allCoord[0]
    # get vector between first and last point - this is the line
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))

    # find the distance from each point to the line:
    # vector between all points and first point
    vecFromFirst = allCoord - firstPoint

    # To calculate the distance to the line, we split vecFromFirst into two
    # components, one that is parallel to the line and one that is perpendicular
    # Then, we take the norm of the part that is perpendicular to the line and
    # get the distance.
    # We find the vector parallel to the line by projecting vecFromFirst onto
    # the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
    # We project vecFromFirst by taking the scalar product of the vector with
    # the unit vector that points in the direction of the line (this gives us
    # the length of the projection of vecFromFirst onto the line). If we
    # multiply the scalar product by the unit vector, we have vecFromFirstParallel

    scalarProduct = np.sum(
        vecFromFirst * npm.repmat(lineVecNorm, nPoints, 1), axis=1)
    
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel

    # distance to line is the norm of vecToLine
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))

    # knee/elbow is the point with max distance value
    idxOfBestPoint = np.argmax(distToLine)

    return(distToLine, idxOfBestPoint)


def init_knee_method(raw_mtx:countMatrix):

    raw_bcs_dec = np.sort(raw_mtx.counts_per_bc)[::-1]
    raw_bcs_dec_ind = np.argsort(raw_mtx.counts_per_bc)[::-1]
    values = list(np.cumsum(raw_bcs_dec))
    
    previous_idxOfBestPoint = 0
    distToLine, idxOfBestPoint = getKneeDistance(values)
    if idxOfBestPoint == 0:
        raise ValueError("Something's gone wrong here!!")

    max_iterations = 100
    iterations = 0
    while idxOfBestPoint - previous_idxOfBestPoint != 0:
        previous_idxOfBestPoint = idxOfBestPoint
        iterations += 1
        if iterations > max_iterations:
            break
        distToLine, idxOfBestPoint = getKneeDistance(values[:idxOfBestPoint*3])

    init_bcs_inds = np.array(raw_bcs_dec_ind[:idxOfBestPoint+1])
    return init_bcs_inds


# ---------------------------------------------------------- non ambient calling ----------------------------------------------------------- #
class NonAmbientBarcodeResult(NamedTuple):
    eval_bcs: np.ndarray  # Candidate barcode indices (n)
    log_likelihood: np.ndarray  # Ambient log likelihoods (n)
    pvalues: np.ndarray  # pvalues (n)
    pvalues_adj: np.ndarray  # B-H adjusted pvalues (n)
    is_nonambient: np.ndarray  # Boolean nonambient calls (n)


def adjust_pvalue_bh(p):
    """Multiple testing correction of p-values using the Benjamini-Hochberg procedure."""
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    # pylint: disable=no-member
    qual = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return qual[np.argsort(descending)]


class SimpleGoodTuringError(Exception):
    pass


def _averaging_transform(r, nr):
    d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
    dr = np.concatenate((0.5 * (d[1:] + d[0:-1]), np.array((d[-1],), dtype=float)))
    return nr.astype(float) / dr


def _rstest(r, coef):
    return r * np.power(1 + 1.0 / r, 1 + coef)


def simple_good_turing(xr: np.ndarray, xnr: np.ndarray):
    """Make a Simple Good-Turing estimate of the frequencies.

    Args:
      xr (np.array(int)): Non-zero item frequencies
      xnr (np.array(int)): Non-zero frequencies of frequencies
    Returns:
      (rstar (np.array(float)), p0 (float)):
        rstar: The adjusted non-zero frequencies
        p0: The total probability of unobserved items
    """
    xr = xr.astype(float)
    xnr = xnr.astype(float)

    xN = np.sum(xr * xnr)

    # Get Linear Good-Turing estimate
    xnrz = _averaging_transform(xr, xnr)
    slope, _, _, _, _ = scipy.stats.linregress(np.log(xr), np.log(xnrz))

    if slope > -1:
        raise SimpleGoodTuringError(
            "The log-log slope is > -1 (%d); the SGT estimator is not applicable to these data."
            % slope
        )

    xrst = _rstest(xr, slope)
    xrstrel = xrst / xr

    # Get traditional Good-Turing estimate
    xrtry = xr == np.concatenate((xr[1:] - 1, np.zeros(1)))
    xrstarel = np.zeros(len(xr))
    xrstarel[xrtry] = (
        (xr[xrtry] + 1) / xr[xrtry] * np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]
    )

    # Determine when to switch from GT to LGT estimates
    tursd = np.ones(len(xr))
    for i in range(len(xr)):
        if xrtry[i]:
            tursd[i] = float(i + 2) / xnr[i] * np.sqrt(xnr[i + 1] * (1 + xnr[i + 1] / xnr[i]))

    xrstcmbrel = np.zeros(len(xr))
    useturing = True
    for r in range(len(xr)):
        if not useturing:
            xrstcmbrel[r] = xrstrel[r]
        else:
            if np.abs(xrstrel[r] - xrstarel[r]) * (1 + r) / tursd[r] > 1.65:
                xrstcmbrel[r] = xrstarel[r]
            else:
                useturing = False
                xrstcmbrel[r] = xrstrel[r]

    # Renormalize the probabilities for observed objects
    sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)

    xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
    p0 = xnr[0] / xN

    return (xr * xrstcmbrel, p0)


def sgt_proportions(frequencies: np.ndarray):
    """Use Simple Good-Turing estimate to adjust for unobserved items.

    Args:
      frequencies (np.array(int)): Nonzero frequencies of items

    Returns:
        pstar (np.array[float]): The adjusted non-zero proportions
        p0 (float): The total probability of unobserved items
    """
    if len(frequencies) == 0:
        raise ValueError("Input frequency vector is empty")
    if np.count_nonzero(frequencies) != len(frequencies):
        raise ValueError("Frequencies must be greater than zero")

    freqfreqs = np.bincount(frequencies)
    assert freqfreqs[0] == 0
    use_freqs = np.flatnonzero(freqfreqs)

    if len(use_freqs) < 10:
        raise SimpleGoodTuringError(
            "Too few non-zero frequency items (%d). Aborting SGT." % len(use_freqs)
        )

    rstar, p0 = simple_good_turing(use_freqs, freqfreqs[use_freqs])

    # rstar contains the smoothed frequencies.
    # Map each original frequency r to its smoothed rstar.
    rstar_dict = dict(zip(use_freqs, rstar))

    rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
    rstar_i = np.fromiter((rstar_dict[f] for f in frequencies), dtype=float, count=len(frequencies))
    pstar = (1 - p0) * (rstar_i / rstar_sum)

    assert np.isclose(p0 + np.sum(pstar), 1)
    return (pstar, p0)


def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """Estimate a gene expression profile by Simple Good Turing.

    Args:
      raw_mat (sparse matrix): Sparse matrix of all counts
      barcode_indices (np.array(int)): Barcode indices to use
      nz_feat (np.array(int)): Indices of features that are non-zero at least once

    Returns:
      profile (np.array(float)): Estimated probabilities of length len(nz_feat).
    """
    # Initial profile estimate
    prof_mat = matrix[:, barcode_indices]

    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    # Simple Good Turing estimate
    p_smoothed, p0 = sgt_proportions(profile[np.flatnonzero(profile)])

    n0 = len(zero_feat)
    if n0 == 0:
        # Renormalize in absence of 0 class
        p_smoothed = p_smoothed / p_smoothed.sum()
        p0_i = -1.0  # unused
    else:
        # Distribute p0 equally among the zero elements.
        p0_i = p0 / n0

    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    assert np.isclose(profile_p.sum(), 1.0)
    return profile_p

# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
    """Estimate a gene expression profile on a given subset of barcodes.

    Use Good-Turing to smooth the estimated profile.

    Args:
      matrix (scipy.sparse.csc_matrix): Sparse matrix of all counts
      use_bcs (np.array(int)): Indices of barcodes to use (col indices into matrix)

    Returns:
      profile (use_features, np.array(float)): Estimated probabilities of length use_features.
    """
    # Use features that are nonzero anywhere in the data
    use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))

    # Estimate background profile
    bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)

    return (use_feats, bg_profile_p)


MIN_UMI_FRAC_OF_MEDIAN = 0.01
MIN_UMIS = 500
MAX_ADJ_PVALUE = 0.01                                   # p-value threshold
NUM_SIMS = 10000


def get_empty_drops_range(num_probe_bcs: int = None):
    """Gets the range of values to use for empty drops background given a chemistry description.

    Args:
        chemistry_description: A string describing the chemistry

    Returns:
        low_index:
        high_index:
    """
    
    N_PARTITIONS = 45000 * num_probe_bcs if num_probe_bcs and num_probe_bcs > 1 else 90000
    return (N_PARTITIONS // 2, N_PARTITIONS)


def eval_multinomial_loglikelihoods(matrix, profile_p, max_mem_gb: float = 0.1):
    """Compute the multinomial log PMF for many barcodes.

    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      profile_p (np.ndarray(float)): Multinomial probability vector
      max_mem_gb (float): Try to bound memory usage.

    Returns:
      log_likelihoods (np.ndarray(float)): Log-likelihood for each barcode
    """
    gb_per_bc = float(matrix.shape[0] * matrix.dtype.itemsize) / (1024**3)
    bcs_per_chunk = max(1, int(np.round(max_mem_gb / gb_per_bc)))
    num_bcs = matrix.shape[1]

    loglk = np.zeros(num_bcs)

    for chunk_start in range(0, num_bcs, bcs_per_chunk):
        chunk = slice(chunk_start, chunk_start + bcs_per_chunk)
        matrix_chunk = matrix[:, chunk].transpose().toarray()
        n = matrix_chunk.sum(1)
        loglk[chunk] = scipy.stats.multinomial.logpmf(matrix_chunk, n, p=profile_p)
    return loglk


def _get_sample(
    profile_p: np.array,
    sample_size: int,
    *,
    existing_sample,
    current_index,
):
    """Gets single samples from a multinomial distribution.

    Code can take an existing array that has len(existing_sample) - current_index "unused" samples in it
    and append that to the start of a new sample.  Done it this way to make sure the sample order of the code
    in this file is always the same, even as we add more vectors.

    Args:
        profile_p:
        sample_size:
        existing_sample:
        current_index:

    Returns:
        An array of sampled indices

    """
    if existing_sample is not None:
        assert current_index is not None
        old_sample = existing_sample[current_index : len(existing_sample)]
        sample_size = sample_size - len(old_sample)
        new_sample = np.random.choice(len(profile_p), size=sample_size, p=profile_p, replace=True)
        return np.concatenate((old_sample, new_sample))
    return np.random.choice(len(profile_p), size=sample_size, p=profile_p, replace=True)


def simulate_multinomial_loglikelihoods(
    profile_p,
    umis_per_bc,
    num_sims: int = 1000,
    jump: int = 1000,
    n_sample_feature_block: int = 1000000,
    verbose: bool = False,
):
    """Simulate draws from a multinomial distribution for various values of N.

    Note, the samples within a given simulation are not independent.  Given two barcodes with UMIs equal to 10 and
    another equal to 20, we simulate the 10 UMI barcode, than draw another 10 samples to get the likelihood for the
    20 UMI barcode.

    Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )

    Args:
      profile_p (np.ndarray(float)): Probability of observing each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode (multinomial N).
      num_sims (int): Number of simulations per distinct N value.
      jump (int): Vectorize the sampling if the gap between two distinct Ns exceeds this.
      n_sample_feature_block (int): Vectorize this many feature samplings at a time.

    Returns:
      distinct_ns (np.ndarray(int)): an array containing the distinct N values
          that were simulated.
      log_likelihoods (np.ndarray(float)): a len(distinct_ns) x num_sims matrix
          containing the simulated log likelihoods.
    """
    np.random.seed(0)

    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))

    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)
    num_all_n = np.max(distinct_n) - np.min(distinct_n)
    if verbose:
        print("Number of distinct N supplied: %d" % len(distinct_n))
        print("Range of N: %d" % num_all_n)
        print("Number of features: %d" % len(profile_p))

    sampled_features = _get_sample(
        profile_p, n_sample_feature_block, existing_sample=None, current_index=None
    )
    k = 0

    log_profile_p = np.log(profile_p)

    for sim_idx in range(num_sims):
        if verbose and sim_idx % 100 == 99:
            sys.stdout.write(".")
            sys.stdout.flush()
        curr_counts = np.ravel(scipy.stats.multinomial.rvs(distinct_n[0], profile_p, size=1))

        curr_loglk = scipy.stats.multinomial.logpmf(curr_counts, distinct_n[0], p=profile_p)

        loglk[0, sim_idx] = curr_loglk

        # There are three strategies here for incrementing the likelihood to go from a Barcode that had X UMIs observed
        # to one with X+Y UMIs observed.  We either get a whole new Y sample, add it to the old one, and recompute the
        # full likelihood, we add a single new observation Y times, or we add a group of Y new observations together at
        # once.  This design comes from trying to maintain compatibility with the original version of the code while
        # improving performance, would not be the way to do this from scratch.
        for i in range(1, len(distinct_n)):
            step = distinct_n[i] - distinct_n[i - 1]
            if step >= jump:
                # Instead of iterating for each n, sample the intermediate ns all at once
                curr_counts += np.ravel(scipy.stats.multinomial.rvs(step, profile_p, size=1))
                curr_loglk = scipy.stats.multinomial.logpmf(curr_counts, distinct_n[i], p=profile_p)
                assert not np.isnan(curr_loglk)
            elif (
                step < 20
            ):  # With small sizes calling np.unique and doing a vector version not worth the overhead
                # Profiling test_filter_barcodes showed 20 was faster than 10 but slower than 100,
                # further optimization of the cutoff might be possible.
                # In this case we'll iteratively sample between the two distinct values of n
                for n in range(distinct_n[i - 1] + 1, distinct_n[i] + 1):
                    j = sampled_features[k]
                    k += 1
                    if k >= n_sample_feature_block:
                        # Amortize this operation
                        sampled_features = _get_sample(
                            profile_p,
                            n_sample_feature_block,
                            existing_sample=None,
                            current_index=None,
                        )
                        k = 0
                    curr_counts[j] += 1
                    curr_loglk += log_profile_p[j] + np.log(float(n) / curr_counts[j])
            else:
                # Vectorized version of an incremental new sampling, grab all new samples at once
                end = k + step
                if end >= n_sample_feature_block:
                    if step >= n_sample_feature_block:
                        n_sample_feature_block = step + 1
                    sampled_features = _get_sample(
                        profile_p,
                        n_sample_feature_block,
                        existing_sample=sampled_features,
                        current_index=k,
                    )
                    k = 0
                    end = step
                feats, cnts = np.unique(sampled_features[k:end], return_counts=True)
                k = end
                old_cnts = curr_counts[feats]
                curr_counts[feats] += cnts
                # Log(p_i^k)
                curr_loglk += np.sum(np.multiply(log_profile_p[feats], cnts))
                # calc N_new! - (N_new-N_old)! to incrementally update likelihood
                factorial_increase = scipy.special.loggamma((distinct_n[i] + 1, distinct_n[i - 1] + 1))
                curr_loglk += factorial_increase[0] - factorial_increase[1]
                # calc X_new! - (X_new - X_old)! and subtract as that's in the denominator
                curr_loglk -= np.sum(scipy.special.loggamma(old_cnts + cnts + 1) - scipy.special.loggamma(old_cnts + 1))
                del old_cnts

            loglk[i, sim_idx] = curr_loglk

    if verbose:
        sys.stdout.write("\n")

    return distinct_n, loglk


def compute_ambient_pvalues(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """Compute p-values for observed multinomial log-likelihoods.

    Args:
      umis_per_bc (nd.array(int)): UMI counts per barcode
      obs_loglk (nd.array(float)): Observed log-likelihoods of each barcode deriving from an ambient profile
      sim_n (nd.array(int)): Multinomial N for simulated log-likelihoods
      sim_loglk (nd.array(float)): Simulated log-likelihoods of shape (len(sim_n), num_simulations)

    Returns:
      pvalues (nd.array(float)): p-values
    """
    assert len(umis_per_bc) == len(obs_loglk)
    assert sim_loglk.shape[0] == len(sim_n)

    # Find the index of the simulated N for each barcode
    sim_n_idx = np.searchsorted(sim_n, umis_per_bc)
    num_sims = sim_loglk.shape[1]

    num_barcodes = len(umis_per_bc)

    pvalues = np.zeros(num_barcodes)

    for i in range(num_barcodes):
        num_lower_loglk = np.sum(sim_loglk[sim_n_idx[i], :] < obs_loglk[i])
        pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
    return pvalues


def find_nonambient_barcodes(
    matrix:countMatrix,
    orig_cell_bcs,
    chemistry_description=None,
    num_probe_bcs=None,
    *,
    min_umi_frac_of_median=MIN_UMI_FRAC_OF_MEDIAN,
    emptydrops_minimum_umis=MIN_UMIS,
    max_adj_pvalue=MAX_ADJ_PVALUE,
    num_sims=NUM_SIMS,
):
    """Call barcodes as being sufficiently distinct from the ambient profile.

    Args:
      matrix (CountMatrix): Full expression matrix.
      orig_cell_bcs (iterable of str): Strings of initially-called cell barcodes.
      chemistry_description: Change ambient RNA estimation for LT chemistry

    Returns:
      NonAmbientBarcodeResult
    """
    # Estimate an ambient RNA profile
    umis_per_bc = matrix.getCountsPerBc()
    bc_order = np.argsort(umis_per_bc)

    low, high = get_empty_drops_range(num_probe_bcs)

    # Take what we expect to be the barcodes associated w/ empty partitions.
    print(f"Range empty barcodes: {low} - {high}")
    empty_bcs = bc_order[::-1][low:high]
    empty_bcs.sort()

    # Require non-zero barcodes
    nz_bcs = np.flatnonzero(umis_per_bc)
    nz_bcs.sort()

    use_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)

    if len(use_bcs) > 0:
        try:
            eval_features, ambient_profile_p = est_background_profile_sgt(matrix.mtx, use_bcs)
        except SimpleGoodTuringError as e:
            print(str(e))
            return None
    else:
        eval_features = np.zeros(0, dtype=np.int64)
        ambient_profile_p = np.zeros(0)

    # Choose candidate cell barcodes
    orig_cell_bc_set = set(orig_cell_bcs)
    print(f"len(matrix.bcs): {len(matrix.bcs)}")
    np_fromiternp = np.fromiter(
            (bc in orig_cell_bc_set for bc in matrix.bcs), count=len(matrix.bcs), dtype=bool
        )
    np_fromiternp_set = set(np_fromiternp)
    
    orig_cells = np.flatnonzero(
        np.fromiter(
            (bc in orig_cell_bc_set for bc in matrix.bcs), count=len(matrix.bcs), dtype=bool
        )
    )

    # No good incoming cell calls
    if orig_cells.sum() == 0:
        return None

    # Look at non-cell barcodes above a minimum UMI count
    eval_bcs = np.ma.array(np.arange(matrix.bcs_dim))
    eval_bcs[orig_cells] = np.ma.masked

    median_initial_umis = np.median(umis_per_bc[orig_cells])
    min_umis = int(
        max(emptydrops_minimum_umis, np.ceil(median_initial_umis * min_umi_frac_of_median))
    )
    print(f"Median UMIs of initial cell calls: {median_initial_umis}")
    print(f"Min UMIs: {min_umis}")

    eval_bcs[umis_per_bc < min_umis] = np.ma.masked
    n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()

    eval_bcs = np.argsort(np.ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[0:n_unmasked_bcs]
    # SORT the barcodes. This is a critical step; eval_bcs is a list of integers, and these indices are used
    # to get the counts via matrix.select_features_by_genome(genome).select_barcodes(eval_bcs).get_counts_per_bc(),
    # which sorts the input, and also to get the string barcode sequences via np.array(matrix.bcs)[eval_bcs],
    # which doesn't. Without sorting here, the matching between the sequences and their counts from both
    # genomes is essentially random, thus the species assignments will be wrong.
    eval_bcs.sort()

    if len(eval_bcs) == 0:
        return None

    assert not np.any(np.isin(eval_bcs, orig_cells))
    print(f"Number of candidate bcs: {len(eval_bcs)}")
    print(f"Range candidate bc umis: {umis_per_bc[eval_bcs].min()}, {umis_per_bc[eval_bcs].max()}")

    eval_mat = matrix.mtx[eval_features, :][:, eval_bcs]

    if len(ambient_profile_p) == 0:
        obs_loglk = np.repeat(np.nan, len(eval_bcs))
        pvalues = np.repeat(1, len(eval_bcs))
        sim_loglk = np.repeat(np.nan, len(eval_bcs))
        return None

    # Compute observed log-likelihood of barcodes being generated from ambient RNA
    obs_loglk = eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

    # Simulate log likelihoods
    distinct_ns, sim_loglk = simulate_multinomial_loglikelihoods(
        ambient_profile_p, umis_per_bc[eval_bcs], num_sims=num_sims, verbose=True
    )

    # Compute p-values
    pvalues = compute_ambient_pvalues(
        umis_per_bc[eval_bcs], obs_loglk, distinct_ns, sim_loglk
    )

    pvalues_adj = adjust_pvalue_bh(pvalues)

    is_nonambient = pvalues_adj <= max_adj_pvalue

    return NonAmbientBarcodeResult(
        eval_bcs=eval_bcs,
        log_likelihood=obs_loglk,
        pvalues=pvalues,
        pvalues_adj=pvalues_adj,
        is_nonambient=is_nonambient,
    )


def emptyDrops_drop_calling(raw_mtx:countMatrix, init_method:str='ordmag', min_umis=MIN_UMIS, max_adj_pvalue=MAX_ADJ_PVALUE):

    if init_method == 'ordmag':
        top_bc_idx, metrics = emptyDrops_init_calling(raw_mtx=raw_mtx)
    elif init_method == 'knee':
        top_bc_idx = init_knee_method(raw_mtx)
    else:
        raise ValueError("init_method must be one of 'ordmag', 'knee'")
    nonambi_calling_res = find_nonambient_barcodes(matrix=raw_mtx, orig_cell_bcs=raw_mtx.bcs[top_bc_idx], emptydrops_minimum_umis=min_umis, max_adj_pvalue=max_adj_pvalue)

    nonambi_bcs = [bc for bc, non_ambi in zip(nonambi_calling_res.eval_bcs, nonambi_calling_res.is_nonambient) if non_ambi]

    return top_bc_idx, nonambi_bcs

