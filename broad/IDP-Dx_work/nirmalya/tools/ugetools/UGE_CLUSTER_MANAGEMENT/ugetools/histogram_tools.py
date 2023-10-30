def make_histogram(lst, bin_size):
    """
    bin the input values and return a count in each bin
    :param lst: list of numeric values
    :param bin_size: width of each bin
    :return: list of tuple(bin_lower_bound, count, count_of_this_bin_plus_larger_bins)
    """
    dct_ret = {}
    for v in lst:
        bin_lower_bound = int(v/bin_size) * bin_size
        dct_ret[bin_lower_bound] = dct_ret.get(bin_lower_bound, 0) + 1
    lst_bins = dct_ret.items()
    lst_bins.sort()
    # sort in reverse numeric order so count of this bin plus larger bins can be computed
    lst_bins.reverse()
    accumulator = 0
    lst_ret = []
    for bin_lower_bound, count in lst_bins:
        lst_ret.append((bin_lower_bound, count, count + accumulator))
        accumulator += count

    lst_ret.reverse()
    return lst_ret