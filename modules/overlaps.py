def overlap(str1, str2, min_len=3):
    """ Locates the longest suffix of str1 matching a prefix of str2, such
        that the overlap is at least min_len characters long. If no overlap
        is found, return 0
    Args:
        str1 (str): The 'suffix' string
        str2 (str): The 'prefix' string
        min_len (int): The minimum length of an overlap, default value = 3

    Returns:
        int : The length of the overlap, or 0 if no overlap found
    """
    # start at the left
    pos = 0
    while True:
        # find str2's prefix in str1
        pos = str1.find(str2[:min_len], pos)
        # If prefix isn't in str1 to the right of pos:
        if pos == -1:
            return 0
        # found prefix; check for full suffix/prefix match
        if str2.startswith(str1[pos:]):
            return len(str1)-pos
        # move past previous match
        pos += 1


def build_dictionary(reads, k_len=3):
    """ Builds a dictionary that associates each k-mer of length 'min_len' with


    Args:
        reads (list( str )): A list of reads from a genome
        k_len (int, optional): The length of k-mers to generate. Defaults to 3.

    Returns:
        dictionary: key: k-mer of length k_len
                    value: a set consisting of the reads that contain k-mer
    """
    dict = {}
    # Get each read from the list of reads
    for read in reads:
        # Get each k-mer from the read, do not run over the end of read
        for offset in range(len(read) - k_len + 1):
            kmer = str(read[offset: offset + k_len])
            # If we haven't seen this kmer, initialize its set
            if kmer not in dict:
                dict[kmer] = set()
            # Add the read to the set, if it isn't already a member
            dict[kmer] |= {read}
    return dict


def overlap_all_pairs(reads, min_len=3):
    """ Generate a set of all pairs of reads in a list of reads that have an
        overlap of at least min_len

    Args:
        reads (list( str )): A list of reads from a genome
        min_len (int, optional): The minimum length of overlap. Defaults to 3.

    Returns:
        set( tuple ): A set of tuples. Each tuple consists of a pair of reads
                        that overlap by at least min_len.
    """
    pairs = set()
    compared = set()

    # build a dictionary of k-mers and the reads that contain them
    dict = build_dictionary(reads, min_len)

    # Get each read from the list
    for read1 in reads:

        # Get each read that contains the k-length suffix of read1
        for read2 in list(dict[read1[-min_len::]]):

            # Don't call overlap if read1 doesn't contain read2's prefix
            if read1 not in dict[read2[:min_len:]]:
                continue

            # Don't call overlap if read1 == read2 or we've already compared
            elif read1 == read2 or (read1, read2) in compared:
                continue

            # Call overlap and add the reads to compared
            else:
                olap = overlap(read1, read2, min_len)
                compared.add((read1, read2))

                # Add the reads to pairs if there is an overlap
                if olap > 0:
                    pairs.add((read1, read2))
    return pairs
