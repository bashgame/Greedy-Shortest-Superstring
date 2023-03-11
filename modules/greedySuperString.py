from modules.overlaps import overlap_all_pairs
from random import randrange, seed


def greedySuperString(reads, max_len=-1):
    # If max_len wasn't provided, assume all reads are the same length
    if max_len == -1:
        max_len = len(reads[0])
    # Loop until we break
    length = max_len
    sup_string = ''
    while True:
        longest_pairs = []
        longest = 0
        pairs = overlap_all_pairs(reads, length)
        for read1, read2, olap in pairs:
            if olap >= longest:
                longest = olap
                longest_pairs += [(read1, read2)]
        if longest == 0 and length <= 1:
            break
        elif longest == 0 and length > 1:
            length -= 1
            continue
        else:
            seed()
            read1, read2 = longest_pairs[randrange(len(longest_pairs))]
            reads.remove(read1)
            reads.remove(read2)
            reads.append(read2 + read1[longest:])
    return sup_string.join(reads)
