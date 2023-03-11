from modules.overlaps import getLargestOverlap


def greedy_scs(reads, k_len=3):
    while True:
        read1, read2, olap = getLargestOverlap(reads, k_len)
        if olap == 0:
            break
        reads.remove(read1)
        reads.remove(read2)
        reads.append(read1 + read2[olap:])

    return ''.join(reads)
