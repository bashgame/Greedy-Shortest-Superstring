def seq_from_fasta(filename):
    """ Strips the initial line from a fasta file containing a single sequence
        for a single organism. Reads the remaining lines and concatenates them
        onto an str.
    Args:
        filename (str): The name of the file to read from in fasta format
    Returns:
        str: The genome sequence contained in the file
    """
    sequence = ''
    try:
        file = open(filename)
    except FileNotFoundError:
        print(f"{filename} not found!")
        return "Error, file not found"
    for line in file:
        line = line.rstrip()
        if not line.startswith('>'):
            sequence += line
    return sequence


def reads_from_fastq(filename):
    """ Reads from a fastq file containing reads from a single organism.
    Args:
        filename (str): The name of the file to read from in fastq format
    Returns:
         str  The name of the read set
        [str] The reads from the file
        [str] The quality values encoded with phred33
    """
    name = ''
    reads = []
    quals = []
    try:
        file = open(filename)
    except FileNotFoundError:
        print(f"{filename} not found!")
        return "Error, file not found", None, None
    while True:
        name = file.readline().rstrip()
        read = file.readline().rstrip()
        file.readline()
        qual = file.readline().rstrip()
        if len(read) == 0:
            break
        reads.append(read)
        quals.append(qual)
    return name, reads, quals
