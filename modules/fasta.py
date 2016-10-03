def write_to_fasta(handle, chrom, seq_str, line_len=60):
    """
    Write sequence string as fasta to file handle.
    """
    handle.write(">{}\n".format(chrom))
    for idx in range(line_len,len(seq_str),line_len):
        handle.write(seq_str[idx-line_len:idx]+"\n")
    if seq_str[idx:]:
        handle.write(seq_str[idx:]+"\n")
