#!/usr/bin/env python3
import pysam
import argparse

def filter_bam(input_bam, ambiguous_bam, te_txp_only_bam, purpl_known_bam):
    # Open the input BAM and prepare output BAMs with the same header.
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    ambiguous_out = pysam.AlignmentFile(ambiguous_bam, "wb", template=in_bam)
    te_txp_out = pysam.AlignmentFile(te_txp_only_bam, "wb", template=in_bam)
    purpl_known_out = pysam.AlignmentFile(purpl_known_bam, "wb", template=in_bam)

    # Process each read in the BAM.
    for read in in_bam.fetch(until_eof=True):
        try:
            xt_tag = read.get_tag("XT")
        except KeyError:
            # Skip reads without an XT tag.
            continue

        tokens = xt_tag.split(',')
        # Determine if there is at least one token from each category.
        has_enst = any(token.startswith("ENST") for token in tokens)
        has_tu = any(token.startswith("TU") for token in tokens)

        if has_enst and has_tu:
            ambiguous_out.write(read)
        elif has_tu and not has_enst:
            te_txp_out.write(read)
        elif has_enst and not has_tu:
            purpl_known_out.write(read)
        # If tokens contain neither, the read is skipped.

    in_bam.close()
    ambiguous_out.close()
    te_txp_out.close()
    purpl_known_out.close()

def main():
    parser = argparse.ArgumentParser(
        description="Split BAM reads based on the XT tag into three categories:\n"
                    "1) ambiguous: both ENST and TU tokens present,\n"
                    "2) te_txp_only: exclusively TU tokens,\n"
                    "3) purpl_known: exclusively ENST tokens."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-a", "--ambiguous", required=True, help="Output BAM file for ambiguous reads")
    parser.add_argument("-t", "--te_txp_only", required=True, help="Output BAM file for TE/transcript-only reads")
    parser.add_argument("-p", "--purpl_known", required=True, help="Output BAM file for purpl known (ENST only) reads")
    args = parser.parse_args()

    filter_bam(args.input, args.ambiguous, args.te_txp_only, args.purpl_known)

if __name__ == "__main__":
    main()
