#!/usr/bin/env python3
import argparse, gzip, csv, sys
from itertools import islice

def parse_args():
    p = argparse.ArgumentParser(
        description="Demultiplex merged FASTQs by barcode into per-well cell files."
    )
    p.add_argument("--annot", required=True,
                   help="annotated TSV with header: Well,Barcode,R1,R2,CellType")
    p.add_argument("--read1", required=True, help="merged R1.fastq.gz")
    p.add_argument("--read2", required=True, help="merged R2.fastq.gz")
    p.add_argument("--cell-type", nargs="+",
                   help="which cell type(s) to extract (default=all)")
    p.add_argument("--unassigned", action="store_true",
                   help="also write reads with no matching barcode")
    return p.parse_args()


def main():
    args = parse_args()

    # Load annotation and filter by cell type
    barcode_to_well = {}
    well_to_cell = {}
    with open(args.annot, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        hdr = next(reader)
        if hdr[0].lower() != 'well':
            f.seek(0)
            reader = csv.reader(f, delimiter="\t")
        for well, barcode, r1, r2, cell, r1_md5, r2_md5 in reader:
            if args.cell_type and cell not in args.cell_type:
                continue
            if barcode in barcode_to_well:
                sys.exit(f"Error: duplicate barcode {barcode}")
            barcode_to_well[barcode] = well
            well_to_cell[well] = cell

    if not barcode_to_well:
        sys.exit("No barcodes found for given --cell-type filter")

    # Prepare output filehandles & paths
    out_handles = {}
    for well, cell in well_to_cell.items():
        r1_name = f"{well}_{cell}_R1.fastq.gz"
        r2_name = f"{well}_{cell}_R2.fastq.gz"
        out_handles[well] = (
            gzip.open(r1_name, "wt"),
            gzip.open(r2_name, "wt")
        )

    # Optional unassigned
    ua1 = ua2 = None
    if args.unassigned:
        ua1 = gzip.open("unassigned_R1.fastq.gz", "wt")
        ua2 = gzip.open("unassigned_R2.fastq.gz", "wt")

    # Determine barcode length
    bc_len = len(next(iter(barcode_to_well)))

    # Demultiplex reads
    with gzip.open(args.read1, "rt") as inf1, gzip.open(args.read2, "rt") as inf2:
        while True:
            r1 = list(islice(inf1, 4))
            r2 = list(islice(inf2, 4))
            if not r1 or not r2:
                break
            key = r1[1].rstrip()[:bc_len]
            if key in barcode_to_well:
                w = barcode_to_well[key]
                o1, o2 = out_handles[w]
                o1.write("".join(r1))
                o2.write("".join(r2))
            elif args.unassigned:
                ua1.write("".join(r1))
                ua2.write("".join(r2))

    # Close handles
    for o1, o2 in out_handles.values():
        o1.close(); o2.close()
    if ua1:
        ua1.close(); ua2.close()

    print("✓ Demultiplexing complete.", file=sys.stderr)

if __name__ == "__main__":
    main()
