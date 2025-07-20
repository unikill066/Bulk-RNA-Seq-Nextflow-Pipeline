#!/usr/bin/env python3
"""
This script infers the strandedness of a set of STAR ReadsPerGene.out.tab files.
For more on strandedness, see STAR documentation: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

With --quantMode GeneCounts option STAR will count number reads per gene while mapping.
A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. The counts coincide with those produced by htseq-count with
default parameters. This option requires annotations (GTF or GFF with –sjdbGTFfile option) used at the genome generation step, or at the mapping step. STAR outputs read counts per gene into
ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
column 1: gene ID
column 2: counts for unstranded RNA-seq
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
Select the output according to the strandedness of your data. Note, that if you have stranded
data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the
count of antisense reads. With --quantMode TranscriptomeSAM GeneCounts, and get both the
Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.


StringTie: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
HTSeq: https://htseq.readthedocs.io/en/latest/

STAR columns (from your screenshot and the STAR docs):
col 2 = unstranded
col 3 = “counts for the 1st read strand aligned with RNA” → corresponds to htseq-count -s yes and StringTie’s fr‑secondstrand (--fr)
col 4 = “counts for the 2nd read strand aligned with RNA” → corresponds to htseq-count -s reverse and StringTie’s fr‑firststrand (--rf)

"""

# imports
import os
import sys
import csv
import glob
import argparse

def infer_strandedness(fpath: str, threshold: float = 0.8) -> str:
    """
    Inspect a STAR ReadsPerGene.out.tab file and return:
      - "second" if col 3 ≥ threshold*(col3+col4)
      - "first"  if col 4 ≥ threshold*(col3+col4)
      - "none"   otherwise
    Skips the first four “_N_*” summary lines.
    """
    c3 = c4 = 0
    with open(fpath) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for i, row in enumerate(reader):
            if i < 4:
                continue
            c3 += int(row[2])
            c4 += int(row[3])

    total = c3 + c4
    if total == 0:
        return "none"
    if c3 / total >= threshold:
        return "second"
    if c4 / total >= threshold:
        return "first"
    return "none"

def gather_files(path: str, recursive: bool=False):
    """
    Return a sorted list of all *ReadsPerGene.out.tab files
    in `path`.  If path is a file, return just that file.
    """
    if os.path.isfile(path):
        return [path]
    if os.path.isdir(path):
        pattern = "**/*ReadsPerGene.out.tab" if recursive else "*ReadsPerGene.out.tab"
        found = glob.glob(os.path.join(path, pattern), recursive=recursive)
        if not found:
            sys.exit(f"ERROR: no files matching '*ReadsPerGene.out.tab' in {path}")
        return sorted(found)
    sys.exit(f"ERROR: path not found: {path}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Infer strandedness from STAR ReadsPerGene.out.tab files"
    )
    p.add_argument(
        "path",
        help="either a ReadsPerGene.out.tab file or a directory containing them"
    )
    p.add_argument(
        "-r", "--recursive",
        action="store_true",
        help="look in subdirectories, too"
    )
    p.add_argument(
        "-t", "--threshold",
        type=float, default=0.8,
        help="fractional cutoff for calling stranded (default: 0.8)"
    )
    args = p.parse_args()

    files = gather_files(args.path, recursive=args.recursive)
    print("Sample\tStrandedness")
    for f in files:
        call = infer_strandedness(f, args.threshold)
        # strip off everything after your sample ID
        sample = os.path.basename(f).split("ReadsPerGene.out.tab")[0]
        print(f"{sample}\t{call}")
