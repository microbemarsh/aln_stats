#!/usr/bin/env python3
import os
import stat
import sys
import re
import csv
from optparse import OptionParser
import numpy
from collections import Counter
import pysam

# Attempt to set the script file as executable for the owner.
try:
    script_path = os.path.realpath(__file__)
    st = os.stat(script_path)
    os.chmod(script_path, st.st_mode | stat.S_IXUSR)
except Exception:
    pass

"""
Script to calculate alignment statistics from a PAF, BAM, or SAM file.
It reports per-read stats (such as read/align identity, coverage, etc.) and
aggregates per-target depth (number of reads mapped to each target).
"""

class PafAlignmentStats(object):
    """
    Helper class to store and compute per-read alignment statistics.
    Works for alignments parsed from PAF or SAM/BAM files.
    """
    def __init__(self, queryName, queryLen, queryStart, queryEnd, strand,
                 targetName, targetLen, targetStart, targetEnd,
                 numMatches, alnBlockLen, readSeqLen, refSeqLen,
                 csTag=None, NMTag=None, globalAlignment=True):
        self.queryName = queryName
        self.queryLen = queryLen
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.strand = strand
        self.targetName = targetName
        self.targetLen = targetLen
        self.targetStart = targetStart
        self.targetEnd = targetEnd
        self.numMatches = numMatches
        self.alnBlockLen = alnBlockLen
        self.readSeqLen = readSeqLen  # from FASTQ or alignment field
        self.refSeqLen = refSeqLen    # from FASTA or alignment field

        self.csTag = csTag
        self.NMTag = NMTag
        self.globalAlignment = globalAlignment

        # Precompute differences if we have a cs tag or NM tag.
        self._parsed = False
        self._numMismatches = 0
        self._numInsertions = 0
        self._numInsertionBases = 0
        self._numDeletions = 0
        self._numDeletionBases = 0

        self._parseDifferences()

    def _parseDifferences(self):
        """
        Parse cs:Z: or NM:i: tags to compute mismatch/insertion/deletion counts.
        """
        if self.csTag is not None:
            # Regex for tokens in cs tag (supports both cs=short and cs=long)
            pattern = re.compile(r'(=[A-Za-z]+|:\d+|\*[A-Za-z]+:[A-Za-z]+|\+[A-Za-z]+|\-[A-Za-z]+)')
            tokens = pattern.findall(self.csTag)
            for token in tokens:
                if token.startswith('='):
                    pass  # matched block already counted in numMatches
                elif token.startswith(':'):
                    pass  # consecutive matches; nothing extra to do
                elif token.startswith('*'):
                    # mismatch token, e.g. "*ac:gt"
                    mismatch_sub = token[1:].split(':')
                    mismatch_length = max(len(mismatch_sub[0]), len(mismatch_sub[1]))
                    self._numMismatches += mismatch_length
                elif token.startswith('+'):
                    inserted_seq = token[1:]
                    self._numInsertions += 1
                    self._numInsertionBases += len(inserted_seq)
                elif token.startswith('-'):
                    deleted_seq = token[1:]
                    self._numDeletions += 1
                    self._numDeletionBases += len(deleted_seq)
        elif self.NMTag is not None:
            # Use NM tag as sum of mismatches, insertions, and deletions.
            nm_val = int(self.NMTag)
            self._numMismatches = nm_val
        self._parsed = True

    def readLength(self):
        """
        Return the read length used in coverage stats.
        For local alignments, only the aligned portion is counted.
        """
        if not self.globalAlignment:
            return max(0, self.queryEnd - self.queryStart)
        return self.queryLen

    def readCoverage(self):
        """
        Coverage = (aligned_length / read_length).
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        aligned_len = max(0, self.queryEnd - self.queryStart)
        return aligned_len / r_len

    def alignmentIdentity(self):
        """
        Alignment identity = numMatches / alignment_length.
        """
        if self.alnBlockLen == 0:
            return 0.0
        return float(self.numMatches) / float(self.alnBlockLen)

    def readIdentity(self):
        """
        Read identity = numMatches / read_length.
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self.numMatches) / r_len

    def mismatchesPerAlignedBase(self):
        """
        mismatches / alignment_length.
        """
        if self.alnBlockLen == 0:
            return 0.0
        return float(self._numMismatches) / float(self.alnBlockLen)

    def insertionsPerReadBase(self):
        """
        insertion bases / read_length.
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self._numInsertionBases) / r_len

    def deletionsPerReadBase(self):
        """
        deletion bases / read_length.
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self._numDeletionBases) / r_len

    def referenceCoverage(self):
        """
        Fraction of the target covered by this alignment.
        """
        if self.targetLen == 0:
            return 0.0
        aligned_len = max(0, self.targetEnd - self.targetStart)
        return aligned_len / float(self.targetLen)


def parseFastaLengths(fastaPath):
    """
    Parse the reference FASTA file to obtain sequence lengths by name.
    """
    lengths = {}
    with open(fastaPath, 'r') as f:
        seqName = None
        seqLen = 0
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seqName is not None:
                    lengths[seqName] = seqLen
                seqName = line[1:].split()[0]
                seqLen = 0
            else:
                seqLen += len(line)
        if seqName is not None:
            lengths[seqName] = seqLen
    return lengths

def parseFastqLengths(fastqPath):
    """
    Parse a FASTQ file to get read lengths by name.
    """
    lengths = {}
    with open(fastqPath, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            readName = header.strip().split()[0][1:]  # remove '@'
            lengths[readName] = len(seq.strip())
    return lengths

def parsePaf(pafPath, readLengths=None, refLengths=None, globalAlignment=True):
    """
    Parse a PAF file and return a list of PafAlignmentStats objects.
    """
    alignments = []
    with open(pafPath, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            (queryName, queryLen, queryStart, queryEnd, strand,
             targetName, targetLen, targetStart, targetEnd,
             numMatches, alnBlockLen, mapq) = parts[:12]

            queryLen = int(queryLen)
            queryStart = int(queryStart)
            queryEnd = int(queryEnd)
            targetLen = int(targetLen)
            targetStart = int(targetStart)
            targetEnd = int(targetEnd)
            numMatches = int(numMatches)
            alnBlockLen = int(alnBlockLen)

            # Use provided FASTQ/FASTA lengths if available.
            readLen = readLengths.get(queryName, queryLen) if readLengths else queryLen
            refLen = refLengths.get(targetName, targetLen) if refLengths else targetLen

            csTag = None
            NMTag = None
            for opt in parts[12:]:
                if opt.startswith("cs:Z:"):
                    csTag = opt[5:]
                elif opt.startswith("NM:i:"):
                    NMTag = opt[5:]
            pafStat = PafAlignmentStats(
                queryName=queryName,
                queryLen=queryLen,
                queryStart=queryStart,
                queryEnd=queryEnd,
                strand=strand,
                targetName=targetName,
                targetLen=targetLen,
                targetStart=targetStart,
                targetEnd=targetEnd,
                numMatches=numMatches,
                alnBlockLen=alnBlockLen,
                readSeqLen=readLen,
                refSeqLen=refLen,
                csTag=csTag,
                NMTag=NMTag,
                globalAlignment=globalAlignment
            )
            alignments.append(pafStat)
    return alignments

def parseBam(bamPath, readLengths=None, refLengths=None, globalAlignment=True, mode="rb"):
    """
    Parse a BAM or SAM file using pysam, returning a list of PafAlignmentStats objects.
    "rb" for BAM files and "r" for SAM files
    """
    alignments = []
    with pysam.AlignmentFile(bamPath, mode) as aln_file:
        for read in aln_file.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            # Query information.
            queryName = read.query_name
            queryLen = (read.query_length if read.query_length is not None
                        else (len(read.query_sequence) if read.query_sequence else 0))
            queryStart = read.query_alignment_start
            queryEnd = read.query_alignment_end
            strand = '-' if read.is_reverse else '+'

            # Reference information.
            targetName = read.reference_name
            targetStart = read.reference_start
            targetEnd = read.reference_end
            targetLen = (refLengths[targetName] if refLengths and targetName in refLengths
                         else read.reference_length)

            # Calculate aligned length (from the query side).
            aligned_length = read.query_alignment_length
            # Sum deletion bases from the CIGAR (operation 2).
            deletion_bases = sum(length for op, length in read.cigartuples if op == 2)
            # Alignment block length = aligned query bases + deletion bases.
            alnBlockLen = aligned_length + deletion_bases

            # Use NM tag if available.
            if read.has_tag("NM"):
                nm_val = int(read.get_tag("NM"))
                NMTag = str(nm_val)
            else:
                nm_val = 0
                NMTag = None

            # Calculate numMatches as (aligned_length - NM + deletion_bases).
            numMatches = aligned_length - nm_val + deletion_bases

            readSeqLen = readLengths.get(queryName, queryLen) if readLengths else queryLen
            refSeqLen = targetLen

            csTag = read.get_tag("cs") if read.has_tag("cs") else None

            bamStat = PafAlignmentStats(
                queryName=queryName,
                queryLen=queryLen,
                queryStart=queryStart,
                queryEnd=queryEnd,
                strand=strand,
                targetName=targetName,
                targetLen=targetLen,
                targetStart=targetStart,
                targetEnd=targetEnd,
                numMatches=numMatches,
                alnBlockLen=alnBlockLen,
                readSeqLen=readSeqLen,
                refSeqLen=refSeqLen,
                csTag=csTag,
                NMTag=NMTag,
                globalAlignment=globalAlignment
            )
            alignments.append(bamStat)
    return alignments

def main():
    usage_str = "usage: %prog alignmentFile referenceFastaFile readFastqFile [options]\n" \
                "   (alignmentFile can be a PAF, BAM, or SAM file)"
    parser = OptionParser(usage=usage_str, version="%prog 0.3")
    
    parser.add_option("--readIdentity", action="store_true", default=False,
                      help="Print readIdentity of alignments.")
    parser.add_option("--alignmentIdentity", action="store_true", default=False,
                      help="Print alignmentIdentity.")
    parser.add_option("--readCoverage", action="store_true", default=False,
                      help="Print read coverage of alignments.")
    parser.add_option("--referenceCoverage", action="store_true", default=False,
                      help="Print reference coverage (targetEnd-targetStart)/targetLen.")
    parser.add_option("--mismatchesPerAlignedBase", action="store_true", default=False,
                      help="Print mismatches per aligned base.")
    parser.add_option("--deletionsPerReadBase", action="store_true", default=False,
                      help="Print deletions per read base.")
    parser.add_option("--insertionsPerReadBase", action="store_true", default=False,
                      help="Print insertions per read base.")
    parser.add_option("--readLength", action="store_true", default=False,
                      help="Print read lengths of aligned reads.")
    parser.add_option("--localAlignment", action="store_true", default=False,
                      help="Ignore unaligned prefix and suffix of each read (local alignment).")
    parser.add_option("--printValuePerReadAlignment", action="store_true", default=False,
                      help="Print the statistic for each read alignment.")
    parser.add_option("--noStats", action="store_true", default=False,
                      help="Do not print summary stats (avg, median, min, max).")
    parser.add_option("--perTargetDepth", action="store_true", default=False,
                      help="Print how many reads aligned to each target (i.e., 'depth').")

    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit("ERROR: Expected three arguments (alignmentFile, referenceFastaFile, readFastqFile).")

    alignmentFile, referenceFastaFile, readFastqFile = args

    # Optionally parse reference & read lengths.
    refLengths  = parseFastaLengths(referenceFastaFile) if os.path.isfile(referenceFastaFile) else {}
    readLengths = parseFastqLengths(readFastqFile) if os.path.isfile(readFastqFile) else {}

    # Determine the file type and parse accordingly.
    aln_lower = alignmentFile.lower()
    if aln_lower.endswith(".bam"):
        mode = "rb"
        alignments = parseBam(
            bamPath=alignmentFile,
            readLengths=readLengths,
            refLengths=refLengths,
            globalAlignment=(not options.localAlignment),
            mode=mode
        )
    elif aln_lower.endswith(".sam"):
        mode = "r"
        alignments = parseBam(
            bamPath=alignmentFile,
            readLengths=readLengths,
            refLengths=refLengths,
            globalAlignment=(not options.localAlignment),
            mode=mode
        )
    else:
        alignments = parsePaf(
            pafPath=alignmentFile,
            readLengths=readLengths,
            refLengths=refLengths,
            globalAlignment=(not options.localAlignment)
        )

    # Build a counter for target depths.
    target_counts = Counter(a.targetName for a in alignments)
    targetDepthDict = dict(target_counts)

    # Collect rows for CSV output.
    all_rows = []
    for a in alignments:
        row = {
            "QueryName": a.queryName,
            "TargetName": a.targetName,
            "ReadIdentity": a.readIdentity(),
            "AlignmentIdentity": a.alignmentIdentity(),
            "ReadCoverage": a.readCoverage(),
            "ReferenceCoverage": a.referenceCoverage(),
            "MismatchesPerAlignedBase": a.mismatchesPerAlignedBase(),
            "InsertionsPerReadBase": a.insertionsPerReadBase(),
            "DeletionsPerReadBase": a.deletionsPerReadBase(),
            "ReadLength": a.readLength(),
            "Depth": targetDepthDict[a.targetName]
        }
        all_rows.append(row)

    # Write CSV output.
    csv_out = "stats_out.csv"
    fieldnames = [
        "QueryName", "TargetName",
        "ReadIdentity", "AlignmentIdentity",
        "ReadCoverage", "ReferenceCoverage",
        "MismatchesPerAlignedBase", "InsertionsPerReadBase",
        "DeletionsPerReadBase", "ReadLength",
        "Depth"
    ]
    with open(csv_out, "w", newline="") as csv_fh:
        writer = csv.DictWriter(csv_fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"[INFO] Wrote stats CSV: {csv_out}")

    # Helper function to print summary statistics.
    def report(values, name):
        if not values:
            return
        if not options.noStats:
            avg_val = float(sum(values)) / len(values)
            med_val = float(numpy.median(values))
            min_val = float(min(values))
            max_val = float(max(values))
            print(f"Average{name}: {avg_val:.4f}")
            print(f"Median{name}: {med_val:.4f}")
            print(f"Min{name}: {min_val:.4f}")
            print(f"Max{name}: {max_val:.4f}")
        if options.printValuePerReadAlignment:
            print("Values{}: {}".format(name, "\t".join(f"{v:.4f}" for v in values)))

    # Report per-alignment stats.
    if options.readIdentity:
        vals = [a.readIdentity() for a in alignments]
        report(vals, "ReadIdentity")
    if options.alignmentIdentity:
        vals = [a.alignmentIdentity() for a in alignments]
        report(vals, "AlignmentIdentity")
    if options.readCoverage:
        vals = [a.readCoverage() for a in alignments]
        report(vals, "ReadCoverage")
    if options.referenceCoverage:
        vals = [a.referenceCoverage() for a in alignments]
        report(vals, "ReferenceCoverage")
    if options.mismatchesPerAlignedBase:
        vals = [a.mismatchesPerAlignedBase() for a in alignments]
        report(vals, "MismatchesPerAlignedBase")
    if options.deletionsPerReadBase:
        vals = [a.deletionsPerReadBase() for a in alignments]
        report(vals, "DeletionsPerReadBase")
    if options.insertionsPerReadBase:
        vals = [a.insertionsPerReadBase() for a in alignments]
        report(vals, "InsertionsPerReadBase")
    if options.readLength:
        vals = [a.readLength() for a in alignments]
        report(vals, "ReadLength")
    if options.perTargetDepth:
        print("\n# per-target depth (number of read alignments):")
        for tname, count in sorted(target_counts.items(), key=lambda x: x[0]):
            print(f"{tname}\t{count}")

if __name__ == '__main__':
    main()