#!/usr/bin/env python3

import os
import sys
import re
import csv
from optparse import OptionParser
import numpy
from collections import Counter

"""
Script to calculate alignment statistics from a PAF file, extended to:
1) Include reference coverage per alignment (fraction of target spanned).
2) Report "depth" = total number of reads that mapped to each target.
"""

class PafAlignmentStats(object):
    """
    A helper class to store and compute per-read alignment statistics
    analogous to what was provided by margin.utils.ReadAlignmentStats for SAM.
    """
    def __init__(self, queryName, queryLen, queryStart, queryEnd, strand,
                 targetName, targetLen, targetStart, targetEnd,
                 numMatches, alnBlockLen, readSeqLen, refSeqLen,
                 csTag=None, NMTag=None, globalAlignment=True):
        """
        Store raw PAF information and, if present, the cs tag or NM tag
        for more detailed computations.
        """
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
        self.readSeqLen = readSeqLen  # from FASTQ or PAF field #2
        self.refSeqLen = refSeqLen    # from FASTA or PAF field #7

        self.csTag = csTag
        self.NMTag = NMTag
        self.globalAlignment = globalAlignment

        # Precompute differences if we have a cs tag or NM tag
        self._parsed = False
        self._numMismatches = 0
        self._numInsertions = 0
        self._numInsertionBases = 0
        self._numDeletions = 0
        self._numDeletionBases = 0

        self._parseDifferences()

    def _parseDifferences(self):
        """
        Parse cs:Z: or NM:i: from PAF to compute mismatch/insertion/deletion counts.
        Handles both 'cs=short' (e.g. :100 *ac +gg -t)
        and 'cs=long' (e.g. =ACGT:10*ac:gt+GG-tt).
        """
        if self.csTag is not None:
            # Regex to capture each token in --cs=short or --cs=long
            #   =ACGT   => a block of matched bases
            #   :10     => short notation for 10 consecutive matches
            #   *ac:gt  => mismatch, 'ac' replaced by 'gt'
            #   +GG     => insertion
            #   -tt     => deletion
            pattern = re.compile(r'(=[A-Za-z]+|:\d+|\*[A-Za-z]+:[A-Za-z]+|\+[A-Za-z]+|\-[A-Za-z]+)')
            tokens = pattern.findall(self.csTag)

            for token in tokens:
                if token.startswith('='):
                    # e.g. "=ACGT" => matched block
                    # This is already accounted for in numMatches (PAF col #10).
                    pass
                elif token.startswith(':'):
                    # :N => N consecutive matches
                    pass
                elif token.startswith('*'):
                    # e.g. "*ac:gt" => mismatch from 'ac' to 'gt'
                    # We count each base as a mismatch.
                    mismatch_sub = token[1:].split(':')  # remove leading '*'
                    mismatch_length = max(len(mismatch_sub[0]), len(mismatch_sub[1]))
                    self._numMismatches += mismatch_length
                elif token.startswith('+'):
                    # e.g. "+GG" => insertion of "GG"
                    inserted_seq = token[1:]
                    self._numInsertions += 1
                    self._numInsertionBases += len(inserted_seq)
                elif token.startswith('-'):
                    # e.g. "-tt" => deletion of "tt"
                    deleted_seq = token[1:]
                    self._numDeletions += 1
                    self._numDeletionBases += len(deleted_seq)

        elif self.NMTag is not None:
            # NM = mismatches + inserted + deleted bases
            nm_val = int(self.NMTag)
            self._numMismatches = nm_val
            # Not distinguishing insertion vs deletion if only NM is available

        self._parsed = True

    def readLength(self):
        """
        Return the read length used in coverage stats.
        If globalAlignment=False, we only count the aligned portion of the read.
        """
        if not self.globalAlignment:
            return max(0, self.queryEnd - self.queryStart)
        return self.queryLen

    def readCoverage(self):
        """
        coverage = (aligned_length / read_length).
        For PAF, aligned_length from the read perspective is (queryEnd - queryStart).
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        aligned_len = max(0, self.queryEnd - self.queryStart)
        return aligned_len / r_len

    def alignmentIdentity(self):
        """
        alignmentIdentity = numMatches / alignment_length
        alignment_length is the sum of matched + mismatched + inserted + deleted bases.
        Minimally, this is 'alnBlockLen' in PAF col #11.
        """
        if self.alnBlockLen == 0:
            return 0.0
        return float(self.numMatches) / float(self.alnBlockLen)

    def readIdentity(self):
        """
        readIdentity = numMatches / read_length
        If local alignment, read_length might be just the aligned portion.
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self.numMatches) / r_len

    def mismatchesPerAlignedBase(self):
        """
        mismatches / alignment_length
        """
        if self.alnBlockLen == 0:
            return 0.0
        return float(self._numMismatches) / float(self.alnBlockLen)

    def insertionsPerReadBase(self):
        """
        insertion bases / read_length
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self._numInsertionBases) / r_len

    def deletionsPerReadBase(self):
        """
        deletion bases / read_length
        """
        r_len = float(self.readLength())
        if r_len == 0:
            return 0.0
        return float(self._numDeletionBases) / r_len

    def referenceCoverage(self):
        """
        Fraction of the target sequence covered by this alignment:
        (targetEnd - targetStart) / targetLen
        """
        if self.targetLen == 0:
            return 0.0
        aligned_len = max(0, self.targetEnd - self.targetStart)
        return aligned_len / float(self.targetLen)


def parseFastaLengths(fastaPath):
    """
    If needed, parse the reference FASTA to get lengths by name.
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
    If needed, parse the read FASTQ to get lengths by name.
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
    Parse the PAF, returning a list of PafAlignmentStats objects.
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

            # If read/ref length dictionaries are provided, override
            if readLengths and queryName in readLengths:
                readLen = readLengths[queryName]
            else:
                readLen = queryLen

            if refLengths and targetName in refLengths:
                refLen = refLengths[targetName]
            else:
                refLen = targetLen

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

def main():
    usage_str = "usage: %prog pafFile readFastqFile referenceFastaFile [options]"
    parser = OptionParser(usage=usage_str, version="%prog 0.1")
    
    parser.add_option("--readIdentity", action="store_true", default=False,
                      help="Print readIdentity of alignments.")
    parser.add_option("--alignmentIdentity", action="store_true", default=False,
                      help="Print alignmentIdentity.")
    parser.add_option("--readCoverage", action="store_true", default=False,
                      help="Print read coverage of alignments.")
    parser.add_option("--referenceCoverage", action="store_true", default=False,
                      help="Print reference coverage of alignments (targetEnd-targetStart)/targetLen.")
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
        sys.exit("ERROR: Expected three arguments (pafFile, readFastqFile, referenceFastaFile).")

    pafFile, readFastqFile, referenceFastaFile = args

    # Optionally parse read & reference lengths
    readLengths = parseFastqLengths(readFastqFile) if os.path.isfile(readFastqFile) else {}
    refLengths  = parseFastaLengths(referenceFastaFile) if os.path.isfile(referenceFastaFile) else {}

    # Parse PAF
    alignments = parsePaf(
        pafPath=pafFile,
        readLengths=readLengths,
        refLengths=refLengths,
        globalAlignment=(not options.localAlignment)
    )

    # Build a counter of target -> read count
    target_counts = Counter(a.targetName for a in alignments)

    # Convert to a dict { targetName: number_of_reads }
    targetDepthDict = dict(target_counts)

    # We'll collect data rows in a list for CSV output:
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
            # Depth for this alignment's target
            "Depth": targetDepthDict[a.targetName]
        }
        all_rows.append(row)

    # Now write out CSV
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

    # Helper to print stats
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
            print("Values{}: {}".format(
                name, "\t".join(f"{v:.4f}" for v in values)
            ))

    # Per-alignment stats
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

    # Finally, for the 'depth' = how many reads mapped to each target
    if options.perTargetDepth:
        target_counts = Counter(a.targetName for a in alignments)
        print("\n# per-target depth (number of read alignments):")
        for tname, count in sorted(target_counts.items(), key=lambda x: x[0]):
            print(f"{tname}\t{count}")

if __name__ == '__main__':
    main()