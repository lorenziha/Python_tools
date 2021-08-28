import matplotlib.pyplot as plt
import re
import argparse

# Initialize parser
parser = argparse.ArgumentParser(description="It parses an allele file")

# Adding optional argument
parser.add_argument("-i", "--input", help="allele file", required=True)
parser.add_argument(
    "-a",
    "--alternative_allele",
    type=bool,
    help="print only alt alleles",
    default=False,
    choices=[True, False],
)
parser.add_argument(
    "-c", "--coverage", help="minimum read coverage of allele ", type=int, default=2
)
parser.add_argument("-q", "--quality", help="mean read quality", type=int, default=10)
parser.add_argument(
    "--min_fr_ratio",
    help="min fwd/rev read ratio",
    type=float,
    choices=[0.1, 0.2, 0.3, 0.4, 0.5],
    default=0.1,
)
parser.add_argument(
    "--max_fr_ratio",
    help="max fwd/rev read ratio",
    type=float,
    choices=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    default=0.9,
)

args = parser.parse_args()

QUAL = args.quality
RATIO_BOTTOM = args.min_fr_ratio
RATIO_TOP = args.max_fr_ratio
AA_FLAG = args.alternative_allele
MIN_COV = args.coverage

print("#Chromosome\tPosition\tA\tT\tC\tG\n")
chrom = ""
pos = 0
allele_info = ""
hist = {}
a, t, c, g, fwd, rev, tq = [0, 0, 0, 0, 0, 0, 0]


with open(args.input, "r") as allele_file:
    for line in allele_file.readlines():
        if re.match(r"CM\d+", line):
            chrom = line.rstrip()
            continue
        else:
            pos, snp_data = line.rstrip().split()
            a, t, c, g, fwd, rev, tq = map(int, snp_data.split(";"))

        total = sum([a, t, c, g])
        if total == 0:  # skip positions with zero reads.
            continue

        # print("total = ", total)
        # Filter out allele supported for less than $MIN_COV
        if (
            (a > 0 & a < MIN_COV)
            & (t > 0 & t < MIN_COV)
            & (c > 0 & c < MIN_COV)
            & (g > 0 & g < MIN_COV)
        ):
            continue

        # calculate allele freqs
        fa = int(100 * (a / total))
        ft = int(100 * (t / total))
        fc = int(100 * (c / total))
        fg = int(100 * (g / total))
        max_freq = max(fa, ft, fc, fg)

        # calculate average read qualitymy
        av_qual = tq / total

        # calculate forward/reverse ratio
        fr = fwd / total

        # Calculate histogram of frequences
        for freq in [fa, ft, fc, fg]:
            # ignore frequences for monoallelic positions.
            if freq > 0 and freq < 100 and av_qual < QUAL:
                # print("freq=", freq, ", hist=", hist[freq])
                if freq in hist:
                    hist[freq] += 1
                else:
                    hist[freq] = 1

            if fr >= RATIO_BOTTOM and fr <= RATIO_TOP and av_qual >= QUAL:
                if AA_FLAG:
                    if max_freq < 100:
                        print(f"{chrom}\t{pos}\t{fa}\t{ft}\t{fc}\t{fg}")
                else:
                    print(f"{chrom}\t{pos}\t{fa}\t{ft}\t{fc}\t{fg}")


# Print out histogram of frequences
print("\n# Histogram\n\n")
for freq in sorted(hist.keys()):
    print(f"{freq}\t{hist[freq]}")

