import matplotlib.pyplot as plt
import re
import argparse
import matplotlib.pyplot as plt

# Initialize parser
parser = argparse.ArgumentParser(
    description="It parses an allele file and plots histogram"
)

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

chrom = ""
pos = 0
allele_info = ""
a, t, c, g, fwd, rev, tq = [0, 0, 0, 0, 0, 0, 0]
freq_list = []

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
            if (
                freq > 0
                and freq < 100
                and av_qual >= QUAL
                and fr >= RATIO_BOTTOM
                and fr <= RATIO_TOP
            ):
                freq_list.append(freq)


# Make histogram
# setting the ranges and no. of intervals
range = (0, 100)
bins = 100

# plotting a histogram
plt.hist(freq_list, bins, range, color="green", histtype="bar", rwidth=0.8)

# x-axis label
plt.xlabel("Allele frequency")
# frequency label
plt.ylabel("Counts")
# plot title
plt.title("Distribution of alternative allele frequences")

# function to show the plot
plt.show()

