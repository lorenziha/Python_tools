#!/Users/lorenziha/bin/miniconda3/bin/python

from os import replace
from typing import Optional
import pandas as pd
import re
import argparse

# Initialize parser
parser = argparse.ArgumentParser(description="It parses an allele file")

# Adding optional argument
parser.add_argument("-t", "--table", help="table file", required=True)
parser.add_argument("-o", "--output_prefix", help="output prefix", required=True)
parser.add_argument(
    "-s",
    "--sample_prefix",
    type=str,
    help="list of sample prefixes to include, default = all",
    default="all",
)

parser.add_argument(
    "-S",
    "--SNP_prevalence_percent",
    type=float,
    help="Minimum SNP prevalence across samples to be included, default = 1",
    default=1.0,
)

parser.add_argument(
    "-i",
    "--skip_heterozygote_positions",
    help="skip heterozygote positions in fasta file (it uses ambiguous nucleotides otherwise), default = False",
    action="store_true",
)
parser.add_argument(
    "-m",
    "--include_missing_positions",
    # type=str,
    help="include missing positions (./.) as '-' or 'N' , default = False",
    choices=["-", "N"],
    default=False,
)

args = parser.parse_args()

input_vcf = args.table
prefix_file = args.sample_prefix
skip_hetero = args.skip_heterozygote_positions
use_missing = args.include_missing_positions
output = args.output_prefix
MIN_PRESENT = args.SNP_prevalence_percent

# Dictionary of IUPAC ambiguities for nucleotides
# '*' is a deletion in GATK, deletions are ignored in consensus, lowercase consensus is udes when an
# 'N' or '*' is part of the genotype. Capitalization is used by some software but ignored by Geneious
# for example
ambiguities = {
    ".": "-",
    "*": "-",
    "A": "A",
    "C": "C",
    "G": "G",
    "N": "N",
    "T": "T",
    "*A": "a",
    "*C": "c",
    "*G": "g",
    "*N": "n",
    "*T": "t",
    "AC": "M",
    "AG": "R",
    "AN": "a",
    "AT": "W",
    "CG": "S",
    "CN": "c",
    "CT": "Y",
    "GN": "g",
    "GT": "K",
    "NT": "t",
    "*AC": "m",
    "*AG": "r",
    "*AN": "a",
    "*AT": "w",
    "*CG": "s",
    "*CN": "c",
    "*CT": "y",
    "*GN": "g",
    "*GT": "k",
    "*NT": "t",
    "ACG": "V",
    "ACN": "m",
    "ACT": "H",
    "AGN": "r",
    "AGT": "D",
    "ANT": "w",
    "CGN": "s",
    "CGT": "B",
    "CNT": "y",
    "GNT": "k",
    "*ACG": "v",
    "*ACN": "m",
    "*ACT": "h",
    "*AGN": "r",
    "*AGT": "d",
    "*ANT": "w",
    "*CGN": "s",
    "*CGT": "b",
    "*CNT": "y",
    "*GNT": "k",
    "ACGN": "v",
    "ACGT": "N",
    "ACNT": "h",
    "AGNT": "d",
    "CGNT": "b",
    "*ACGN": "v",
    "*ACGT": "N",
    "*ACNT": "h",
    "*AGNT": "d",
    "*CGNT": "b",
    "*ACGNT": "N",
}


def func_getNuc(genotype):
    a = "".join(sorted(set(genotype.replace("/", ""))))
    my_nucl = ambiguities[a]
    return my_nucl.upper()


def func_clean(my_key):
    return my_key.replace(".GT", "")


my_prefixes = list()
if prefix_file == "all":
    my_prefixes.append("ALL")
else:
    with open(prefix_file, "r") as prefixes:
        for prefix in prefixes.readlines():
            my_prefixes.append(prefix.rstrip())

vcf_chunck = pd.read_csv(
    input_vcf, sep="\t", index_col=["CHROM", "POS"], chunksize=1000
)

index_tmp = list()
seq = dict()

# print("********* skip_hetero =", skip_hetero)
# print("********* use_missing =", use_missing)

for piece in vcf_chunck:
    if my_prefixes[0] != "ALL":
        piece = piece[my_prefixes]

    samples = list(map(func_clean, piece.columns))
    for index, gt in piece.iterrows():
        flag_is_homo = True
        flag_use_missing = True
        count_present_snp = 0
        for i in range(0, len(samples)):
            if gt[i] != "./.":
                count_present_snp += 1
                # print(count_present_snp)

            if (not use_missing) and (gt[i] == "./."):
                flag_use_missing = False

            if (skip_hetero == True) and (len(set(gt[i])) > 2):
                flag_is_homo = False

            # print(count_present_snp, len(samples))

        if (count_present_snp / len(samples)) < MIN_PRESENT:
            flag_use_missing = False

        # print(flag_is_homo, flag_missing)
        if flag_is_homo and flag_use_missing:
            index_tmp.append(index)
            for i in range(0, len(samples)):
                if samples[i] in seq.keys():
                    seq[samples[i]] += func_getNuc(gt[i])
                else:
                    seq[samples[i]] = func_getNuc(gt[i])

        # print(seq[samples[i]])

# Write fasta file
output_fasta = output + ".fasta"
output_index = output + ".index.txt"
my_fasta = open(output_fasta, "w")
for i in range(0, len(samples)):
    if samples[i] in seq.keys():
        if use_missing == "N":
            seq[samples[i]] = seq[samples[i]].replace("-", "N")

        my_seq = ">" + samples[i] + "\n" + seq[samples[i]] + "\n"
        my_fasta.write(my_seq)
    else:
        Warning("Missing samples[i] sample")

my_fasta.close()

# Write index_file
my_index_file = open(output_index, "w")
for i in index_tmp:
    my_snp_pos = i[0] + "\t" + str(i[1]) + "\n"
    my_index_file.write(my_snp_pos)

my_index_file.close()

