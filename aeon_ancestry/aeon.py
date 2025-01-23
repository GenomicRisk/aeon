"""
aeon

Naomi Warren 2023
"""

import argparse
import logging
import os

import colorama
import pandas as pd
import torch
from colorama import Fore, Style

from aeon_ancestry.estimatorModels import PopulationMixtureModelRandom
from aeon_ancestry.genotype_from_vcf import Genotypes
from aeon_ancestry.visualisePCA import saveIndividualPCAplot, saveTrioPCAplot, transformToPCA


def aeon(args):
    print(f"{Fore.BLUE}*** AEON -- Ancestry Estimation ***{Style.RESET_ALL}")

    # Process input parameters
    outfix = args.out
    if not outfix:
        in_split = os.path.basename(args.vcffile).split("_")
        if len(in_split) > 1:
            outfix = in_split[0]
        else:
            outfix = in_split.split(".")[0]

    af_data = pd.read_table(args.allele_freqs)
    loci_list = af_data["VAR_ID"]
    pop_names = af_data.columns[4:]
    af_floats = af_data.drop(columns=["CHROM", "START", "STOP", "VAR_ID"]).astype("float32")
    af_tensor = torch.tensor(af_floats.values)

    # Extract genotypes from file
    print()
    print("Extracting genotypes ...")
    g = Genotypes(args.vcffile, loci_list)
    samples = g.getSamples()
    print()

    # Save PCA visualisation plots
    if not args.no_visualisation:
        print("Making PCA visualisation ... ", end="")
        if args.allele_freqs == "refs/g1k_allele_freqs.txt":
            subset = False
        else:
            subset = True
        pca_df = transformToPCA(g, subset=subset)

        if args.inheritance and len(samples) <= 3:
            saveTrioPCAplot(outfix, pca_df)
        else:
            for sample in samples:
                saveIndividualPCAplot(sample, pca_df)
        log_df = pca_df.T
        log_df.reset_index(inplace=True, names="Sample")
        print()
    else:
        log_df = pd.DataFrame({"Sample": samples})

    # Log processing info
    log_df["FractionLociImputed"] = g.imputationFraction()
    log_df["FractionLociNonCalled"] = log_df.Sample.apply(g.nonCalledFraction)
    log_df["NumberMultiAllelic"] = log_df.Sample.apply(g.multiAllelic)

    # Estimate population memberships for each sample in input VCF
    result_df = pd.read_table(args.population_labels)
    result_df.set_index("Population", inplace=True)

    losses = []
    for sample in samples:
        print(f"Estimating membership for sample {sample} ...")
        dosages = g.dosageForSample(sample)
        model = PopulationMixtureModelRandom(dosages, af_tensor)
        result_mle, loss = model.est_mle(trace=100)

        losses.append(loss / len(dosages))

        new = pd.DataFrame(torch.round(result_mle, decimals=2), index=pop_names, columns=[sample])
        result_df = result_df.merge(new, left_index=True, right_index=True)

    log_df["LossPerLoci"] = losses
    log_df.to_csv(f"{outfix}_ae_stats.csv", index=False)
    result_df.to_csv(f"{outfix}_ae.csv")

    # Report results
    superpop_res = result_df.groupby(["Superpopulation"]).sum(numeric_only=True)
    print()
    print(f"{Fore.GREEN}                ------------ SUMMARY ------------")
    print(f"{superpop_res}{Style.RESET_ALL}")
    print("See outfiles for further details.")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="""Predict population memberships given a VCF with input samples and known population allele frequencies. Note: \
                                     any variants present in AF file not present in the VCF will be imputed for all samples as homozygous reference."""
    )
    parser.add_argument(
        "vcffile", help="Input VCF (or VCF.GZ or BCF. Compressed files must be indexed.)"
    )
    parser.add_argument(
        "-a",
        "--allele_freqs",
        default="refs/g1k_allele_freqs.txt",
        help="Tab-delimited file containing allele frequencies for each population with one variant per line, and the header line CHROM     START   STOP    VAR_ID    POP1   [POP2 ...]",
    )
    parser.add_argument(
        "--population_labels",
        default="refs/pop2super.txt",
        help="Tab-delimited file listing populations in allele_freqs file and their corresponding superpopulation, and the header line: Superpopulation  Population",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=False,
        help="PREFIX for output files. Default is to capture all chars of input filename before the first underscore e.g. 450-100_variants_file.vcf -> 450-100",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Run in verbose mode (prints INFO level logs to stderr as well as WARNINGs)",
    )
    parser.add_argument(
        "--inheritance",
        action="store_true",
        help="Run in inheritance mode - all samples from VCF will be plotted/visualised together. Note: only works for <=3 samples.",
    )
    parser.add_argument(
        "--no_visualisation", action="store_true", help="Do not output any visualisation files."
    )

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)

    colorama.init()

    aeon(args)