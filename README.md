# AEon<sup>1</sup> - Ancestry Estimation #

_aeon_ performs global genetic Ancestry Estimation from genome-wide SNPs. By default, aeon compares variants from your input VCF file
against the per-population allele frequencies (AFs) provided in `refs/g1k_allele_freqs.txt`. This file contains AFs at
128097 ancestry-informative loci, calculated for the 26 populations provided on [1000 Genomes](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

Aeon accepts input as VCF, VCF.GZ or BCF. If you are using a VCF.GZ or BCF, aeon requires the corresponding index file
to exist in the same directory.

The tool provides 2 output CSV files:

- out_ae.csv       -> table of estimated ancestry fractions
- out_ae_stats.csv      -> table with top 3 PC projections and some information about the run,

and one or more visualisation files.

## Installing from PyPi ##

To install:

```bash
pipx install aeon-ancestry
```

Then you can run aeon from anywhere using:

```bash
aeon -h
```

## Running aeon using python poetry ##

To install:

```bash
cd aeon
poetry install
```

Then you can run aeon with:

```bash
poetry run aeon -h
```

## Running aeon using the BitBucket repo ##

Before starting, make sure you have a working version of python (tested using Python 3.10) and pip.

Clone the repo to your local machine and install the requirements in `requirements.txt` using pip (you may want to do this in
a virtual environment e.g. using conda).

```bash
$ git clone git@bitbucket.org:wrensong/aeon.git
$ cd aeon

# Install requirements (recommended using a virtual env):
$ conda create -n my_env python=3.10
$ conda activate my_env
(my_env)$ pip install -r requirements.txt

```

You should then be able to run aeon from the command line like so (make sure you have activated your virtual env, if you made one):

```bash
(my_env)$ python aeon.py sample_variants.bcf -o my_output
```

## Running using the Docker Image ##

A docker image is available [here](https://hub.docker.com/r/naomiwren/aeon). Within the container you can run aeon like so:

```bash
aeon.py sample_variants.bcf -o my_output
```

## Example usage ##

To see all possible input options, run `python aeon.py -h`:

```bash
(my_env)$ python aeon.py -h
usage: aeon.py [-h] [-a ALLELE_FREQS] [--population_labels POPULATION_LABELS] [-o OUT] [-v] [--inheritance] [--no_visualisation] vcffile

Predict population memberships given a VCF with input samples and known population allele frequencies. Note: any variants present in AF file
not present in the VCF will be imputed for all samples as homozygous reference.

positional arguments:
  vcffile               Input VCF (or VCF.GZ or BCF. Compressed files must be indexed.)

options:
  -h, --help            show this help message and exit
  -a ALLELE_FREQS, --allele_freqs ALLELE_FREQS
                        Tab-delimited file containing allele frequencies for each population with one variant per line, and the header line
                        CHROM START STOP VAR_ID POP1 [POP2 ...]
  --population_labels POPULATION_LABELS
                        Tab-delimited file listing populations in allele_freqs file and their corresponding superpopulation, and the header
                        line: Superpopulation Population
  -o OUT, --out OUT     PREFIX for output files. Default is to capture all chars of input filename before the first underscore e.g.
                        450-100_variants_file.vcf -> 450-100
  -v, --verbose         Run in verbose mode (prints INFO level logs to stderr as well as WARNINGs)
  --inheritance         Run in inheritance mode - all samples from VCF will be plotted/visualised together. Note: only works for <=3 samples.
  --no_visualisation    Do not output any visualisation files.
```

**For a generic single- or multi-sample VCF from WGS data:**

`$ python aeon.py sample_variants.bcf -o my_output`

Output:

```txt
my_output_ae.csv
my_output_ae_stats.csv
sample1_PCA_plot.png
[sample2_PCA_plot.png ...]
```

**For a duo/trio multi-sample VCF:**

`$ python aeon.py sample_variants.bcf -o my_output --inheritance`

**Without a specified output file prefix (use with caution):**

It is possible not to specify an output file prefix. All characters before the first `_` in the input VCF filename will be used as a prefix.

`$ python aeon.py prefix_example_sample_variants.bcf`

Output:

```txt
prefix_ae.csv
prefix_ae_stats.csv
sample1_PCA_plot.png
[sample2_PCA_plot.png ...]
```

**For a region-bound VCF:**
_**(e.g. data generated from WES)**_

If your VCF does not cover the whole genome, e.g. the data was generated with Whole Exome Sequencing, you will want to subset the reference AF file to your regions of coverage. Using a bed file containing your regions of interest, you can easily do this using bedtools (tested with bedtools v2.31.0; installation instructions [here](https://bedtools.readthedocs.io/en/latest/content/installation.html)) and the following command:

```
(my_env) $ bedtools intersect -a refs/g1k_allele_freqs.txt -b regions.bed -header > intersected_AFs.txt
(my_env) $ python aeon.py sample_variants.bcf -o my_output -a intersected_AFs.txt
```

**With an alternative reference AF file:**

`$ python aeon.py sample_variants.bcf -a reference_afs.txt --population_labels pops_table.txt -o my_output`

See [Reference Populations](#markdown-header-reference-populations) for information regarding the format of `reference_afs.txt` and `pops_table.txt`.

### Example data ###

The `example` directory contains example trio input data taken from 1000 Genomes in the file `g1k_trio_ASW.bcf`, as well as corresponding
output data. Once you have [set up aeon on your machine](#markdown-header-running-aeon-using-the-bitbucket-repo), you can run the example data yourself and compare your results:

```
(my_env) $ python aeon.py example/g1k_trio_ASW.bcf --inheritance -o ASW_example
```

This should produce 3 files: `ASW_example_ae.csv`, `ASW_example_ae_stats.csv` and `ASW_example_PCA_plot.png`. Looking in `ASW_example_ae.csv`, you will notice that the majority of the estimated ancestry (0.6-0.7) is assigned to the ASW (African Ancestry in Southwest US) population. The remainder is split amongst the other populations within the AFR superpopulation. This reflects the heterogeneity within and overlap between African populations. By examining the plots in `ASW_example_PCA_plot.png`, you can see how the trio is located amongst the mint green ASW samples, with some samples in different shades of green (other African populations) in the vicinity.

**Note:** Since this example data is from 1000 Genomes, the parent samples were used as part of the reference set to calculate the allele frequencies provided in `refs/g1k_allele_freqs.txt`.
However, the child sample was not used for calculating AFs, as only the 2504 _unrelated_ individuals from the 1000 Genomes dataset were selected.

## Interpreting Output ##

PLEASE NOTE that the underlying model assumes that the ancestry of any input sample can be completely explained by the 26 reference populations provided. Ancestry from a non-characterised population will not be captured; instead it will be assigned to the closest known population. If you would like to use your own allele frequency data including a different population of particular interest to your study, see [Reference Populations](#markdown-header-reference-populations).

### ae.csv file ###

Aeon provides ancestry estimation fractions per population to 2 decimal places. From analysis on trio (mother/father/child) data:

- score > 0.1    -> significant
- 0.05 < score < 0.1  -> likely significant
- score < 0.05   -> likely insignificant
- score < 0.02   -> noise

It is important to bear in mind that ancestral populations are not completely isolated from each other - some populations are more closely related to each other than others due to their historical geographic context, and this is reflected in grouping into _superpopulations_. The `ae.csv` file records which superpopulation each population belongs to, to assist interpretation. Some guidelines:

- Superpopulation estimates are more robust, as superpopulations are more easily separated than populations. A superpopulation score > 0.1 is likely significant, even if population-level resolution is ambiguous.
- The EUR superpopulation has a high degree of overlap between populations. It is not unreasonable to see estimated ancestry fractions split across all EUR populations for a EUR sample.
- The AFR and AMR superpopulations show a high degree of variation within populations. As such, low scores in these populations are more ambiguous.
- PUR is the closest AMR population to the EUR superpopulation, reflecting historical admixture. If the _only_ AMR population to receive a score > 0.02 is PUR, and the sample scores highly for EUR populations, treat the AMR estimate with caution. Visualisation can help in these cases.

### ae_sample_stats.csv file ###

This output file provides some basic information about how samples were processed. For each sample, it provides the following information:

- The **PC1**, **PC2** and **PC3** values calculated from the genotype vector, used to plot the sample against reference populations
- **FractionLociImputed:** Number of ancestry-informative loci without a record present in the input VCF, divided by 128097. All such loci are imputed as reference. This fraction should be the same for all samples in a VCF.
- **FractionLociNonCalled:** Number of ancestry-informative alleles with a record present in the input VCF, but non-called genotype, divided by 2\*128097. This often happens if you have merged multiple VCFs where one file contains a variant absent in the other file. All such loci are imputed as reference. This fraction will likely differ across samples in a VCF.
- **NumberMultiAllelic:** Number of ancestry-informative loci with a different variant allele from the one listed in the allele frequency file. All non-reference values are compressed to 1.
- **LossPerLoci:** Measure of model 'loss' (as returned by `pyro.infer.SVI().step()`) divided by number of loci.

### Visualisation ###

To help visualise the distribution of ancestral populations, Principal Component Analysis was performed on the variant allele frequencies provided in
`refs/g1k_allele_freqs.txt`. Aeon then returns one `.png` file per sample in your input VCF, plotting each sample against a backdrop of all the reference samples
projected onto the top 3 principal components (PCs). The figure is composed of 3 plots, showing 2 dimensions at a time - PC1 vs PC2 in the top left, PC3 vs PC2 in the top right, and PC1 vs PC3 in the bottom left.

The output for each sample will be saved to `sampleName_PCA_plot.png`. If your VCF file contains a duo or trio of related samples (e.g. parents and child), you can optioinally add the `--inheritance` flag, which will plot all input samples on the same figure. In this case, the output file will be named according to the specified output file prefix.

This visualisation can assist with interpretation of unusual cases. It is particularly interesting to observe that cases with admixture are very poorly characterised by PCA, often appearing on the plot outside the known population clusters and closest to an unrelated 'in-between' cluster.

Note that if you use a modified allele frequency file that contains different populations from the default file, this will NOT be reflected in the output plots. Each circular point in the output visualisation represents a sample from the 1000 Genomes dataset, whose coordinates have been pre-computed based on their genotype, NOT generated from the input allele frequency file.

## Background ##

### Model ###

Aeon is based on probabilistic model, whereby population membership fractions _p_ are estimated from a sample's genotype _G_ given a data likelihood model _Pr(G|p)_.

Assume an individual _i_ is genotyped at each of _n_ unlinked biallelic autosomal loci to give a vector of alt allele dosages _d<sub>i</sub>_ where _d<sub>i,j</sub>_ ∈ {0, 1, 2}. _d<sub>i,j</sub>_ is sampled on two draws from a binomial distribution with probabilities depending on the genotype  _g<sub>i,j</sub>_ and variant allele frequencies A<sub>k,j</sub> in each population _k_ ∈ 1..._m_. The best estimate for population membership fractions is then calculated using Maximum Likelihood Estimation on this model.

Note that several assumptions are implicit in this model:
- All contributing ancestral populations are represented in the input allele frequency matrix _A_
  - as such, if you expect a different population to significantly contribute to ancestry, use a modified AF file with the flag `-a`
- Loci are biallelic
  - see 'NumberMultiAllelic' column in your out_ae_stats.csv file to see if this assumption holds for your data. A small proportion of input loci that are multi-allelic is not likely to have a significant effect on the results.
- Loci are unlinked
  - The loci in the default AF file are selected to specifically adhere to this assumption. If you supply your own modified AF file with different loci, it is up to you to make sure they are unlinked.
- No missing values in input
  - Due to the compressed nature of VCF files, if a locus/variant from the AF file is missing from your input VCF, the tool assumes that all samples in the VCF are homozygous reference at this locus. The fraction of loci imputed as reference in this manner can be seen in the 'FractionLociImputed' column in your out_ae_stats.csv file. It is not uncommon to se a fraction between 0.4-0.5 if your input does not include rows that are homozygous reference.
  - If you DO NOT WANT this behaviour (e.g. you have only sequenced specific regions of the genome like exons, and don't want all introns to be imputed as reference!), you will need to [subset the AF file](#markdown-header-for-a-region-bound-vcf) for the desired regions using bedtools and run aeon with this subsetted AF file. Note that using a subset of loci provides less information to the tool, resulting in a faster running time but less accurate predictions.

### Reference Populations ###

Allele frequencies from 26 reference populations were obtained using data from [1000 Genomes](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).
As such, this model cannot predict membership in populations _outside_ this dataset - if other ancestry is present in your sample of interest,
it will be captured (erroneously) across the 26 known populations.

However, you can use your own allele frequency reference file if you have access to population AFs of specific interest to your study. Make sure your allele frequency file is tab-delimited with the header `CHROM    START    STOP    VAR_ID   POP1   [POP2 ...]`, and format your variant IDs as `chrN_pos_REF_ALT`. You will also need to supply a tab-delimited file of population labels and their corresponding superpopulations.
See `g1k_allele_freqs.txt` and `pop2super.txt` in the `refs` directory for examples.

### 128097 'ancestry-informative' loci ###

An initial list of (GRCh37) variants of interest was retrieved from [mpinese/mgrb-manuscript](https://github.com/mpinese/mgrb-manuscript/blob/master/manuscript/pca-1000g/01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list). This list contained 133872 variants, selected using several criteria:

- **Easy to sequence:** to reduce the number of missing or poor-quality values in both training and test sets
- **Ancestry informative:** loci that demonstrated Hardy-Weinberg Equilibrium within populations, but not across all populations
- **LD-pruned:** variants were pruned based on Linkage Disequilibrium in order to discard loci correlated through co-inheritance patterns
- **Excluded highly conserved regions:** since homogeneity across all populations would give no discriminative power.
- For further info, see [Pinese et al. (2020)]( https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6978518/).

This initial list was then further refined to suit the reference set from 1000 Genomes. Variants in GRCh37 that had become reference in GRCh38 were excluded, as were variants for which the genotype was missing for some of the 2504 reference samples, resulting in the final list of 128097 variants.

## Who do I talk to? ##

Contact Naomi for more info at <nwarren@ccia.org.au>

<sup>1</sup> _**aeon** /ˈiːən/ (noun):_ an indefinite and very long period of time
