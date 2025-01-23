import argparse, logging
import pandas as pd
import torch

from pysam import VariantFile

class Genotypes ():
    def __init__(self, vcfpath, loci_list):
        self.vcf = VariantFile(vcfpath)
        self.samples = self.vcf.header.samples
        self.imputed_frac = 0.0
        self.noncalled = {}
        self.multi_allelic = {}

        try:
            self.vcf.fetch(region="chr1:1-2")
        except ValueError as e:
            print(e)
            print("Index file required for VCF.GZ or BCF file.")
            exit(2)

        self.geno = self.extractGenotypes(loci_list)

    def getSamples (self):
        return self.samples
    
    def imputationFraction (self):
        return round(self.imputed_frac, 4)
    
    def nonCalledFraction (self, sample):
        return (self.noncalled[sample] / (2*128097) )
    
    def multiAllelic (self, sample):
        return self.multi_allelic[sample]

    def extractGenotypes (self, var_list):

        genos = {}
        var_ids = []
        for sample in self.samples:
            genos[sample] = []

        imputed = 0
        for variant in var_list:
            rid = variant.strip()
            var_ids.append(rid)
            (chr, pos, ref, var)  = rid.split("_")

            vars = self.vcf.fetch(region=f"{chr}:{pos}-{pos}")
            rec = next(vars, False)
            if rec:
                for s_name, s_values in rec.samples.items():
                    genos[s_name].append(s_values['GT'])
            else:
                imputed += 1
                for s_name in self.samples:
                    genos[s_name].append((0,0))

        for sample, geno in genos.items():
            geno = pd.Series(geno)
            self.multi_allelic[sample] = geno.apply(self.isMultiAllelic).sum()
            self.noncalled[sample] = geno.apply(self.notCalled).sum()
            genos[sample] = torch.as_tensor(geno.apply(self.genoTupleToInt))

        genos['var_ids'] = var_ids

        self.imputed_frac = imputed/len(var_ids)
        logging.info(f"    {sum(self.noncalled.values())/2} genotypes (across all samples) were not called and will be imputed as reference.")
        logging.info(f"    {imputed} loci had no entry in VCF and were imputed as homozygous reference ({round(self.imputed_frac*100, 2)}%).")
        return genos
    
    def dosageForSample(self, sample):
        return self.geno[sample]
    
    def variantList(self):
        return self.geno['var_ids']

    @staticmethod
    def isMultiAllelic (gtuple):
        if gtuple[0] is not None and gtuple[0] > 1 :
            return True
        elif gtuple[1] is not None and gtuple[1] > 1:
            return True
        else:
            return False
        
    @staticmethod
    def notCalled (gtuple):
        nc = 0
        for g in gtuple:
            if g is None:
                nc += 1
        return nc

    # HELPER FUNCTIONS
    @staticmethod
    def genoTupleToInt (gtuple):
        # Note: compresses multi-allelic sites
        HOMREF, HET, HOMALT = range(3)

        #Note: imputes missing (None) as reference
        if gtuple[0] == 0 or gtuple[0] is None:
            if gtuple[1] == 0 or gtuple[1] is None:
                return HOMREF
            else:
                return HET
        else:
            if gtuple[1] == 0 or gtuple[1] is None:
                return HET
            else:
                return HOMALT

def main (args):
    with open(args.regions) as vars:
        loci = vars.readlines()

    g = Genotypes(args.vcffile, loci)
    df = g.geno

    print(len(df), len(df.columns))
    print(df.head())
    for sample in g.getSamples():
        print (g.dosageForSample(sample)[0:5])


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='Get genotypes of samples from VCF file as a dataframe.')
    parser.add_argument('vcffile',
                        help='Input VCF')
    parser.add_argument('-r', '--regions', required=False,
                        help="Text file containing variants of interest (one variant per line), e.g. chr1_1000000_A_T. Positions are 1-indexed, as in a bcftools query region.")

    args = parser.parse_args()
    main(args)