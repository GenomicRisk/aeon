import logging, argparse

import pandas as pd
import torch
import torch.distributions.constraints as constraints
import pyro
import pyro.distributions as dist
import pyro.infer
import pyro.optim

import time

import colorama
from colorama import Fore
from colorama import Style

from aeon.genotype_from_vcf import Genotypes

class PopulationMixtureModelRandom:
    def __init__(self, dosages: torch.Tensor, allele_freqs: torch.Tensor):
        '''
        dosages: n-vector of observed alternate allele dosages at n loci, values in {0, 1, 2}
        allele_freqs: n x m matrix of alternate allele frequencies at n loci across m populations
        '''
        pyro.clear_param_store()
        assert dosages.size(0) == allele_freqs.size(0)
        self.dosages = dosages
        self.allele_freqs = allele_freqs


    @staticmethod
    def model_stochastic(d, f, val_args=True):
        num_pops = f.size(1)
        alpha0 = torch.ones(num_pops)
        p = pyro.sample('p', dist.Dirichlet(alpha0))
        r = torch.matmul(f, p)
        pr = torch.column_stack(((1-r)**2, 2*r*(1-r), r**2))
        with pyro.plate('locus', d.size(0)) as i:
            pyro.sample('d_{}'.format(i), dist.Categorical(pr[i,:], validate_args=val_args), obs=d[i])


    @staticmethod
    def model_fixed(d, f, val_args=True):
        num_pops = f.size(1)
        p = pyro.param('p', torch.ones(num_pops)/num_pops, constraint=constraints.simplex)
        r = torch.matmul(f, p)
        pr = torch.column_stack(((1-r)*(1-r), 2*r*(1-r), r*r))
        with pyro.plate('locus', d.size(0)) as i:
            pyro.sample('d_{}'.format(i), dist.Categorical(pr[i,:], validate_args=val_args), obs=d[i])


    @staticmethod
    def guide_svi(d, f):
        num_pops = f.size(1)
        alpha_q = pyro.param('alpha_q', torch.ones(num_pops), constraint=constraints.positive)
        pyro.sample('p', dist.Dirichlet(alpha_q))


    def est_mle(self, lr=0.05, min_steps=200, max_steps=1000, termination_stable=200, tol_loss=2, lrd=0.1**(1/1000), trace=0, val_args=True):
        optimizer = pyro.optim.ClippedAdam({'lr': lr, 'lrd': lrd})

        svi = pyro.infer.SVI(lambda d, f: self.model_fixed(d, f, val_args=val_args), lambda d, f, : None, optimizer, loss=pyro.infer.Trace_ELBO())

        best_loss = None
        for step in range(max_steps):
            loss = svi.step(self.dosages, self.allele_freqs)
            if best_loss == None or loss < best_loss-tol_loss:
                best_loss = loss
                best_loss_step = step
            if step > min_steps and step-best_loss_step >= termination_stable:
                break
            if trace > 0 and (step+1) % trace == 0:
                logging.info('MLE step {}: loss {:.2f} (best loss {:.2f} at step {})'.format(step+1, loss, best_loss, best_loss_step+1))

        return pyro.param('p').detach(), loss

    def est_svi(self, lr=0.1, min_steps=1000, max_steps=10000, termination_stable=1000, tol_loss=2, tol_prec=0.01, lrd=0.1**(1/5000), trace=0, val_args=True):
        optimizer = pyro.optim.ClippedAdam({'lr': lr, 'lrd': lrd})

        svi = pyro.infer.SVI(lambda d, f: self.model_stochastic(d, f, val_args=val_args), self.guide_svi, optimizer, loss=pyro.infer.Trace_ELBO(), val_args=val_args)

        best_loss = None
        best_precision = None
        for step in range(max_steps):
            loss = svi.step(self.dosages, self.allele_freqs)
            precision = torch.sum(pyro.param('alpha_q')).item()
            if best_loss == None or loss < best_loss-tol_loss:
                best_loss = loss
                best_loss_step = step
            if best_precision == None or precision > best_precision*(1+tol_prec):
                best_precision = precision
                best_precision_step = step
            if step > min_steps and step-best_loss_step >= termination_stable and step-best_precision_step >= termination_stable:
                break
            if trace > 0 and (step+1) % trace == 0:
                logging.info('SVI step {}: loss {:.2f}, precision {:.2f} (best loss {:.2f} at step {}, best precision {:.2f} at step {})'.format(step+1, loss, precision, best_loss, best_loss_step+1, best_precision, best_precision_step+1))

        return pyro.param('alpha_q').detach()

    def est_mcmc(self, warmup=100, n_steps=200, trace=0, val_args=True):
        kernel = pyro.infer.mcmc.NUTS(lambda d, f: self.model_stochastic(d, f, val_args=val_args))
        mcmc = pyro.infer.MCMC(kernel, num_samples=n_steps, warmup_steps=warmup, disable_progbar=trace==0)
        mcmc.run(self.dosages, self.allele_freqs)
        return mcmc


def simulateData (true_popmix, nloci):
    # Generate some simulated data
    torch.manual_seed(314159)

    true_popmix = true_popmix / torch.sum(true_popmix)
    #print('True population mix:', true_popmix)
    npops = len(true_popmix)
    #print('Num populations:', npops)
    pop_afs = torch.rand([nloci, npops])
    
    indiv_afs = torch.matmul(pop_afs, true_popmix)
    indiv_dosage = 1*(torch.rand(nloci) < indiv_afs) + 1*(torch.rand(nloci) < indiv_afs)

    return [pop_afs, indiv_dosage]

def profile():
    print('pops\tloci\trep\tmethod\ttime')
    num_reps = 10
    for num_pops in [5, 10, 20]:
        for num_loci in [1000, 10000, 100000]:
            for rep in range(num_reps):
                # Generate some simulated data
                true_popmix = torch.rand([num_pops])
                (pop_afs, indiv_dosage) = simulateData(true_popmix, num_loci)

                model = PopulationMixtureModelRandom(indiv_dosage, pop_afs)

                pyro.clear_param_store()
                time_start = time.time()
                result = model.est_mle(trace=0)
                time_compute = time.time() - time_start
                print('{}\t{}\t{}\tMLE\t{}'.format(num_pops, num_loci, rep+1, time_compute))

                pyro.clear_param_store()
                time_start = time.time()
                result = model.est_svi(trace=0)
                time_compute = time.time() - time_start
                print('{}\t{}\t{}\tSVI\t{}'.format(num_pops, num_loci, rep+1, time_compute))

                pyro.clear_param_store()
                time_start = time.time()
                result = model.est_mcmc(trace=0)
                time_compute = time.time() - time_start
                print('{}\t{}\t{}\tMCMC\t{}'.format(num_pops, num_loci, rep+1, time_compute))


def testSimulatedMissingPopulation():
    npops = 5
    true_popmix = torch.tensor([0.5, 0.25, 0, 0, 0.25])
    nloci = 100000
    (pop_afs, ind_dosage) = simulateData(true_popmix, nloci)

    # Test missing population -- we deliberately leave out a minor population 
    # (25% of the individual) to test performance in this case.
    model = PopulationMixtureModelRandom(ind_dosage, pop_afs[:,:(npops-1)])

    print('\nMaximum likelihood estimation...')
    result_mle = model.est_mle(trace=100)
    print('MLE estimates:\n', result_mle)

    print('\nStochastic variational inference...')
    result_svi = model.est_svi(trace=1000)
    print('SVI estimates:\n', result_svi)
    print('SVI mean:', result_svi / torch.sum(result_svi))
    # Torch is missing a Beta iCDF so just approximate by sampling
    print('SVI PCI:\n', torch.quantile(torch.distributions.Dirichlet(result_svi).sample(torch.tensor([10000])), torch.tensor([0.05, 0.95]), dim=0))

    print('\nMarkov chain Monte Carlo...')
    result_mcmc = model.est_mcmc()
    print('MCMC estimate:')
    print(f"{Fore.GREEN}{result_mcmc.summary}{Style.RESET_ALL}")

def predictPopMemberships(ind_dosage, pop_afs, pop_order, name):
    model = PopulationMixtureModelRandom(ind_dosage, pop_afs)

    print('\nMaximum likelihood estimation...')
    result_mle = model.est_mle(trace=100)
    print('MLE estimates:\n', torch.round(result_mle, decimals=3))
    result_all = pd.DataFrame(torch.round(result_mle, decimals=3), index=pop_order, columns=["MLE"])

    print('\nStochastic variational inference...')
    result_svi = model.est_svi(trace=1000)
    print('SVI estimates:\n', result_svi)
    print(f"SVI mean: {result_svi / torch.sum(result_svi)}")
    # Torch is missing a Beta iCDF so just approximate by sampling
    print('SVI PCI:\n', torch.quantile(torch.distributions.Dirichlet(result_svi).sample(torch.tensor([10000])), torch.tensor([0.05, 0.95]), dim=0))
    result_all["SVI (mean)"] = result_svi / torch.sum(result_svi)


    print('\nMarkov chain Monte Carlo...')
    result_mcmc = model.est_mcmc()
    print('MCMC estimate:')
    result_mcmc.summary()
    print("Summary table:")
    print(f"{Fore.GREEN}{result_all}{Style.RESET_ALL}")
    result_all.to_csv(f"{name}_estimates.csv", index=False)

def main(args):
    af_data = pd.read_table(args.allele_freqs)
    loci_list = af_data["VAR_ID"]
    pop_names = af_data.columns[1:]
    af_floats = af_data.drop(columns="VAR_ID").astype('float32')
    af_tensor = torch.tensor(af_floats.values)
    g = Genotypes(args.vcffile, loci_list)
    print("Populations:", end=" ")
    i = 0
    for pop in pop_names:
        print (f"{Fore.LIGHTGREEN_EX}{i}: {pop}{Style.RESET_ALL}", end="  ")
        i += 1
    print()
    for sample in g.getSamples():
        print (f"{Fore.BLUE}---------- {sample} ----------{Style.RESET_ALL}")
        predictPopMemberships(g.dosageForSample(sample), af_tensor, pop_names, sample)
        print()
    

if __name__ == '__main__':
    import time
    logging.basicConfig(level=logging.WARN)
    colorama.init()

    # profile()
    # testSimulatedMissingPopulation()

    parser = argparse.ArgumentParser(description='Predict population memberships given a VCF with input samples and known population allele frequencies. Note: \
                                     any variants present in AF file not present in the VCF will be imputed for all samples as homozygous reference.')
    parser.add_argument('vcffile',
                        help='Input VCF')
    parser.add_argument('-a', '--allele_freqs', required=True,
                        help="Tab-delimited file containing allele frequencies for each population with one variant per line, and the header line VAR_ID    POP1   [POP2 ...]")

    args = parser.parse_args()
    main(args)
    