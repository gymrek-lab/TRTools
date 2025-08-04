"""
Functions for SISTR simulations
"""

import msprime
import numpy as np
from scipy.stats import geom
import stdpopsim

def GetMuPrime(baseline_mu, baseline_mu_allele, L, target_allele, min_mu, max_mu):
    log_mu_prime = np.log10(baseline_mu)+L*(target_allele-baseline_mu_allele)
    mu_prime = 10**log_mu_prime
    if mu_prime < min_mu:
        return min_mu
    if mu_prime > max_mu:
        return max_mu    
    return mu_prime

def GetStepSizeProb(a1, a2, beta, p):
    step_size = (a2-a1)
    up_prob = max([0.01, 0.5*(1-beta*p*a1)])
    up_prob = min(up_prob, 0.99)
    down_prob = 1-up_prob
    if step_size>0: dir_prob = up_prob
    else: dir_prob = down_prob
    step_prob = geom.pmf(abs(step_size), p)
    return dir_prob*step_prob

def GetEffectivePopSizes(n_effective, max_iter, demog_model, popid):
    """
    Get list of effective pop sizes to use for each generation
    starting at t=0...t=max_iter
    """
    # If n_effective is set, use constant value
    if n_effective is not None:
        return [n_effective]*max_iter
    # Otherwise, get list of n_effs to use based on demographic model
    times = np.linspace(0, max_iter, num=max_iter)
    target_popindex = [item.name for item in demog_model.populations].index(popid)
    debugger = demog_model.model.debug()
    popsizes = debugger.population_size_trajectory(np.linspace(0, max_iter, num=max_iter))[:,target_popindex]
    # Update based on if lineage was present or not
    pll = debugger.possible_lineage_locations()
    for epoch, values in pll.items():
        # If the lineage wasn't present
        # Figure out time of the epoch and
        # which population we split from
        if not values[target_popindex]:
            start_ind = int(epoch[0])
            end_ind = epoch[1]
            if np.isinf(end_ind):
                end_ind = len(popsizes)
            else: end_ind = int(end_ind)
            relevant_events = [item for item in demog_model.model.events \
                if int(item.time)==start_ind and \
                type(item)==msprime.demography.MassMigration and \
                item.source==target_popindex
            ]
            if len(relevant_events) != 1:
                common.WARNING("Error: problem parsing stdpopsim model")
                return None
            # Replace with trajectory of the ancestral pop at that time
            # TODO!!! what if there were additional splits?
            # Need to further refine and check we only include pops at
            # times they were active
            # the check above should catch this for now
            anc_popind = relevant_events[0].dest
            anc_popsizes = debugger.population_size_trajectory(times)[:,anc_popind]
            popsizes[start_ind:end_ind] = anc_popsizes[start_ind:end_ind+1]    
    popsizes = popsizes[::-1] # since we are doing forward simulations
    return popsizes

def GetTransitionMatrix(num_alleles, mu, beta, p, L, min_mu, max_mu):
    # Initialize matrix (optimal=0)
    transition_matrix = np.zeros((num_alleles, num_alleles))

    # Fill in probability to transition from a1 to a2
    for i in range(num_alleles):
        for j in range(num_alleles):
            a1 = -1*int(num_alleles/2)+i
            a2 = -1*int(num_alleles/2)+j
            mu_prime = GetMuPrime(mu, 0, L, a1, min_mu, max_mu)
            if a1==a2: transition_matrix[i,j] = 1-mu_prime
            else:
                prob = GetStepSizeProb(a1, a2, beta, p)
                transition_matrix[i,j] = mu_prime*prob
    # Rescale each row to sum to 1 
    for i in range(num_alleles):
        rowsum = np.sum(transition_matrix[i,:])
        transition_matrix[i,:] = transition_matrix[i,:]/rowsum
    return transition_matrix

def GetMarginalFitnessVector(allele_freqs, fitness_matrix):
    return np.matmul(fitness_matrix, allele_freqs)

def GetGradient(marginal_fitness_vector):
    return 2*marginal_fitness_vector

def GetMeanFitness(marginal_fitness_vector, allele_freqs):
    return np.dot(marginal_fitness_vector, allele_freqs)

def GetCovarianceMatrix(allele_freqs):
    
    allele_freqs_matrix = np.reshape(allele_freqs, (1,len(allele_freqs)))
    allele_freqs_matrix_trans = np.reshape(allele_freqs, (len(allele_freqs),1))
    
    # Fill in off-diagonal elements C[i,j] = -p[i]p[j]/2 (divison by 2 not included here)
    covariance_matrix = np.matmul(allele_freqs_matrix_trans, -1*allele_freqs_matrix)
    
    # Fill in diagonal elements C[i,i] = p[i](1-p[i])/2 (division by 2 not included here)
    for i in range(0, len(allele_freqs)):
        covariance_matrix.itemset((i,i),allele_freqs[i]*(1-allele_freqs[i]))
        
    return covariance_matrix

def GetFitnessMatrix(num_alleles, s):
    fitness_matrix = np.zeros((num_alleles, num_alleles))
    for i in range (0, num_alleles):
        for j in range (0, num_alleles):
            a1 = -1*int(num_alleles/2)+i
            a2 = -1*int(num_alleles/2)+j
           
            # Get fitness of each allele
            w_a1 = 1-abs(a1)*s #np.exp(-1*abs(a1)*s) Chance allele won't die, higher w = higher fitness
            w_a2 = 1-abs(a2)*s #np.exp(-1*abs(a2)*s) Chance allele won't die, higher w = higher fitness
            
            # Genotype fitness cannot be less than 0
            fitness_matrix[i,j] = max(0, w_a1/2 + w_a2/2)
    return fitness_matrix

def RunSimulation(
        sval,
        transition_matrix_transpose,
        set_start_equal=False,
        max_iter=1,
        n_effective=10000,
        demog_model=None,
        popid=None,
        use_drift=True,
        end_samp_n=0,
        ancestral_allele=None
    ):

    # Set the starting vector of allele frequencies
    num_alleles = transition_matrix_transpose.shape[0]
    if set_start_equal:
        allele_freqs = np.full((num_alleles), 1.0/num_alleles)
    else:
        if ancestral_allele is None:
            ancestral_allele = 0
        allele_freqs = np.zeros(num_alleles)
        allele_freqs[int(num_alleles/2)+ancestral_allele] = 1

    # Compute fitness matrix
    fitness_matrix = GetFitnessMatrix(num_alleles, sval)

    # Simulate the desired number of generations
    effective_pop_sizes_by_generation = GetEffectivePopSizes(n_effective, max_iter, demog_model, popid)
    if effective_pop_sizes_by_generation is None:
        return None
    t = 0
    while t < max_iter:
        # Determine N_e based on demographic model
        N_e = effective_pop_sizes_by_generation[t]

        # Calculate marginal fitness w*(a[i]) for each allele
        marginal_fitness_vector = GetMarginalFitnessVector(allele_freqs, fitness_matrix)

        # Calculate gradient vector of partial derivatives of mean fitness
        gradient = GetGradient(marginal_fitness_vector)

        # Calculate mean fitness
        mean_fitness = GetMeanFitness(marginal_fitness_vector, allele_freqs)

        # Calculate covariance matrix
        covariance_matrix = GetCovarianceMatrix(allele_freqs)

        # Calculate new allele_freqs
        # Applying selection
        allele_freqs = allele_freqs + (np.matmul(covariance_matrix, gradient))/(2*mean_fitness)

        # Applying mutation
        allele_freqs = np.matmul(transition_matrix_transpose, allele_freqs)

        if use_drift == True:
            # Use multinomial sampling
            allele_counts = np.random.multinomial(2*N_e, allele_freqs)

            # Rescale allele_freqs to sum to 1
            rowsum = np.sum(allele_counts)
            allele_freqs = allele_counts/rowsum
        t += 1

    # End sampling step
    # Use multinomial sampling on smaller sample size
    if end_samp_n > 0:
        allele_counts = np.random.multinomial(end_samp_n, allele_freqs)
        # Rescale allele_freqs to sum to 1
        rowsum = np.sum(allele_counts)
        allele_freqs = allele_counts/rowsum

    # Set up and return results dictionary
    res = {}
    res["afreqs"] = allele_freqs
    res["afreqs_string"] = ",".join([str(item) for item in allele_freqs])
    return res