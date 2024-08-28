"""
Functions for SISTR simulations
"""

import numpy as np
from scipy.stats import geom

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

def GetEffectivePopSize(N_e, t, max_iter, model=None):
    return N_e # TODO see https://github.com/BonnieCSE/SISTR/blob/master/simulations/Simulation_functions.py#L254

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
        use_drift=True,
        end_samp_n=0
    ):

    # Set the starting vector of allele frequencies
    num_alleles = transition_matrix_transpose.shape[0]
    if set_start_equal:
        allele_freqs = np.full((num_alleles), 1.0/num_alleles)
    else:
        allele_freqs = np.zeros(num_alleles)
        allele_freqs[int(num_alleles/2)] = 1

    # Compute fitness matrix
    fitness_matrix = GetFitnessMatrix(num_alleles, sval)

    # Simulate the desired number of generations
    t = 0
    N_e = n_effective
    while t < max_iter:
        # Determine N_e based on demographic model
        N_e = GetEffectivePopSize(N_e, t, max_iter, model=None)

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