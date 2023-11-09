import data
import numpy as np
import pandas as pd


        
def refine_data(hdxms_data_list):
    '''
    Remove peptides with large deviation across replicates
    std > 80% * max difference between states
    '''
    

def find_overlapped_peptides(protein_state):
    '''
    find overlapped peptides in the hdxms_data object that could be used for substation 

    Parameters
    ----------
    protein_state : hdxms_data.states[0]
        the protein state of interest
    '''

    # check the input
    if not isinstance(protein_state, data.ProteinState):
        raise TypeError("The input should be a protein state of hdxms_data object.")
    
    # create a dictionary to store the subgroups
    start_subgroups = {} # subgroups with the same start position
    end_subgroups = {} # subgroups with the same end position

    # iterate over all the peptides
    for peptide in protein_state.peptides:
        
        # check if the start position is already in the dictionary
        if peptide.start in start_subgroups:
            start_subgroups[peptide.start].append(f"{peptide.identifier}")
        else:
            start_subgroups[peptide.start] = [f"{peptide.identifier}"]
        
        # check if the end position is already in the dictionary
        if peptide.end in end_subgroups:
            end_subgroups[peptide.end].append(f"{peptide.identifier}")
        else:
            end_subgroups[peptide.end] = [f"{peptide.identifier}"]
        
    # combine two subgroups
    combined = {**start_subgroups, **end_subgroups}
    subgroups = {}
    for key in combined.keys():
        value = list(set(combined[key]))
        if len(value) > 1:
            subgroups[key] = [protein_state.get_peptide(idf) for idf in value]

    return subgroups

def subtract_peptides(peptide_1, peptide_2):
    """
    Subtract two peptides and create a new peptide.
    """
    # Check if the two peptides have the same length
    if peptide_1.sequence == peptide_2.sequence:
        raise ValueError("Cannot subtract the same two peptides.")
    
    if peptide_1.start != peptide_2.start and peptide_1.end != peptide_2.end:
        raise ValueError("start or end need to be the different.")
    
    if peptide_1.protein_state.state_name != peptide_2.protein_state.state_name:
        raise  ValueError("Cannot subtract peptides from different states.")
    
    # get longer peptide and shorter peptide
    if len(peptide_1.sequence) > len(peptide_2.sequence):
        longer_peptide = peptide_1
        shorter_peptide = peptide_2
    else:
        longer_peptide = peptide_2
        shorter_peptide = peptide_1

    timepoints1 = set([tp.deut_time for tp in longer_peptide.timepoints])
    timepoints2 = set([tp.deut_time for tp in shorter_peptide.timepoints])
    common_timepoints = list(timepoints1.intersection(timepoints2))
    common_timepoints.sort()

    if shorter_peptide.start == longer_peptide.start:
        start = shorter_peptide.end + 1
        end = longer_peptide.end 
        new_sequence = longer_peptide.sequence[shorter_peptide.end - shorter_peptide.start +1 : longer_peptide.end - longer_peptide.start +1]
    else:
        start = longer_peptide.start
        end = shorter_peptide.start - 1
        new_sequence = longer_peptide.sequence[0: shorter_peptide.start - longer_peptide.start]

    # Create a new peptide (n_fastamides=0)
    new_peptide = data.Peptide(sequence=new_sequence, start=start, end=end, protein_state=peptide_1.protein_state, n_fastamides=0)

    # iterate over all the timepoints
    for tp in common_timepoints:
        std = np.sqrt(longer_peptide.get_timepoint(tp).stddev**2 + shorter_peptide.get_timepoint(tp).stddev**2)
        timepoints = data.Timepoint(new_peptide, tp, longer_peptide.get_deut(tp) - shorter_peptide.get_deut(tp), std)
        new_peptide.add_timepoint(timepoints)
    new_peptide.note = f"Subtracted from {longer_peptide.identifier} to {shorter_peptide.identifier}"
    return new_peptide



def exponential_sum_decorator(n):
    def decorator(func):
        def wrapper(*args):
            if len(args) != 2*n + 1:
                raise ValueError(f"Expected {2*n + 1} arguments, got {len(args)}")
            
            t = args[0]  # assuming the first argument is the time variable 't'
            k_values = args[1:n+1]  # extracting the rate constants
            max_d_values = args[n+1:]  # extracting the max deuteration levels
            
            result = 0
            for k, max_d in zip(k_values, max_d_values):
                result += func(t, k, max_d)  # call the original function with the ith rate constant, time, and max_d
            return result
        return wrapper
    return decorator


def fit_func(n=1):
    @exponential_sum_decorator(n=n)
    def single_amide_ex(t, k, max_d):
        return max_d * (1 - np.exp(-k*t))
    
    return single_amide_ex


def exchange_fit(x, a, b, c, d, e, f, g, max_d):
    return max_d - a * np.exp(-d * x) - b * np.exp(-e * x) - c * np.exp(-f * x) - g

def exchange_fit_low(x, b, c, e, f, g, max_d):
    return max_d - b * np.exp(-e * x) - c * np.exp(-f * x) - g



from numba import njit

@njit
def custom_kl_divergence(p, q):
    """Calculate the Kullback-Leibler divergence between two probability distributions."""
        
    # Ensure the probabilities are normalized
    x1 = p / np.sum(p)
    x2 = q / np.sum(q)

    divergence = np.sum(x1 * np.log(x1 / x2))
    return divergence

@njit
def jensen_shannon_divergence(p, q):
    x1 = p / np.sum(p)
    x2 = q / np.sum(q)
    m = (x1 + x2) / 2
    return (custom_kl_divergence(x1, m) + custom_kl_divergence(x2, m)) / 2

@njit
def custom_pad(array, target_length, pad_value=1e-10):
    #replace 0s with a small value
    array[array == 0] = pad_value
    padded_array = np.full(target_length, pad_value)  
    padded_array[:len(array)] = array 
    return padded_array


@njit
def norm_pad_two_arrays(p, q):
    '''
    normalize two arrays and pad them to the same size and 
    '''
    x1 = p / np.sum(p)
    x2 = q / np.sum(q)

    size = max(len(x1), len(x2))
    x1 = custom_pad(x1, size)
    x2 = custom_pad(x2, size)
    return x1, x2


@njit
def get_divergence(p, q, method='KL'):
    
    x1, x2 = norm_pad_two_arrays(p, q)
    
    if method == 'KL':
        return custom_kl_divergence(x1, x2)
    elif method == 'JS':
        return jensen_shannon_divergence(x1, x2)
    

@njit
def get_mse(p, q):

    x1, x2 = norm_pad_two_arrays(p, q)
    
    return np.mean(np.square(x1 - x2))

@njit
def get_mae(p, q):

    x1, x2 = norm_pad_two_arrays(p, q)
    
    return np.mean(np.abs(x1 - x2))


@njit
def get_sum_ae(p, q):

    x1, x2 = norm_pad_two_arrays(p, q)
    return np.abs(x1 - x2).sum()

@njit
def event_probabilities(p_array):
    n = p_array.shape[0]
    # dp[i][j] stores the probability of j successes in the first i events
    dp = np.zeros((n+1, n+1))
    dp[0][0] = 1  # Base case: probability of 0 successes in 0 events is 1
    
    for i in range(1, n+1):
        p = p_array[i-1]
        q = 1 - p
        dp[i][0] = dp[i-1][0] * q  # Probability of 0 successes in i events
        
        for j in range(1, i+1):
            # Probability of j successes in i events
            dp[i][j] = dp[i-1][j-1] * p + dp[i-1][j] * q
    
    # The probabilities of 0, 1, 2, ..., n successes in n events
    probabilities = dp[n]
    
    return probabilities



def remove_tps_from_state(removing_tps, state):
    
    #n_tp = len([tp for pep in state.peptides for tp in pep.timepoints])
    #print(f'Number of timepoints before removing: {n_tp}')

    # remove the replicates from the state
    for pep in state.peptides:
        for tp in pep.timepoints:
            if tp in removing_tps:
                pep.timepoints.remove(tp)
                #print(f'{pep.identifier} {tp.deut_time} {tp.charge_state} removed')

        # Remove peptides without timepoints or with no time 0
        if pep.timepoints == []:
            state.peptides.remove(pep)
            print(f'{pep.sequence} removed')

    for pep in state.peptides:
        tp0_tps = [tp for tp in pep.timepoints if tp.deut_time == 0 ]
        if len(tp0_tps) == 0:
            state.peptides.remove(pep)
            print(f'{pep.sequence} removed')
    
    #n_tp = len([tp for pep in state.peptides for tp in pep.timepoints])
    #print(f'Number of timepoints after removing: {n_tp}')
    
