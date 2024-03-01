import data
import numpy as np
import pandas as pd
from spectra import get_theo_ms


        
def refine_data(hdxms_data_list, std_threshold=1.0):
    '''
    Remove peptides with large deviation across replicates
    all tps has std > 1.0 will be removed
    '''
    all_tps = [tp for data in hdxms_data_list for state in data.states for peptide in state.peptides for tp in peptide.timepoints if tp.deut_time != np.inf and tp.deut_time != 0.0]
    grouped_tps = group_by_attributes(all_tps, ['peptide.protein_state.state_name', 'peptide.identifier','deut_time'])
    
    high_std_tps = {}
    all_std = []
    for k, v in grouped_tps.items():
        if len(v) > 1:
            std = np.std([i.num_d for i in v])
            if std > std_threshold:
                high_std_tps[k] = v
            all_std.append(std)
    for k,v in high_std_tps.items():
        std = np.std([i.num_d for i in v])
        #print(k, std)     

    all_std = np.array(all_std)  
    print(f'Average timepoint std: {np.mean(all_std)}')

    for k,v in high_std_tps.items():
        std = np.std([i.num_d for i in v])
        states = [data.get_state(k[0]) for data in hdxms_data_list]
        for state in states:
            remove_tps_from_state(v, state)
        #print(k, std)    
    


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


    # return None if:
    #skip if new_peptide has less than 3 timepoints
    if len(common_timepoints) <= 3:
        return None
    #skip if new_peptide is single Proline
    if new_peptide.sequence == 'P':
        return None  
    if new_peptide.max_d == 0:
        return None 

    # iterate over all the timepoints
    for tp in common_timepoints:
        
        longer_tp = longer_peptide.get_timepoint(tp)
        shorter_tp = shorter_peptide.get_timepoint(tp)

        std = np.sqrt(longer_peptide.get_timepoint(tp).stddev**2 + shorter_peptide.get_timepoint(tp).stddev**2)
        timepoint = data.Timepoint(new_peptide, tp, longer_tp.num_d - shorter_tp.num_d, std)
        #timepoint = data.Timepoint(new_peptide, tp, longer_peptide.get_deut(tp) - shorter_peptide.get_deut(tp), std)

        if hasattr(longer_tp, 'isotope_envelope') and hasattr(shorter_tp, 'isotope_envelope') and tp != np.inf:

            #timepoint.isotope_envelope = deconvolute(longer_peptide.get_timepoint(tp).isotope_envelope, shorter_peptide.get_timepoint(tp).isotope_envelope)
            p1 = longer_tp.isotope_envelope.copy()
            p2 = shorter_tp.isotope_envelope.copy()
            #p3 = get_theo_ms(timepoint)
            p3 = force_to_envelope(deconvolute_fft(p1, p2))
            best_model = deconvolute_iso(p1, p2, p3, steps=2000)
            if best_model[1] is None:
                continue 
            else:
                timepoint.isotope_envelope = best_model[1]
                if abs(get_centroid_iso(timepoint) - timepoint.num_d) > 0.5:
                    continue
                new_peptide.add_timepoint(timepoint)

        else:
            new_peptide.add_timepoint(timepoint)
    
    if len(new_peptide.timepoints) < 3:
        return None        
    
    #skip if new_peptide is negative
    if np.average([tp.num_d for tp in new_peptide.timepoints]) < -0.3 or new_peptide.max_d < 0.0:
        return None  

    new_peptide.note = f"Subtraction: {longer_peptide.identifier} - {shorter_peptide.identifier}"

    # if 0.0 not in [tp.deut_time for tp in new_peptide.timepoints]:
    #     print(longer_peptide.identifier, shorter_peptide.identifier)
    #     raise ValueError("Cannot subtract peptides without time 0.")

    return new_peptide


def average_timepoints(tps_list):

    # check if the same state, peptide, and timepoint
    if len(set([(tp.peptide.protein_state.state_name, tp.peptide.identifier, tp.deut_time) for tp in tps_list])) != 1:
        raise ValueError('Timepoints must be from the same state, peptide, and timepoint')

    deut_time = tps_list[0].deut_time
    num_d = np.average([tp.num_d for tp in tps_list])
    std_d = np.std([tp.num_d for tp in tps_list])
    
    avg_timepoint = data.Timepoint(tps_list[0].peptide, deut_time, num_d, std_d)
    try:
        
        size = max([tp.isotope_envelope.size for tp in tps_list])
        iso_list = [custom_pad(tp.isotope_envelope, size, pad_value=0.0) for tp in tps_list]
        avg_iso = np.average(iso_list, axis=0)
        avg_timepoint.isotope_envelope = avg_iso/np.sum(avg_iso)
                
    except:
        #print('No isotope envelope found for timepoint')
        pass
    
    return avg_timepoint


def deconvolute_fft(p, p1):
    # Compute the Fourier Transforms
    P = np.fft.fft(p)
    P1 = np.fft.fft(p1, n=len(p))  # Make sure the arrays have the same length
    
    # Division in the Frequency Domain
    P2 = P / P1
    
    # Compute the Inverse Fourier Transform
    p2 = np.fft.ifft(P2).real  # Take the real part to ignore small numerical errors

    #wipe the negative values to 1e-10
    p2[p2 < 0] = 1e-10
    p2 = p2 / np.sum(p2)

    return p2



def metropolis_criterion(old_loss, new_loss, temperature):
    if new_loss < old_loss:
        return True
    else:
        delta_loss = new_loss - old_loss
        probability = np.exp(-delta_loss / temperature)
        return np.random.rand() < probability
    

def force_to_envelope(p3_array):
    
    max_index = np.argmax(p3_array)

    before_peak = np.sort(p3_array[0:max_index])
    after_peak = np.sort(p3_array[max_index+1:])[::-1]
    
    #make a new envelope
    new_envelope = np.concatenate((before_peak, np.array(p3_array[max_index]),after_peak), axis=None)
    
    return new_envelope


def deconvolute_iso(p1, p2, p3, temperature=0.05, steps=2000, keep_best=True):
    '''
    Deconvolution function that updates p3 all at once.
    p3: initial guess
    '''

    max_length = max(p1.size, p2.size)
    p1 = custom_pad(p1, max_length)
    p2 = custom_pad(p2, max_length)
    p3 = p3[:max_length] / np.sum(p3[:max_length])
    
    p1_estimated = convolve(p3, p2)
    p1_estimated = normlize(p1_estimated)

    num_peaks = np.count_nonzero(p1>0.1)  #number of major peaks
    previous_loss = get_sum_ae(p1, p1_estimated)/num_peaks
    
    best_models = []
    
    continus_rejects = 0

    for i in range(steps):
        for j in range(len(p3)):
            original_p3_j = p3[j]
            change = np.random.normal(0, 0.1)  # Generate change for the current peak
            p3[j] += change
            p3[j] = max(p3[j], 1e-10)  # Prevent negative values
            
            p3 = normlize(p3)  # Normalize p3 after all individual changes
            p3 = force_to_envelope(p3)  # Force to envelope

            # Calculate new loss 
            p1_estimated = convolve(p3, p2)
            p1_estimated = normlize(p1_estimated)
            new_loss = get_sum_ae(p1, p1_estimated) / num_peaks

            # Decide whether to accept the change
            if metropolis_criterion(previous_loss, new_loss, temperature):
                previous_loss = new_loss
                continus_rejects = 0
            else:
                # Revert the change for the current peak
                p3[j] = original_p3_j
                continus_rejects += 1

        best_models.append((previous_loss, p3.copy()))
        temperature *= 0.99  # Decrease the temperature

        if new_loss < 0.1 or (0.1 < new_loss < 0.2 and continus_rejects > 20):
            best_models.sort(key=lambda x: x[0])
            return best_models[0]
        elif new_loss > 0.2 and continus_rejects > 500:
            return (previous_loss, None)

        #Optional: Uncomment to see progress every 100 steps
        # if i % 100 == 0:
        #     print(i, previous_loss)

    

def get_centroid_iso(tp):
    mz = np.arange(0, len(tp.isotope_envelope))
    intensity = tp.isotope_envelope
    return mz @ intensity


def get_t0_centroid_iso(tp):
    t0_centroid_iso = [get_centroid_iso(tp) for tp in tp.peptide.timepoints if tp.deut_time == 0]
    return np.average(t0_centroid_iso)
    
    
def get_num_d_from_iso(tp):
    return (get_centroid_iso(tp) - get_t0_centroid_iso(tp))/1.00866491588



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
def normlize(p):
    return p / np.sum(p)

@njit
def convolve(p1, p2):
    return np.convolve(p1, p2)


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
        pep.timepoints = [tp for tp in pep.timepoints if tp not in removing_tps]
        #print(f'{pep.identifier} {tp.deut_time} {tp.charge_state} removed')

        # Remove peptides without timepoints or with no time 0
        if pep.timepoints == []:
            state.peptides.remove(pep)
            print(f'{pep.identifier} removed')

    for pep in state.peptides:
        tp0_tps = [tp for tp in pep.timepoints if tp.deut_time == 0 ]
        if len(tp0_tps) == 0:
            state.peptides.remove(pep)
            print(f'{pep.identifier} removed')

    for pep in state.peptides:   
        # with less than 3 timepoints not counting time 0, inf
        real_tps_num = len(set(tp.deut_time for tp in pep.timepoints if tp.deut_time != np.inf and tp.deut_time != 0.0))
        if real_tps_num < 3:
            state.peptides.remove(pep)
            print(f'{pep.identifier} removed')
    
    #n_tp = len([tp for pep in state.peptides for tp in pep.timepoints])
    #print(f'Number of timepoints after removing: {n_tp}')
    


from collections import defaultdict

def group_by_attributes(objects, attributes):
    grouped_data = defaultdict(list)
    for obj in objects:
        key = []
        for attr in attributes:
            if '.' in attr:
                nested_attrs = attr.split('.')
                value = obj
                for na in nested_attrs:
                    value = getattr(value, na)
                key.append(value)
            else:
                key.append(getattr(obj, attr))
        grouped_data[tuple(key)].append(obj)
    return grouped_data



def backexchange_for_peps_no_data(hdxms_data_list):
    '''
    Calculate the average backexchange for all peptides in the list of hdxms_data objects and correct the peptides has no data
    '''
    for hdxms_data in hdxms_data_list:

        all_peps = [pep for state in hdxms_data.states for pep in state.peptides]
        pep_with_exp_max_d = [pep for pep in all_peps if pep.get_timepoint(np.inf) is not None]
        pep_with_no_exp_mad_d = [pep for pep in all_peps if not pep in pep_with_exp_max_d]
        avg_backexchange = np.mean([pep.max_d/pep.theo_max_d for pep in pep_with_exp_max_d])
        
        for pep in pep_with_no_exp_mad_d:
            max_d = pep.theo_max_d * avg_backexchange
            inf_tp = data.Timepoint(pep, np.inf, max_d, np.nan)
            pep.add_timepoint(inf_tp)
        
        print('Number of peptides with no data: {}'.format(len(pep_with_no_exp_mad_d)))
        print('Average backexchange for peptides with no data: {}'.format(avg_backexchange))



def calculate_coverages(hdxms_datas, state_name):
    states = [data.get_state(state_name) for data in hdxms_datas]
    coverage = np.zeros(len(states[0].hdxms_data.protein_sequence))
    for state in states:
        for pep in state.peptides:
            coverage[pep.start-1:pep.end] += 1 
    return coverage



def generate_bayesian_hdx_script(install_dir,exp_names, protein_sequence, protein_state, base_directory, making_chunks=False, if_original_bayesian_hdx=False):
    """
    Generates a script for the Bayesian HDX analysis based on the provided parameters and a template file.

    :param raw_spectra_paths: List of paths to raw spectra.
    :param protein_sequence: The protein sequence.
    :param protein_state: The state of the protein.
    :param base_directory: The base directory for the data.
    :param template_file: Path to the template file.
    :return: A string containing the generated script.
    """
    # Read the template from the file
    if making_chunks:
        with open('../../lib/run_bayesian_hdx_template_chunks.txt', 'r') as file:
            script_template = file.read()
    else:
        if if_original_bayesian_hdx:
            with open('../../lib/run_original_bayesian_hdx_template.txt', 'r') as file:
                script_template = file.read()
        else:
            with open('../../lib/run_bayesian_hdx_template.txt', 'r') as file:
                script_template = file.read()

    # Generate the script using the template
    generated_script = script_template.format(
        install_dir=install_dir,
        protein_state=protein_state,
        sequence=protein_sequence,
        exp_names=exp_names,
        base_directory=base_directory
    )

    return generated_script

