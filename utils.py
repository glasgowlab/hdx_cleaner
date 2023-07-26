import math
from scipy.optimize import curve_fit
import numpy as np

def compile_exchange_info(cleaned, states, states_dict):
    peptide_exchange_dict = {}
    stdev_dict_dict = {}

    for state in states:
        for peptide in states_dict[state]['Sequence']:
            if peptide not in peptide_exchange_dict:
                peptide_exchange_dict[peptide] = {}
            peptide_exchange = (cleaned.loc[(cleaned['Sequence'] == peptide) & 
                                            (cleaned['Protein State'] == state)][['#D']]).values[:]#[:9]
            peptide_exchange_dict[peptide][state] = peptide_exchange

        stdev_dict_dict[state] = {peptide: (cleaned.loc[(cleaned['Sequence'] == peptide) & 
                                                        (cleaned['Protein State'] == state)][['Stddev']]).values
                                  for peptide in states_dict[state]['Sequence']}
    return peptide_exchange_dict, stdev_dict_dict

def fit_functions(peptides, peptide_exchange_dict, timepoints):
    trialT = np.logspace(1.5, np.log10(max(timepoints)*2), 1000)
    peptide_fit_dict = {}
    peptide_params_dict = {}
    peptide_err_dict = {}

    for peptide in peptides:
        num_prolines = peptide[2:].count('P')
        max_protons = len(peptide) - 2 - num_prolines

        peptide_fit_dict[peptide] = []
        peptide_params_dict[peptide] = []
        peptide_err_dict[peptide] = []

        element = peptide_exchange_dict.get(peptide) # element is state
        if element is not None:
            for state in element.keys():
                sub_fit_dict, sub_params_dict, sub_err_dict = fit_single_peptide(peptide, element, state, max_protons, trialT, timepoints)
                peptide_fit_dict[peptide].append(sub_fit_dict)
                peptide_params_dict[peptide].append(sub_params_dict)
                peptide_err_dict[peptide].append(sub_err_dict)
    return trialT, peptide_fit_dict, peptide_params_dict, peptide_err_dict

def fit_single_peptide(peptide, element, state, max_protons, trialT, timepoints):
    def exchange_fit(x, a, b, c, d, e, f, g):
        return max_protons - a * np.exp(-d * x) - b * np.exp(-e * x) - c * np.exp(-f * x) - g

    def exchange_fit_low(x, b, c, e, f, g):
        return max_protons - b * np.exp(-e * x) - c * np.exp(-f * x) - g

    exchange_list = element[state] # these are the avg tp measurements for the state

    sub_fit_dict = {}
    sub_params_dict = {}
    sub_err_dict = {}

    p1_index = 0
    peptide1_tps = []
    peptide1_ex = []
    for tp in exchange_list:
        if not math.isnan(float(tp)):
            peptide1_tps.append(timepoints[p1_index])
            peptide1_ex.append(float(tp))
        p1_index = p1_index + 1
    try:
        if peptide1_ex[-1] > .5:
            popt, pcov = curve_fit(f = exchange_fit, xdata = peptide1_tps, ydata = peptide1_ex,
                                    bounds = (0, [max_protons, max_protons, max_protons, 1, .1, .01, max_protons]),
                                    maxfev = 1000)
            exchange_peptide1 = exchange_fit(trialT, *popt)
            perr = np.sqrt(np.diag(pcov))

        else:
            popt, pcov = curve_fit(f = exchange_fit_low, xdata = peptide1_tps, ydata = peptide1_ex,
                                    bounds = (0, [max_protons, max_protons, .1, .01, max_protons]),
                                    maxfev = 1000)
            exchange_peptide1 = exchange_fit_low(trialT, *popt)
            perr = np.sqrt(np.diag(pcov))
    except:
        #popt, pcov = np.zeros(popt.shape), np.zeros(pcov.shape)
        popt, pcov = np.zeros(7), np.zeros(7)
        exchange_peptide1 = np.zeros(len(trialT))
        perr = np.zeros(7)
    sub_fit_dict[state] = exchange_peptide1
    sub_params_dict[state] = popt
    sub_err_dict[state] = perr
    return sub_fit_dict, sub_params_dict, sub_err_dict


# Create a function to generate pairs
def generate_pairs(states, args):
    pairs = []
    if len(states) > 1:
        for state1 in states:
            for state2 in states:
                if [state2, state1] not in pairs and state1 != state2:
                    pairs.append([state1, state2])
        if (args.s1 is not None) and (args.s2 is not None):
            pairs.append([args.s1, args.s2])
    return pairs