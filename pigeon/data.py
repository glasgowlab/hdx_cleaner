import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import itertools
from tools import find_overlapped_peptides, subtract_peptides, exchange_fit, exchange_fit_low, fit_func, average_timepoints
import spectra
from hdxrate import k_int_from_sequence
import random
import math
import re
import warnings


class RangeList:
    def __init__(self, range_list_file=None, range_df=None):

        if range_list_file is not None:
            df = pd.read_csv(range_list_file)
            if len(df.columns) != 2:
                df = pd.read_csv(range_list_file,  skiprows=1)
                self.range_list = df[['Start', 'End']].drop_duplicates().sort_values(by='Start').reset_index(drop=True)
            else:   
                self.range_list = df[["Start", "End"]]
        else:
            self.range_list = range_df

    def to_set(self):
        return set(tuple(row) for _, row in self.range_list.iterrows())
    
    def to_dataframe(self):
        return self.range_list
    
    def to_csv(self, path):
        self.range_list.to_csv(path, index=False)
    
    def union(self, other):
        if isinstance(other, RangeList):
            other_set = other.to_set()
        elif isinstance(other, set):
            other_set = other
        else:
            print("Invalid input. Please provide a RangeList object or a set.")
            return None
        
        result_set = self.to_set().union(other_set)
        result_df = pd.DataFrame(list(result_set), columns=["Start", "End"])
        result_df = result_df.astype(int)
        result_df = result_df.sort_values(by="Start").reset_index(drop=True)
        return RangeList(range_df=result_df)
    
    def intersection(self, other):
        if isinstance(other, RangeList):
            other_set = other.to_set()
        elif isinstance(other, set):
            other_set = other
        else:
            print("Invalid input. Please provide a RangeList object or a set.")
            return None
        
        result_set = self.to_set().intersection(other_set)
        result_df = pd.DataFrame(list(result_set), columns=["Start", "End"])
        result_df = result_df.astype(int)
        result_df = result_df.sort_values(by="Start").reset_index(drop=True)
        return RangeList(range_df=result_df)
    
    def difference(self, other):
        if isinstance(other, RangeList):
            other_set = other.to_set()
        elif isinstance(other, set):
            other_set = other
        else:
            print("Invalid input. Please provide a RangeList object or a set.")
            return None
        
        result_set = self.to_set().difference(other_set)
        result_df = pd.DataFrame(list(result_set), columns=["Start", "End"])
        result_df = result_df.astype(int)
        result_df = result_df.sort_values(by="Start").reset_index(drop=True)
        return RangeList(range_df=result_df)


class HDXMSDataCollection:
    def __init__(self, hdxms_data_list):
        self.hdxms_data_list = hdxms_data_list


class HDXMSData:
    def __init__(self, protein_name, n_fastamides=2, protein_sequence=None):
        self.protein_name = protein_name
        self.states = []
        self.n_fastamides = n_fastamides
        if protein_sequence is None:
            raise ValueError("Protein sequence is required")
        self.protein_sequence = protein_sequence

    def add_state(self, state):
        # Check if state already exists
        for existing_state in self.states:
            if existing_state.state_name == state.state_name:
                raise ValueError("State already exists")
        self.states.append(state)
        return state
    
    def load_protein_sequence(self, sequence):
        self.protein_sequence = sequence

    @property
    def num_states(self):
        return len(self.states)
    
    def get_state(self, state_name):
        for state in self.states:
            if state.state_name == state_name:
                return state
        return None
     
    def plot_res_coverage(self):
        res_cov = ResidueCoverage(self)
        return res_cov.plot()


    def reindex_peptide_from_pdb(self, pdb_file, first_residue_index=1):

        def pdb2seq(pdb_file):

            import warnings
            from Bio import SeqIO
            from Bio import BiopythonWarning

            # Suppress all Biopython-specific warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', BiopythonWarning)
                records = list(SeqIO.parse(pdb_file, 'pdb-atom'))
                return str(records[0].seq)
    
        def find_peptide(seq, peptide):
            start_index = seq.find(peptide)
            if start_index == -1:
                return (-1, -1)
            end_index = start_index + len(peptide) - 1
            return (start_index, end_index)

        pdb_sequence = pdb2seq(pdb_file)
        a_middle_pep = self.states[0].peptides[int(len(self.states[0].peptides)/2)]
        pdb_start, pdb_end= find_peptide(pdb_sequence, a_middle_pep.sequence)
        index_offset = a_middle_pep.start - pdb_start - first_residue_index

        for state in self.states:
            for peptide in state.peptides:
                peptide.start -= index_offset
                peptide.end -= index_offset
                
                # identifier
                old_idf = peptide.identifier
                idf_start, idf_end = re.match(r'(-?\d+)-(-?\d+)', old_idf).group(1,2)
                idf_seq = old_idf.split(' ')[1]
                new_idf = f"{int(idf_start) - index_offset}-{int(idf_end) - index_offset} {idf_seq}"
                peptide.identifier = new_idf

        print(f"Peptide reindexed with offset {-1*index_offset}")

    def to_dataframe(self, if_percent=False):
        return revert_hdxmsdata_to_dataframe(self, if_percent=if_percent)
    
    def to_bayesianhdx_format(self, OUTPATH=None):
        convert_dataframe_to_bayesianhdx_format(self.to_dataframe(), self.protein_name, OUTPATH)

    def _drop_peptides(self, drop_list):
        for state in self.states:
            for peptide in drop_list:
                state.peptides.remove(peptide)
          

class ProteinState:
    def __init__(self, state_name, hdxms_data=None):
        self.peptides = []
        self.state_name = state_name
        self.if_subtracted = False
        self.num_subtracted_added = 0
        self.hdxms_data = hdxms_data

    def add_peptide(self, peptide):
        # Check if peptide already exists
        for existing_peptide in self.peptides:
            if existing_peptide.identifier == peptide.identifier:
                raise ValueError(f"{ peptide.identifier} Peptide already exists")
        self.peptides.append(peptide)
        return peptide

    @property
    def num_peptides(self):
        return len(self.peptides)

    def get_peptide(self, identifier):
        for peptide in self.peptides:
            if peptide.identifier == identifier:
                return peptide
        return None
    
    def add_new_peptides_by_subtract(self):
        '''
        add new peptides to the protein state by subtracting the overlapped peptides
        '''

        # check the input
        if not isinstance(self, ProteinState):
            raise TypeError("The input should be a protein state of hdxms_data object.")

        subgroups = find_overlapped_peptides(self)

        new_peptide_added = []
        for key in sorted(subgroups.keys()):
            subgroup = sorted(subgroups[key], key=lambda x: x.start)
            if len(subgroup) >= 2:
                # create all possible pairs of items
                pairs = list(itertools.combinations([i for i in subgroup], 2))
                
                for pair in pairs:
                    #print(pair[0].identifier, pair[1].identifier)

                    # skip if the new peptide already exists
                    note1 = f"Subtraction: {pair[0].identifier} - {pair[1].identifier}"
                    note2 = f"Subtraction: {pair[1].identifier} - {pair[0].identifier}"
                    notes = [pep.note for pep in self.peptides if pep.note is not None]
                    if note1 in notes or note2 in notes:
                        continue

                    try:
                        new_peptide = subtract_peptides(pair[0],pair[1])
                    except:
                        continue

                    # skil if a None peptide is returned
                    if new_peptide is None:
                        continue
                    
                    #skip if new_peptide has high recaculated num_d error
                    if hasattr(new_peptide.timepoints[0], 'isotope_envelope'):
                        from tools import get_num_d_from_iso
                        avg_error = np.abs(np.average([tp.num_d - get_num_d_from_iso(tp) for tp in new_peptide.timepoints if tp.deut_time != np.inf]))
                        if avg_error > 0.15:
                            continue
                    
                    # add the new peptide to the protein state object
                    try:
                        self.add_peptide(new_peptide)
                        new_peptide_added.append(new_peptide)
                    except:
                        pass
                        #print(f"Peptide {new_peptide.sequence} already exists in the protein state.")

        print(f"{len(new_peptide_added)} new peptides added to the protein state.")   

        self.subtracted_peptides = new_peptide_added
        self.num_subtracted_added += len(new_peptide_added)
        return new_peptide_added

    def add_all_subtract(self):
        if self.if_subtracted:
            print(f"{self.num_subtracted_added} subtracted peptides have already been subtracted.")

        while True:
            new_peptides = self.add_new_peptides_by_subtract()
            if len(new_peptides) == 0:
                break
        
        self.if_subtracted = True


class Peptide:
    def __init__(self, sequence, start, end, protein_state=None, n_fastamides=0):
        self.identifier = f"{start}-{end} {sequence}" # raw sequence without any modification
        self.sequence = sequence[n_fastamides:]
        self.start = start + n_fastamides
        self.end = end
        self.timepoints = []
        self.note = None
        self.n_fastamides = n_fastamides

        if protein_state is not None:
            self.protein_state = protein_state
        
        #self.max_d = self.get_max_d()

    def add_timepoint(self, timepoint):
        # Check if timepoint already exists
        for existing_timepoint in self.timepoints:
            if existing_timepoint.deut_time == timepoint.deut_time and existing_timepoint.charge_state == timepoint.charge_state:
                raise ValueError(f"{self.start}-{self.end} {self.sequence}: {timepoint.deut_time} (charge: {timepoint.charge_state})Timepoint already exists")
        
        self.timepoints.append(timepoint)
        return  timepoint
    
    @property
    def num_timepoints(self):
        return len(self.timepoints)

    @property
    def max_d(self):

        # if no inf timepoint, return the theoretical max_d
        inf_tp = self.get_timepoint(np.inf)

        if inf_tp is None:
            max_d = self.theo_max_d
        else:
            max_d = inf_tp.num_d

        return max_d
    
    
    @property
    def theo_max_d(self):
        num_prolines = self.sequence.count('P')
        theo_max_d = len(self.sequence) - num_prolines
        return theo_max_d
    
    @property
    def fit_results(self):

        try:
            max_timepoint = max([tp.deut_time for tp in self.timepoints])
            trialT = np.logspace(1.5, np.log10(max_timepoint*2), 100)
            
            if self.timepoints[-1].num_d > 0.5:
                x = [tp.deut_time for tp in self.timepoints]
                y = [tp.num_d for tp in self.timepoints]
                popt, pcov = curve_fit(f = exchange_fit, xdata = x, ydata = y,
                        #bounds = (0, [self.max_d, self.max_d, self.max_d, 1, .1, .01, self.max_d, self.max_d]),
                        bounds = (0, [self.max_d, self.max_d, self.max_d, 1, .1, .01, self.max_d, self.max_d]),
                        #bounds = (0, [np.inf, np.inf, self.max_d]),
                        maxfev = 100000)
                y_pred = exchange_fit(trialT, *popt)
                perr = np.sqrt(np.diag(pcov))
            else:
                x = [tp.deut_time for tp in self.timepoints]
                y = [tp.num_d for tp in self.timepoints]
                popt, pcov = curve_fit(f = exchange_fit_low, xdata = x, ydata = y,
                        bounds = (0, [self.max_d, self.max_d, .1, .01, self.max_d, self.max_d]),
                        maxfev = 100000)
                y_pred = exchange_fit_low(trialT, *popt)
                perr = np.sqrt(np.diag(pcov))
            return trialT, y_pred, popt, perr
        except:
            raise ValueError(f"Error in fitting peptide: {self.start}-{self.end} {self.sequence}")

    def fit_hdx_stats(self, start = {'a': None, 'b': 0.001, 'd':0, 'p': 1}):

        if start['a'] is None:
            start['a'] = self.max_d
    
        max_timepoint = max([tp.deut_time for tp in self.timepoints])
        trialT = np.logspace(1.5, np.log10(max_timepoint*2), 100)

        x = [tp.deut_time for tp in self.timepoints]
        y = [tp.num_d for tp in self.timepoints]

        def model_func(timepoint, a, b, p, d):
            return a * (1 - np.exp(-b * (timepoint ** p))) + d
        
        try:
            popt, pcov = curve_fit(model_func, x, y, p0=list(start.values()), maxfev=100000, method='lm')
            y_pred = model_func(trialT, *popt)
            perr = np.sqrt(np.diag(pcov))
        except Exception as e:
            print(e)
            print(f"Error in fitting peptide: {self.start}-{self.end} {self.sequence}")
            y_pred = np.zeros(len(trialT))
            popt = np.zeros(4)
            perr = np.zeros(4)

        return trialT, y_pred, popt, perr


    def new_fit(self):
        from sklearn.metrics import mean_squared_error

        x = np.array([tp.deut_time for tp in self.timepoints if tp.deut_time != np.inf])
        y = np.array([tp.num_d for tp in self.timepoints if tp.deut_time != np.inf])
        max_timepoint = max([tp.deut_time for tp in self.timepoints if tp.deut_time != np.inf])
        min_timepoint = min([tp.deut_time for tp in self.timepoints if tp.deut_time != 0])
        trialT = np.logspace(np.log10(min_timepoint), np.log10(max_timepoint), 100)

        fit_resluts = {}
        for exp_num in range(2, 4):

            n = exp_num
            
            try:
                popt, pcov = curve_fit(fit_func(n=n), x, y, p0=[0.01]*n + [1]*n, bounds=(0, [np.inf]*n + [self.max_d]*n), maxfev=1000)
                y_pred = fit_func(n=n)(trialT, *popt)
                perr = np.sqrt(np.diag(pcov))

                mse = mean_squared_error(y, fit_func(n=n)(x, *popt))
                loss = mse + np.sqrt(np.sum(perr)) * 0.1 

                fit_resluts[exp_num] = {'popt': popt, 'pcov': pcov, 'perr': perr, 'mse': mse, 'loss': loss, 'y_pred': y_pred, 'trialT': trialT}
            except Exception as e:
                #print(f"Error in fitting peptide: exp_num={exp_num}")
                fit_resluts[exp_num] = {'loss': np.inf}
                warnings.warn(f"Error in fitting peptide: exp_num={exp_num}")

        best_fit = min(fit_resluts, key=lambda x: fit_resluts[x]['loss'])
        best_model = fit_resluts[best_fit]

        try:
            return best_model['trialT'], best_model['y_pred'], best_model['popt'], best_model['perr']
        except:
            return None, None, None, None


    def get_deut(self, deut_time):
        for timepoint in self.timepoints:
            if timepoint.deut_time == deut_time:
                return timepoint.num_d
        return None
    
    def get_deut_percent(self, deut_time):
        for timepoint in self.timepoints:
            if timepoint.deut_time == deut_time:
                return timepoint.d_percent
        return None
    
    def get_timepoint(self, deut_time, charge_state=None):

        if charge_state is None:
            timepoints = [tp for tp in self.timepoints if tp.deut_time == deut_time]

            #if not empty return average timepoint
            if len(timepoints) == 1:
                return timepoints[0]
            
            elif len(timepoints) > 1:
                #avg_timepoint = Timepoint(self, deut_time, np.average([tp.num_d for tp in timepoints]), np.std([tp.num_d for tp in timepoints]))
                avg_timepoint = average_timepoints(timepoints)
                return avg_timepoint
            else:
                return None
        else:
            timepoints = [tp for tp in self.timepoints if tp.deut_time == deut_time and tp.charge_state == charge_state]
            if len(timepoints) == 1:
                return timepoints[0]
            else:
                return None

class Timepoint:
    def __init__(self, peptide, deut_time, num_d, stddev, charge_state=None):
        self.peptide = peptide
        self.deut_time = deut_time
        self.num_d = num_d
        self.stddev = stddev
        #self.d_percent = num_d / peptide.max_d
        self.charge_state = charge_state

    def load_raw_ms_csv(self, csv_file):
        df = pd.read_csv(csv_file, names=['m/z', 'Intensity'])
        # normalize intensity to sum to 1
        #df['Intensity'] = df['Intensity'] / df['Intensity'].sum()
        self.raw_ms = df
        
        iso = spectra.get_isotope_envelope(self, add_sn_ratio_to_tp=True)
        if iso is not None:
            self.isotope_envelope = iso['Intensity'].values
        else:
            self.isotope_envelope = None

    @property
    def d_percent(self):
        return round(self.num_d / self.peptide.max_d * 100, 2)


class HDXStatePeptideCompares:
    def __init__(self, state1_list, state2_list):
        self.state1_list = state1_list
        self.state2_list = state2_list
        self.peptide_compares = []
    
    @property
    def common_idfs(self):
        
        #indetifer f"{pep.start}-{pep.end} {pep.sequence}"
        peptides1 = set([pep.identifier for state1 in self.state1_list for pep in state1.peptides])
        peptides2 = set([pep.identifier for state2 in self.state2_list for pep in state2.peptides])
        common_sequences = peptides1.intersection(peptides2)                
        return common_sequences

    def add_all_compare(self):
        import re
        peptide_compares = []
        for idf in self.common_idfs:
            peptide1_list = [state1.get_peptide(idf) for state1 in self.state1_list]
            peptide2_list = [state2.get_peptide(idf) for state2 in self.state2_list]
            peptide_compare = PeptideCompare(peptide1_list, peptide2_list)
            peptide_compares.append(peptide_compare)
        self.peptide_compares = sorted(peptide_compares, key=lambda x: int(re.search(r"(-?\d+)--?\d+ \w+", x.compare_info).group(1)))
    
    def to_dataframe(self):
        df = pd.DataFrame()
        for peptide_compare in self.peptide_compares:
            peptide_compare_df = pd.DataFrame([{'Sequence': peptide_compare.compare_info.split(': ')[1],'deut_diff_avg': peptide_compare.deut_diff_avg}])
            df = pd.concat([df, peptide_compare_df], ignore_index=True)
        df = df.pivot_table(index='Sequence', values='deut_diff_avg')
        df.index = pd.CategoricalIndex(df.index, categories=[i.compare_info.split(': ')[1] for i in self.peptide_compares])
        df = df.sort_index()
        return df

class PeptideCompare:
    def __init__(self, peptide1_list, peptide2_list):
        self.peptide1_list = [peptide for peptide in peptide1_list if peptide is not None]
        self.peptide2_list = [peptide for peptide in peptide2_list if peptide is not None]
        
        #if not same sequence, raise error
        set_1 = set([peptide.sequence for peptide in self.peptide1_list])
        set_2 = set([peptide.sequence for peptide in self.peptide1_list])
        if set_1 != set_2:
            raise ValueError('Cannot compare peptides with different sequences')
    
    @property
    def compare_info(self):
        peptide1 = self.peptide1_list[0]
        peptide2 = self.peptide2_list[0]

        return f'{peptide1.protein_state.state_name}-{peptide2.protein_state.state_name}: {peptide1.start}-{peptide1.end} {peptide1.sequence}'
    
    
    @property
    def common_timepoints(self):

        timepoints1 = set([tp.deut_time for peptide1 in self.peptide1_list for tp in peptide1.timepoints if tp is not None])
        timepoints2 = set([tp.deut_time for peptide2 in self.peptide2_list for tp in peptide2.timepoints if tp is not None])
        common_timepoints = list(timepoints1.intersection(timepoints2))
        common_timepoints.sort()

        return np.array(common_timepoints)

    @property
    def deut_diff(self):

        deut_diff = []
        for timepoint in self.common_timepoints:
            deut_diff.append(self.get_deut_diff(timepoint))

        deut_diff = np.array(deut_diff)
        return deut_diff

    @property
    def deut_diff_sum(self):
        return np.sum(self.deut_diff)

    @property
    def deut_diff_avg(self):
        return np.average(self.deut_diff)

    def get_deut_diff(self, timepoint):

        deut1_array = np.array([pep1.get_deut_percent(timepoint) for pep1 in self.peptide1_list if pep1.get_deut_percent(timepoint) is not None])
        deut2_array = np.array([pep2.get_deut_percent(timepoint) for pep2 in self.peptide2_list if pep2.get_deut_percent(timepoint) is not None])
        
        result = deut1_array.mean() - deut2_array.mean()
        return result


class HDXStateResidueCompares:
    def __init__(self, resids, state1_list, state2_list):
        self.resids = resids
        self.residue_compares = []
        self.state1_list = state1_list
        self.state2_list = state2_list
        
    def add_all_compare(self):

        for resid in self.resids:
            res_compare = ResidueCompare(resid, self.state1_list, self.state2_list)
            if not res_compare.if_empty:
                self.residue_compares.append(res_compare)   

    def get_residue_compare(self, resid):
        return self.residue_compares[self.resids.index(resid)]
    

class ResidueCompare:
    def __init__(self, resid, state1_list, state2_list):
        self.resid = resid
        self.state1_list = state1_list
        self.state2_list = state2_list

    
    #@property
    #def resname(self):
    #   a_peptide = self.containing_peptides1[0]
    #    return a_peptide.sequence[a_peptide.start - self.resid]

    @property
    def compare_info(self):
        return f'{self.state1_list[0].state_name}-{self.state2_list[0].state_name}: {self.resid}'
    
    def find_peptides_containing_res(self, state_list):
        res_containing_peptides = []
        for state in state_list:
            for pep in state.peptides:
                if self.resid >= pep.start and self.resid <= pep.end:
                #if self.resid - pep.start < 5 and pep.end - self.resid < 5: 
                    res_containing_peptides.append(pep)
        return res_containing_peptides

    @property
    def containing_peptides1(self):
        return self.find_peptides_containing_res(self.state1_list)
    
    @property
    def containing_peptides2(self):
        return self.find_peptides_containing_res(self.state2_list)


    @property
    def if_empty(self):
        if len(self.containing_peptides1) == 0 or len(self.containing_peptides2) == 0:
            return True
        else:
            return False

    @property
    def common_timepoints(self):
        tp1 = set([tp.deut_time for pep1 in self.containing_peptides1 for tp in pep1.timepoints if tp is not None])
        tp2 = set([tp.deut_time for pep2 in self.containing_peptides2 for tp in pep2.timepoints if tp is not None])
        common_timepoints = list(tp1.intersection(tp2))
        common_timepoints.sort()

        return np.array(common_timepoints)

    @property
    def deut_diff(self):

        if len(self.common_timepoints) == 0:
            return np.nan
        else:
            deut_diff = []
            for timepoint in self.common_timepoints:
                    deut_diff.append(self.get_deut_diff(timepoint))
            deut_diff = np.array(deut_diff)

            return deut_diff

    @property
    def deut_diff_sum(self):
        return np.sum(self.deut_diff)

    @property
    def deut_diff_avg(self):
        return np.average(self.deut_diff)

    def get_deut_diff(self, timepoint):
        
        
        deut1_array = np.array([pep1.get_deut_percent(timepoint) for pep1 in self.containing_peptides1 if pep1.get_deut_percent(timepoint) is not None])
        deut2_array = np.array([pep2.get_deut_percent(timepoint) for pep2 in self.containing_peptides2 if pep2.get_deut_percent(timepoint) is not None])

        result = deut1_array.mean() - deut2_array.mean()
        return result


class SimulatedData:
    def __init__(self, length=100, seed=42, noise_level=0):

        random.seed(seed)
        self.length = length
        self.gen_seq()
        self.gen_logP()
        self.cal_k_init()
        self.timepoints = np.insert(np.logspace(1, 12, 20), 0, 0)
        #self.timepoints = np.insert(self.timepoints, -1, 1e12)
        self.noise_level = noise_level



    def gen_seq(self,):

        AA_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.sequence = ''.join(random.choices(AA_list, k=self.length))

    
    def gen_logP(self,):
        logP = np.array([random.uniform(2., 10.) for i in range(self.length)])
        # Proline residues are not exchangeable
        for i in range(self.length):
            if self.sequence[i] == 'P':
                logP[i] = 0.
        self.logP = logP

    def cal_k_init(self):
        self.log_k_init = np.log10(k_int_from_sequence(self.sequence, 293., 7.))


    def calculate_incorporation(self):
        def calculate_simple_deuterium_incorporation(rate, time):
            # Calculates the deuterium incorporation for a log(kex)
            # at time = t (seconds) assuming full saturation
            return 1 - np.exp(-1*(10**rate)*time)
        incorporations = []
        for res_i in range(self.length):
            log_kex = self.log_k_init[res_i] - self.logP[res_i]
            incorporations.append(calculate_simple_deuterium_incorporation(log_kex, self.timepoints))

        self.incorporations = np.array(incorporations)



    def _gen_peptides(self, min_len=5, max_len=12, max_overlap=5, num_peptides=30):
            """
            Chunks a given sequence into peptides of random length.
            :param sequence: The sequence to be chunked.
            :param min_len: Minimum length of a chunk.
            :param max_len: Maximum length of a chunk.
            :param max_overlap: Maximum overlap between consecutive chunks.
            :param num_peptides: The desired number of peptides.
            :return: A list of chunked peptides.
            """
            chunks = []
            sequence_length = len(self.sequence)
            covered = [False] * sequence_length

            while len(chunks) < num_peptides and not all(covered):
                start = random.randint(0, sequence_length - 1)
                end = min(start + random.randint(min_len, max_len), sequence_length)

                # Check if the current chunk overlaps significantly with already covered parts
                if not all(covered[start:end]):
                    chunks.append(self.sequence[start:end])
                    for i in range(start, end):
                        covered[i] = True

                # Introduce random overlaps
                if max_overlap and end < sequence_length:
                    overlap_end = min(end + random.randint(0, max_overlap), sequence_length)
                    for i in range(end, overlap_end):
                        covered[i] = False

            self.peptides = [i for i in sorted(chunks, key=lambda x: self.sequence.find(x)) if len(i) > 3]
        

    def gen_peptides(self, min_len=5, max_len=12, max_overlap=5, num_peptides=30):

        chunks = []
        avg_peptide_length = math.ceil(min_len + max_len) / 2
        sequence_length = len(self.sequence)    
         
        avg_coverage = math.ceil(num_peptides/sequence_length)

        for i in range(sequence_length):
            count = 0
            while count < avg_coverage:
                start = max(0, random.randint(int(i-avg_peptide_length-2), i)) # n_fastamides = 2
                end = min(random.randint(i, int(i+avg_peptide_length)), sequence_length)
                pep_len = len([i for i in self.sequence[start:end] if i != 'P'])
                if pep_len > 3:
                    chunks.append(self.sequence[start:end])
                    count += 1
                else:
                    continue
        
        covered = [False] * sequence_length
        while not all(covered):
            reduced_chunks = random.sample(chunks, k=num_peptides,)
            
            for pep in reduced_chunks:
                start = self.sequence.find(pep)
                end = start + len(pep)
                for i in range(start, end):
                    covered[i] = True
                    
        print(len(reduced_chunks))

        self.peptides = [i for i in sorted(reduced_chunks, key=lambda x: self.sequence.find(x)) if len(i) > 3]


    def convert_to_hdxms_data(self):
        hdxms_data = HDXMSData('simulated_data', protein_sequence=self.sequence)
        protein_state = ProteinState('SIM', hdxms_data=hdxms_data)
        hdxms_data.add_state(protein_state)

        #calculate incorporation
        self.calculate_incorporation()


        for peptide in self.peptides:
            start = self.sequence.find(peptide) + 1
            end = start + len(peptide) - 1

            peptide_obj = Peptide(peptide, start, end, protein_state, n_fastamides=2) 
            try:
                protein_state.add_peptide(peptide_obj)
                for tp_i, tp in enumerate(self.timepoints):
                    pep_incorp = sum(self.incorporations[peptide_obj.start-1:peptide_obj.end][:,tp_i])
                    #random_stddev = peptide_obj.theo_max_d * self.noise_level * random.uniform(-1, 1)
                    random_stddev = 0
                    tp_obj = Timepoint(peptide_obj, tp, pep_incorp + random_stddev, random_stddev,)
                    
                    isotope_envelope = spectra.get_theo_ms(tp_obj)
                    isotope_noise = np.array([random.uniform(-1, 1)* self.noise_level *peak for peak in isotope_envelope])
                    tp_obj.isotope_envelope = isotope_envelope + isotope_noise
                    #tp_obj.isotope_envelope = isotope_envelope
                
                    peptide_obj.add_timepoint(tp_obj)

            except:
                continue
            
        self.hdxms_data = hdxms_data

    
                #peptide_obj.add_timepoint(tp, self.incorporations[start:end])

