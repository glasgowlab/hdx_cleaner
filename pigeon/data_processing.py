import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from utils import compile_exchange_info, fit_functions


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

def process_table(table):
    """Process a single table"""
    bigdf = pd.read_csv(table)
    header = bigdf.iloc[:, 0:8].copy()
    header.columns = header.iloc[0]
    header.drop(0, inplace=True)

    new_df = pd.DataFrame()
    tp_frames = bigdf.columns[8::8]

    for i, tp_frame in enumerate(tp_frames):
        new_df = process_tp_frame(header, bigdf, tp_frame, i, new_df)

    new_df['State'] = new_df['State'].str.upper().str[:3]

    return new_df


def process_tp_frame(header, bigdf, tp_frame, i, new_df):
    """Process a single tp_frame"""
    tpdict = bigdf.iloc[:, 8+8*i:16+8*i].copy()
    tpdict.columns = tpdict.iloc[0]
    tpdict.drop(0, inplace=True)
    tpdict = pd.concat([header, tpdict], axis=1)
    tpdict['Deut Time (sec)'] = float(tp_frame.split('s')[0])
    new_df = pd.concat([new_df, tpdict])

    return new_df

def process_cleaned_data(df):
    df = df.groupby(['Sequence', 'Deut Time (sec)', 'State', 'Start', 'End'])['#D'].agg(['mean', 'std', 'count'])
    df.reset_index(inplace=True)
    df.rename(columns={'mean': '#D', 'std': 'Stddev', 'count': '#Rep', 'State': 'Protein State'}, inplace=True)

    if 0 not in df['Deut Time (sec)'].unique():
        tp0s = df.groupby(['Sequence', 'Protein State', 'Start', 'End']).sum()
        tp0s[['Deut Time (sec)', '#D', 'Stddev']] = 0
        tp0s['#Rep'] = 1
        df = pd.concat([df, tp0s.reset_index()]).sort_values(by=['Start', 'End', 'Deut Time (sec)', 'Protein State'])

    return df

def load_ranges_file(ranges_file, newbigdf, exclude=False):
    """Load the ranges file"""
    rangeslist = RangeList(ranges_file).to_dataframe()

    def exclude_ranges(newbigdf, rangeslist):
        exc = pd.merge(newbigdf, rangeslist, how='left', indicator=True)
        cleaned = exc[exc['_merge'] == 'left_only'].drop('_merge', axis=1)
        return cleaned

    def include_ranges(newbigdf, rangeslist):
        cleaned = pd.merge(rangeslist, newbigdf)
        return cleaned

    if exclude:
        cleaned = exclude_ranges(newbigdf, rangeslist)
        print('rangeslist excluded !')
    elif rangeslist is not None:
        cleaned = include_ranges(newbigdf, rangeslist)
        print('rangeslist included !')
    else:
        cleaned = newbigdf.drop_duplicates()
        print('no rangeslist !')

    return cleaned


class HDXMSDataCollection:
    def __init__(self, hdxms_data_list):
        self.hdxms_data_list = hdxms_data_list


class HDXMSData:
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.states = []

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

        print(f"Peptide reindexed with offset {-1*index_offset}")


class ProteinState:
    def __init__(self, state_name):
        self.peptides = []
        self.state_name = state_name

    def add_peptide(self, peptide):
        # Check if peptide already exists
        for existing_peptide in self.peptides:
            if existing_peptide.sequence == peptide.sequence:
                raise ValueError("Peptide already exists")
        self.peptides.append(peptide)
        return peptide

    @property
    def num_peptides(self):
        return len(self.peptides)

    def get_peptide(self, sequence):
        for peptide in self.peptides:
            if peptide.sequence == sequence:
                return peptide
        return None

class Peptide:
    def __init__(self, sequence, start, end, protein_state=None):
        self.sequence = sequence
        self.start = start
        self.end = end
        self.timepoints = []

        if protein_state is not None:
            self.protein_state = protein_state

    def add_timepoint(self, timepoint):
        # Check if timepoint already exists
        for existing_timepoint in self.timepoints:
            if existing_timepoint.deut_time == timepoint.deut_time:
                raise ValueError("Timepoint already exists")
        
        self.timepoints.append(timepoint)
        return  timepoint
    
    @property
    def num_timepoints(self):
        return len(self.timepoints)

    @property
    def max_d(self):
        num_prolines = self.sequence[2:].count('P')
        max_d = len(self.sequence) - 2 - num_prolines
        return max_d
    
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
    
    def get_timepoint(self, deut_time):
        for timepoint in self.timepoints:
            if timepoint.deut_time == deut_time:
                return timepoint
        return None

def exchange_fit(x, a, b, c, d, e, f, g, max_d):
    return max_d - a * np.exp(-d * x) - b * np.exp(-e * x) - c * np.exp(-f * x) - g

def exchange_fit_low(x, b, c, e, f, g, max_d):
    return max_d - b * np.exp(-e * x) - c * np.exp(-f * x) - g



class Timepoint:
    def __init__(self, peptide, deut_time, num_d, stddev):
        self.peptide = peptide
        self.deut_time = deut_time
        self.num_d = num_d
        self.stddev = stddev
        self.d_percent = num_d / peptide.max_d


def load_data_to_hdxmsdata(df, protein_name="LacI"):
    ''' 
    Load data from dataframe to HDXMSData object

    example usage:
    hdxms_data = load_data_to_hdxmsdata(cleaned)

    cleaned: dataframe containing cleaned data
    
    '''
    hdxms_data = HDXMSData(protein_name)
    
    # Iterate over rows in the dataframe
    for _, row in df.iterrows():
        # Check if protein state exists
        protein_state = None
        for state in hdxms_data.states:
            if state.state_name == row['Protein State']:
                protein_state = state
                break
        
        # If protein state does not exist, create and add to HDXMSData
        if not protein_state:
            protein_state = ProteinState(row['Protein State'])
            hdxms_data.add_state(protein_state)
        
        # Check if peptide exists in current state
        peptide = None
        for pep in protein_state.peptides:
            if pep.sequence == row['Sequence']:
                peptide = pep
                break
        
        # If peptide does not exist, create and add to ProteinState
        if not peptide:
            # skip if peptide is less than 4 residues
            if len(row['Sequence']) < 4:
                continue
            peptide = Peptide(row['Sequence'], row['Start'], row['End'], protein_state)
            protein_state.add_peptide(peptide)
        
        # Add timepoint data to peptide
        timepoint = Timepoint(peptide, row['Deut Time (sec)'], row['#D'], row['Stddev'])
        
        # if timepoint has no data, skip
        if np.isnan(timepoint.num_d):
            continue
        else:
            peptide.add_timepoint(timepoint)
    
    return hdxms_data


class HDXStatePeptideCompares:
    def __init__(self, state1_list, state2_list):
        self.state1_list = state1_list
        self.state2_list = state2_list
        self.peptide_compares = []
    
    @property
    def common_sequences(self):

        peptides1 = set([peptide.sequence for state1 in self.state1_list for peptide in state1.peptides])
        peptides2 = set([peptide.sequence for state2 in self.state2_list for peptide in state2.peptides])
        common_sequences = peptides1.intersection(peptides2)                
        return common_sequences

    def add_all_compare(self):
        import re
        peptide_compares = []
        for sequence in self.common_sequences:
            peptide1_list = [state1.get_peptide(sequence) for state1 in self.state1_list]
            peptide2_list = [state2.get_peptide(sequence) for state2 in self.state2_list]
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
                if self.resid > pep.start and self.resid < pep.end:
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


