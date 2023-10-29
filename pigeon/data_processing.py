import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import itertools
from utils import compile_exchange_info, fit_functions
import os
from glob import glob
#from plot_functions import ResidueCoverage


def read_hdx_tables(tables, ranges, exclude=False):
    newbigdf = pd.DataFrame()
    cleaned_list = []
    
    # Process all tables
    for table, range_file in zip(tables, ranges):
        newbigdf = process_table(table)

        # Convert columns to the appropriate data types
        newbigdf['Start'] = newbigdf['Start'].apply(np.int64)
        newbigdf['End'] = newbigdf['End'].apply(np.int64)
        newbigdf['#D'] = newbigdf['#D'].apply(float)

        cleaned = load_ranges_file(range_file, newbigdf, exclude)
        cleaned_list.append(cleaned)

    cleaned = pd.concat(cleaned_list, ignore_index=True)
    cleaned = clean_data_df(cleaned)
    cleaned.reset_index(inplace=True,drop=True)
    return cleaned


# Function to find chunks between "Start RT" and "Conf"
def find_chunks(series):
    chunks = []
    start_index = None
    
    for i, label in enumerate(series):
        if label == 'Start RT':
            start_index = i  # Store the index when "Start RT" is found
        elif label == 'Conf' and start_index is not None:
            chunks.append((start_index, i))  # Store the tuple of start and end indices
            start_index = None  # Reset the start index
    
    return chunks


def process_table(table):
    """Process a single table"""
    bigdf = pd.read_csv(table)
    header = bigdf.iloc[:, 0:8].copy()
    header.columns = header.iloc[0]
    header.drop(0, inplace=True)

    new_df = pd.DataFrame()
    chunks = find_chunks(bigdf.iloc[0])

    for chunk in chunks:
        new_df = process_tp_frame(header, bigdf, chunk, new_df)

    new_df['State'] = new_df['State'].str.upper().str[:3]

    return new_df


def process_tp_frame(header, bigdf, tp_chunk, new_df):
    """
    Process a single tp_frame, tp_chunk is the index of the tp_frame
    """
    tp_df = bigdf.iloc[:, tp_chunk[0]:tp_chunk[1]].copy()
    tp_df.columns = tp_df.iloc[0]
    tp_time = bigdf.columns[tp_chunk[0]]
    tp_df.drop(0, inplace=True)
    tp_df = pd.concat([header, tp_df], axis=1)
    if tp_time.startswith('Full-D'):
        tp_df['Deut Time (sec)'] = np.Inf
    else:
        tp_df['Deut Time (sec)'] = float(tp_time.split('s')[0])
    
    new_df = pd.concat([new_df, tp_df])

    return new_df


def clean_data_df(df):
    df = df.groupby(['Sequence', 'Deut Time (sec)', 'State', 'Start', 'End', 'Charge'])['#D'].agg(['mean', 'std', 'count'])
    df.reset_index(inplace=True)
    df.rename(columns={'mean': '#D', 'std': 'Stddev', 'count': '#Rep', 'State': 'Protein State'}, inplace=True)

    if 0 not in df['Deut Time (sec)'].unique():
        tp0s = df.groupby(['Sequence', 'Protein State', 'Start', 'End', 'Charge']).sum()
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
    def __init__(self, protein_name, n_fastamides=2):
        self.protein_name = protein_name
        self.states = []
        self.n_fastamides = n_fastamides

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
                peptide.identifier = f"{peptide.start}-{peptide.end} {peptide.sequence}"

        print(f"Peptide reindexed with offset {-1*index_offset}")

    def to_dataframe(self):
        return revert_hdxmsdata_to_dataframe(self)
    
    def to_bayesianhdx_format(self, OUTPATH=None):
        convert_dataframe_to_bayesianhdx_format(self.to_dataframe(), self.protein_name, OUTPATH)

    def _drop_peptides(self, drop_list):
        for state in self.states:
            for peptide in drop_list:
                state.peptides.remove(peptide)
          

class ProteinState:
    def __init__(self, state_name):
        self.peptides = []
        self.state_name = state_name
        self.if_subtracted = False
        self.num_subtracted_added = 0

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
        for subgroup in subgroups.values():
            if len(subgroup) >= 2:
                # create all possible pairs of items
                pairs = list(itertools.combinations([i for i in subgroup], 2))
                
                for pair in pairs:
                    new_peptide = subtract_peptides(pair[0],pair[1])

                    #skip if new_peptide has less than 3 timepoints
                    if new_peptide.num_timepoints <= 3:
                        continue
                    #skip if new_peptide is single Proline
                    if new_peptide.sequence == 'P':
                        continue  

                    #skip if new_peptide is negative
                    if np.average([tp.num_d for tp in new_peptide.timepoints]) < -0.5:
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

    def  add_all_subtract(self):
        if self.if_subtracted:
            print(f"{self.num_subtracted_added} subtracted peptides have already been subtracted.")

        while True:
            new_peptides = self.add_new_peptides_by_subtract()
            if len(new_peptides) == 0:
                break
        
        self.if_subtracted = True


class Peptide:
    def __init__(self, sequence, start, end, protein_state=None, n_fastamides=2):
        self.identifier = f"{start}-{end} {sequence}" # raw sequence without any modification
        self.sequence = sequence[n_fastamides:]
        self.start = start + n_fastamides
        self.end = end
        self.timepoints = []
        self.note = None

        if protein_state is not None:
            self.protein_state = protein_state

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
            num_prolines = self.sequence[2:].count('P')
            max_d = len(self.sequence) - num_prolines # the peptide already skipped the first two residues
            max_d = len(self.sequence)
        else:
            max_d = inf_tp.num_d

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


    def new_fit(self):
        from sklearn.metrics import mean_squared_error

        x = np.array([tp.deut_time for tp in self.timepoints])
        y = np.array([tp.num_d for tp in self.timepoints])
        max_timepoint = max([tp.deut_time for tp in self.timepoints])
        trialT = np.logspace(1.5, np.log10(max_timepoint*2), 100)

        fit_resluts = {}
        for exp_num in range(1, 5):



            n = exp_num
            
            try:
                popt, pcov = curve_fit(fit_func(n=n), x, y, p0=[0.01]*n + [1]*n, bounds=(0, [np.inf]*n + [self.max_d]*n), maxfev=1000)
                y_pred = fit_func(n=n)(trialT, *popt)
                perr = np.sqrt(np.diag(pcov))

                mse = mean_squared_error(y, fit_func(n=n)(x, *popt))
                loss = mse + np.sqrt(np.sum(perr)) * 0.1 

                fit_resluts[exp_num] = {'popt': popt, 'pcov': pcov, 'perr': perr, 'mse': mse, 'loss': loss, 'y_pred': y_pred, 'trialT': trialT}
            except Exception as e:
                print(f"Error in fitting peptide: exp_num={exp_num}")
                fit_resluts[exp_num] = {'loss': np.inf}

        best_fit = min(fit_resluts, key=lambda x: fit_resluts[x]['loss'])
        best_model = fit_resluts[best_fit]
        return best_model['trialT'], best_model['y_pred'], best_model['popt'], best_model['perr']


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
                avg_timepoint = Timepoint(self, deut_time, np.average([tp.num_d for tp in timepoints]), np.std([tp.num_d for tp in timepoints]))
                return avg_timepoint
            else:
                return None
        else:
            timepoints = [tp for tp in self.timepoints if tp.deut_time == deut_time and tp.charge_state == str(charge_state)]
            if len(timepoints) == 1:
                return timepoints[0]
            else:
                return None

def exchange_fit(x, a, b, c, d, e, f, g, max_d):
    return max_d - a * np.exp(-d * x) - b * np.exp(-e * x) - c * np.exp(-f * x) - g

def exchange_fit_low(x, b, c, e, f, g, max_d):
    return max_d - b * np.exp(-e * x) - c * np.exp(-f * x) - g



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
        self.raw_ms = df\
        
    @property
    def d_percent(self):
        return round(self.num_d / self.peptide.max_d * 100, 2)
    


def load_dataframe_to_hdxmsdata(df, protein_name="Test", n_fastamides=2):
    ''' 
    Load data from dataframe to HDXMSData object

    example usage:
    hdxms_data = load_data_to_hdxmsdata(cleaned)

    cleaned: dataframe containing cleaned data
    
    '''
    hdxms_data = HDXMSData(protein_name, n_fastamides)
    
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
            #identifier = f"{row['Start']+n_fastamides}-{row['End']} {row['Sequence'][n_fastamides:]}"
            identifier = f"{row['Start']}-{row['End']} {row['Sequence']}"
            if pep.identifier == identifier:
                peptide = pep
                break
        
        # If peptide does not exist, create and add to ProteinState
        if not peptide:
            # skip if peptide is less than 4 residues
            if len(row['Sequence']) < 4:
                continue
            peptide = Peptide(row['Sequence'], row['Start'], row['End'], protein_state, n_fastamides=n_fastamides) 
            #peptide = Peptide(row['Sequence'], row['Start'], row['End'], protein_state) 
            protein_state.add_peptide(peptide)
        
        # Add timepoint data to peptide
        timepoint = Timepoint(peptide, row['Deut Time (sec)'], row['#D'], row['Stddev'],row['Charge'])
        
        # if timepoint has no data, skip
        if np.isnan(timepoint.num_d):
            continue
        else:
            peptide.add_timepoint(timepoint)
    
    return hdxms_data


import pandas as pd

def revert_hdxmsdata_to_dataframe(hdxms_data):
    '''
    Convert HDXMSData object to DataFrame
    '''
    # List to hold data
    data_list = []
    
    # Iterate over states in HDXMSData
    for state in hdxms_data.states:
        # Iterate over peptides in ProteinState
        for pep in state.peptides:
            # Iterate over timepoints in Peptide
            for timepoint in pep.timepoints:
                # Dictionary to hold data for this timepoint
                data_dict = {
                    'Protein State': state.state_name,
                    'Sequence': pep.sequence,
                    'Start': pep.start,
                    'End': pep.end,
                    'Deut Time (sec)': timepoint.deut_time,
                    '#D': timepoint.num_d,
                    'Stddev': timepoint.stddev,
                    'Charge': timepoint.charge_state
                }
                data_list.append(data_dict)
    
    # Create DataFrame from data_list
    df = pd.DataFrame(data_list)
    return df


def convert_dataframe_to_bayesianhdx_format(data_df, protein_name='protein', OUTPATH=None):
    if OUTPATH is None:
        OUTPATH = './'
        print('No output path specified, using current directory')
    
    data_df_renamed = data_df.rename(columns={'Sequence':'peptide_seq', 'Start':'start_res', 'End':'end_res', 'Deut Time (sec)':'time', '#D':'D_inc',
                                              'Charge':'charge_state'})

    # filter out peptides with negative start or end residue (due to the His tag)
    data_df_renamed = data_df_renamed[data_df_renamed['start_res'] > 0]
    data_df_renamed = data_df_renamed[data_df_renamed['end_res'] > 0]

    # filter out the single proline
    data_df_renamed = data_df_renamed[data_df_renamed['peptide_seq'] != 'P']
    
    # Save data for each state
    states = set(data_df_renamed['Protein State'])
    for state in states:
        state_df = data_df_renamed[data_df_renamed['Protein State'] == state]
        state_df.to_csv(f'{OUTPATH}/bayesian_hdx_{protein_name}_'+state+'.dat', index=False)

    print(f'Data saved to {OUTPATH}')


        
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
    if not isinstance(protein_state, ProteinState):
        raise TypeError("The input should be a protein state of hdxms_data object.")
    
    # create a dictionary to store the subgroups
    start_subgroups = {} # subgroups with the same start position
    end_subgroups = {} # subgroups with the same end position

    # iterate over all the peptides
    for peptide in protein_state.peptides:
        
        # check if the start position is already in the dictionary
        if peptide.start in start_subgroups:
            start_subgroups[peptide.start].append(f"{peptide.start}-{peptide.end} {peptide.sequence}")
        else:
            start_subgroups[peptide.start] = [f"{peptide.start}-{peptide.end} {peptide.sequence}"]
        
        # check if the end position is already in the dictionary
        if peptide.end in end_subgroups:
            end_subgroups[peptide.end].append(f"{peptide.start}-{peptide.end} {peptide.sequence}")
        else:
            end_subgroups[peptide.end] = [f"{peptide.start}-{peptide.end} {peptide.sequence}"]
        
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
    new_peptide = Peptide(sequence=new_sequence, start=start, end=end, protein_state=peptide_1.protein_state, n_fastamides=0)

    # iterate over all the timepoints
    for tp in common_timepoints:
        std = np.sqrt(longer_peptide.get_timepoint(tp).stddev**2 + shorter_peptide.get_timepoint(tp).stddev**2)
        timepoints = Timepoint(new_peptide, tp, longer_peptide.get_deut(tp) - shorter_peptide.get_deut(tp), std)
        new_peptide.add_timepoint(timepoints)
    new_peptide.note = f"Subtracted from {longer_peptide.identifier} to {shorter_peptide.identifier}"
    return new_peptide


def load_raw_ms_to_hdxms_data(hdxms_data, raw_spectra_path):
    '''
    Load raw MS data from csv files to hdxms_data object. 
    !!! use it before reindex_peptide_from_pdb
    '''
    for state in hdxms_data.states:

        state_raw_spectra_path = os.path.join(raw_spectra_path, state.state_name)
        
        # glob all the folders
        path_dict = {}
        folders = sorted(glob(state_raw_spectra_path + '/*'))

        for folder in folders:
            start, end, seq = folder.split('/')[-1].split('-')
            #start, end, seq = int(start)+2, int(end), seq[2:]     # skip first two res
            start, end, seq = int(start), int(end), seq
            pep_idf = f'{start}-{end} {seq}'
            #pep_idf = f'{start+hdxms_data.n_fastamides}-{end}'
            path_dict[pep_idf] = folder

        
        # iterate through all peptides
        for peptide in state.peptides:
            pep_sub_folder =  path_dict[peptide.identifier]


            for tp in peptide.timepoints:
                if tp.deut_time == 0:
                    csv_name = f'Non-D-1-z{tp.charge_state}.csv'
                else:
                    csv_name = f'{int(tp.deut_time)}s-1-z{tp.charge_state}.csv'
                
                csv_file_path = os.path.join(pep_sub_folder, csv_name)
                df = tp.load_raw_ms_csv(csv_file_path)
                #print(csv_file_path)

    print('Done loading raw MS data.')
     

class HDXStatePeptideCompares:
    def __init__(self, state1_list, state2_list):
        self.state1_list = state1_list
        self.state2_list = state2_list
        self.peptide_compares = []
    
    @property
    def common_sequences(self):
        
        #indetifer f"{pep.start}-{pep.end} {pep.sequence}"
        peptides1 = set([f"{pep.start}-{pep.end} {pep.sequence}" for state1 in self.state1_list for pep in state1.peptides])
        peptides2 = set([f"{pep.start}-{pep.end} {pep.sequence}" for state2 in self.state2_list for pep in state2.peptides])
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


