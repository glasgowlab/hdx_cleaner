import pandas as pd
import numpy as np
from utils import compile_exchange_info, fit_functions


class RangeList:
    def __init__(self, range_list_path):
        df = pd.read_csv(range_list_path)
        if len(df.columns) != 2:
            self.range_list = df[['Start', 'End']].drop_duplicates().sort_values(by='Start').reset_index(drop=True)
        else:   
            self.range_list = df[["Start", "End"]]
    
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
        return result_df
    
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
        return result_df
    
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
        return result_df


def process_table(table):
    """Process a single table"""
    bigdf = pd.read_csv(table)
    header = bigdf.iloc[:, 0:8].copy()  # create a copy
    header.columns = header.iloc[0]
    header.drop(0, inplace=True)
 

    new_df = pd.DataFrame()
    tp_frames = bigdf.columns[8::8]

    for i, tp_frame in enumerate(tp_frames):
        # pass a copy of header and bigdf to the function
        new_df = process_tp_frame(header.copy(), bigdf.copy(), tp_frame, i, new_df)

    #capitalise the state names and keep only the first 3 letters
    new_df['State'] = new_df['State'].str.upper().str[:3]
    
    return new_df


def process_tp_frame(header, bigdf, tp_frame, i, new_df):
    """Process a single tp_frame"""
    tpdict = bigdf.iloc[:, 8+8*i:16+8*i].copy()
    tpdict.columns = tpdict.iloc[0]
    tpdict.drop(0, inplace=True)
    tpdict = pd.concat([header, tpdict], axis=1)
    tpdict.loc[:, 'Deut Time (sec)'] = float(tp_frame.split('s')[0])
    new_df = pd.concat([new_df, tpdict])
    


    return new_df


def load_ranges_file(ranges_file, newbigdf, exclude=False):
    """Load the ranges file"""

    #rangeslist = pd.read_csv(ranges_file)
    rangeslist = RangeList(ranges_file).to_dataframe()

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



def exclude_ranges(newbigdf, rangeslist):
    """Exclude ranges from the data"""
    exc = pd.merge(newbigdf, rangeslist, how='left', indicator=True)
    cleaned = exc[exc['_merge']=='left_only'].drop('_merge', axis=1)

    return cleaned

def include_ranges(newbigdf, rangeslist):
    """Include only the specified ranges in the data"""
    cleaned = pd.merge(rangeslist, newbigdf)

    return cleaned


def get_unique_sorted(series):
    return sorted(set(series))

# Instead of converting to dict then back to list, directly use set to remove duplicates and sorted to sort
#states = get_unique_sorted(cleaned['State'])
#peptides = get_unique_sorted(cleaned['Sequence'])
#timepoints = get_unique_sorted(cleaned['Deut Time (sec)'])

def process_cleaned_data(df):
    df = df.groupby(['Sequence', 'Deut Time (sec)', 'State', 'Start', 'End'])['#D'].agg(['mean','std', 'count'])
    df.reset_index(inplace=True)
    df.rename(columns={'mean':'#D', 'std': 'Stddev', 'count':'#Rep', 'State':'Protein State'}, inplace=True)

    # Only add 0 timepoint if it doesn't exist
    if 0 not in df['Deut Time (sec)'].unique():
        tp0s = df.groupby(['Sequence', 'Protein State', 'Start', 'End']).sum()
        tp0s[['Deut Time (sec)', '#D', 'Stddev']] = 0
        tp0s['#Rep'] = 1
        df = pd.concat([df, tp0s.reset_index()]).sort_values(by=['Start', 'End', 'Deut Time (sec)', 'Protein State'])

    return df


# Timepoints also needs to include 0 if it wasn't already included
#timepoints = get_unique_sorted(cleaned['Deut Time (sec)'])

def create_sequence_dict(df, states):
    sequence_df = df.loc[(df['Protein State'] == states[0]) & (df['Deut Time (sec)'] == 0)][['Start', 'End', 'Sequence']].drop_duplicates().reset_index(drop=True)
    states_dict = {state: sequence_df.copy() for state in states}
    peptides_2 = create_peptides_list(sequence_df)

    return states_dict, sequence_df['Start'], peptides_2

#states_dict, first_res = create_sequence_dict(cleaned, states)

def create_peptides_list(sequence_df):
    return [f"{sequence_df['Start'][i]}-{sequence_df['End'][i]}- {sequence_df['Sequence'][i]}" for i in range(len(sequence_df))]

#peptides_2 = create_peptides_list(sequence_df)



def load_data(args):
    newbigdf = pd.DataFrame()

    # Process all tables
    for table in args.table:
        newbigdf = pd.concat([newbigdf, process_table(table)])

    # Convert columns to the appropriate data types
    newbigdf['Start'] = newbigdf['Start'].apply(np.int64)
    newbigdf['End'] = newbigdf['End'].apply(np.int64)
    newbigdf['#D'] = newbigdf['#D'].apply(float)

    return newbigdf


def clean_data(args, newbigdf):
    cleaned = load_ranges_file(args.ranges, newbigdf, args.exclude)
    cleaned = process_cleaned_data(cleaned)

    states = list(dict.fromkeys(cleaned['Protein State']))
    peptides = list(dict.fromkeys(cleaned['Sequence']))
    timepoints = list(dict.fromkeys(cleaned['Deut Time (sec)']))

    states_dict, first_res, peptides_2 = create_sequence_dict(cleaned, states)

    return cleaned, states, peptides, timepoints, states_dict, first_res, peptides_2


def fit_data(cleaned, states, states_dict, peptides, timepoints):
    # Compile the exchange information
    peptide_exchange_dict, stdev_dict_dict = compile_exchange_info(
        cleaned, states, states_dict)

    # Fit the exchange functions
    trialT, peptide_fit_dict, peptide_params_dict, peptide_err_dict = fit_functions(
        peptides, peptide_exchange_dict, timepoints)

    return peptide_exchange_dict, stdev_dict_dict, trialT, peptide_fit_dict, peptide_params_dict, peptide_err_dict


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


class Peptide:
    def __init__(self, sequence, start, end,):
        self.sequence = sequence
        self.start = start
        self.end = end
        #self.charge = charge
        #self.search_rt = search_rt
        #self.max_d = max_d

        self.timepoints = []

    def add_timepoint(self, timepoint):
        # Check if timepoint already exists
        for existing_timepoint in self.timepoints:
            if existing_timepoint.deut_time == timepoint.deut_time:
                print("Timepoint already exists")
                return
        self.timepoints.append(timepoint)
        return  timepoint
    
    @property
    def num_timepoints(self):
        return len(self.timepoints)


class Timepoint:
    def __init__(self, deut_time, num_d, stddev):
        self.deut_time = deut_time
        self.num_d = num_d
        self.stddev = stddev


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
            peptide = Peptide(row['Sequence'], row['Start'], row['End'])
            protein_state.add_peptide(peptide)
        
        # Add timepoint data to peptide
        timepoint = Timepoint(row['Deut Time (sec)'], row['#D'], row['Stddev'])
        peptide.add_timepoint(timepoint)
    
    return hdxms_data


