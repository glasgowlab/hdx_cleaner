from data import HDXMSData, ProteinState, Peptide, Timepoint, RangeList
from spectra import refine_large_error_reps
import tools
import numpy as np
import pandas as pd
import os
from glob import glob
import time

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


def load_dataframe_to_hdxmsdata(df, protein_name="Test", n_fastamides=2, protein_sequence=None, fulld_approx=False, saturation=1.0):
    ''' 
    Load data from dataframe to HDXMSData object

    example usage:
    hdxms_data = load_data_to_hdxmsdata(cleaned)

    cleaned: dataframe containing cleaned data
    
    '''
    hdxms_data = HDXMSData(protein_name, n_fastamides, protein_sequence=protein_sequence, saturation=saturation)
    
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
            protein_state = ProteinState(row['Protein State'], hdxms_data=hdxms_data)
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
        timepoint = Timepoint(peptide, row['Deut Time (sec)'], row['#D'], row['Stddev'],int(row['Charge']))
        
        # if timepoint has no data, skip
        if np.isnan(timepoint.num_d):
            continue
        else:
            if timepoint.deut_time == np.inf and fulld_approx:
                # add inf timepoint as a real timepoint, 1e8s
                real_inf_timepoint = Timepoint(peptide, 100000000, row['#D'], row['Stddev'],int(row['Charge']))
                peptide.add_timepoint(real_inf_timepoint)
            peptide.add_timepoint(timepoint)
    
    return hdxms_data


import pandas as pd

def revert_hdxmsdata_to_dataframe(hdxms_data, if_percent=False):
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
                t_inf_same_charge = pep.get_timepoint(np.inf, timepoint.charge_state)
                # Dictionary to hold data for this timepoint

                #raw idf with fastamides
                if pep.n_fastamides == 2:
                    start, end, = pep.identifier.split(' ')[0].split('-')
                    start, end = int(start), int(end)
                    seq = pep.identifier.split(' ')[1]
                elif pep.n_fastamides == 0:
                    start = pep.start - 2
                    end = pep.end
                    seq = hdxms_data.protein_sequence[start-1:end] # -1 to account for python indexing

                data_dict = {
                    'Protein State': state.state_name,
                    'Sequence': seq,
                    'Start': start,
                    'End': end,
                    'Deut Time (sec)': timepoint.deut_time,
                    '#D': timepoint.num_d,
                    'Stddev': timepoint.stddev,
                    'Charge': str(timepoint.charge_state),
                    'Max #D': pep.max_d,
                }
                if if_percent:
                    data_dict['#D'] = timepoint.d_percent

                if t_inf_same_charge is not None:
                    data_dict['Max #D'] = t_inf_same_charge.num_d

                data_list.append(data_dict)
    
    # Create DataFrame from data_list
    df = pd.DataFrame(data_list)
    print('Reminder: sequence contains fastamides !!!')
    return df


def convert_dataframe_to_bayesianhdx_format(data_df, protein_name='protein', OUTPATH=None):
    if OUTPATH is None:
        OUTPATH = './'
        print('No output path specified, using current directory')

    # filter out the inf timepoint
    data_df = data_df[data_df['Deut Time (sec)'] != np.inf]
    
    data_df_renamed = data_df.rename(columns={'Sequence':'peptide_seq', 'Start':'start_res', 'End':'end_res', 'Deut Time (sec)':'time', '#D':'D_inc',
                                              'Charge':'charge_state', 'Max #D':'max_D'})

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
                if tp.deut_time == np.inf:
                    tp.isotope_envelope = None
                    continue
                elif tp.deut_time == 0:
                    csv_name = f'Non-D-1-z{tp.charge_state}.csv'
                    csv_file_path = os.path.join(pep_sub_folder, csv_name)
                elif tp.deut_time == 100000000:
                    csv_file_path = glob(f'{pep_sub_folder}/Full-D-*-z{tp.charge_state}.csv')[0]

                else:
                    csv_name = f'{int(tp.deut_time)}s-1-z{tp.charge_state}.csv'
                    csv_file_path = os.path.join(pep_sub_folder, csv_name)
                try:
                    df = tp.load_raw_ms_csv(csv_file_path)
                except:
                    print(peptide.identifier, tp.deut_time, tp.charge_state)
                    print(csv_file_path)
                
        bad_timepoints = [tp for peptide in state.peptides for tp in peptide.timepoints if tp.isotope_envelope is None and tp.deut_time != np.inf]
        high_back_ex_tps = [tp for pep in state.peptides for tp in pep.timepoints if pep.max_d/pep.theo_max_d < 0.5]
        tools.remove_tps_from_state(bad_timepoints, state)
        tools.remove_tps_from_state(high_back_ex_tps, state)

    print('Done loading raw MS data.')
     
    for state in hdxms_data.states:
        refine_large_error_reps(state)


def export_iso_files(hdxms_data, outdir, overwrite=True):
    import shutil
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        if overwrite:
            shutil.rmtree(outdir)
            os.mkdir(outdir)

    all_tps = [tp for state in hdxms_data.states for peptide in state.peptides for tp in peptide.timepoints if tp.deut_time != np.inf]
    
    for tp in all_tps:
        state = tp.peptide.protein_state.state_name
        
        #raw idf with fastamides
        if tp.peptide.n_fastamides == 2:
            start, end, = tp.peptide.identifier.split(' ')[0].split('-')
            seq = tp.peptide.identifier.split(' ')[1]
        elif tp.peptide.n_fastamides == 0:
            start = tp.peptide.start - 2
            end = tp.peptide.end
            seq = hdxms_data.protein_sequence[start-1:end] # -1 to account for python indexing

        idf = f'{start}-{end}-{seq}'
        npy_file_name = f'{state}_{idf}_tp{int(tp.deut_time)}_ch{tp.charge_state}.npy'
        np.save(os.path.join(outdir, npy_file_name), tp.isotope_envelope)
    print(f'Isotope files saved to {outdir}')
    print('Reminder: sequence contains fastamides !!!')