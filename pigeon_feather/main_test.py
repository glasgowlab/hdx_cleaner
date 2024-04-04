from data_processing import *
from plot_functions import *
import argparse
import numpy as np
import pandas as pd


OUTDIR = '/Users/chenlin/Downloads/rbsb_output/'

parser = argparse.ArgumentParser(description='Statistical analysis of HDX/MS data for curated RbsR peptides in one or more functional states.')
parser.add_argument('--pm', dest='pm', help="path to pymol structure", required=True)
parser.add_argument('--t', '--table', dest='table', help="path to uptake table", nargs='+', required=True)
parser.add_argument('--r', '--ranges', dest='ranges', help='path to ranges list csv')
parser.add_argument('--e', '--exclude', dest='exclude', action='store_true', help='exclude rather than include rangeslist')
parser.add_argument('--s1', dest='s1', help='first state to compare')
parser.add_argument('--s2', dest='s2', help='second state to compare')
parser.add_argument('--compare', dest='compare', help='dna, ligand, both', nargs=3)
parser.add_argument('--cbarmax', dest='cbarmax', type=float, help='max value for colorbar axis for dDbar')
parser.add_argument('--ldmin', dest='ldmin', type=float, help='in dDbar, minimum difference threshold between ligand/dna states')



#args = parser.parse_args()
args = parser.parse_args(args=['--t', './EcPFK_20230707_BioRep1/peptide_pool_20230724.csv','./EcPFK_20230424/peptide_pool_20230502.csv'
                               '--r', './EcPFK_20230707_BioRep1/ranges_list_e.csv', 
                               '--pm', '../example/2dri_protein.pdb',])

args.table = ['../example/10-25_rbsb_pool_results-CLEANED.csv']
args.ranges = ['../example/rangeslist-10-25-select-cleaned.csv']



global colorbar_max, delta_LD_threshold
colorbar_max = 0.05 if args.cbarmax is None else args.cbarmax
delta_LD_threshold = 0.075 if args.ldmin is None else args.ldmin



newbigdf = pd.DataFrame()

cleaned_list = []
# Process all tables
for table, range_file in zip(args.table[:2], args.ranges[:2]):
    newbigdf = process_table(table)

    # Convert columns to the appropriate data types
    newbigdf['Start'] = newbigdf['Start'].apply(np.int64)
    newbigdf['End'] = newbigdf['End'].apply(np.int64)
    newbigdf['#D'] = newbigdf['#D'].apply(float)



    cleaned  = load_ranges_file(range_file, newbigdf, args.exclude)
    cleaned_list.append(cleaned)

cleaned = pd.concat(cleaned_list, ignore_index=True)
cleaned = process_cleaned_data(cleaned)


hdxms_data = load_data_to_hdxmsdata(cleaned)

# make a uptake plot for all the peptides in hdms_data
uptakes = UptakePlotsCollection(if_plot_fit=True, if_reindex_title_from_pdb=False, pdb_file='../example/2dri_protein.pdb')
uptakes.add_plot_all(hdxms_data)
uptakes.save_plots(OUTDIR)


# make all the comparison plots and save them

from itertools import product

items =[state.state_name for state in hdxms_data.states]
combinations = product(['APO'], [x for x in items if x != 'APO'])

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

for state1, state2 in combinations:
    compare = HDXStateCompare(hdxms_data.get_state(state1), hdxms_data.get_state(state2))
    compare.add_all_compare()

    heatmap_compare_tp = create_heatmap_compare_tp(compare, colorbar_max)
    heatmap_compare_tp.savefig(f'{OUTDIR}/{state1}-{state2}-heatmap-tp.png')

    heatmap_compare = create_heatmap_compare(compare, colorbar_max)
    heatmap_compare.savefig(f'{OUTDIR}/{state1}-{state2}-heatmap.png')

    create_compare_pymol_plot(compare, colorbar_max=colorbar_max, pdb_file='../example/2dri_protein.pdb', path=OUTDIR)



