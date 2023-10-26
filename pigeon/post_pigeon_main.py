from data_processing import *
from plot_functions import *
import argparse
import numpy as np
import pandas as pd
import glob
import os
from pymol import cmd

parser = argparse.ArgumentParser(description='Statistical analysis of HDX/MS data for curated RbsR peptides in one or more functional states.')
parser.add_argument('--pm', dest='pm', help="path to pymol structure")
parser.add_argument('--t', '--tables', dest='tables', help="path to uptake table", nargs='+') #No commas between files
parser.add_argument('--r', '--ranges', dest='ranges', help='path to ranges list csv', nargs='+') #No commas between files
parser.add_argument('--e', '--exclude', dest='exclude', action='store_true', help='exclude rather than include rangeslist')
parser.add_argument('--s1', dest='s1', help='first state to compare')
parser.add_argument('--s2', dest='s2', help='second state to compare')
parser.add_argument('--compare', dest='compare', help='dna, ligand, both', nargs=3)
parser.add_argument('--cbarmax', dest='cbarmax', type=float, help='max value for colorbar axis for dDbar')
parser.add_argument('--ldmin', dest='ldmin', type=float, help='in dDbar, minimum difference threshold between ligand/dna states')
parser.add_argument('--o', dest='outdir', help='path to output directory')
parser.add_argument('--sub', dest='subtract', action='store_true', help='subtract peptides from each other to get more peptides and higher resolution') 

args = parser.parse_args()
print(args)
print(args.tables)
cmd.load(args.pm)

#Example usage:
#python post_pigeon_main.py --t TechRep1_peptide_pool_result.csv TechRep2_peptide_pool_result.csv  --r TechRep1_peptide_rangelist.csv TechRep2_peptide_rangelist.csv --pm 1pfk_Xray.pdb --o outputdir --cbarmax 0.5 --sub

OUTDIR = args.outdir 
if args.outdir is None:
    OUTDIR = os.getcwd()

# load the data
hdxms_data_list = []
for i in range(len(args.tables)):

    if args.exclude:
        cleaned = read_hdx_tables([args.tables[i]], [args.ranges[i]], exclude=True)
        hdxms_data = load_dataframe_to_hdxmsdata(cleaned)
        hdxms_data.reindex_peptide_from_pdb(args.pm, first_residue_index=1)
    else:
        cleaned = read_hdx_tables([args.tables[i]], [args.ranges[i]])
        hdxms_data = load_dataframe_to_hdxmsdata(cleaned)
        hdxms_data.reindex_peptide_from_pdb(args.pm, first_residue_index=1)
    
    hdxms_data_list.append(hdxms_data)

# if args.subtract == True then do the peptide subtraction
if args.subtract:
    print('Peptides will be subtracted for increased resolution')
    # Add new peptides to protein state after subtraction 
    [state.add_all_subtract() for data in hdxms_data_list for state in data.states]
else:
    print('Peptides will not be subtracted')

# make a uptake plot for all the peptides in hdms_data
#uptakes = UptakePlotsCollection(if_plot_fit=False) #If False, no fitting (just data points)
#uptakes.add_plot_all(hdxms_data_list)
#uptakes.save_plots(OUTDIR)


from itertools import product

items =[state.state_name for state in hdxms_data.states]
combinations = product(['APO'], [x for x in items if x != 'APO'])

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

for state1_name, state2_name in combinations:

    state1_list = [i.get_state(state1_name) for i in hdxms_data_list]
    state2_list = [i.get_state(state2_name) for i in hdxms_data_list]

    compare = HDXStatePeptideCompares(state1_list, state2_list)
    compare.add_all_compare()

    heatmap_compare_tp = create_heatmap_compare_tp(compare, 0.5)
    heatmap_compare_tp.savefig(f'{OUTDIR}/{state1_name}-{state2_name}-heatmap-tp.png')

    heatmap_compare = create_heatmap_compare(compare, 0.5)
    heatmap_compare.savefig(f'{OUTDIR}/{state1_name}-{state2_name}-heatmap.png')

    heatmap_compare_separated = create_heatmap_with_dotted_line(compare,0.5)
    heatmap_compare_separated.savefig(f'{OUTDIR}/{state1_name}-{state2_name}-heatmap-sep.png')

    create_compare_pymol_plot(compare, colorbar_max=0.2, pdb_file=args.pm, path=OUTDIR)

    res_compares = HDXStateResidueCompares([i for i in range(1, 320)], state1_list, state2_list)
    res_compares.add_all_compare()

    create_compare_pymol_plot(res_compares, 0.2, pdb_file=args.pm, path=OUTDIR, save_pdb=True)