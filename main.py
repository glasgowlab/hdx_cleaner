# %% [markdown]
# ### load the data

# %%
from data_processing import process_table, load_ranges_file, process_cleaned_data
import argparse
import numpy as np
import pandas as pd
from pymol import cmd



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



args = parser.parse_args(args=['--t', './example/10-25_rbsb_pool_results-CLEANED.csv',
                               '--r', './example/rangeslist-10-25.csv', 
                               '--pm', './example/2dri_protein.pdb',
                               '--ldmin','0.3'])


print(args)
cmd.load(args.pm)
colorbar_max = 0.2 if args.cbarmax is None else args.cbarmax
delta_LD_threshold = 0.075 if args.ldmin is None else args.ldmin



newbigdf = pd.DataFrame()

# Process all tables
for table in args.table:
    newbigdf = pd.concat([newbigdf, process_table(table)])

# Convert columns to the appropriate data types
newbigdf['Start'] = newbigdf['Start'].apply(np.int64)
newbigdf['End'] = newbigdf['End'].apply(np.int64)
newbigdf['#D'] = newbigdf['#D'].apply(float)

# %% [markdown]
# ### clean the data

# %%

from data_processing import get_unique_sorted, create_sequence_dict

cleaned  = load_ranges_file(args.ranges, newbigdf, args.exclude)
cleaned = process_cleaned_data(cleaned)


states = list(dict.fromkeys(cleaned['Protein State']))
#states.sort()
peptides = list(dict.fromkeys(cleaned['Sequence']))
#peptides.sort()
timepoints = list(dict.fromkeys(cleaned['Deut Time (sec)']))
#timepoints.sort()


states_dict, first_res, peptides_2 = create_sequence_dict(cleaned, states)

# %% [markdown]
# ### non linear fitting of the data

# %%
from utils import compile_exchange_info, fit_functions


# Compile the exchange information
peptide_exchange_dict, stdev_dict_dict = compile_exchange_info(cleaned, states, states_dict)

# Fit the exchange functions
trialT, peptide_fit_dict, peptide_params_dict, peptide_err_dict = fit_functions(peptides, peptide_exchange_dict, timepoints)


# %% [markdown]
# ### plot the uptake plots

# %%
import matplotlib.pyplot as plt
from plot_functions import assign_colors, create_plot, handle_element, handle_list_item, adjust_plot, save_plot

color_dict = assign_colors(states)
ax1 = create_plot()

for peptide in peptide_exchange_dict:
    element = peptide_exchange_dict.get(peptide)
    if element is not None:
        handle_element(peptide, element, timepoints, stdev_dict_dict, color_dict, peptide_fit_dict, ax1, trialT)
        for list_item in peptides_2:
            plot_title, y_lim = handle_list_item(list_item, peptide, ax1, peptides_2)
            if plot_title and y_lim is not None:
                adjust_plot(ax1, plot_title, y_lim)
                save_plot(plot_title)
        ax1.clear()

plt.close('all')

# %%



