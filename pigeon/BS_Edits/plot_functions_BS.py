import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import matplotlib.colors as col
from matplotlib.patches import Rectangle
from matplotlib import cm
import pandas as pd
import seaborn as sns
import numpy as np
from data_processing import *


font = {'family' : 'Arial', 'weight' : 'normal', 'size' : 36}
axes = {'titlesize' : 36, 'titleweight' : 'bold', 'labelsize' : 36}
plt.rc('font', **font)
plt.rc('axes', **axes)
plt.rc('lines', lw = 3)
colors = ['k', 'red', 'blue', 'purple', 'gray', 'orange', 'yellow', 'green', 'brown']

def assign_colors(states):
    color_dict = {}
    for i in range(len(states)):
        color_dict[states[i]] = colors[i]
    return color_dict

def shift_seq_index(peptides_2_item, shift):
    start = int(peptides_2_item.split('- ')[0].split('-')[0]) + shift
    end = int(peptides_2_item.split('- ')[0].split('-')[1]) + shift
    return str(start) + '-' + str(end) + '- ' + peptides_2_item.split('- ')[1]

def create_plot():
    figure, (ax1) = plt.subplots(1, 1, figsize=(9,8))
    return ax1

def handle_element(peptide, element, timepoints, stdev_dict_dict, color_dict, peptide_fit_dict, ax1, trialT):
    for state in element:
        for key in stdev_dict_dict:
            if state == key:
                current_stdev_dict = stdev_dict_dict.get(key)
                current_stdev = current_stdev_dict.get(peptide).flatten()[0:len(timepoints)]
        if element.get(state).size != 0:
            #print(f"Plotting {peptide} in {state} state")
            reduced_timepoints = timepoints[0:len(element.get(state))]
            ax1.plot(reduced_timepoints, element.get(state), 'o', label = state, markersize = 18, alpha = 0.5,
                color = color_dict.get(state))
        for element_2 in peptide_fit_dict.get(peptide):
            for e2_key in element_2:
                if e2_key == state:

                    ax1.plot(trialT, element_2.get(state), '-', color = color_dict.get(state))

def handle_list_item(list_item, peptide, ax1, peptides_2):
    title_string = "- " + peptide
    y_lim = len(peptide) - 2 - peptide[2:].count('P') + 0.25
    if list_item.endswith(title_string):
        plot_title = shift_seq_index(list_item, -14)
        return plot_title, y_lim
    return None, None

def adjust_plot(ax1, plot_title, y_lim):
    ax1.set_ylabel('# Deuterons')
    ax1.set_xlabel('Time (seconds)')
    ax1.set_title(plot_title)
    ax1.set_xscale('log')
    ax1.set_ylim(0, y_lim)
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

def save_plot(plot_title):
    plt.legend(frameon = False, bbox_to_anchor=(1.03, 1.03))
    if not os.path.exists('plots/'):
        os.mkdir('plots/')
    plt.savefig('plots/' + plot_title + '.png', bbox_inches='tight')


def plot_uptake_plots(peptide_exchange_dict, timepoints, stdev_dict_dict, color_dict, peptide_fit_dict, ax1, trialT, peptides_2):
    for peptide in peptide_exchange_dict:
        element = peptide_exchange_dict.get(peptide)
        if element is not None:
            handle_element(peptide, element, timepoints, stdev_dict_dict,
                           color_dict, peptide_fit_dict, ax1, trialT)
            for list_item in peptides_2:
                plot_title, y_lim = handle_list_item(
                    list_item, peptide, ax1, peptides_2)
                if plot_title and y_lim is not None:
                    adjust_plot(ax1, plot_title, y_lim)
                    save_plot(plot_title)
            ax1.clear()

    plt.close('all')

class UptakePlot:
    def __init__(self, hdxms_datas, identifier:str, color_dict=None, if_plot_fit=True,  figure=None, ax=None):
        '''
        hdxms_datas: list of class HDXMSData objects
        '''
        self.hdxms_datas = hdxms_datas
        self.identifier = identifier
        self.sequence = identifier.split(' ')[1]
        self.color_dict = self.make_color_dict(color_dict)
        self.if_plot_fit = if_plot_fit
        self.figure = figure
        self.ax = ax
        #self.title = self.make_title()
        #self.title = identifier
        self.uptakeplot = self.make_uptakeplot()


    def make_uptakeplot(self):

        plt.rcParams['legend.fontsize'] = 22

        if self.figure is None and self.ax is None:
            figure, ax = plt.subplots(1, 1, figsize=(9,8))

        scatter_shapes = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']

        sns.lineplot(data=self.hdxms_datas_df, 
                     x='time', y='deut', hue='state', errorbar='sd', 
                     err_style='bars', marker='o', linestyle=' ', markersize=18, alpha=0.5,
                     ax=ax, palette=self.color_dict)
        
        # Plot the fit
        if self.if_plot_fit:
            for state_name in self.hdxms_datas_df.state.unique():
                avg_peptide = self.get_average_peptide(state_name)
                #trialT, y_pred, popt, perr =  avg_peptide.fit_results
                trialT, y_pred, popt, perr =  avg_peptide.new_fit()
                #trialT, y_pred, popt, perr =  avg_peptide.fit_hdx_stats()
                ax.plot(trialT, y_pred, '-', color=self.color_dict[state_name])

        # set up the plot
        ax.set_ylabel('# Deuterons')
        ax.set_xlabel('Time (seconds)')
        ax.set_xscale('log')
        avg_peptide = self.get_average_peptide(self.hdxms_datas_df.state.unique()[0])
        
        ax.set_ylim(min(self.hdxms_datas_df['deut'])-1, max(self.hdxms_datas_df['deut'])+1)
        #ax.set_ylim(-0.3, avg_peptide.max_d*1.1)
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.legend()

        ax.set_title(self.identifier)
        plt.close()

        return figure
    
    @property
    def hdxms_datas_df(self):
        
        hdxms_datas_df = pd.DataFrame()

        for hdxms_data_index, hdxms_data in enumerate(self.hdxms_datas):
            hdxms_data_df = pd.DataFrame()
            
            for state in hdxms_data.states:
                peptide = state.get_peptide(self.identifier)
                
                if peptide is not None:
                    peptide_df_i = pd.DataFrame({
                        'time': [tp.deut_time for tp in peptide.timepoints],
                        'deut': [tp.num_d for tp in peptide.timepoints],
                        'state': state.state_name
                    })
                    
                    hdxms_data_df = pd.concat([hdxms_data_df, peptide_df_i], ignore_index=True)
            
            hdxms_data_df['data_set_index'] = hdxms_data_index
            hdxms_datas_df = pd.concat([hdxms_datas_df, hdxms_data_df], ignore_index=True)
        return hdxms_datas_df

    def get_average_peptide(self, state_name):
        grouped_df = self.hdxms_datas_df.groupby('state')
        group = grouped_df.get_group(state_name)

        final_group = group.groupby('time')
        group_mean =final_group.mean(numeric_only=True).reset_index()
        group_std =final_group.std(numeric_only=True).reset_index()


        start = self.identifier.split(' ')[0].split('-')[0]
        end = self.identifier.split(' ')[0].split('-')[1]
        peptide = Peptide(self.sequence, start, end, f"averaged peptide: {state_name}")
        for i in range(len(group_mean)):
            timepoint = Timepoint(peptide, group_mean['time'][i], group_mean['deut'][i], group_std['deut'][i])
            peptide.add_timepoint(timepoint)
        
        return peptide


    def make_title(self):
        for hdxms_data in self.hdxms_datas:
            for state in hdxms_data.states:
                try:
                    start = state.get_peptide(self.sequence).start
                    end = state.get_peptide(self.sequence).end
                    return f'{start}-{end} {self.sequence}'
                except AttributeError:
                    pass
        return f'Missing data for {self.sequence}'
    
    def make_color_dict(self, color_dict=None):
        if color_dict is None:
            colors = ['k', 'red', 'blue', 'purple', 'gray', 'orange', 'yellow', 'green', 'brown']
            
            color_dict = {}
            state_names = list(set([state.state_name for hdxms_data in self.hdxms_datas for state in hdxms_data.states]))
            state_names.sort()
            for i, state_name in enumerate(state_names):
                color_dict[state_name] = colors[i]
        return color_dict

#new residue coverage plotting script heat map thing - sav
class ResidueCoverage:
    def __init__(self, hdxms_data):
        self.hdxms_data = hdxms_data

    def calculate_coverages(self, state_name):
        state = self.hdxms_data.get_state(state_name)
        max_end = max(pep.end for pep in state.peptides)
        coverage = np.zeros(max_end + 1)
        for pep in state.peptides:
            coverage[pep.start:pep.end + 1] += 1
        return coverage

    def plot(self):

        row_num = len(self.hdxms_data.states)
        fig, axes = plt.subplots(row_num, 1, figsize=(20, 3*row_num), sharex=True, sharey=True)
        
        coverage_list = []
        for state in self.hdxms_data.states:
            coverage = self.calculate_coverages(state.state_name)
            coverage_list.append(coverage)
        
        max_coverage = max([max(coverage) for coverage in coverage_list])
        min_coverage = min([min(coverage) for coverage in coverage_list])

        for i, state in enumerate(self.hdxms_data.states):   
            ax = axes.flatten()[i]
            im = ax.imshow(coverage_list[i].reshape(1, -1), aspect='auto', cmap='Blues', extent=[0, 300, 0, 1], vmin=min_coverage, vmax=max_coverage)
            ax.set_yticks([])
            ax.set_ylabel(state.state_name)   
        
        ax.set_xlabel('Resid')  

        cbar_ax = fig.add_axes([0.92, 0.4, 0.02, 0.4])  # [left, bottom, width, height]
        fig.colorbar(im, cax=cbar_ax)
        
        plt.tight_layout(rect=[0, 0, 0.95, 1])  # Adjust the layout to make room for the colorbar
        plt.show()
        


class UptakePlotsCollection:
    def __init__(self, color_dict=None, if_plot_fit=True, pdb_file=None):
        #self.hdxms_datas = hdxms_datas
        self.plots = []
        self.color_dict = color_dict
        self.if_plot_fit = if_plot_fit
        self.pdb_file = pdb_file

        
    def add_plot(self, hdxms_datas, sequence:str):
        plot = UptakePlot(hdxms_datas, sequence, if_plot_fit=self.if_plot_fit, color_dict=self.color_dict)
        self.plots.append(plot)
        
    def add_plot_all(self, hdxms_datas):

        def get_unique_sequences(hdxms_datas):
            sequences = []
            for hdxms_data in hdxms_datas:
                for state in hdxms_data.states:
                    for peptide in state.peptides:
                        idf = f'{peptide.start}-{peptide.end} {peptide.sequence}'
                        sequences.append(idf)
            sequences = list(set(sequences))
            sequences.sort()
            return sequences
        
        unique_sequences = get_unique_sequences(hdxms_datas)
        for sequence in unique_sequences:
            self.add_plot(hdxms_datas, sequence)

    def save_plots(self, path):
        folder_name = os.path.join(path, 'uptake_plots')
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        for plot in self.plots:
            fig = plot.uptakeplot
            fig.savefig(f'{folder_name}/{plot.identifier}.png', bbox_inches='tight')


def create_heatmap_compare(compare, colorbar_max, colormap="RdBu"):

    import matplotlib.colors as col
    from matplotlib.patches import Rectangle
    from matplotlib import cm

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14
        }
    axes = {'titlesize' : 18,
            'titleweight' : 'bold',
            'labelsize' : 16
        }

    plt.rc('font', **font)
    plt.rc('axes', **axes)

    plt.rcParams['figure.figsize'] = (4, 25)
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    colormap = cm.get_cmap(colormap)

    fig, ax = plt.subplots(figsize=(20,10))

    leftbound = compare.peptide_compares[0].peptide1_list[0].start-10
    rightbound = compare.peptide_compares[-1].peptide1_list[0].end+10
    ax.set_xlim(leftbound,rightbound)
    ax.xaxis.set_ticks(np.arange(round(leftbound,-1), round(rightbound,-1), 10))
    ax.set_ylim(-5,110)
    ax.grid(axis='x')
    ax.yaxis.set_ticks([])

    #sns.heatmap(compare_df, cmap="RdBu", linewidths=1, vmin=-colorbar_max, vmax=colorbar_max, ax=ax)
    norm = col.Normalize(vmin=-colorbar_max,vmax=colorbar_max)

    for i,peptide_compare in enumerate(compare.peptide_compares):
        for peptide in peptide_compare.peptide1_list:
            rect = Rectangle((peptide.start, (i % 20) * 5 + ((i // 20) % 2) * 2.5), 
                             peptide.end - peptide.start,
                             4,
                             fc=colormap(norm(peptide_compare.deut_diff_avg)))
            ax.add_patch(rect)


    fig.colorbar(cm.ScalarMappable(cmap=colormap, norm=norm))

    ax.set_title(compare.state1_list[0].state_name + '-' + compare.state2_list[0].state_name)
    fig.tight_layout()
    plt.close()

    return fig

from matplotlib import cm

# Create a function to create a heatmap
def create_heatmap_compare_tp(compare, colorbar_max, colormap="RdBu"):

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14
        }
    axes = {'titlesize' : 18,
            'titleweight' : 'bold',
            'labelsize' : 16
        }

    plt.rc('font', **font)
    plt.rc('axes', **axes)

    plt.rcParams['figure.figsize'] = (10, 25)
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'

    fig, ax = plt.subplots()

    df = pd.DataFrame()
    for pep_compare in compare.peptide_compares:
        pep_dict = {}
        pep_dict['title'] = pep_compare.compare_info.split(': ')[1]
        for tp, uptate_v in zip(pep_compare.common_timepoints, pep_compare.deut_diff):
            pep_dict[int(tp)] = uptate_v
        df = pd.concat([df, pd.DataFrame(pep_dict, index=[0])], ignore_index=True)
    df = df.set_index('title')
    #sort the columns by timepoint
    df = df.reindex(sorted(df.columns), axis=1)
    #sort the rows by sequence
    #import re
    #df['Start'] = df.index.map(lambda x: int(re.search(r"(-?\d+)--?\d+ \w+", x).group(1)))
    #df = df.sort_values(by=['Start']).drop('Start', axis=1)

    ax = sns.heatmap(df, cmap=colormap, linewidths=.75, vmin=-colorbar_max, vmax=colorbar_max)
    ax.set_title(compare.state1_list[0].state_name + '-' + compare.state2_list[0].state_name)
    ax.set_ylabel('')
    fig.tight_layout()
    plt.close()
    
    return fig
    
#Function to make a heatmap with positive and negative deuterium uptake values separated by a dotted line
def create_heatmap_with_dotted_line(compare, colorbar_max, colormap="RdBu"):
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': 14
        }
    axes = {'titlesize': 18,
            'titleweight': 'bold',
            'labelsize': 16
        }

    plt.rc('font', **font)
    plt.rc('axes', **axes)

    plt.rcParams['figure.figsize'] = (4, 25)
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    colormap = cm.get_cmap(colormap)

    fig, ax = plt.subplots(figsize=(20, 10))

    leftbound = compare.peptide_compares[0].peptide1_list[0].start - 10
    rightbound = compare.peptide_compares[-1].peptide1_list[0].end + 10
    ax.set_xlim(leftbound, rightbound)
    ax.xaxis.set_ticks(np.arange(round(leftbound, -1), round(rightbound, -1), 10))
    y_min = -105
    y_max = 105
    ax.set_ylim(y_min, y_max)
    y_middle=(y_min+y_max)/2
    ax.grid(axis='x')
    ax.yaxis.set_ticks([])
    #add a horizontal dotted line that matches with the 0.0 on the colorbar
    ax.axhline(y=y_middle, color='k', linestyle='--', linewidth=1)

    norm = col.Normalize(vmin=-colorbar_max, vmax=colorbar_max)

    fig.colorbar(cm.ScalarMappable(cmap=colormap, norm=norm))

    for i,peptide_compare in enumerate(compare.peptide_compares):
            for peptide in peptide_compare.peptide1_list:
                deut_diff_avg = peptide_compare.deut_diff_avg
                #print(deut_diff_avg)
                
                if deut_diff_avg > 0:
                    y_position = (i % 20) * 5 + ((i // 20) % 2) * 2.5 + y_middle + 2
                    
                else:
                    y_position = y_middle - ((i % 20) * 5 + ((i // 20) % 2) * 2.5) - 5  # Below the line

                rect = Rectangle((peptide.start, y_position),
                                    peptide.end - peptide.start,
                                    3,
                                    fc=colormap(norm(deut_diff_avg)))
                ax.add_patch(rect)


    ax.set_title(compare.state1_list[0].state_name + '-' + compare.state2_list[0].state_name)
    fig.tight_layout()
            
    return fig

def create_compare_pymol_plot(compares, colorbar_max, colormap="RdBu", pdb_file=None, path=None,
                              save_pdb=False):

    rgb_df = gen_rgb_df(compares, colorbar_max, colormap)

    from pymol import cmd
    cmd.delete('all')
    cmd.load(pdb_file)
    cmd.color("gray")

    if isinstance(compares, HDXStatePeptideCompares):
        for i, seq in enumerate(rgb_df['title']):
            parts = seq.split()
            if len(parts) > 0:
                resi = parts[0]  # Extract the first part (e.g., "ABC")
                seq = parts[-1]
                seq_n = "res_" + seq
                #print(seq_n)
                #print(f"Selecting: {seq}, 'resi {resi}'")
                cmd.select(seq_n, 'resi ' + resi)
                cmd.set_color(f'{seq_n}', [rgb_df['r'][i], rgb_df['g'][i], rgb_df['b'][i]])
                cmd.color(f'{seq_n}', seq_n)
                #print(f'{seq_n}', seq)

    elif isinstance(compares, HDXStateResidueCompares):
            for i, seq in enumerate(rgb_df['title']):
                if np.isnan(rgb_df['s'].values[i]):
                    continue
                parts = seq.split()
                if len(parts) > 0:
                    resi = parts[0]  # Extract the first part (e.g., "ABC")
                    cmd.select(seq, 'resi ' + resi)
                    cmd.set_color(f'res_{seq}', [rgb_df['r'][i], rgb_df['g'][i], rgb_df['b'][i]])
                    cmd.color(f'res_{seq}', seq)
                    cmd.delete(seq)
    
    cmd.ray(1000,1000)
    
    if path is not None:
        if isinstance(compares, HDXStatePeptideCompares):
            full_path = os.path.join(path, f'{compares.state1_list[0].state_name}-{compares.state2_list[0].state_name}_{colorbar_max}_pepcompare-pm.pse')
        elif isinstance(compares, HDXStateResidueCompares):
            full_path = os.path.join(path, f'{compares.state1_list[0].state_name}-{compares.state2_list[0].state_name}_{colorbar_max}_rescompare-pm.pse')
        cmd.save(full_path)
    else:
        raise ValueError('Please provide a path to save the pymol session')
    
    if isinstance(compares, HDXStateResidueCompares):
        if save_pdb:
            pdb_full_path = os.path.join(path, f'{compares.state1_list[0].state_name}-{compares.state2_list[0].state_name}_rescompare-pm.pdb')
            save_pdb = plot_on_pdb(pdb_file, compares, pdb_full_path)



def plot_on_pdb(pdb, residue_compares, path=None):

    import MDAnalysis

    u = MDAnalysis.Universe(pdb)
    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms
    
    for res_compare in residue_compares.residue_compares:
        res_id = res_compare.resid
        res = protein.select_atoms('resnum {}'.format(res_id))
        res.tempfactors = res_compare.deut_diff_avg

    if path is None:
        raise ValueError('Please provide a path to save the pymol session')
    else:
        u.atoms.write(path)

    

def gen_rgb_df(compare, colorbar_max, colormap="RdBu"):
    colormap = plt.get_cmap(colormap)
    df = pd.DataFrame()

    compare_list = compare.residue_compares if isinstance(compare, HDXStateResidueCompares) else compare.peptide_compares

    for compare_i in compare_list:
        pep_dict = {}
        pep_dict['title'] = compare_i.compare_info.split(': ')[1]
        s_i = (compare_i.deut_diff_avg + colorbar_max) / (2 * colorbar_max)
        pep_dict['s'] = s_i
        df = pd.concat([df, pd.DataFrame(pep_dict, index=[0])], ignore_index=True)

    df['s'].clip(upper=1, lower=0, inplace=True)
    rgb = colormap(df['s']) * 255
    df['r'] = rgb[:, 0]
    df['g'] = rgb[:, 1]
    df['b'] = rgb[:, 2]

    return df