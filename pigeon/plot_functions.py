import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import matplotlib.colors as col
from matplotlib.patches import Rectangle
from matplotlib import cm
import pandas as pd
import seaborn as sns
import numpy as np


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
    def __init__(self, hdxms_datas, sequence:str, color_dict=None, if_plot_fit=True,  figure=None, ax=None):
        '''
        hdxms_datas: list of class HDXMSData objects
        '''
        self.hdxms_datas = hdxms_datas
        self.sequence = sequence
        self.color_dict = self.make_color_dict(color_dict)
        self.if_plot_fit = if_plot_fit
        self.figure = figure
        self.ax = ax
        self.uptakeplot = self.make_uptakeplot()

    def make_uptakeplot(self):
        if self.figure is None and self.ax is None:
            figure, ax = plt.subplots(1, 1, figsize=(9,8))

        scatter_shapes = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']

        for hdxms_data_index, hdxms_data in enumerate(self.hdxms_datas):
            for state in hdxms_data.states:
                
                peptide = state.get_peptide(self.sequence)

                if peptide is not None:

                    ax.plot([tp.deut_time for tp in peptide.timepoints], 
                            [tp.num_d for tp in peptide.timepoints], 
                            scatter_shapes[hdxms_data_index],
                            markersize = 18, alpha = 0.5, label =state.state_name, color=self.color_dict[state.state_name])

                    # Plot the fit
                    if self.if_plot_fit:
                        trialT, y_pred, popt, perr =  peptide.fit_results
                        ax.plot(trialT, y_pred, '-', color=self.color_dict[state.state_name])

                    # set up the plot
                    ax.set_ylabel('# Deuterons')
                    ax.set_xlabel('Time (seconds)')
                    ax.set_xscale('log')
                    ax.set_ylim(0, peptide.max_d*1.1)
                    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.legend()

                    self.title = f'{peptide.start}-{peptide.end} {peptide.sequence}'
                    ax.set_title(self.title)
                    plt.close()

        return figure
    
    def make_title(self):
        start = self.hdxms_datas[0].states[0].get_peptide(self.sequence).start
        end = self.hdxms_datas[0].states[0].get_peptide(self.sequence).end
        return f'{start}-{end} {self.sequence}'
    
    def make_color_dict(self, color_dict=None):
        if color_dict is None:
            colors = ['k', 'red', 'blue', 'purple', 'gray', 'orange', 'yellow', 'green', 'brown']
            
            color_dict = {}
            state_names = list(set([state.state_name for hdxms_data in self.hdxms_datas for state in hdxms_data.states]))
            state_names.sort()
            for i, state_name in enumerate(state_names):
                color_dict[state_name] = colors[i]
        return color_dict
    
    def reindex_title_from_pdb(self, pdb_file, first_residue_index=1):
        
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


        pdb_seq = pdb2seq(pdb_file)
        start, end = find_peptide(pdb_seq, self.sequence)
        if start == -1:
            raise ValueError(f'Peptide {self.sequence} not found in pdb file {pdb_file}')
        
        self.title = f'{start+first_residue_index}-{end+first_residue_index} {self.sequence}'

        self.uptakeplot.axes[0].set_title(self.title)
        
    
    
class UptakePlotsCollection:
    def __init__(self, color_dict=None, if_plot_fit=True, if_reindex_title_from_pdb=False, pdb_file=None):
        #self.hdxms_data = hdxms_data
        self.plots = []
        self.color_dict = color_dict
        self.if_plot_fit = if_plot_fit
        self.if_reindex_title_from_pdb = if_reindex_title_from_pdb
        self.pdb_file = pdb_file

        
    def add_plot(self, hdxms_data, sequence:str):
        if self.if_reindex_title_from_pdb:
            plot = UptakePlot([hdxms_data], sequence, if_plot_fit=self.if_plot_fit, color_dict=self.color_dict)
            plot.reindex_title_from_pdb(self.pdb_file)
        else:
            plot = UptakePlot([hdxms_data], sequence)
        self.plots.append(plot)
        
    def add_plot_all(self, hdxms_data):

        def get_unique_sequences(hdxms_data):
            sequences = []
            for state in hdxms_data.states:
                for peptide in state.peptides:
                    sequences.append((peptide.start, peptide.end, peptide.sequence))
            sequences = list(set(sequences))
            sequences.sort()
            return sequences
        
        unique_sequences = get_unique_sequences(hdxms_data)
        for sequence in unique_sequences:
            self.add_plot(hdxms_data, sequence[2])

    def save_plots(self, path):
        folder_name = os.path.join(path, 'uptake_plots')
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        for plot in self.plots:
            fig = plot.uptakeplot
            fig.savefig(f'{folder_name}/{plot.title}.png', bbox_inches='tight')


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

    leftbound = compare.peptide_compares[0].peptide1.start-10
    rightbound = compare.peptide_compares[-1].peptide1.end+10
    ax.set_xlim(leftbound,rightbound)
    ax.xaxis.set_ticks(np.arange(round(leftbound,-1), round(rightbound,-1), 10))
    ax.set_ylim(-5,110)
    ax.grid(axis='x')
    ax.yaxis.set_ticks([])

    #sns.heatmap(compare_df, cmap="RdBu", linewidths=1, vmin=-colorbar_max, vmax=colorbar_max, ax=ax)
    norm = col.Normalize(vmin=-colorbar_max,vmax=colorbar_max)

    for i,peptide_compare in enumerate(compare.peptide_compares):
        
        rect = Rectangle((peptide_compare.peptide1.start, (i % 20) * 5 + ((i // 20) % 2) * 2.5), 
                         peptide_compare.peptide1.end - peptide_compare.peptide1.start,
                         4,
                         fc=colormap(norm(peptide_compare.deut_diff_avg)))
        
        ax.add_patch(rect)

    fig.colorbar(cm.ScalarMappable(cmap=colormap, norm=norm))

    ax.set_title(compare.state1.state_name + '-' + compare.state2.state_name)
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
    
    ax = sns.heatmap(df, cmap=colormap, linewidths=.75, vmin=-colorbar_max, vmax=colorbar_max)
    ax.set_title(compare.state1.state_name + '-' + compare.state2.state_name)
    ax.set_ylabel('')
    plt.close()

    return fig
