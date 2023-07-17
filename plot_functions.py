import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os

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
            ax1.plot(timepoints, element.get(state), 'o', label = state, markersize = 18, alpha = 0.5,
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
