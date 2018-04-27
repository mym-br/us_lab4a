# This file is in the public domain.

import matplotlib.pyplot as plt
from matplotlib import rc

# 1: gui
# 2: article
# 3: thesis
current_plot_type = 1

class PlotConfig:
    def __init__(self, plot_dir):
        if type(plot_dir) != str:
            raise TypeError('plot_dir is not string.')
        self.plot_dir = plot_dir

        if current_plot_type == 1:
            pass
        elif current_plot_type == 2:
            #font_properties = {'family' : 'sans-serif',
            #                   'sans-serif' : ['Helvetica'],
            #                   'size' : 12}
            font_properties = {'family' : 'serif',
                               'serif' : ['Times'],
                               'size' : 12}
            rc('font', **font_properties)
            ## for Palatino and other serif fonts use:
            #rc('font',**{'family':'serif','serif':['Palatino']})
            #plt.rc('font', family='serif')

            plt.rc('text', usetex=True)
            plt.rcParams['text.latex.unicode'] = True
            #plt.rcParams['text.latex.preamble'] = r'\usepackage{cmbright}'
            # Use Helvetica for math (when using LaTeX).
            #plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'
            plt.rcParams['text.latex.preamble'] = r'\usepackage{siunitx}'

            # Use scientific notation if log10 of the axis range is smaller than the
            # first or larger than the second.
            plt.rcParams['axes.formatter.limits'] = (-3, 4)
            #plt.rcParams['axes.formatter.limits'] = (-7, 7)
        elif current_plot_type == 3:
            font_properties = {'family' : 'serif',
                               'serif' : ['Times'],
                               'size' : 12}
            rc('font', **font_properties)

            plt.rc('text', usetex=True)
            plt.rcParams['text.latex.unicode'] = True
            plt.rcParams['text.latex.preamble'] = r'\usepackage{siunitx}'

            # Use scientific notation if log10 of the axis range is smaller than the
            # first or larger than the second.
            plt.rcParams['axes.formatter.limits'] = (-3, 4)
            #plt.rcParams['axes.formatter.limits'] = (-7, 7)
        else:
            raise ValueError('Invalid plot type: ' + str(current_plot_type))

    def set_tick_labels_padding(self, padding):
        # Padding between the axis and the tick labels.
        plt.rcParams['xtick.major.pad'] = padding
        plt.rcParams['ytick.major.pad'] = padding

    def t1(self):
        return current_plot_type == 1
    def t2(self):
        return current_plot_type == 2
    def t3(self):
        return current_plot_type == 3
