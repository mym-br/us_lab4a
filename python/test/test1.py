#!/bin/env python
# -*- coding: utf-8 -*-
# This file is in the public domain.
"""
Demo of TeX rendering.

You can use TeX to render all of your matplotlib text if the rc
parameter text.usetex is set.  This works currently on the agg and ps
backends, and requires that you have tex and the other dependencies
described at http://matplotlib.sf.net/matplotlib.texmanager.html
properly installed on your system.  The first time you run a script
you will see a lot of output from tex and associated tools.  The next
time, the run may be silent, as a lot of the information is cached in
~/.tex.cache

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

font_properties = {'family' : 'sans-serif',
                   'sans-serif' : ['Helvetica'],
                   'size' : 12}
rc('font', **font_properties)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('font', family='serif')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.unicode'] = True
#plt.rcParams['text.latex.preamble'] = r'\usepackage{cmbright}'
# Use Helvetica for math (when using LaTeX).
plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'

# Padding between the axis and the tick labels.
plt.rcParams['xtick.major.pad'] = 8.0
plt.rcParams['ytick.major.pad'] = 8.0

# Use scientific notation if log10 of the axis range is smaller than the
# first or larger than the second.
plt.rcParams['axes.formatter.limits'] = (-3, 4)
#plt.rcParams['axes.formatter.limits'] = (-7, 7)

#plt.rcParams['axes.formatter.use_mathtext'] = True
#plt.rcParams['mathtext.default'] = 'regular'

# Changes the font of the tick labels when using LaTeX.
#def prepare_figure():
#    ax = plt.gca()
#    ax.set_xticklabels(ax.get_xticks(), font_properties)
#    ax.set_yticklabels(ax.get_yticks(), font_properties)





# Example data
t = np.arange(-1.5, 1.0 + 0.01, 0.01)
s = (np.cos(4.2 * np.pi * t) + 0.5) * 1e12

plt.figure(figsize=(6, 3))
plt.plot(t * 1e5, s)
#plt.plot(t, s * 1e-12)
#prepare_figure()

#plt.margins(0.3, 0.3)

#plt.xlabel(r'\textbf{time} (s)')
plt.xlabel(r'time ($\mu \mathrm{s}$) ($\mathrm{10^5} \mu$s)')
#plt.ylabel(r'\textit{voltage} (mV)',fontsize=16)
plt.ylabel(r'voltage ($\mathrm{10^{12}}$ mV)')
#plt.title(ur'\TeX\ is Number áùmôçü'
#          r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!",
#          fontsize=16)
plt.title(r'\TeX\ is Number áùmôçü', fontsize=16)


##          fontsize=16, color='gray')
# Make room for the ridiculously large title.
#plt.subplots_adjust(top=0.8, bottom=0.2)
#plt.subplots_adjust(top=0.8)

# Padding between the x and y labels and the figure borders.
plt.tight_layout(pad=0.2)

plt.savefig('tex_demo.pdf')


plt.figure(figsize=(8, 4))
#prepare_figure()
plt.plot(t, s)
plt.xlabel(r'time ($\alpha \beta$) time ($\mu \mathrm{s}$) ($\mathrm{10^5} \mu$s)')
plt.ylabel(r'Foobar ($\zeta$)')

# Padding between the x and y labels and the figure borders.
plt.tight_layout(pad=0.2)
plt.savefig('tex_demo2.pdf')
plt.savefig('tex_demo2.png', dpi=300)

plt.show()
