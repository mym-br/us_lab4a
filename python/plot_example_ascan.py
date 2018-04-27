#!/usr/bin/env python3
# This file is in the public domain.

#from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from scipy.signal import resample

from util.plotutil import PlotConfig
pcfg = PlotConfig('fig/plot_example_ascan')
pcfg.set_tick_labels_padding(6.0)



# Soma de dois fasores, para defasagem de zero a 2 pi.
a = 10
b = 9
th = np.arange(0.0, 2.0, 0.01)
if   pcfg.t1(): plt.figure()
elif pcfg.t3(): plt.figure(figsize=(6.7, 5))
plt.plot(th, np.abs(a + b * np.exp(1j * th * np.pi)), 'k')
plt.grid(True)
#plt.xlabel("$\\Delta\\phi / \\pi$")
if pcfg.t1():
    pass
elif pcfg.t3():
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + "/sum_phasors.pdf")

# A-scan (40 Msamples/s).
ascan = np.loadtxt("../data/example_ascan.txt")
t = np.linspace(0.0, (len(ascan) - 1) / 40.0, len(ascan))
if   pcfg.t1(): plt.figure()
elif pcfg.t3(): plt.figure(figsize=(5, 3))
plt.plot(t, ascan / 2048.0, 'k')
plt.grid(True)
if pcfg.t1():
    plt.xlabel("t ($\\mathrm{{\\mu}s}$)")
    plt.ylabel("Valor normalizado")
elif pcfg.t3():
    plt.xlabel(r'$t$ ($\si{\micro\second}$)')
    plt.ylabel("Valor normalizado")
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + "/example_ascan.pdf")

# Detalhe do A-scan (40 Msamples/s).
if   pcfg.t1(): plt.figure()
elif pcfg.t3(): plt.figure(figsize=(5, 3))
ascanDetail = ascan[1250:1400]
plt.plot(t[1250:1400], ascanDetail / 2048.0, 'k')
plt.grid(True)
if pcfg.t1():
    plt.xlabel("t ($\\mathrm{{\\mu}s}$)")
    plt.ylabel("Valor normalizado")
elif pcfg.t3():
    plt.xlabel(r'$t$ ($\si{\micro\second}$)')
    plt.ylabel("Valor normalizado")
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + "/example_ascan_detail.pdf")

# Detalhe do A-scan.
if   pcfg.t1(): plt.figure()
elif pcfg.t3(): plt.figure(figsize=(5, 3))
ascanDetail2 = resample(ascanDetail, len(ascanDetail) * 4)
tAscanDetail2 = np.linspace(1250.0 / 40.0, 1399.0 / 40.0, len(ascanDetail) * 4)
ascanEnv = abs(hilbert(ascanDetail2))
plt.plot(tAscanDetail2, ascanDetail2 / 2048.0, 'k--', tAscanDetail2, ascanEnv / 2048.0, 'k')
plt.grid(True)
plt.legend( ('sinal RF', 'envelope') )
if pcfg.t1():
    plt.xlabel("t ($\\mathrm{{\\mu}s}$)")
    plt.ylabel("Valor normalizado")
elif pcfg.t3():
    plt.xlabel(r'$t$ ($\si{\micro\second}$)')
    plt.ylabel("Valor normalizado")
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + "/example_ascan_envelope.pdf")

plt.show()
