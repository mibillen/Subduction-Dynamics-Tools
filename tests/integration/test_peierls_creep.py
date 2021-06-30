# -*- coding: utf-8 -*-
#!/opt/anaconda3/bin/python3
r"""Test for peierls rheology

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest 

descriptions:
    every function is a separate test, combined usage with pytest module
""" 

# makes  a deformation mechanism map in P-T space for constant grain size and strainrate

import os
import shutil
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# import module from python scripts folder
import python_scripts.flow_law_functions as flf

test_dir = "test_output"
if not os.path.isdir(test_dir):
    os.mkdir(test_dir)
# source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_run_time')

def test_plot_peierls_hansen19():
    # test using implementation from Hansen etal 2019
    # this function reproduces fig 8 from the paper
    T = 298  # room temperature
    edot = 1e-5  # guessed
    ds = 10**np.linspace(-2, 5, 100) * 1e-6
    viscs = flf.peierls_hansen19_visc(T,edot,ds, sigma_b = 1.8e9)
    stresses = 2 * viscs * edot
    std = np.array([[-0.9861750909741067, 0.9732704402515723],
     [-0.8179781662144259, 0.9150943396226414],
     [-0.6329799541448438, 0.8522012578616351],
     [-0.38082917902442137, 0.7735849056603773],
     [-0.19598872552112967, 0.7201257861635219],
     [-0.01120085820660055, 0.6698113207547169],
     [0.22385940556572148, 0.6132075471698113],
     [0.475589491176037, 0.5597484276729559],
     [0.6768368355735053, 0.5251572327044025],
     [0.8947014156202009, 0.4968553459119496],
     [1.0622935992091032, 0.4748427672955975],
     [1.246555604635999, 0.4559748427672956],
     [1.4977598283586788, 0.4339622641509433],
     [1.8492985002418951, 0.4119496855345912],
     [2.1170673734250434, 0.39937106918238996],
     [2.4015060684461824, 0.3899371069182389],
     [2.736059401358826, 0.3836477987421385],
     [2.987053280326454, 0.37421383647798745],
     [3.2880566248080596, 0.37421383647798745],
     [3.62255737153194, 0.371069182389937],
     [3.890168486148797, 0.3679245283018868]])
    # plot fig 8
    fig = plt.figure(figsize=(5, 10))
    ax1 = fig.add_subplot(2,1,1)
    ax1.loglog(ds/1e-6, stresses/1e9, 'b', label='test, edot=%.4e' % edot)
    ax1.loglog(10.0**std[:, 0], 10.0**std[:, 1], 'k--', label='standard')
    ax1.grid()
    ax1.set_xlim([1e-2, 1e5])
    ax1.set_ylim([1.0, 30.0])
    ax1.set_xlabel('Length scale(um)')
    ax1.set_ylabel('Differential stress(GPa)')
    ax1.legend()
    ax1.set_title('fig 8')
    pdffile = os.path.join(test_dir, 'peierls_hansen19.pdf')
    # fig 10,c (temperature dependence)
    Ts = np.linspace(200.0, 1000.0, 100)
    d = 0.01  # 1 cm
    viscs = flf.peierls_hansen19_visc(Ts,edot,d, sigma_b=0.0)  # sigma_b = 0, corresponding to epsl = 0.0
    stresses = 2 * viscs * edot
    std = np.array([[344.146218960934, 2.2154696132596685],
    [576.6130514142911, 1.541436464088398],
    [786.6619299040804, 0.9558011049723758],
    [1045.8894292609425, 0.22651933701657434]])
    # plot fig 10, c
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(Ts, stresses/1e9, 'b', label='test, edot=%.4e' % edot)
    ax2.plot(std[:, 0], std[:, 1], 'k--', label='standard')
    ax2.set_xlim([0, 1200.0])
    ax2.set_ylim([0, 3.0])
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Differential stress(GPa)')
    ax2.legend()
    ax2.set_title('fig 9')
    fig.savefig(pdffile, bbox_inches='tight')
    assert(os.path.isfile(pdffile))