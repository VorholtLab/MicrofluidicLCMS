# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:59:36 2024

@author: azamuner
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 17:27:31 2024

@author: azamuner
"""
import re
import numpy as np
import emzed
from matplotlib import pyplot as plt


def get_concentrations(t):
    standards = get_standards(t)
    standards.add_or_replace_column('conc_std', 
                                    standards.apply(_get_conc, 
                                                    standards.filename), 
                                    float, 
                                    insert_after='filename')
    insert_inj_n(standards)
    coeffs = {}
    for c in standards.compound.unique_values():
        stds_c = standards.filter(standards.compound == c)
        coeffs[c] = get_calibration_curve(stds_c.conc_std.to_list(), 
                                          stds_c.area_chromatogram.to_list(),
                                          stds_c.inj_n.to_list(),
                                          c)

        
def get_standards(t):
    t_ = t.copy()
    t_.add_or_replace_column('check', t_.apply(_check, t_.filename), bool)
    stds = t_.filter(t_.check)
    stds.drop_columns('check')
    return stds

def _check(filename):
    if '_std_' in filename:
        return True
    else:
        return False
    
def get_samples(t):
    t_ = t.copy()
    t_.add_or_replace_column('check', t_.apply(_check_for_sample, t_.filename), bool)
    samples = t_.filter(t_.check)
    samples.drop_columns('check')
    return samples

def _check_for_sample(filename):
    if '_std_' not in filename:
        return True
    else:
        return False
    
def _get_conc(filename):
    nameToConc = get_NameToConc()
    match = re.search("_std_(.+).mzML", filename).group(1)
    return nameToConc.get(match)

def get_NameToConc():
    return {'500nM': 500,
            '250nM': 250,
            '100nM': 100,
            '75nM' : 75,
            '50nM' : 50,
            '25nM' : 25,
            '10nM' : 10,
            '7_5nM': 7.5,
            '5nM'  : 5,
            '2_5nM': 2.5,
            '1nM'  : 1,
            '0_75nM':0.75,
            '0_5nM':0.5,
            '0_25nM': 0.25,
            '0_1nM': 0.1,
            '75pM': 0.075,
            '50pM': 0.05,
            '25pM':0.025,
            '10pM':0.01,
            '0nM': 0}

def insert_inj_n(pt, inj_n_start_bottle = -1):
    pt.add_or_replace_column('inj_n_absolut', 
                             pt.apply(_extract_inj_n, pt.filename), 
                             int, 
                             insert_after = 'filename'
                             )
    if inj_n_start_bottle == -1:    
        inj_n_start_bottle = min(pt.inj_n_absolut.to_list())-1
    pt.add_or_replace_column('inj_n', 
                             pt.apply(lambda x: x - inj_n_start_bottle, pt.inj_n_absolut),
                             int, 
                             insert_after= 'inj_n_absolut'
                             )
    pt.drop_columns('inj_n_absolut')
    
def _extract_inj_n(source):
    match = re.search("_(.+)_AZ", source)
    return match.group(1)

def get_calibration_curve(conc, area, inj_n, compound):
    conc = np.array(conc)
    area = np.array(area)
    inj_n = np.array(inj_n)

    plt.figure()
    
    plt.scatter(conc, area, c = inj_n, cmap='viridis', s=20, edgecolor='k')
    for i in range(len(conc)):
        plt.text(conc[i], area[i], f'{inj_n[i]:.2f}', fontsize=8, ha='right', va='bottom')    
    # plt.plot(conc, area, 'o', markersize = 3)
    plt.colorbar(label='Value of third array')
    plt.title(f'{compound}', fontsize=10)




