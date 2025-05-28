# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 14:39:37 2024

@author: azamuner
"""

import re
import numpy as np
import uncertainties
from uncertainties import ufloat
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

script_dir = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(script_dir, '..', 'Data')


def add_sample_information(t, dataset):
    add_sample_name(t, dataset)
    insert_inj_n(t)
    clean_and_add_condition_and_replicate(t, dataset)
    return t


def add_sample_name(t, dataset):
    t.add_or_replace_column('sample_name', 
                            t.apply(_get_sample_name, 
                                    t.filename,
                                    dataset), 
                            str, 
                            insert_after='compound')


def _get_sample_name(filename, dataset):
    if dataset == '1to1':
        return re.search("_1to1_(.+).mzML", filename).group(1)
    else:
        return re.search("_1to10_(.+).mzML", filename).group(1)


def clean_and_add_condition_and_replicate(t, dataset):
    t.replace_column('sample_name', 
                     t.apply(_clean_sample_name, 
                             t.sample_name,
                             dataset),
                     str)
    t.add_or_replace_column('condition', 
                            t.apply(_add_condition, 
                                    t.sample_name),
                         str,
                         insert_after = 'sample_name')
    t.add_or_replace_column('replicate_number', 
                            t.apply(_add_replicate_number, 
                                    t.sample_name,
                                    t.inj_n),
                            str,
                            insert_after = 'condition')


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


def _clean_sample_name(filename, dataset):
    if dataset == '1to1':
        return re.search("(.+)_1to1", filename).group(1)
    else:
        return re.search("(.+)_1to10", filename).group(1)


def _add_condition(filename):
    return filename.split('_')[0]


def _add_replicate_number(filename, inj_n):
    if 'MM' not in filename:
        return filename.split('_')[1]
    else:
        return str(inj_n)


def plot_data(df):
    condition_order = ['MM', 'neg', 'L68', 'L203', 'L265']
    df['condition'] = pd.Categorical(df['condition'], categories=condition_order, ordered=True)
    compounds = df['compound'].unique()
    
    for compound in compounds:
        plt.figure(figsize=(8, 6))
        
        subset = df[df['compound'] == compound]
        subset = subset.sort_values('condition')
        conditions = subset['condition'].values
        nominal_values = [x.nominal_value for x in subset['average_conc']]
        std_devs = [x.std_dev for x in subset['average_conc']]
        loq = get_limit_of_quantification(compound)
        not_reached_loq = subset.groupby('condition', observed=False)['reached_LOQ'].apply(lambda x: (x == 'below LOQ').sum())
        plt.bar(conditions, nominal_values, yerr=std_devs, capsize=5)
        
        subset['predicted_c[nM]'] = subset.apply(lambda row: check_and_set_loq(row, loq), axis=1)
        sns.stripplot(x=subset['condition'], y=[x.nominal_value for x in subset['predicted_c[nM]']], 
                      hue=subset['reached_LOQ'], jitter=0.1, dodge=False, marker="o", palette={'above LOQ': 'red', 'below LOQ': 'black'}, zorder=2)
        plt.axhline(y=loq, color='blue', linestyle='--', label='LOQ', zorder=1)
        y_min, y_max = plt.gca().get_ylim()
        for condition in conditions:
            if not_reached_loq[condition] > 0:
                plt.text(condition, y_min + (y_max - y_min) * 0.025, f'n<LOQ = {int(not_reached_loq[condition])}', color='black', ha='center')
        plt.ylabel('predicted concentration [nM]')
        plt.title(f'{compound}')
        plt.xticks(rotation=45)
        plt.ylim(0 - (y_max) * 0.01, y_max)
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.savefig(os.path.join(data_folder,'results', compound+'_results.svg' ), 
                    dpi = 1000, bbox_inches='tight', format = 'svg')


def check_and_set_loq(row, LOQ):
    if row['predicted_c[nM]'].nominal_value < LOQ:
        return ufloat(0, 0)
    else:
        return row['predicted_c[nM]']


def get_limit_of_quantification(c):
    cTolimit =  {'Nicotinamide': 0.25,
            'Pantothenate': 0.25,
            'Thiamine' : 0.01,
            'Biotin' : 0.25,
            'Dethiobiotin' : 0.01}
    return cTolimit.get(c)


def add_mean_concentration(t):
    t.add_or_replace_column('average_conc', 
                            t.group_by(t.compound, 
                                       t.condition).aggregate(_get_average_concentration, 
                                                              t['predicted_c[nM]'], 
                                                              t.compound), 
                            object)
                                                              
def _get_average_concentration(conc, compound):
    LOQ = get_limit_of_quantification(list(set(compound))[0])
    conc_filt = [c for c in conc if c.nominal_value >= LOQ]
    if len(conc_filt)==1:
        return conc_filt[0]
    elif len(conc_filt)==0:
        return ufloat(0,0)
    nominal_values = np.vectorize(uncertainties.nominal_value)(conc_filt)
    st_error_replicates = np.vectorize(uncertainties.std_dev)(conc_filt)
    average = np.average(nominal_values, weights=1/ (st_error_replicates**2))
    variance_between_replicates = np.var(nominal_values, ddof=1)
    # weighted_variance = np.average(st_error_replicates**2, weights=1/(st_error_replicates**2))
    # combined_standard_error = np.sqrt((weighted_variance + variance_between_replicates) / len(conc_filt))
    combined_standard_error = np.sqrt(variance_between_replicates / len(conc_filt))
    return ufloat(average, combined_standard_error)


def reached_LOQ(t):
    t.add_or_replace_column('reached_LOQ', 
                            t.apply(_reached_LOQ,t['predicted_c[nM]'], t.compound),
                            str)
    
def _reached_LOQ(conc, compound):
    LOQ = get_limit_of_quantification(compound)
    if conc.nominal_value >= LOQ:
        return 'above LOQ'
    else:
        return 'below LOQ'