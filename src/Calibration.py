# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 17:27:31 2024

@author: azamuner
"""
import re
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from uncertainties import ufloat
from post_processing import insert_inj_n

script_dir = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(script_dir, '..', 'Data')


def extract_samples_and_predict_concentrations(t, dataset):
    if dataset == '1to10':
        t = t.filter(t.compound == 'Pantothenate')

    standards = get_standards(t)
    standards.add_or_replace_column('std_concs_std', 
                                    standards.apply(_get_std_concs, 
                                                    standards.filename), 
                                    float, 
                                    insert_after='filename')
    standards.save(os.path.join(data_folder,'Calibration_results', dataset+'_standards.table' ), overwrite=True)
    insert_inj_n(standards)
    coeffs = {}
    for c in standards.compound.unique_values():
        stds_c = standards.filter(standards.compound == c)
        coeffs[c] = get_calibration_curve(stds_c.std_concs_std.to_list(), 
                                          stds_c.area_chromatogram.to_list(), 
                                          c,
                                          dataset)
    samples = get_samples(t)
    samples.add_or_replace_column('predicted_c[nM]', 
                                  samples.apply(predict_concentration, 
                                                samples.area_chromatogram, 
                                                samples.compound,
                                                coeffs,
                                                dataset),
                                  object,
                                  insert_after='compound')
    return samples
        

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
    
    
def _get_std_concs(filename):
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


def get_NameToConc_panto():
    return {'500nM': 500,
            '100nM': 100,
            '75nM' : 75,
            '50nM' : 50,
            '10nM' : 10,
            '7_5nM': 7.5,
            '5nM'  : 5,
            '2_5nM': 2.5,
            '1nM'  : 1 }


def get_calibration_curve(std_concs, area, compound, dataset):
    std_concs = np.array(std_concs)
    area = np.array(area)
    
    loq = get_limit_of_quantification(compound, dataset)
    std_concs_filt = std_concs[np.logical_or(std_concs >= loq, std_concs ==0)]
    area_filt = area[np.logical_or(std_concs >= loq, std_concs ==0)]

    if dataset=='1to1':
        weights = 1/ np.sqrt(area_filt)
        weights_sorted = np.sort(weights)
        if weights_sorted[-1] == np.inf:
            weights[weights == weights_sorted[-1]] = 10*weights[weights == weights_sorted[-2]]
        coeffs, cov = np.polyfit(std_concs_filt, area_filt, deg=1, w=weights, cov=True)
    else:
        area_replicates = []
        for std_conc in np.unique(std_concs):
            if std_conc >= loq or (std_conc == 0):
                area_replicates.append(area[std_concs == std_conc])
        mean_measured = np.array([np.mean(replicate) for replicate in area_replicates])
        weights = 1/(np.array([np.std(replicate) for replicate in area_replicates]))
        coeffs, cov = np.polyfit(np.unique(std_concs), mean_measured, deg=1, w =weights , cov=True)
       
    slope, intercept = coeffs
    st_error_slope, st_error_intercept = np.sqrt(np.diag(cov))  # Diagonal of cov gives variance
    slope_u = ufloat(slope, st_error_slope)
    intercept_u = ufloat(intercept, st_error_intercept)
    
    
    ### plot and save calibration curve ###
    area_predicted = np.polyval(coeffs, std_concs_filt)
    rsquared = get_rsquared(area_predicted, area_filt)
    plt.figure(figsize=(18,18))
    plt.plot(std_concs_filt, area_predicted)
    plt.plot(std_concs, area, 'o', color='red',  markersize = 3, label = 'Below LOQ')
    plt.plot(std_concs_filt, area_filt, 'o', color = 'black', markersize = 3, label = 'Above LOQ')
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Below LOQ',
                          markerfacecolor='red', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='Above LOQ',
                          markerfacecolor='black', markersize=10)]
    plt.legend(handles=legend_elements)
    plt.xlabel('Concentration [nM]')
    plt.ylabel('Measured peak area')
    plt.title(f'{compound}: slope= {slope_u}, intercept: {intercept_u}, R-squared: {round(rsquared, 3)}', fontsize=10)
    plt.savefig(os.path.join(data_folder,'Calibration_results', dataset+'_'+compound+'_calibration_curve.svg' ), 
                dpi = 10000, bbox_inches='tight', format = 'svg')
    
    return slope_u, intercept_u


def get_rsquared(area_pred, area):
    ss_res = np.sum((area - area_pred) ** 2)  # Residual sum of squares
    ss_tot = np.sum((area - np.mean(area)) ** 2)  # Total sum of squares
    return 1 - (ss_res / ss_tot)


def predict_concentration(area, compound, coeffs, dataset):
    slope_u,intercept_u = coeffs.get(compound)
    if dataset == '1to1':
        dilution_factor = 1
    else:
        dilution_factor = 10
    return (area - intercept_u)*dilution_factor/ slope_u


def get_limit_of_quantification(c, dataset):
    if dataset== '1to1':
        cTolimit =  {'Nicotinamide': 0.25,
                'Pantothenate': 0.25,
                'Thiamine' : 0.01,
                'Biotin' : 0.25,
                'Dethiobiotin' : 0.01}
    else:
        cTolimit =  {
                'Pantothenate': 0.1}
    return cTolimit.get(c)

    
