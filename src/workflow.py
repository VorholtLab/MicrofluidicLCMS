# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:43:06 2024

@author: azamuner
"""
import os
import emzed
here = os.path.abspath(os.path.dirname(__file__))
os.chdir(here)
from extract_peaks import extract_peaks
from Calibration import extract_samples_and_predict_concentrations
from post_processing import add_sample_information
from post_processing import reached_LOQ, add_mean_concentration, plot_data
script_dir = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(script_dir, '..', 'Data')


def workflow():
    samples_1to1 = get_table('mzML_files_1to1')
    samples_1to10 = get_table('mzML_files_1to10')
    t = emzed.Table.stack_tables([samples_1to1, samples_1to10])
    t.save(os.path.join(data_folder, 'results', 'results_samples.table'), overwrite = True)
    final_data = t.extract_columns('compound', 'condition', 'average_conc', 'predicted_c[nM]', 'reached_LOQ').to_pandas()
    final_data.to_excel(os.path.join(data_folder,'results', 'results_samples.xlsx'), index = False)
    plot_data(final_data)
    
    
def get_table(files):
    paths = os.path.join(data_folder,files, "*.mzML")
    pms = emzed.io.load_peak_maps(paths) 
    ref_peak_table = emzed.io.load_excel(os.path.join(data_folder,'peaks_table.xlsx'))
    if '1to10' in paths:
        dataset = '1to10'
    else:
        dataset = '1to1'
    peak_table = extract_peaks(ref_peak_table, pms)
    samples = extract_samples_and_predict_concentrations(peak_table, dataset)
    add_sample_information(samples, dataset)
    if dataset == '1to1':
        samples = samples.filter((samples.compound != 'Pantothenate') | ((samples.compound == 'Pantothenate') & (samples.condition.is_in(['MM', 'L68', 'L265']))))
    if dataset == '1to10':
        samples = samples.filter(samples.condition.is_in(['neg', 'L203']))
    reached_LOQ(samples)
    add_mean_concentration(samples)
    return samples


workflow()





    






