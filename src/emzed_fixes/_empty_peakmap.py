# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 12:54:44 2023

@author: pkiefer
"""

from emzed import Spectrum, PeakMap
from emzed.utils.sqlite import Connection
from emzed.ms_data.peak_map import create_table
import numpy as np


def create_empty_peak_map(polarity, ms_level):
    spectra = _empty_spectrums(polarity, ms_level)
    return _create_peakmap(spectra, polarity, ms_level)


def _empty_spectrums(polarity, ms_level):
    peaks = np.array(
        [
            [0.0, 0.0],
        ]
    )
    return [Spectrum(0, 0.0, ms_level, polarity, 0.0, peaks)]


def _create_peakmap(spectra, polarity, ms_level):
    conn = Connection(None)
    conn, access_name = create_table(conn=conn)
    info = {
        "rt_ranges": {
            ms_level: (spectra[0].rt, spectra[-1].rt),
            None: (spectra[0].rt, spectra[-1].rt),
        },
        "mz_ranges": {ms_level: (0.0, 0.0), None: (0.0, 0.0)},
        "ms_levels": {ms_level},
        "polarities": {polarity},
    }
    meta_data = {"full_source": "", "source": ""}
    pm = PeakMap(conn, access_name, meta_data, info)
    with pm.spectra_for_modification() as sp:
        [sp.add_spectrum(spec) for spec in spectra]
    return pm
