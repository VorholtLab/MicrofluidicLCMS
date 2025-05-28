# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 15:37:01 2024

@author: pkiefer
"""
from pybaselines import Baseline
from emzed.ms_data.peak_map import MSChromatogram, Chromatogram, ChromatogramType
import numpy as np


def subtract_baselines(t, subtract_baseline=True, **kwargs):
    update = t.add_or_replace_column
    if subtract_baseline:
        pstfx = "__0"
        names = [name for name in t.col_names if name.endswith("chromatogram")]
        # names = ['mz', 'rtmin_chromatogram', 'rtmax_chromatogram', 'chromatogram']
        name2type = dict(zip(t.col_names, t.col_types))
        for name in names:
            update(name + pstfx, t[name], name2type[name], format_=None)
        t.replace_column(
            "chromatogram",
            t.apply(_subtract_baseline, t.chromatogram, ignore_nones=False),
            MSChromatogram,
        )


def _subtract_baseline(ms_chromatogram):
    if ms_chromatogram is not None:
        bline = _get_baseline_from_mixture_model(ms_chromatogram)
        ints_corr = ms_chromatogram.intensities - bline
        # to avouid negative trends
        ints_corr[ints_corr < 0] = 0
        return build_chromatogram(ms_chromatogram, ints_corr)


def _get_baseline_from_mixture_model(chromatogram):
    rts = chromatogram.rts
    ints = chromatogram.intensities
    if len(ints) <4:
        return np.array([0]*len(ints))
    baseline_fitter = Baseline(x_data=rts)
    lam = len(ints) / 1000
    num_knots = len(ints) // 100
    if lam < 1:
        lam = 1
    if num_knots < 4:
        num_knots = 4
    baseline = baseline_fitter.mixture_model(ints, lam=lam, num_knots=num_knots)[0]
    # to avoid increasing values:
    baseline[baseline < 0] = 0
    # we check for succesfull fitting
    if any(np.isnan(baseline)):
        return np.array([0] * len(ints))
    return baseline


def build_chromatogram(chromatogram, intensities):
    if type(chromatogram) == MSChromatogram:
        type_ = _get_type(chromatogram)
        return MSChromatogram(
            chromatogram.mz,
            chromatogram.precursor_mz,
            chromatogram.rts,
            intensities,
            type_,
        )
    if type(chromatogram) == Chromatogram:
        return Chromatogram(chromatogram.rts, intensities)


def _get_type(ms_chromatogram):
    srm = ChromatogramType.SELECTED_REACTION_MONITORING_CHROMATOGRAM
    sim = ChromatogramType.SELECTED_ION_MONITORING_CHROMATOGRAM
    d = {srm.name: srm, sim.name: sim}
    return d.get(ms_chromatogram.type)
