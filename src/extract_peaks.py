# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:24:32 2024

@author: azamuner
"""

import os
import numpy as np
from emzed import Table, PeakMap, MzType, RtType
from emzed.quantification import integrate_chromatograms
from emzed import peak_picking as pp
from in_out import load_config
from emzed_fixes._empty_peakmap import create_empty_peak_map
from adapt_retention_time import local_adjust
from subtract_baselines import subtract_baselines

script_dir = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(script_dir, '..', 'Data')

def extract_peaks(peaks_table, samples):
    config_path = os.path.join(data_folder, "default_config.txt")
    config = load_config(config_path)
    kwargs = config.get('extract_peaks')
    _update_id_column(peaks_table)
    tables = []
    for sample in samples:
        t = add_sample_to_peaks_table(peaks_table, sample, kwargs)
        subtract_baselines(t, **kwargs)
        t = _extract(t, **kwargs)
        tables.append(t)
    return Table.stack_tables(tables)



def add_sample_to_peaks_table(pt, peakmap, kwargs):

    pt = add_peakmap_from_spectra(pt, peakmap, kwargs)
    pt.add_column_with_constant_value("filename", peakmap.meta_data["source"], str)

    return pt


def _extract(t, integration_algorithm, **kwargs):
    local_adjust(t)
    integrate_chromatograms(t, integration_algorithm, in_place=True, show_progress=False)
    t = _fall_back_integration(t)
    return t


def add_peakmap_from_spectra(pt, peakmap, kwargs):
    pt = pt.copy()
    _update_mz_range(pt, **kwargs)
    t = create_ms2_peakmap_table(peakmap)
    pt = assign_peakmaps_by_precursor(pt, t, **kwargs)
    pt = convert_peakmap_to_chromatogram(pt, **kwargs)
    return pt


def create_ms2_peakmap_table(peakmap):
    prec2pm = peakmap.split_by_precursors()
    precs, pms = zip(*list(prec2pm.items()))
    precursors = np.array(precs).reshape(1, -1)[0].tolist()
    rows = list(zip(precursors, pms))
    columns = ["precursor_mz", "peakmap"]
    types = [MzType, PeakMap]
    t = Table.create_table(columns, types, rows=rows)
    return t


def _fall_back_integration(t):
    subs = []
    for sub in t.split_by_iter("valid_model_chromatogram"):
        valids = list(set(sub.valid_model_chromatogram))
        if all(valids):
            subs.append(sub)
        else:
            integrate_chromatograms(sub, "linear", in_place=True)
            subs.append(sub)
    return Table.stack_tables(subs)


def assign_peakmaps_by_precursor(
    pt, t, precursor_mz_tol, mz_tol_abs=0, mz_tol_rel=0, **kwargs
):
    # since mz_tol_rel is given in ppm
    mz_tol_rel = mz_tol_rel / 1e6
    _update_pm_ranges(t, mz_tol_abs, mz_tol_rel)
    # precursor_mzs match
    condition1 = pt.precursor_mz.approx_equal(t.precursor_mz, precursor_mz_tol, 0)
    # peak mz range is covered by peakmap mz range we sel
    condition2 = pt.mz.in_range(t.mzmin, t.mzmax)
    # the peak range is covered by measure
    condition3 = (pt.rtmax >= t.rtmin) & (pt.rtmin <= t.rtmax)

    t = pt.left_join(t, condition1 & condition2 & condition3)

    drop_cols = [
        cn
        for cn in t.col_names
        if cn.endswith("__0") and cn.startswith("peakmap") == False
    ]
    t.drop_columns(*drop_cols)
    t.rename_postfixes(__0="")
    _replace_none_by_empty_peakmap(t)
    return t.copy()


def _update_mz_range(pt, mz_tol_abs=0, mz_tol_rel=0, **kwargs):
    mz_tol_rel *= 1e-6  # since mz_tol_rel is given in ppm
    pt.add_or_replace_column(
        "mzmax", pt.mz * (1 + mz_tol_rel) + mz_tol_abs, MzType, insert_after="mz"
    )
    pt.add_or_replace_column(
        "mzmin", pt.mz * (1 - mz_tol_rel) - mz_tol_abs, MzType, insert_after="mz"
    )


def _update_pm_ranges(t, mz_tol_abs, mz_tol_rel):
    t.add_column(
        "rtmin", t.apply(_get_rt_value, t.peakmap, 0, ignore_nones=False), RtType
    )
    t.add_column(
        "rtmax", t.apply(_get_rt_value, t.peakmap, 1, ignore_nones=False), RtType
    )
    t.add_column(
        "mzmin", t.apply(_get_mz_value, t.peakmap, 0, ignore_nones=False), MzType
    )
    t.replace_column("mzmin", t.mzmin * (1 - mz_tol_rel) - mz_tol_abs, MzType)
    t.add_column(
        "mzmax", t.apply(_get_mz_value, t.peakmap, 1, ignore_nones=False), MzType
    )
    t.replace_column("mzmax", t.mzmax * (1 + mz_tol_rel) + mz_tol_abs, MzType)


def _get_rt_value(peakmap, i):
    return peakmap.rt_range()[i]


def _get_mz_value(peakmap, i):
    return peakmap.mz_range()[i]


def _replace_none_by_empty_peakmap(t, mslevel=2):
    pm_empty = create_empty_peak_map(None, mslevel)
    # pm_empty_neg = create_empty_peak_map('-', mslevel)
    pms = [_replace(pm, pm_empty) for pm in t.peakmap]
    t.replace_column("peakmap", pms, PeakMap)


def _replace(pm, empty):
    return empty if pm is None else pm
    

def convert_peakmap_to_chromatogram(
    t, chromatogram_boundary_factor=3, ms_level=None, **kwargs
):
    # _replace_bounderies_by_peakmap_range(t)
    _enlarge_bouderies(t, chromatogram_boundary_factor)
    t = pp.extract_chromatograms(t, ms_level=ms_level)
    _update_original_rts(t)
    drop_cols = ("mzmin", "mzmax", "peakmap", "precursor_mz", "rtmin", "rtmax", 'rt_range')
    cols = [col for col in t.col_names if col in drop_cols]
    t.drop_columns(*cols)
    return t


def _enlarge_bouderies(t, chromatogram_boundary_factor=3):
    t.add_column("rtmin_", t.rtmin, RtType)
    t.add_column("rtmax_", t.rtmax, RtType)
    t.add_column("rt_range", t.apply(_get_rt_range, t.peakmap), object)
    # rtmin, rtmax = t.peakmap.unique_value().rt_range()
    t.replace_column(
        "rtmin", t.rtmin - chromatogram_boundary_factor * t.rt_max_shift, RtType
    )
    t.replace_column(
        "rtmin",
        (t.rtmin < t.apply(lambda v: v[0], t.rt_range)).then_else(
            t.apply(lambda v: v[0], t.rt_range), t.rtmin
        ),
        RtType,
    )
    t.replace_column(
        "rtmax", t.rtmax + chromatogram_boundary_factor * t.rt_max_shift, RtType
    )
    t.replace_column(
        "rtmax",
        (t.rtmax > t.apply(lambda v: v[1], t.rt_range)).then_else(
            t.apply(lambda v: v[1], t.rt_range), t.rtmax
        ),
        RtType,
    )
    
    
def _get_rt_range(pm):
    return pm.rt_range()


def _update_original_rts(t):
    t.replace_column("rtmin_chromatogram", t.rtmin_, RtType)
    t.replace_column("rtmax_chromatogram", t.rtmax_, RtType)
    t.rename_columns(rt="rt_chromatogram")
    t.drop_columns("rtmin_", "rtmax_")
    
    
def _update_id_column(t):
    if not "id" in t.col_names:
        t.add_enumeration()