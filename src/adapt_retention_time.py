# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 12:56:10 2023

@author: pkiefer
Concept
There are 2 approaches 
1) via fit (in that case the alignment peaks_table originates from a measurement)
2) from internal standard: The simplest approach: we select most important 
peak based on internal standard: Problems of mass isomers;
-> solutions 
The risk in case of several in

"""
import numpy as np
from emzed import RtType
from emzed.quantification import integrate_chromatograms


######################################################
# Local adaptation


from scipy.signal import savgol_filter


def local_adjust(t, min_fwhm=2.0, integration_algorithm="linear"):
    _determine_shift(t)
    _set_windows(t, min_fwhm, "rt_chromatogram")
    integrate_chromatograms(t, integration_algorithm, in_place=True)
    t.add_or_replace_column(
        "rt_chromatogram",
        t.apply(determine_rt, t.model_chromatogram, t.rt_chromatogram),
        RtType,
    )


def _determine_shift(t):
    t.add_column("lboundary", t.rtmin_chromatogram - t.rt_chromatogram, RtType)
    t.add_column("rboundary", t.rtmax_chromatogram - t.rt_chromatogram, RtType)
    _update_spectral_distance(t)
    t.replace_column(
        "rtmin_chromatogram",
        (t.rt_max_shift > 0).then_else(
            t.rt_chromatogram - t.rt_max_shift, t.rt_chromatogram - t.min_tol
        ),
        RtType,
    )
    t.replace_column(
        "rtmax_chromatogram",
        (t.rt_max_shift > 0).then_else(
            t.rt_chromatogram + t.rt_max_shift, t.rt_chromatogram + t.min_tol
        ),
        RtType,
    )
    integrate_chromatograms(t, "linear", in_place=True)
    t.add_column(
        "temp", t.apply(determine_rt, t.model_chromatogram, t.rt_chromatogram), RtType
    )
    t.replace_column("temp", t.temp.if_not_none_else(t.rt_chromatogram), RtType)
    t.replace_column("rtmin_chromatogram", t.lboundary + t.temp, RtType)
    t.replace_column("rtmax_chromatogram", t.rboundary + t.temp, RtType)
    integrate_chromatograms(t, "linear", in_place=True)
    t.drop_columns("lboundary", "rboundary", "temp", "min_tol")


def _update_spectral_distance(t):
    t.add_column(
        "min_tol",
        t.apply(
            _max_delta_spec,
            t.chromatogram,
            t.rtmin_chromatogram,
            t.rtmax_chromatogram,
            t.rt_max_shift,
            ignore_nones=False,
        ),
        RtType,
    )


def _max_delta_spec(*kwargs):
    epsilon = 1e-2  # value to increase numerical stability but with no significant influence on
    # rt shift
    if not _check_params(kwargs):
        return epsilon
    chrom, rtmin, rtmax, rttol = kwargs
    rts = chrom.rts
    peak_range = rts[(rts >= rtmin - rttol) & (rts <= rtmax + rttol)]
    deltas = np.diff(peak_range)
    return max(deltas) + epsilon if len(deltas) else epsilon


def _check_params(kwargs):
    chrom, rtmin, rtmax, rttol = kwargs
    cond1 = len(chrom) if chrom is not None else False
    cond2 = [v is not None for v in [rtmin, rtmax, rttol]]
    return cond1 and cond2


def _set_windows(t, min_fwhm=2.0, rt_value_column=None):
    params = [t.model_chromatogram, min_fwhm]
    params.append(t[rt_value_column]) if rt_value_column else params.append(None)
    t.add_column("temp", t.apply(adapt_window, *params), object)
    t.add_column("_rtmin", t.apply(lambda v: v[0], t.temp), RtType)
    t.add_column("_rtmax", t.apply(lambda v: v[1], t.temp), RtType)
    t.replace_column(
        "rtmin_chromatogram",
        t._rtmin.is_not_none().then_else(t._rtmin, t.rtmin_chromatogram),
        RtType,
    )
    t.replace_column(
        "rtmax_chromatogram",
        t._rtmax.is_not_none().then_else(t._rtmax, t.rtmax_chromatogram),
        RtType,
    )

    t.drop_columns("temp", "_rtmin", "_rtmax")


def adapt_window(model, min_fwhm, rt=None):
    rts, ints = model.graph()
    if not len(rts):
        return None, None
    smoothed = _get_smoothed(ints)
    rt_ = _determine_rt(rts, smoothed)
    rt = rt if rt_ is None else rt_
    if rt is None:
        return None, None
    fwhm = _determine_fwhm(rts, smoothed, rt, min_fwhm)
    f_asym = _determine_asymmetric_factor(rts, smoothed, rt)
    left_pos = np.argwhere(np.logical_and(rts >= rt - 2 * f_asym * fwhm, rts < rt))
    right_pos = np.argwhere(np.logical_and(rts <= rt + 2 / f_asym * fwhm, rts > rt))
    start = _min_value_pos(left_pos, smoothed, False)
    start = start if start else 0
    stop = _min_value_pos(right_pos, smoothed)
    stop = stop if stop else len(rts) - 1
    return rts[start], rts[stop]


def _min_value_pos(positions, intensities, rightside=True):

    if len(positions):
        i = np.where(intensities[positions] == min(intensities[positions]))
        index = 0 if rightside else -1
        return positions[i][index]


def determine_rt(model, rt):
    # rtmin, rtmax = limits
    rts, ints = model.graph()
    if len(ints):
        smoothed = _get_smoothed(ints)
        rt = _determine_rt(rts, smoothed)
    return rt


def _get_smoothed(values):
    try:
        smoothed = savgol_filter(values, 7, 3)
        # we exclude negative intensity values
        smoothed[smoothed < 0] = 0
        return smoothed
    except:
        return values


def _determine_rt(rts, smoothed):
    i = np.where(smoothed == max(smoothed))[0]
    i = int((min(i) + max(i)) / 2)
    return rts[i]


def _determine_asymmetric_factor(rts, smoothed, rt):
    max_int = _determine_max_int(rt, rts, smoothed)
    pos = np.where(np.logical_and(smoothed >= 0.2 * max_int, smoothed <= max_int))[0]
    diff_pos = np.diff(pos)
    check = np.where(diff_pos > 1)[0] + 1
    blocks = np.split(pos, check)
    i = np.where(abs(rts - rt) == np.min(abs(rts - rt)))[0]
    for block in blocks:
        if i in block:
            left = len(block[block <= i])
            right = len(block[block >= i])
            return left / right
    return 1.0


def _determine_fwhm(rts, smoothed, rt, min_fwhm=2.0):
    fwhm = min_fwhm
    max_int = _determine_max_int(rt, rts, smoothed)
    pos = np.where(np.logical_and(smoothed >= 0.5 * max_int, smoothed <= max_int))[0]
    diff_pos = np.diff(pos)
    check = np.where(diff_pos > 1)[0] + 1
    blocks = np.split(pos, check)
    i = np.where(abs(rts - rt) == np.min(abs(rts - rt)))[0]
    for block in blocks:
        if i in block:
            fwhm = rts[max(block)] - rts[min(block)]
    return fwhm if fwhm > min_fwhm else min_fwhm


def _determine_max_int(rt, rts, ints):
    diffs = abs(rts - rt)
    pos = np.where(diffs == min(diffs))
    return max(ints[pos])


