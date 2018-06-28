#! /usr/bin/env python

import sys
import os
import re
import math
import glob
import logging

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
_LOG = logging.getLogger(os.path.basename(__file__))

import pycoevolity
import project_util

import matplotlib as mpl

# Use TrueType (42) fonts rather than Type 3 fonts
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
tex_font_settings = {
        "text.usetex": True,
        "font.family": "sans-serif",
        # "font.serif": [
        #         "Computer Modern Roman",
        #         "Times",
        #         ],
        # "font.sans-serif": [
        #         "Computer Modern Sans serif",
        #         "Helvetica",
        #         ],
        # "font.cursive": [
        #         "Zapf Chancery",
        #         ],
        # "font.monospace": [
        #         "Computer Modern Typewriter",
        #         "Courier",
        #         ],
        "text.latex.preamble" : [
                "\\usepackage[T1]{fontenc}",
                "\\usepackage[cm]{sfmath}",
                ]
}

mpl.rcParams.update(tex_font_settings)

import matplotlib.pyplot as plt
from matplotlib import gridspec

def get_nevents_probs(
        results_paths,
        nevents = 1):
    probs = []
    prob_key = "num_events_{0}_p".format(nevents)
    for d in pycoevolity.parsing.spreadsheet_iter(results_paths):
        probs.append((
                float(d[prob_key]),
                int(int(d["true_num_events"]) == nevents)
                ))
    return probs

def bin_prob_correct_tuples(probability_correct_tuples, nbins = 20):
    bin_upper_limits = list(get_sequence_iter(0.0, 1.0, nbins+1))[1:]
    bin_width = (bin_upper_limits[1] - bin_upper_limits[0]) / 2.0
    bins = [[] for i in range(nbins)]
    n = 0
    for (p, t) in probability_correct_tuples:
        n += 1
        binned = False
        for i, l in enumerate(bin_upper_limits):
            if p < l:
                bins[i].append((p, t))
                binned = True
                break
        if not binned:
            bins[i].append((p, t))
    total = 0
    for b in bins:
        total += len(b)
    assert total == n
    assert len(bins) == nbins
    est_true_tups = []
    for i, b in enumerate(bins):
        ests = [p for (p, t) in b]
        est = sum(ests) / float(len(ests))
        correct = [t for (p, t) in b]
        true = sum(correct) / float(len(correct))
        est_true_tups.append((est, true))
    return bins, est_true_tups

def get_nevents_estimated_true_probs(
        results_paths,
        nevents = 1,
        nbins = 20):
    nevent_probs = get_nevents_probs(
            results_paths = results_paths,
            nevents = nevents)
    _LOG.info("\tparsed results for {0} simulations".format(len(nevent_probs)))
    bins, tups = bin_prob_correct_tuples(nevent_probs, nbins = nbins)
    _LOG.info("\tbin sample sizes: {0}".format(
            ", ".join(str(len(b)) for b in bins)
            ))
    return bins, tups

def plot_nevents_estimated_vs_true_probs(
        results_paths,
        nevents = 1,
        nbins = 20,
        plot_file_prefix = ""):
    bins, est_true_probs = get_nevents_estimated_true_probs(
            results_paths = results_paths,
            nevents = nevents,
            nbins = nbins)

    plt.close('all')
    fig = plt.figure(figsize = (4.0, 3.5))
    ncols = 1
    nrows = 1
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    ax = plt.subplot(gs[0, 0])
    x = [e for (e, t) in est_true_probs]
    y = [t for (e, t) in est_true_probs]
    sample_sizes = [len(b) for b in bins]
    line, = ax.plot(x, y)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 3.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    for i, (label, lx, ly) in enumerate(zip(sample_sizes, x, y)):
        if i == 0:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (1, 1),
                    textcoords = "offset points",
                    horizontalalignment = "left",
                    verticalalignment = "bottom")
        elif i == len(x) - 1:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (-1, -1),
                    textcoords = "offset points",
                    horizontalalignment = "right",
                    verticalalignment = "top")
        else:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (-1, 1),
                    textcoords = "offset points",
                    horizontalalignment = "right",
                    verticalalignment = "bottom")
    ylabel_text = ax.set_ylabel("True probability", size = 14.0)
    ax.text(0.5, -0.14,
            "Posterior probability of one divergence",
            horizontalalignment = "center",
            verticalalignment = "top",
            size = 14.0)
    identity_line, = ax.plot(
            [0.0, 1.0],
            [0.0, 1.0])
    plt.setp(identity_line,
            color = '0.8',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)

    gs.update(left = 0.10, right = 0.995, bottom = 0.18, top = 0.91)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_file_name = "est-vs-true-prob-nevent-1.pdf"
    if plot_file_prefix:
        plot_file_name = plot_file_prefix + "-" + plot_file_name
    plot_path = os.path.join(plot_dir,
            plot_file_name)
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def get_sequence_iter(start = 0.0, stop = 1.0, n = 10):
    assert(stop > start)
    step = (stop - start) / float(n - 1)
    return ((start + (i * step)) for i in range(n))

def truncate_color_map(cmap, min_val = 0.0, max_val = 10, n = 100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(
                    n = cmap.name,
                    a = min_val,
                    b = max_val),
            cmap(list(get_sequence_iter(min_val, max_val, n))))
    return new_cmap

def get_errors(values, lowers, uppers):
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def ci_width_iter(results, parameter_str):
    n = len(results["eti_95_upper_{0}".format(parameter_str)])
    for i in range(n):
        upper = float(results["eti_95_upper_{0}".format(parameter_str)][i])
        lower = float(results["eti_95_lower_{0}".format(parameter_str)][i])
        yield upper - lower

def absolute_error_iter(results, parameter_str):
    n = len(results["true_{0}".format(parameter_str)])
    for i in range(n):
        t = float(results["true_{0}".format(parameter_str)][i])
        e = float(results["mean_{0}".format(parameter_str)][i])
        yield math.fabs(t - e)


def plot_ess_versus_error(
        parameters,
        results_paths,
        column_labels = None,
        row_labels = None,
        parameter_label = "divergence time",
        plot_file_prefix = None):
    _LOG.info("Generating ESS vs CI scatter plots for {0}...".format(parameter_label))

    assert(len(parameters) == len(set(parameters)))
    if row_labels:
        assert len(row_labels) ==  len(results_paths)
    if column_labels:
        assert len(column_labels) == len(results_paths[0])

    nrows = len(results_paths)
    ncols = len(results_paths[0])

    if not plot_file_prefix:
        plot_file_prefix = parameters[0] 
    plot_file_prefix_ci = plot_file_prefix + "-ess-vs-ci-width"
    plot_file_prefix_error = plot_file_prefix + "-ess-vs-error"

    # Very inefficient, but parsing all results to get min/max for parameter
    ess_min = float('inf')
    ess_max = float('-inf')
    ci_width_min = float('inf')
    ci_width_max = float('-inf')
    error_min = float('inf')
    error_max = float('-inf')
    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)
            for parameter_str in parameters:
                ci_widths = tuple(ci_width_iter(results, parameter_str))
                errors = tuple(absolute_error_iter(results, parameter_str))
                ess_min = min(ess_min,
                        min(float(x) for x in results["ess_sum_{0}".format(parameter_str)]))
                ess_max = max(ess_max,
                        max(float(x) for x in results["ess_sum_{0}".format(parameter_str)]))
                ci_width_min = min(ci_width_min, min(ci_widths))
                ci_width_max = max(ci_width_max, max(ci_widths))
                error_min = min(error_min, min(errors))
                error_max = max(error_max, max(errors))
    ess_axis_buffer = math.fabs(ess_max - ess_min) * 0.05
    ess_axis_min = ess_min - ess_axis_buffer
    ess_axis_max = ess_max + ess_axis_buffer
    ci_width_axis_buffer = math.fabs(ci_width_max - ci_width_min) * 0.05
    ci_width_axis_min = ci_width_min - ci_width_axis_buffer
    ci_width_axis_max = ci_width_max + ci_width_axis_buffer
    error_axis_buffer = math.fabs(error_max - error_min) * 0.05
    error_axis_min = error_min - error_axis_buffer
    error_axis_max = error_max + error_axis_buffer

    plt.close('all')
    w = 1.6
    h = 1.5
    fig_width = (ncols * w) + 1.0
    fig_height = (nrows * h) + 0.7
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            x = []
            y = []
            for parameter_str in parameters:
                x.extend(float(x) for x in results["ess_sum_{0}".format(parameter_str)])
                y.extend(ci_width_iter(results, parameter_str))

            assert(len(x) == len(y))
            ax = plt.subplot(gs[row_index, column_index])
            line, = ax.plot(x, y)
            plt.setp(line,
                    marker = 'o',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    linestyle = '',
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(ess_axis_min, ess_axis_max)
            ax.set_ylim(ci_width_axis_min, ci_width_axis_max)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "Effective sample size of {0}".format(parameter_label),
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            "CI width {0}".format(parameter_label),
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = 0.08, right = 0.98, bottom = 0.08, top = 0.97)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix_ci))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


    _LOG.info("Generating ESS vs error scatter plots for {0}...".format(parameter_label))
    plt.close('all')
    w = 1.6
    h = 1.5
    fig_width = (ncols * w) + 1.0
    fig_height = (nrows * h) + 0.7
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            x = []
            y = []
            for parameter_str in parameters:
                x.extend(float(x) for x in results["ess_sum_{0}".format(parameter_str)])
                y.extend(absolute_error_iter(results, parameter_str))
                

            assert(len(x) == len(y))
            ax = plt.subplot(gs[row_index, column_index])
            line, = ax.plot(x, y)
            plt.setp(line,
                    marker = 'o',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    linestyle = '',
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(ess_axis_min, ess_axis_max)
            ax.set_ylim(error_axis_min, error_axis_max)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "Effective sample size of {0}".format(parameter_label),
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            "Absolute error of {0}".format(parameter_label),
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = 0.08, right = 0.98, bottom = 0.08, top = 0.97)

    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix_error))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))



def generate_scatter_plots(
        parameters,
        results_paths,
        column_labels = None,
        row_labels = None,
        parameter_label = "divergence time",
        parameter_symbol = "\\tau",
        plot_width = 1.9,
        plot_height = 1.8,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        x_label_size = 18.0,
        y_label_size = 18.0,
        y_label = None,
        plot_file_prefix = None):
    _LOG.info("Generating scatter plots for {0}...".format(parameter_label))

    assert(len(parameters) == len(set(parameters)))
    if row_labels:
        assert len(row_labels) ==  len(results_paths)
    if column_labels:
        assert len(column_labels) == len(results_paths[0])

    nrows = len(results_paths)
    ncols = len(results_paths[0])

    if not plot_file_prefix:
        plot_file_prefix = parameters[0] 

    # Very inefficient, but parsing all results to get min/max for parameter
    parameter_min = float('inf')
    parameter_max = float('-inf')
    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)
            for parameter_str in parameters:
                parameter_min = min(parameter_min,
                        min(float(x) for x in results["true_{0}".format(parameter_str)]))
                parameter_max = max(parameter_max,
                        max(float(x) for x in results["true_{0}".format(parameter_str)]))
                parameter_min = min(parameter_min,
                        min(float(x) for x in results["mean_{0}".format(parameter_str)]))
                parameter_max = max(parameter_max,
                        max(float(x) for x in results["mean_{0}".format(parameter_str)]))
    axis_buffer = math.fabs(parameter_max - parameter_min) * 0.05
    axis_min = parameter_min - axis_buffer
    axis_max = parameter_max + axis_buffer

    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            x = []
            y = []
            y_upper = []
            y_lower = []
            for parameter_str in parameters:
                x.extend(float(x) for x in results["true_{0}".format(parameter_str)])
                y.extend(float(x) for x in results["mean_{0}".format(parameter_str)])
                y_lower.extend(float(x) for x in results["eti_95_lower_{0}".format(parameter_str)])
                y_upper.extend(float(x) for x in results["eti_95_upper_{0}".format(parameter_str)])

            assert(len(x) == len(y))
            assert(len(x) == len(y_lower))
            assert(len(x) == len(y_upper))
            proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
                    x,
                    y_lower,
                    y_upper)
            rmse = pycoevolity.stats.root_mean_square_error(x, y)
            ax = plt.subplot(gs[row_index, column_index])
            line = ax.errorbar(
                    x = x,
                    y = y,
                    yerr = get_errors(y, y_lower, y_upper),
                    ecolor = '0.65',
                    elinewidth = 0.5,
                    capsize = 0.8,
                    barsabove = False,
                    marker = 'o',
                    linestyle = '',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(axis_min, axis_max)
            ax.set_ylim(axis_min, axis_max)
            identity_line, = ax.plot(
                    [axis_min, axis_max],
                    [axis_min, axis_max])
            plt.setp(identity_line,
                    color = '0.7',
                    linestyle = '-',
                    linewidth = 1.0,
                    marker = '',
                    zorder = 0)
            ax.text(0.02, 0.97,
                    "\\scriptsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                            parameter_symbol,
                            proportion_within_ci),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    size = 6.0,
                    zorder = 200)
            ax.text(0.02, 0.87,
                    # "\\scriptsize\\noindent$\\textrm{{\\sffamily RMSE}} = {0:.2e}$".format(
                    "\\scriptsize\\noindent RMSE = {0:.2e}".format(
                            rmse),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    size = 6.0,
                    zorder = 200)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "True {0} (${1}$)".format(parameter_label, parameter_symbol),
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = x_label_size)
    if y_label is None:
        y_label = "Estimated {0} ($\\hat{{{1}}}$)".format(parameter_label, parameter_symbol)
    fig.text(0.005, 0.5,
            y_label,
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = y_label_size)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def generate_specific_scatter_plot(
        parameters,
        results_paths,
        parameter_label = "divergence time",
        parameter_symbol = "\\tau",
        plot_title = None,
        include_x_label = True,
        include_y_label = True,
        include_rmse = True,
        include_ci = True,
        plot_width = 3.5,
        plot_height = 3.0,
        xy_label_size = 16.0,
        title_size = 16.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        plot_file_prefix = None,
        variable_only = False):
    _LOG.info("Generating scatter plots for {0}...".format(parameter_label))

    assert(len(parameters) == len(set(parameters)))
    if not plot_file_prefix:
        plot_file_prefix = parameters[0] 

    # Very inefficient, but parsing all results to get min/max for parameter
    parameter_min = float('inf')
    parameter_max = float('-inf')
    results = pycoevolity.parsing.get_dict_from_spreadsheets(
            results_paths,
            sep = "\t",
            offset = 0)
    for parameter_str in parameters:
        parameter_min = min(parameter_min,
                min(float(x) for x in results["true_{0}".format(parameter_str)]))
        parameter_max = max(parameter_max,
                max(float(x) for x in results["true_{0}".format(parameter_str)]))
        parameter_min = min(parameter_min,
                min(float(x) for x in results["mean_{0}".format(parameter_str)]))
        parameter_max = max(parameter_max,
                max(float(x) for x in results["mean_{0}".format(parameter_str)]))
    axis_buffer = math.fabs(parameter_max - parameter_min) * 0.05
    axis_min = parameter_min - axis_buffer
    axis_max = parameter_max + axis_buffer

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    x = []
    y = []
    y_upper = []
    y_lower = []
    for parameter_str in parameters:
        x.extend(float(x) for x in results["true_{0}".format(parameter_str)])
        y.extend(float(x) for x in results["mean_{0}".format(parameter_str)])
        y_lower.extend(float(x) for x in results["eti_95_lower_{0}".format(parameter_str)])
        y_upper.extend(float(x) for x in results["eti_95_upper_{0}".format(parameter_str)])

    assert(len(x) == len(y))
    assert(len(x) == len(y_lower))
    assert(len(x) == len(y_upper))
    proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
            x,
            y_lower,
            y_upper)
    rmse = pycoevolity.stats.root_mean_square_error(x, y)
    ax = plt.subplot(gs[0, 0])
    line = ax.errorbar(
            x = x,
            y = y,
            yerr = get_errors(y, y_lower, y_upper),
            ecolor = '0.65',
            elinewidth = 0.5,
            capsize = 0.8,
            barsabove = False,
            marker = 'o',
            linestyle = '',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 2.5,
            zorder = 100,
            rasterized = True)
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    identity_line, = ax.plot(
            [axis_min, axis_max],
            [axis_min, axis_max])
    plt.setp(identity_line,
            color = '0.7',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    if include_ci:
        ax.text(0.02, 0.97,
                "\\normalsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                        parameter_symbol,
                        proportion_within_ci),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 200)
    if include_rmse:
        ax.text(0.02, 0.87,
                "\\normalsize\\noindent RMSE = {0:.2e}".format(
                        rmse),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 200)
    if include_x_label:
        ax.set_xlabel(
                "True {0} (${1}$)".format(parameter_label, parameter_symbol),
                fontsize = xy_label_size)
    if include_y_label:
        ax.set_ylabel(
                "Estimated {0} ($\\hat{{{1}}}$)".format(parameter_label, parameter_symbol),
                fontsize = xy_label_size)
    root_shape, root_scale = get_root_gamma_parameters(root_shape_setting, root_scale_setting)
    if plot_title:
        ax.set_title(plot_title,
                fontsize = title_size)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_histograms(
        parameters,
        results_paths,
        column_labels = None,
        row_labels = None,
        parameter_label = "Number of variable sites",
        parameter_discrete = True,
        range_key = "range",
        number_of_digits = 0,
        plot_width = 1.9,
        plot_height = 1.8,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        plot_file_prefix = None):
    _LOG.info("Generating histograms for {0}...".format(parameter_label))

    assert(len(parameters) == len(set(parameters)))
    if row_labels:
        assert len(row_labels) ==  len(results_paths)
    if column_labels:
        assert len(column_labels) == len(results_paths[0])

    nrows = len(results_paths)
    ncols = len(results_paths[0])

    if not plot_file_prefix:
        plot_file_prefix = parameters[0] 
    root_alpha_pattern = re.compile(r'root-(?P<alpha_setting>\S+)-\d00k')
    locus_size_pattern = re.compile(r'root-\d+-\d00k-(?P<locus_size>\d+)u?l')

    # Very inefficient, but parsing all results to get min/max for parameter
    parameter_min = float('inf')
    parameter_max = float('-inf')
    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)
            for parameter_str in parameters:
                parameter_min = min(parameter_min,
                        min(float(x) for x in results["{0}".format(parameter_str)]))
                parameter_max = max(parameter_max,
                        max(float(x) for x in results["{0}".format(parameter_str)]))

    axis_buffer = math.fabs(parameter_max - parameter_min) * 0.05
    axis_min = parameter_min - axis_buffer
    axis_max = parameter_max + axis_buffer

    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    hist_bins = None
    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            x = []
            for parameter_str in parameters:
                if parameter_discrete:
                    x.extend(int(x) for x in results["{0}".format(parameter_str)])
                else:
                    x.extend(float(x) for x in results["{0}".format(parameter_str)])

            summary = pycoevolity.stats.get_summary(x)
            _LOG.info("0.025, 0.975 quantiles: {0:.2f}, {1:.2f}".format(
                    summary["qi_95"][0],
                    summary["qi_95"][1]))

            x_range = (parameter_min, parameter_max)
            if parameter_discrete:
                x_range = (int(parameter_min), int(parameter_max))
            ax = plt.subplot(gs[row_index, column_index])
            n, bins, patches = ax.hist(x,
                    # normed = True,
                    weights = [1.0 / float(len(x))] * len(x),
                    bins = hist_bins,
                    range = x_range,
                    cumulative = False,
                    histtype = 'bar',
                    align = 'mid',
                    orientation = 'vertical',
                    rwidth = None,
                    log = False,
                    color = None,
                    edgecolor = '0.5',
                    facecolor = '0.5',
                    fill = True,
                    hatch = None,
                    label = None,
                    linestyle = None,
                    linewidth = None,
                    zorder = 10,
                    )
            if hist_bins is None:
                hist_bins = bins
            ax.text(0.98, 0.98,
                    "\\scriptsize {mean:,.{ndigits}f} ({lower:,.{ndigits}f}--{upper:,.{ndigits}f})".format(
                            # int(round(summary["mean"])),
                            # int(round(summary[range_key][0])),
                            # int(round(summary[range_key][1]))),
                            mean = summary["mean"],
                            lower = summary[range_key][0],
                            upper = summary[range_key][1],
                            ndigits = number_of_digits),
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    zorder = 200)

            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # make sure y-axis is the same
    y_max = float('-inf')
    all_axes = fig.get_axes()
    for ax in all_axes:
        ymn, ymx = ax.get_ylim()
        y_max = max(y_max, ymx)
    for ax in all_axes:
        ax.set_ylim(0.0, y_max)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            parameter_label,
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            # "Density",
            "Frequency",
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-histograms.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_model_plots(
        results_paths,
        column_labels = None,
        row_labels = None,
        number_of_comparisons = 5,
        plot_width = 1.6,
        plot_height = 1.5,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        y_label_size = 18.0,
        y_label = None,
        number_font_size = 12.0,
        plot_file_prefix = None):
    _LOG.info("Generating model plots...")

    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)

    if row_labels:
        assert len(row_labels) ==  len(results_paths)
    if column_labels:
        assert len(column_labels) == len(results_paths[0])

    nrows = len(results_paths)
    ncols = len(results_paths[0])

    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_paths_row in enumerate(results_paths):
        for column_index, paths in enumerate(results_paths_row):
            results = pycoevolity.parsing.get_dict_from_spreadsheets(
                    paths,
                    sep = "\t",
                    offset = 0)

            true_map_nevents = []
            true_map_nevents_probs = []
            for i in range(number_of_comparisons):
                true_map_nevents.append([0 for i in range(number_of_comparisons)])
                true_map_nevents_probs.append([[] for i in range(number_of_comparisons)])
            true_nevents = tuple(int(x) for x in results["true_num_events"])
            map_nevents = tuple(int(x) for x in results["map_num_events"])
            true_nevents_cred_levels = tuple(float(x) for x in results["true_num_events_cred_level"])
            true_model_cred_levels = tuple(float(x) for x in results["true_model_cred_level"])
            assert(len(true_nevents) == len(map_nevents))
            assert(len(true_nevents) == len(true_nevents_cred_levels))
            assert(len(true_nevents) == len(true_model_cred_levels))

            true_nevents_probs = []
            map_nevents_probs = []
            for i in range(len(true_nevents)):
                true_nevents_probs.append(float(
                    results["num_events_{0}_p".format(true_nevents[i])][i]))
                map_nevents_probs.append(float(
                    results["num_events_{0}_p".format(map_nevents[i])][i]))
            assert(len(true_nevents) == len(true_nevents_probs))
            assert(len(true_nevents) == len(map_nevents_probs))

            mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
            median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

            nevents_within_95_cred = 0
            model_within_95_cred = 0
            ncorrect = 0
            for i in range(len(true_nevents)):
                true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
                true_map_nevents_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
                if true_nevents_cred_levels[i] <= 0.95:
                    nevents_within_95_cred += 1
                if true_model_cred_levels[i] <= 0.95:
                    model_within_95_cred += 1
                if true_nevents[i] == map_nevents[i]:
                    ncorrect += 1
            p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
            p_model_within_95_cred = model_within_95_cred / float(len(true_nevents))
            p_correct = ncorrect / float(len(true_nevents))

            _LOG.info("p(nevents within CS) = {0:.4f}".format(p_nevents_within_95_cred))
            _LOG.info("p(model within CS) = {0:.4f}".format(p_model_within_95_cred))
            ax = plt.subplot(gs[row_index, column_index])

            ax.imshow(true_map_nevents,
                    origin = 'lower',
                    cmap = cmap,
                    interpolation = 'none',
                    aspect = 'auto'
                    # extent = [0.5, 3.5, 0.5, 3.5]
                    )
            for i, row_list in enumerate(true_map_nevents):
                for j, num_events in enumerate(row_list):
                    ax.text(j, i,
                            str(num_events),
                            horizontalalignment = "center",
                            verticalalignment = "center",
                            size = number_font_size)
            ax.text(0.98, 0.02,
                    "\\scriptsize$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                            p_nevents_within_95_cred),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
            ax.text(0.02, 0.98,
                    "\\scriptsize$p(\\hat{{k}} = k) = {0:.3f}$".format(
                            p_correct),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes)
            ax.text(0.98, 0.98,
                    "\\scriptsize$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
                            median_true_nevents_prob),
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        # Make sure ticks correspond only with number of events
        ax.xaxis.set_ticks(range(number_of_comparisons))
        ax.yaxis.set_ticks(range(number_of_comparisons))
        if ax.is_last_row() and ax.is_first_col():
            xtick_labels = [item for item in ax.get_xticklabels()]
            for i in range(len(xtick_labels)):
                xtick_labels[i].set_text(str(i + 1))
            ytick_labels = [item for item in ax.get_yticklabels()]
            for i in range(len(ytick_labels)):
                ytick_labels[i].set_text(str(i + 1))
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)
        else:
            xtick_labels = ["" for item in ax.get_xticklabels()]
            ytick_labels = ["" for item in ax.get_yticklabels()]
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "True number of events ($k$)",
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    if y_label is None:
        y_label = "Estimated number of events ($\\hat{{k}}$)"
    fig.text(0.005, 0.5,
            y_label,
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    if plot_file_prefix:
        plot_path = os.path.join(plot_dir,
                "{0}-nevents.pdf".format(plot_file_prefix))
    else:
        plot_path = os.path.join(plot_dir,
                "nevents.pdf")
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def generate_specific_model_plots(
        results_paths,
        number_of_comparisons = 5,
        plot_title = None,
        include_x_label = True,
        include_y_label = True,
        include_median = True,
        include_cs = True,
        include_prop_correct = True,
        plot_width = 3.5,
        plot_height = 3.0,
        xy_label_size = 16.0,
        title_size = 16.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        lower_annotation_y = 0.02,
        upper_annotation_y = 0.92,
        plot_file_prefix = None,
        variable_only = False):
    _LOG.info("Generating model plots...")

    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    results = pycoevolity.parsing.get_dict_from_spreadsheets(
            results_paths,
            sep = "\t",
            offset = 0)

    true_map_nevents = []
    true_map_nevents_probs = []
    for i in range(number_of_comparisons):
        true_map_nevents.append([0 for i in range(number_of_comparisons)])
        true_map_nevents_probs.append([[] for i in range(number_of_comparisons)])
    true_nevents = tuple(int(x) for x in results["true_num_events"])
    map_nevents = tuple(int(x) for x in results["map_num_events"])
    true_nevents_cred_levels = tuple(float(x) for x in results["true_num_events_cred_level"])
    true_model_cred_levels = tuple(float(x) for x in results["true_model_cred_level"])
    assert(len(true_nevents) == len(map_nevents))
    assert(len(true_nevents) == len(true_nevents_cred_levels))
    assert(len(true_nevents) == len(true_model_cred_levels))

    true_nevents_probs = []
    map_nevents_probs = []
    for i in range(len(true_nevents)):
        true_nevents_probs.append(float(
            results["num_events_{0}_p".format(true_nevents[i])][i]))
        map_nevents_probs.append(float(
            results["num_events_{0}_p".format(map_nevents[i])][i]))
    assert(len(true_nevents) == len(true_nevents_probs))
    assert(len(true_nevents) == len(map_nevents_probs))

    mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
    median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

    nevents_within_95_cred = 0
    model_within_95_cred = 0
    ncorrect = 0
    for i in range(len(true_nevents)):
        true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
        true_map_nevents_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
        if true_nevents_cred_levels[i] <= 0.95:
            nevents_within_95_cred += 1
        if true_model_cred_levels[i] <= 0.95:
            model_within_95_cred += 1
        if true_nevents[i] == map_nevents[i]:
            ncorrect += 1
    p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
    p_model_within_95_cred = model_within_95_cred / float(len(true_nevents))
    p_correct = ncorrect / float(len(true_nevents))

    _LOG.info("p(nevents within CS) = {0:.4f}".format(p_nevents_within_95_cred))
    _LOG.info("p(model within CS) = {0:.4f}".format(p_model_within_95_cred))
    ax = plt.subplot(gs[0, 0])

    ax.imshow(true_map_nevents,
            origin = 'lower',
            cmap = cmap,
            interpolation = 'none',
            aspect = 'auto'
            )
    for i, row_list in enumerate(true_map_nevents):
        for j, num_events in enumerate(row_list):
            ax.text(j, i,
                    str(num_events),
                    horizontalalignment = "center",
                    verticalalignment = "center")
    if include_cs:
        ax.text(0.98, lower_annotation_y,
                "$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                        p_nevents_within_95_cred),
                horizontalalignment = "right",
                verticalalignment = "bottom",
                transform = ax.transAxes)
    if include_prop_correct:
        ax.text(0.02, upper_annotation_y,
                "$p(\\hat{{k}} = k) = {0:.3f}$".format(
                        p_correct),
                horizontalalignment = "left",
                verticalalignment = "bottom",
                transform = ax.transAxes)
    if include_median:
        ax.text(0.98, upper_annotation_y,
                "$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
                        median_true_nevents_prob),
                horizontalalignment = "right",
                verticalalignment = "bottom",
                transform = ax.transAxes)
    if include_x_label:
        ax.set_xlabel("True \\# of events ($k$)",
                # labelpad = 8.0,
                fontsize = xy_label_size)
    if include_y_label:
        ax.set_ylabel("Estimated \\# of events ($\\hat{{k}}$)",
                labelpad = 8.0,
                fontsize = xy_label_size)
    root_shape, root_scale = get_root_gamma_parameters(root_shape_setting, root_scale_setting)
    if plot_title:
        ax.set_title(plot_title,
                fontsize = title_size)

    # Make sure ticks correspond only with number of events
    ax.xaxis.set_ticks(range(number_of_comparisons))
    ax.yaxis.set_ticks(range(number_of_comparisons))
    xtick_labels = [item for item in ax.get_xticklabels()]
    for i in range(len(xtick_labels)):
        xtick_labels[i].set_text(str(i + 1))
    ytick_labels = [item for item in ax.get_yticklabels()]
    for i in range(len(ytick_labels)):
        ytick_labels[i].set_text(str(i + 1))
    ax.set_xticklabels(xtick_labels)
    ax.set_yticklabels(ytick_labels)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    plot_dir = os.path.join(project_util.VAL_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-nevents.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def get_msbayes_results(
        true_path,
        sim_dir,
        number_of_pairs = 5,
        number_of_sims = 500,
        posterior_sample_size = 2000,
        prior_sample_size = 600000):
    results = {
            'divergence time': {'true': [], 'mean': [], 'lower': [], 'upper': []},
            'root population size': {'true': [], 'mean': [], 'lower': [], 'upper': []},
            'leaf population size': {'true': [], 'mean': [], 'lower': [], 'upper': []},
            'nevents': {'true': [], 'map': [], 'cred_level': [], 'true_prob': [], 'map_prob': []},
            }
    column_header_prefixes = {
            'divergence time': ("PRI.t.",),
            'root population size': ("PRI.aTheta.",),
            'leaf population size': ("PRI.d1Theta.", "PRI.d2Theta."),
            }
    true_values = pycoevolity.parsing.get_dict_from_spreadsheets(
            [true_path],
            sep = '\t')
    nsims = len(true_values.values()[0])
    assert (nsims == number_of_sims)
    for sim_idx in range(nsims):
        posterior_path = os.path.join(sim_dir,
                "d1-m1-s{sim_num}-{prior_sample_size}-posterior-sample.txt.gz".format(
                        sim_num = sim_idx + 1,
                        prior_sample_size = prior_sample_size))
        _LOG.info("Parsing {0}".format(posterior_path))
        posterior = pycoevolity.parsing.get_dict_from_spreadsheets(
                [posterior_path],
                sep = '\t')
        for pair_idx in range(number_of_pairs):
            for parameter_key, header_prefixes in column_header_prefixes.items():
                for header_prefix in header_prefixes:
                    header = "{0}{1}".format(header_prefix, pair_idx + 1)
                    true_val = float(true_values[header][sim_idx])
                    post_sum = pycoevolity.stats.get_summary(float(x) for x in posterior[header])
                    results[parameter_key]['true'].append(true_val)
                    results[parameter_key]['mean'].append(post_sum['mean'])
                    results[parameter_key]['lower'].append(post_sum['qi_95'][0])
                    results[parameter_key]['upper'].append(post_sum['qi_95'][1])
        true_nevents = int(true_values["PRI.Psi"][sim_idx])
        nevent_freqs = pycoevolity.stats.get_freqs(int(x) for x in posterior["PRI.Psi"])
        sorted_nevent_freqs = sorted(nevent_freqs.items(), reverse = True,
                key = lambda x: x[1])
        map_nevents = sorted_nevent_freqs[0][0]
        true_nevents_prob = nevent_freqs[true_nevents]
        map_nevents_prob = nevent_freqs[map_nevents]
        results['nevents']['true'].append(true_nevents)
        results['nevents']['map'].append(map_nevents)
        results['nevents']['true_prob'].append(true_nevents_prob)
        results['nevents']['map_prob'].append(map_nevents_prob)
        cred_level = 0.0
        for n, p in sorted_nevent_freqs:
            if n == true_nevents:
                break
            cred_level += p
        results['nevents']['cred_level'].append(cred_level)
    assert (len(results['divergence time']['true']) == number_of_sims * number_of_pairs)
    assert (len(results['root population size']['true']) == number_of_sims * number_of_pairs)
    assert (len(results['leaf population size']['true']) == number_of_sims * number_of_pairs * 2)
    assert (len(results['nevents']['true']) == number_of_sims)
    return results

def generate_bake_off_plots(
        number_of_pairs = 5,
        number_of_sims = 500,
        posterior_sample_size = 2000,
        prior_sample_size = 500000,
        include_msbayes = False):
    _LOG.info("Generating plots of bake-off results...")
    plot_dir = os.path.join(project_util.BAKE_OFF_DIR, "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    msbayes_results = None
    if include_msbayes:
        msbayes_results = get_msbayes_results(
                true_path = os.path.join(project_util.BAKE_OFF_DIR,
                        "results",
                        "msbayes",
                        "pymsbayes-results",
                        "observed-summary-stats",
                        "observed-1.txt"),
                sim_dir = os.path.join(project_util.BAKE_OFF_DIR,
                        "results",
                        "msbayes",
                        "pymsbayes-results",
                        "pymsbayes-output",
                        "d1",
                        "m1"),
                number_of_pairs = number_of_pairs,
                number_of_sims = number_of_sims,
                posterior_sample_size = posterior_sample_size,
                prior_sample_size = prior_sample_size)
    dpp_msbayes_results = get_msbayes_results(
            true_path = os.path.join(project_util.BAKE_OFF_DIR,
                    "results",
                    "dpp-msbayes",
                    "pymsbayes-results",
                    "observed-summary-stats",
                    "observed-1.txt"),
            sim_dir = os.path.join(project_util.BAKE_OFF_DIR,
                    "results",
                    "dpp-msbayes",
                    "pymsbayes-results",
                    "pymsbayes-output",
                    "d1",
                    "m1"),
            number_of_pairs = number_of_pairs,
            number_of_sims = number_of_sims,
            posterior_sample_size = posterior_sample_size,
            prior_sample_size = prior_sample_size)
    results_paths = glob.glob(os.path.join(project_util.VAL_DIR,
            "03pairs-dpp-root-0100-040k-0200l",
            "batch00[12345]",
            "results.csv.gz"))
    results = pycoevolity.parsing.get_dict_from_spreadsheets(
            results_paths,
            sep = "\t",
            offset = 0)

    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)

    parameter_dict = {
            "divergence time": [
                    "root_height_c1sp1",
                    "root_height_c2sp1",
                    "root_height_c3sp1",
                    ],
            "root population size": [
                    "pop_size_root_c1sp1",
                    "pop_size_root_c2sp1",
                    "pop_size_root_c3sp1",
                    ],
            "leaf population size": [
                    "pop_size_c1sp1",
                    "pop_size_c2sp1",
                    "pop_size_c3sp1",
                    "pop_size_c1sp2",
                    "pop_size_c2sp2",
                    "pop_size_c3sp2",
                    ],
            "nevents": [],
    }
    col_headers = [
            'ecoevolity',
            'dpp-msbayes',
            ]
    if include_msbayes:
        col_headers.append('msbayes')
    for parameter_key, parameters in parameter_dict.items():
        parameter_label = parameter_key 
        parameter_symbol = "\\tau"
        if parameter_key.endswith("size"):
            parameter_symbol = "N_e\\mu"

        plt.close('all')
        if include_msbayes:
            fig = plt.figure(figsize = (9, 3))
        else:
            fig = plt.figure(figsize = (6, 3))
        nrows = 1
        ncols = len(col_headers)
        gs = gridspec.GridSpec(nrows, ncols,
                wspace = 0.0,
                hspace = 0.0)

        if parameter_key == "nevents":
            ecoevolity_results = {'true': [], 'map': [], 'cred_level': [],
                    'true_prob': [], 'map_prob': []}
            ecoevolity_results['true'].extend(int(x) for x in results["true_num_events"])
            ecoevolity_results['map'].extend(int(x) for x in results["map_num_events"])
            ecoevolity_results['cred_level'].extend(float(x) for x in results["true_num_events_cred_level"])
            for i in range(len(results["map_num_events"])):
                map_n = results["map_num_events"][i]
                map_p = float(results["num_events_{0}_p".format(map_n)][i])
                true_n = results["true_num_events"][i]
                true_p = float(results["num_events_{0}_p".format(true_n)][i])
                ecoevolity_results["map_prob"].append(map_p)
                ecoevolity_results["true_prob"].append(true_p)
            all_results = {
                'ecoevolity' : ecoevolity_results,
                'dpp-msbayes': dpp_msbayes_results[parameter_key],
                }
            if include_msbayes:
                all_results['msbayes'] = msbayes_results[parameter_key]

            row_idx = 0
            for col_idx, col_header in enumerate(col_headers):
                r = all_results[col_header]
                true_nevents = r['true']
                map_nevents = r['map']
                map_nevents_probs = r['map_prob']
                true_nevents_probs = r['true_prob']
                true_nevents_cred_levels = r['cred_level']

                assert(len(true_nevents) == len(map_nevents))
                assert(len(true_nevents) == len(true_nevents_cred_levels))
                assert(len(true_nevents) == len(map_nevents_probs))
                assert(len(true_nevents) == len(true_nevents_probs))

                mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
                median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

                nevents_within_95_cred = 0
                ncorrect = 0
                true_map_nevents = []
                true_map_probs = []
                for i in range(number_of_pairs):
                    true_map_nevents.append([0 for i in range(number_of_pairs)])
                    true_map_probs.append([[] for i in range(number_of_pairs)])
                for i in range(len(true_nevents)):
                    true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
                    true_map_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
                    if true_nevents_cred_levels[i] <= 0.95:
                        nevents_within_95_cred += 1
                    if true_nevents[i] == map_nevents[i]:
                        ncorrect += 1
                p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
                p_correct = ncorrect / float(len(true_nevents))

                _LOG.info("p(nevents within CS) = {0:.4f}".format(p_nevents_within_95_cred))
                ax = plt.subplot(gs[row_idx, col_idx])

                ax.imshow(true_map_nevents,
                        origin = 'lower',
                        cmap = cmap,
                        interpolation = 'none',
                        aspect = 'auto'
                        # extent = [0.5, 3.5, 0.5, 3.5]
                        )
                for i, row_list in enumerate(true_map_nevents):
                    for j, num_events in enumerate(row_list):
                        # if num_events > 0:
                        #     map_probs = true_map_probs[i][j]
                        #     assert len(map_probs) == num_events
                        #     mean_prob = sum(map_probs) / len(map_probs)
                        #     median_prob = pycoevolity.stats.median(map_probs)
                        #     ax.text(j, i,
                        #             "{0:d}\n{1:.3f}".format(num_events, median_prob),
                        #             horizontalalignment = "center",
                        #             verticalalignment = "center")
                        # else:
                        ax.text(j, i,
                                str(num_events),
                                horizontalalignment = "center",
                                verticalalignment = "center")
                ax.text(0.98, 0.02,
                        "\\scriptsize$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                                p_nevents_within_95_cred),
                        horizontalalignment = "right",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
                ax.text(0.02, 0.98,
                        "\\scriptsize$p(\\hat{{k}} = k) = {0:.3f}$".format(
                                p_correct),
                        horizontalalignment = "left",
                        verticalalignment = "top",
                        transform = ax.transAxes)
                # ax.text(0.99, 0.99,
                #         "\\scriptsize$\\overline{{pp(k)}} = {0:.3f}$".format(
                #                 mean_true_nevents_prob),
                #         horizontalalignment = "right",
                #         verticalalignment = "top",
                #         transform = ax.transAxes)
                ax.text(0.98, 0.98,
                        "\\scriptsize$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
                                median_true_nevents_prob),
                        horizontalalignment = "right",
                        verticalalignment = "top",
                        transform = ax.transAxes)
                if row_idx == 0:
                    ax.text(0.5, 1.015,
                            col_header,
                            horizontalalignment = "center",
                            verticalalignment = "bottom",
                            size = 16.0,
                            transform = ax.transAxes)

            # show only the outside ticks
            all_axes = fig.get_axes()
            for ax in all_axes:
                if not ax.is_last_row():
                    ax.set_xticks([])
                if not ax.is_first_col():
                    ax.set_yticks([])

            # show tick labels only for lower-left plot 
            all_axes = fig.get_axes()
            for ax in all_axes:
                ax.xaxis.set_ticks(range(number_of_pairs))
                ax.yaxis.set_ticks(range(number_of_pairs))
                if ax.is_last_row() and ax.is_first_col():
                    xtick_labels = [item for item in ax.get_xticklabels()]
                    for i in range(len(xtick_labels)):
                        xtick_labels[i].set_text(str(i + 1))
                    ytick_labels = [item for item in ax.get_yticklabels()]
                    for i in range(len(ytick_labels)):
                        ytick_labels[i].set_text(str(i + 1))
                    ax.set_xticklabels(xtick_labels)
                    ax.set_yticklabels(ytick_labels)
                else:
                    xtick_labels = ["" for item in ax.get_xticklabels()]
                    ytick_labels = ["" for item in ax.get_yticklabels()]
                    ax.set_xticklabels(xtick_labels)
                    ax.set_yticklabels(ytick_labels)

            # avoid doubled spines
            all_axes = fig.get_axes()
            for ax in all_axes:
                for sp in ax.spines.values():
                    sp.set_visible(False)
                    sp.set_linewidth(2)
                if ax.is_first_row():
                    ax.spines['top'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                else:
                    ax.spines['bottom'].set_visible(True)
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    ax.spines['right'].set_visible(True)
                else:
                    ax.spines['right'].set_visible(True)

            fig.text(0.5, 0.002,
                    "True number of events ($k$)",
                    horizontalalignment = "center",
                    verticalalignment = "bottom",
                    size = 12.0)
            fig.text(0.005, 0.5,
                    "Estimated number of events ($\\hat{{k}}$)",
                    horizontalalignment = "left",
                    verticalalignment = "center",
                    rotation = "vertical",
                    size = 12.0)

            gs.update(left = 0.08, right = 0.995, bottom = 0.14, top = 0.91)

            plot_path = os.path.join(plot_dir,
                    "nevents.pdf")
            plt.savefig(plot_path)
            _LOG.info("Plots written to {0!r}\n".format(plot_path))
            continue

        parameter_min = float('inf')
        parameter_max = float('-inf')
        for parameter_str in parameters:
            parameter_min = min(parameter_min,
                    min(float(x) for x in results["true_{0}".format(parameter_str)]))
            parameter_max = max(parameter_max,
                    max(float(x) for x in results["true_{0}".format(parameter_str)]))
            parameter_min = min(parameter_min,
                    min(float(x) for x in results["mean_{0}".format(parameter_str)]))
            parameter_max = max(parameter_max,
                    max(float(x) for x in results["mean_{0}".format(parameter_str)]))
        parameter_min = min(parameter_min,
                min(dpp_msbayes_results[parameter_key]['true']))
        parameter_min = min(parameter_min,
                min(dpp_msbayes_results[parameter_key]['mean']))
        parameter_max = max(parameter_max,
                max(dpp_msbayes_results[parameter_key]['true']))
        parameter_max = max(parameter_max,
                max(dpp_msbayes_results[parameter_key]['mean']))
        if include_msbayes:
            parameter_min = min(parameter_min,
                    min(msbayes_results[parameter_key]['true']))
            parameter_min = min(parameter_min,
                    min(msbayes_results[parameter_key]['mean']))
            parameter_max = max(parameter_max,
                    max(msbayes_results[parameter_key]['true']))
            parameter_max = max(parameter_max,
                    max(msbayes_results[parameter_key]['mean']))

        axis_buffer = math.fabs(parameter_max - parameter_min) * 0.05
        axis_min = parameter_min - axis_buffer
        axis_max = parameter_max + axis_buffer

        ecoevolity_results = {'true': [], 'mean': [], 'lower': [], 'upper': []}
        for parameter_str in parameters:
            ecoevolity_results['true'].extend(float(x) for x in results["true_{0}".format(parameter_str)])
            ecoevolity_results['mean'].extend(float(x) for x in results["mean_{0}".format(parameter_str)])
            ecoevolity_results['lower'].extend(float(x) for x in results["eti_95_lower_{0}".format(parameter_str)])
            ecoevolity_results['upper'].extend(float(x) for x in results["eti_95_upper_{0}".format(parameter_str)])

        all_results = {
            'ecoevolity' : ecoevolity_results,
            'dpp-msbayes': dpp_msbayes_results[parameter_key],
            }
        if include_msbayes:
            all_results['msbayes'] = msbayes_results[parameter_key]

        row_idx = 0
        for col_idx, col_header in enumerate(col_headers):
            r = all_results[col_header]
            x = r['true']
            y = r['mean']
            y_lower = r['lower']
            y_upper = r['upper']

            assert(len(x) == len(y))
            assert(len(x) == len(y_lower))
            assert(len(x) == len(y_upper))
            proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
                    x,
                    y_lower,
                    y_upper)
            rmse = pycoevolity.stats.root_mean_square_error(x, y)
            _LOG.info("p(within CI) = {0:.4f}".format(proportion_within_ci))
            _LOG.info("RMSE = {0:.2e}".format(rmse))
            ax = plt.subplot(gs[row_idx, col_idx])
            line = ax.errorbar(
                    x = x,
                    y = y,
                    yerr = get_errors(y, y_lower, y_upper),
                    ecolor = '0.65',
                    elinewidth = 0.5,
                    capsize = 0.8,
                    barsabove = False,
                    marker = 'o',
                    linestyle = '',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(axis_min, axis_max)
            ax.set_ylim(axis_min, axis_max)
            identity_line, = ax.plot(
                    [axis_min, axis_max],
                    [axis_min, axis_max])
            plt.setp(identity_line,
                    color = '0.7',
                    linestyle = '-',
                    linewidth = 1.0,
                    marker = '',
                    zorder = 0)
            ax.text(0.02, 0.97,
                    "\\scriptsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                            parameter_symbol,
                            proportion_within_ci),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    size = 6.0,
                    zorder = 200)
            ax.text(0.02, 0.87,
                    # "\\scriptsize\\noindent$\\textrm{{\\sffamily RMSE}} = {0:.2e}$".format(
                    "\\scriptsize\\noindent RMSE = {0:.2e}".format(
                            rmse),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    size = 6.0,
                    zorder = 200)
            if row_idx == 0:
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        size = 16.0,
                        transform = ax.transAxes)

        # show only the outside ticks
        all_axes = fig.get_axes()
        for ax in all_axes:
            if not ax.is_last_row():
                ax.set_xticks([])
            if not ax.is_first_col():
                ax.set_yticks([])

        # show tick labels only for lower-left plot 
        all_axes = fig.get_axes()
        for ax in all_axes:
            if ax.is_last_row() and ax.is_first_col():
                continue
            xtick_labels = ["" for item in ax.get_xticklabels()]
            ytick_labels = ["" for item in ax.get_yticklabels()]
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)

        # avoid doubled spines
        all_axes = fig.get_axes()
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                sp.set_linewidth(2)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
            else:
                ax.spines['bottom'].set_visible(True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
            else:
                ax.spines['right'].set_visible(True)

        fig.text(0.5, 0.001,
                "True {0} (${1}$)".format(parameter_label, parameter_symbol),
                horizontalalignment = "center",
                verticalalignment = "bottom",
                size = 14.0)
        fig.text(0.005, 0.5,
                "Estimated {0} ($\\hat{{{1}}}$)".format(parameter_label, parameter_symbol),
                horizontalalignment = "left",
                verticalalignment = "center",
                rotation = "vertical",
                size = 12.0)

        gs.update(left = 0.11, right = 0.995, bottom = 0.17, top = 0.92)

        plot_file_prefix = parameter_label.replace(" ", "-")
        plot_path = os.path.join(plot_dir,
                "{0}-scatter.pdf".format(plot_file_prefix))
        plt.savefig(plot_path, dpi=600)
        _LOG.info("Plots written to {0!r}\n".format(plot_path))


def main_cli(argv = sys.argv):
    results_paths = [
            [
                glob.glob(os.path.join(project_util.VAL_DIR,
                        "05pairs-00500sites",
                        "batch00?",
                        "results.csv.gz")),
                glob.glob(os.path.join(project_util.VAL_DIR,
                        "05pairs-01000sites",
                        "batch00?",
                        "results.csv.gz")),
                glob.glob(os.path.join(project_util.VAL_DIR,
                        "05pairs-02000sites",
                        "batch00?",
                        "results.csv.gz")),
                glob.glob(os.path.join(project_util.VAL_DIR,
                        "05pairs-10000sites",
                        "batch00?",
                        "results.csv.gz")),
            ]
    ]
    column_labels = [
            "500 sites",
            "1000 sites",
            "2000 sites",
            "10000 sites",
    ]

    parameter_symbol = "t"
    generate_scatter_plots(
            parameters = [
                    "root_height_c1sp1",
                    "root_height_c2sp1",
                    "root_height_c3sp1",
                    "root_height_c4sp1",
                    "root_height_c5sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "divergence time",
            parameter_symbol = parameter_symbol,
            plot_width = 1.9,
            plot_height = 2.2,
            pad_left = 0.1,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            y_label_size = 18.0,
            y_label = "Estimated time ($\\hat{{{0}}}$)".format(parameter_symbol),
            plot_file_prefix = "div-time")
    parameter_symbol = "N_e\\mu"
    generate_scatter_plots(
            parameters = [
                    "pop_size_root_c1sp1",
                    "pop_size_root_c2sp1",
                    "pop_size_root_c3sp1",
                    "pop_size_root_c4sp1",
                    "pop_size_root_c5sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "root population size",
            parameter_symbol = parameter_symbol,
            plot_width = 1.9,
            plot_height = 2.2,
            pad_left = 0.11,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            y_label_size = 18.0,
            y_label = "Estimated size ($\\hat{{{0}}}$)".format(parameter_symbol),
            plot_file_prefix = "root-pop-size")
    generate_scatter_plots(
            parameters = [
                    "pop_size_c1sp1",
                    "pop_size_c2sp1",
                    "pop_size_c3sp1",
                    "pop_size_c4sp1",
                    "pop_size_c5sp1",
                    "pop_size_c1sp2",
                    "pop_size_c2sp2",
                    "pop_size_c3sp2",
                    "pop_size_c4sp2",
                    "pop_size_c5sp2",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "leaf population size",
            parameter_symbol = parameter_symbol,
            plot_width = 1.9,
            plot_height = 2.2,
            pad_left = 0.11,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            y_label_size = 18.0,
            y_label = "Estimated size ($\\hat{{{0}}}$)".format(parameter_symbol),
            plot_file_prefix = "leaf-pop-size")

    generate_model_plots(
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            number_of_comparisons = 5,
            plot_width = 2.5,
            plot_height = 3.0,
            pad_left = 0.065,
            pad_right = 0.995,
            pad_bottom = 0.18,
            pad_top = 0.92,
            y_label_size = 18.0,
            y_label = "Estimated number ($\\hat{{k}}$)",
            number_font_size = 6.0,
            plot_file_prefix = None)

    generate_histograms(
            parameters = [
                    "n_var_sites_c1",
                    "n_var_sites_c2",
                    "n_var_sites_c3",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "Number of variable sites",
            parameter_discrete = True,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "number-of-variable-sites")
    generate_histograms(
            parameters = [
                    "ess_sum_ln_likelihood",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "Effective sample size of log likelihood",
            parameter_discrete = False,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "ess-ln-likelihood")
    generate_histograms(
            parameters = [
                    "ess_sum_root_height_c1sp1",
                    "ess_sum_root_height_c2sp1",
                    "ess_sum_root_height_c3sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "Effective sample size of divergence time",
            parameter_discrete = False,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "ess-div-time")
    generate_histograms(
            parameters = [
                    "ess_sum_pop_size_root_c1sp1",
                    "ess_sum_pop_size_root_c2sp1",
                    "ess_sum_pop_size_root_c3sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "Effective sample size of root population size",
            parameter_discrete = False,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "ess-root-pop-size")
    generate_histograms(
            parameters = [
                    "psrf_ln_likelihood",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "PSRF of log likelihood",
            parameter_discrete = False,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "psrf-ln-likelihood")
    generate_histograms(
            parameters = [
                    "psrf_root_height_c1sp1",
                    "psrf_root_height_c2sp1",
                    "psrf_root_height_c3sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "PSRF of divergence time",
            parameter_discrete = False,
            range_key = "range",
            number_of_digits = 0,
            plot_width = 2.0,
            plot_height = 2.3,
            pad_left = 0.09,
            pad_right = 0.995,
            pad_bottom = 0.22,
            pad_top = 0.91,
            plot_file_prefix = "psrf-div-time")

    plot_ess_versus_error(
            parameters = [
                    "root_height_c1sp1",
                    "root_height_c2sp1",
                    "root_height_c3sp1",
                    ],
            results_paths = results_paths,
            column_labels = column_labels,
            row_labels = None,
            parameter_label = "divergence time",
            plot_file_prefix = "div-time")

    # generate_bake_off_plots(
    #         number_of_pairs = 5,
    #         number_of_sims = 500,
    #         posterior_sample_size = 2000,
    #         prior_sample_size = 500000)

    # plot_nevents_estimated_vs_true_probs(
    #         nevents = 1,
    #         sim_dir = "03pairs-dpp-root-0100-100k",
    #         nbins = 5,
    #         plot_file_prefix = "100k-sites")


if __name__ == "__main__":
    main_cli()
