#! /usr/bin/env python2.7

import argparse
import os
import sys

import matplotlib

sys.path.append('{}/..'.format(os.path.abspath(os.path.dirname(__file__))))
from src.data.xrdata import open_ds, apply_ds_expr, check_da


def parse_args():
    """"""

    opt_list = []

    parser = argparse.ArgumentParser(description='')

    # Positional
    parser.add_argument('files', help='One or more file glob strings, each defining a dataset', type=str, nargs='+')
    parser.add_argument('var', help='Variable expression', type=str)

    # Optional
    args_data = parser.add_argument_group('Data options')
    args_data.add_argument('--preprocess', help='Variable expression (applied before "var")', type=str)
    opt_list += ['preprocess']

    args_axes = parser.add_argument_group('Axes options')
    args_axes.add_argument('--title', help='Figure title', type=str)
    args_axes.add_argument('--sp_x', help='Divide data by "var" (variable expression), "ds" (dataset), '
                                          'dimension or integer to generate subplots (x axis)', type=str)
    args_axes.add_argument('--sp_y', help='Divide data by "var" (variable expression), "ds" (dataset), '
                                          'dimension or integer to generate subplots (y axis)', type=str)
    args_axes.add_argument('--st_fmt', help='Template for subplot title (using "{x}" and "{y}")', type=str)
    args_axes.add_argument('--st_x', help='Values for {x} in --st_fmt', type=str, nargs='*')
    args_axes.add_argument('--st_y', help='Values for {y} in --st_fmt', type=str, nargs='*')
    args_axes.add_argument('--xlims', help='X axis range', type=float, nargs=2)
    args_axes.add_argument('--ylims', help='Y axis range (alternatively --clims)', type=float, nargs=2, dest='clims')
    args_axes.add_argument('--xname', help='Name of coordinate to use for x axis', type=str)
    args_axes.add_argument('--xlabel', help='X axis label', type=str)
    args_axes.add_argument('--ylabel', help='Y axis label', type=str)
    opt_list += ['title', 'sp_x', 'sp_y', 'st_fmt', 'st_x', 'st_y', 'xlims', 'xname', 'xlabel', 'ylabel']

    args_line = parser.add_argument_group('Line plot options (1D data)')
    args_line.add_argument('--labels', help='Legend labels', type=str, nargs='*')
    args_line.add_argument('--color', help='Line colours (or colour maps) to cycle through', type=str, nargs='*')
    args_line.add_argument('--linestyle', help='Line styles to cycle through', type=str, nargs='*', default=['solid'])
    args_line.add_argument('--linewidth', help='Line widths to cycle through', type=float, nargs='*', default=[1.5])
    args_line.add_argument('--style_cycle', help='Which line styles to cycle', type=str, nargs='*',
                           choices=['color', 'linestyle', 'linewidth'])
    opt_list += ['labels', 'color', 'linestyle', 'linewidth', 'style_cycle']

    args_contour = parser.add_argument_group('Contouring options (2D data)')
    args_contour.add_argument('--yname', help='Name of coordinate to use for y axis', type=str)
    args_contour.add_argument('--levels', help='Contour levels', type=float, nargs='*')
    args_contour.add_argument('--clims', help='Contour range (alternatively --ylims)', type=float, nargs=2)
    args_contour.add_argument('--nlev', help='Number of contour levels', type=int, default=25)
    args_contour.add_argument('--log_decs', help='Number of logarithmic decades (ignores --nlev)', type=int, default=0)
    args_contour.add_argument('--cmap', help='Colour map', type=str)
    args_contour.add_argument('--cb_fmt', help='Colour bar string format', type=str)
    opt_list += ['yname', 'levels', 'clims', 'nlev', 'log_decs', 'cmap', 'cb_fmt']

    args_output = parser.add_argument_group('Output options')
    args_output.add_argument('--savefig', help='Save figure as filename', type=str)

    return parser.parse_args(), opt_list


def check_field(da):
    """"""

    check_da(da)
    if da.ndim not in (1, 2):
        raise ValueError('Shape of field to be plotted must be 1D or 2D:\n{}'.format(da))

    ### Variables and datasets will occupy their own dimensions, so da must always be 2D and we pass xname to 
    ### identify the x axis dimension


def read_field(fnames, expr):
    """"""

    ds = open_ds(fnames)
    da = apply_ds_expr(ds, expr).squeeze()
    check_field(da)

    ### Need axis collapsing methods here

    return da


if __name__ == '__main__':
    args, opts = parse_args()
    kwargs = {i: getattr(args, i) for i in opts}

    ### Overlap between plot_field and 2D sections; generalise field.py?

    da_plot = read_field(args.files, args.var)

    if args.savefig:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from src.plot.line import figure

    figure(da_plot, **kwargs)

    if args.savefig:
        plt.savefig(args.savefig)
    else:
        plt.show()
