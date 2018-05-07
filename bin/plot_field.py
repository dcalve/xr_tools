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
    args_data.add_argument('--xname', help='Name of x coordinate', type=str, default='nav_lon')
    args_data.add_argument('--yname', help='Name of y coordinate', type=str, default='nav_lat')
    opt_list += ['xname', 'yname']

    args_axes = parser.add_argument_group('Axes options')
    args_axes.add_argument('--title', help='Figure title', type=str)
    opt_list += ['title']

    args_contour = parser.add_argument_group('Contouring options')
    args_contour.add_argument('--levels',   help='Contour levels', type=float, nargs='*')
    args_contour.add_argument('--clims',    help='Contour range', type=float, nargs=2)
    args_contour.add_argument('--nlev',     help='Number of contour levels', type=int, default=25)
    args_contour.add_argument('--log_decs', help='Number of logarithmic decades (ignores --nlev)', type=int, default=0)
    args_contour.add_argument('--cmap',     help='Colour map', type=str)
    args_contour.add_argument('--cb_fmt',   help='Colour bar string format', type=str)
    opt_list += ['levels', 'clims', 'nlev', 'log_decs', 'cmap', 'cb_fmt']

    args_output = parser.add_argument_group('Output options')
    args_output.add_argument('--savefig', help='Save figure as filename', type=str)

    return parser.parse_args(), opt_list


def check_field(da):
    """"""

    check_da(da)
    if da.ndim != 2:
        raise ValueError('Shape of field to be plotted is not 2D:\n{}'.format(da))


def read_field(fnames, expr):
    """"""

    ds = open_ds(fnames)
    da = apply_ds_expr(ds, expr).squeeze()
    check_field(da)

    return da


if __name__ == '__main__':
    args, opts = parse_args()
    kwargs = {i: getattr(args, i) for i in opts}

    da_plot = read_field(args.files, args.var)

    if args.savefig:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from src.plot.field import figure

    figure(da_plot, **kwargs)

    if args.savefig:
        plt.savefig(args.savefig)
    else:
        plt.show()
