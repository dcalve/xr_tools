import re

import xarray as xr

from ..util.general import safe_eval


def open_ds(fnames):
    """"""

    ds_all = []

    for i_file, fname in enumerate(fnames):
        ds = xr.open_mfdataset(fname).assign_coords(n_file=i_file)

        # NOTE: Tidier with swap_dims?
        for dim in ds.dims:
            i_dim = 'i_{}'.format(dim)
            ds = ds.rename({dim: i_dim}).assign_coords(
                **{dim: xr.DataArray(ds[dim].values, dims=i_dim), i_dim: range(ds[dim].size)})

        ds_all += [ds]

    ds = xr.concat(ds_all, dim='n_file')

    return ds


def apply_ds_expr(ds, expr):
    """"""

    # Apply expressions in sequence (separated by '=>') starting from ds

    ds_wrk = ds
    expr_parsed = []

    for e in expr.split('=>'):
        e = e.strip()

        # '{}' == dataset resulting from previous expression
        e = re.sub('{}', 'ds_wrk', e)
        # '{n}' == nth dataset
        e = re.sub('{(?P<n_file>[0-9]+)}', 'ds.isel(n_file=\g<n_file>)', e)

        # '{var}' == variable from '{}'
        e = re.sub('{(?P<var>[\w]+)}', 'ds_wrk["\g<var>"]', e)
        # '{var=n}' == variable from '{n}'
        e = re.sub('{(?P<var>[\w]+)=(?P<n_file>[0-9]+)}', 'ds["\g<var>"].isel(n_file=\g<n_file>)', e)
        # 'var' == '{var}'
        e = re.sub('^(?P<var>[\w]+)$', 'ds_wrk["\g<var>"]', e)

        # '{[dim=ind]}' == index ind of dimension dim of '{}'
        e = re.sub('{\[(?P<dim>[^=]+)=(?P<ind>[^=]+)]}', 'ds_wrk.isel(i_\g<dim>=\g<ind>)', e)
        # '{n[dim=ind]}' == index ind of dimension dim of '{n}'
        e = re.sub('{(?P<n_file>[0-9]+)\[(?P<dim>[^=]+)=(?P<ind>[^=]+)]}',
                   'ds_wrk.isel(n_file=\g<n_file>, i_\g<dim>=\g<ind>)', e)

        # '{var[dim=ind]}' == index ind of dimension dim of '{var}'
        e = re.sub('{(?P<var>[\w]+)\[(?P<dim>[^=]+)=(?P<ind>[^=]+)]}', 'ds_wrk["\g<var>"].isel(i_\g<dim>=\g<ind>)', e)
        # '{var=n[dim=ind]}' == index ind of dimension dim of '{var=n}'
        e = re.sub('{(?P<var>[\w]+)=(?P<n_file>[0-9]+)\[(?P<dim>[^=]+)=(?P<ind>[^=]+)]}',
                   'ds_wrk["\g<var>"].isel(n_file=\g<n_file>, i_\g<dim>=\g<ind>)', e)

        ds_wrk = safe_eval(e, ds=ds, ds_wrk=ds_wrk, allow_calls=None, allow_methods='isel')
        expr_parsed += [e]

    print ' => '.join(expr_parsed)

    return ds_wrk


def check_da(da):
    """"""

    if not isinstance(da, xr.DataArray):
        raise TypeError('Data must be an xarray.DataArray object, found:\n{}'.format(da))
