import cartopy
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cartopy import img_transform
from cartopy.mpl.geoaxes import GeoAxesSubplot
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, SymLogNorm
from matplotlib.ticker import LogLocator, MaxNLocator

from ..util.general import wrap_longitude

# Use scientific notation for 1e-3 <= n <= 1e3 in *default* tick labels
matplotlib.rcParams.update({'font.size': 10, 'axes.formatter.limits': (-3, 3)})


def figure(da, figsize=(8, 5), **kwargs):
    """"""

    plt.figure(figsize=figsize)
    FieldPlot(da, **kwargs).panel()


class FieldPlot(object):
    """"""

    def __init__(self, da, xname='nav_lon', yname='nav_lat', projection='PlateCarree',
                 extent=None, central_lon=None, levels=None, title=None,
                 cmap=None, nlev=25, clims=None, log_decs=0, cb_fmt=None,
                 subplot=111, **kwargs):
        """"""

        self.ref_crs = cartopy.crs.PlateCarree()

        # Other plotting args
        self.subplot = subplot
        self.title = title
        self.kwargs = kwargs

        # Axes / projection args
        self.extent = extent
        self.c_lon = central_lon if central_lon else self.default_central_lon
        self.nx = int(plt.gcf().get_dpi() * plt.gcf().get_figwidth())
        self.ny = int(plt.gcf().get_dpi() * plt.gcf().get_figheight())

        # Data args
        self.xname = xname
        self.yname = yname
        self.x, self.y, self.field = self.data(da, projection)

        # Contouring args
        self.vmin = clims[0] if clims else self.field.min()
        self.vmax = clims[1] if clims else self.field.max()
        self.pm = (self.vmax > 0 > self.vmin)
        self.levels = levels if levels else self.default_levels(nlev)
        self.norm = self.bnorm(log_decs)
        self.cmap = self.read_cmap(cmap) if cmap else self.default_cmap
        self.cb_fmt = cb_fmt

    def panel(self):
        """"""

        # Plotting
        clev = self.plot(self.x, self.y, self.field)

        # Adjust axes before adding further annotation
        plt.gcf().tight_layout(rect=(.05, 0, .95, .95))

        # Colorbar
        self.colorbar(clev)

        # Further annotation
        self.annotate()

    def plot(self, x, y, field, method='mesh'):
        """"""

        # Do the plotting
        if method == 'filled':
            clev = plt.contourf(x, y, field,
                                extend='both',
                                norm=self.norm, cmap=self.cmap)
        elif method == 'mesh':
            clev = plt.pcolormesh(x, y, field,
                                  norm=self.norm, cmap=self.cmap)
        # TODO: This is not yet working properly (e.g. doesn't use a colorbar, needs clabels, needs white axisbg,
        #       needs to know self.norm boundaries for all norm classes)
        #        elif method == 'line':
        #            clev = plt.contour(x, y, field, extend='both', colors='k')
        else:
            raise ValueError('Plot method "{}" not recognised'.format(method))

        return clev

    def colorbar(self, clev):
        """"""

        cb = plt.colorbar(clev, orientation='horizontal', extend='both', pad=.1, format=self.cb_fmt)

        # Adjust colorbar ticks for log scale
        if isinstance(clev.norm, (LogNorm, SymLogNorm)):
            clim = list(clev.get_clim())

            # For symmetric log scales, the central tick is the linear threshold
            if self.pm:
                clim[0] = clev.norm.linthresh

            # Determine major (every decade) and minor (every tenth of a decade) ticks
            minor_ticks = LogLocator(subs=range(2, 10)).tick_values(*clim)
            major_ticks = LogLocator().tick_values(*clim)
            minor_ticks = minor_ticks[(clim[0] <= minor_ticks) & (clim[1] >= minor_ticks)]
            major_ticks = major_ticks[(clim[0] <= major_ticks) & (clim[1] >= major_ticks)]

            # For symmetric log scales, ticks are mirrored around the central tick
            if self.pm:
                minor_ticks = np.concatenate([-minor_ticks, minor_ticks])
                major_ticks = np.concatenate([-major_ticks[1:], major_ticks])

            cb.set_ticks(major_ticks)

            # TODO: When self.cb_fmt = None, this drops the offset label (no longer a ScalarFormatter)- need to add
            # it back manually or use a FuncFormatter For symmetric log scales, update the central tick label to
            # indicate that it represents the linear threshold
            if self.pm:
                major_ticklabels = [i.get_text() for i in cb.ax.xaxis.get_ticklabels()]
                major_ticklabels[len(major_ticks) / 2] = '< |{}|'.format(major_ticklabels[len(major_ticks) / 2])
                cb.set_ticklabels(major_ticklabels)

            cb.ax.xaxis.set_ticks(clev.norm(minor_ticks), minor=True)

    def annotate(self):
        """"""

        is_geoax = isinstance(plt.gca(), GeoAxesSubplot)
        is_rectilinear = (
                is_geoax and isinstance(plt.gca().projection, (cartopy.crs.PlateCarree, cartopy.crs.Mercator)))

        # Set coastlines and land
        if is_geoax:
            plt.gca().coastlines()
            plt.gca().background_patch.set_facecolor('lightgrey')
        else:
            plt.gca().set_facecolor('lightgrey')

        # Set aspect ratio to fill panel
        plt.gca().set_aspect('auto')

        # Title
        if self.title:
            plt.title(self.title, fontsize=12)

        # Set limits
        if self.extent:
            if is_geoax:
                plt.gca().set_extent(self.extent, crs=self.ref_crs)
            else:
                plt.gca().set_xlim([self.extent[0], self.extent[1]])
                plt.gca().set_ylim([self.extent[2], self.extent[2]])
        elif not is_geoax:
            plt.gca().autoscale(tight=True)

        # Gridlines (rectilinear GeoAxes)
        if is_rectilinear:
            gl = plt.gca().gridlines(crs=self.ref_crs, draw_labels=True)
            #            gl.xlines = gl.ylines = False
            gl.xlabels_top = gl.ylabels_right = False
            gl.xlabel_style = gl.ylabel_style = {'size': 9, 'color': 'grey'}
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

            # Prevent duplication of ticks at the dateline
            plt.draw()

            for xlab in gl.xlabel_artists:
                if xlab.get_position()[0] == -180:
                    xlab.set_visible(False)

        # Gridlines (other GeoAxes)
        elif is_geoax:
            ylocs = MaxNLocator(3).tick_values(*self.extent[-2:])
            plt.gca().gridlines(crs=self.ref_crs, color='dimgrey', linewidth=1.5, ylocs=ylocs)

            # Labels have to be drawn by hand
            for yloc in ylocs:
                if self.extent[-2] <= yloc <= self.extent[-1]:
                    plt.text(self.c_lon, yloc + .2 * abs(ylocs[1] - ylocs[0]) * np.sign(ylocs[0]),
                             '{:.0f}$^o${}'.format(abs(yloc), 'N' if yloc > 0 else 'S'),
                             color='grey', fontsize=9, transform=self.ref_crs)

    def data(self, da, projection):
        """"""

        # Reproject the data
        if projection and (da[self.xname].ndim == 2):

            projection = getattr(cartopy.crs, projection)()

            # Set the central longitude (if applicable)
            if self.c_lon:
                try:
                    projection = projection.__class__(central_longitude=self.c_lon)
                except:
                    pass

            # TODO: Need to detect and add cyclic point insertion
            x, y = img_transform.mesh_projection(projection, nx=self.nx, ny=self.ny)[:2]
            field = img_transform.regrid(da.to_masked_array(), da[self.xname].values, da[self.yname].values,
                                         self.ref_crs, projection, x, y)

            plt.gcf().add_subplot(self.subplot, projection=projection)

        # Subsample the data and plot on grid coords
        else:
            # NOTE: Assume here that dims are (y, x); I don't think auto-transpose happens in mpl so should be
            # obvious from figure if this assumption is wrong
            da = da[slice(None, None, max(1, da.shape[0] / self.ny)),
                    slice(None, None, max(1, da.shape[1] / self.nx))]
            x = da[self.xname]
            y = da[self.yname]
            field = da.to_masked_array()

            plt.gcf().add_subplot(self.subplot)

        return x, y, field

    @property
    def default_central_lon(self):
        """"""

        # Determine central longitude; if this does not lie within extent then extent longitudes will not be set
        if self.extent:
            span = wrap_longitude(self.extent[1], base=self.extent[0]) - self.extent[0]

            # Set halfway between extent limits, unless extent is longitudinally global
            c_lon = self.extent[0] + span / 2. if span != 0 else 0
        else:
            c_lon = 0

        # Using a central longitude of +180 causes issues with unconstrained latitudes for some reason
        if c_lon == 180:
            c_lon = -180

        return c_lon

    @property
    def default_cmap(self):
        """"""

        if self.pm:
            return plt.cm.RdBu_r
        else:
            return plt.cm.viridis

    def default_levels(self, nlev):
        """"""

        loc = MaxNLocator(nlev + 1, symmetric=self.pm)
        loc.create_dummy_axis()
        loc.set_bounds(self.vmin, self.vmax)
        loc.set_bounds(*loc.autoscale())
        levs = loc()

        return levs

    def bnorm(self, log_decs=0):
        """"""

        # Logarithmic scale
        if log_decs:
            if self.pm:
                min_range = min(abs(self.vmin), abs(self.vmax))
                linthresh = 10 ** int(math.floor(math.log10(min_range / float(10 ** log_decs))))
                norm = SymLogNorm(linthresh, linscale=0, vmin=self.vmin, vmax=self.vmax)
            else:
                norm = LogNorm(vmin=self.vmin, vmax=self.vmax)

        # Linear scale
        else:
            norm = ExBoundaryNorm(self.levels, ncolors=255, set_middle=self.pm)

        return norm

    # TODO: Make independent function
    def read_cmap(self, cmap):
        """"""

        try:
            with open(cmap) as fid:
                rgb = [np.array(line.strip().lstrip('[').rstrip(']').strip().split(), dtype=float) for line in fid]
                cmap = ListedColormap(rgb)
        except IOError:
            try:
                cmap = getattr(plt.cm, cmap)
            except AttributeError:
                pkg, cmap = cmap.split('.')
                cmap = getattr(__import__(pkg), cmap)
            except:
                raise RuntimeError('Could not read colour map "{}"'.format(cmap))

        return cmap


class ExBoundaryNorm(BoundaryNorm):
    """
    Subclass of BoundaryNorm that reserves some end colours for extrema.
    """

    def __init__(self, boundaries, ncolors, set_middle=False):
        """
        Pick evenly-spaced colours for the given boundaries, reserving space at
        either end of the colour map for extrema.

        Arguments:
            boundaries-
                The boundaries for which to set colours.
            ncolors-
                The number of colours to use.

        Keywords:
            clip-
                Not used.
            set_middle-
                Set the middle boundaries to the central colour value.
        """

        self.set_middle = set_middle
        self.ncolors = ncolors
        self.space = ncolors / len(boundaries)

        super(ExBoundaryNorm, self).__init__(boundaries, ncolors - 2 * self.space)

        # Determine mid points for odd and even numbers of boundaries
        if self.N % 2 != 0:
            self.mid = [self.boundaries[self.N / 2 + i] for i in (-1, 1)]
        else:
            self.mid = [self.boundaries[self.N / 2 + i] for i in (-1, 0)]

    def __call__(self, x, clip=None):
        """
        Return colours for the requested boundaries.

        """

        # Do not modify colours returned for masked data
        if np.ma.is_masked(x):
            not_masked = ~x.mask
        else:
            not_masked = slice(None)

        col = super(ExBoundaryNorm, self).__call__(x)
        x_valid = x[not_masked]

        # Offset non-extrema colours
        col_valid = col[not_masked] + self.space
        col_valid[x_valid < self.vmin] = -1
        col_valid[x_valid >= self.vmax] = self.ncolors
        #         col_valid[x_valid < self.boundaries.min()] = -1
        #         col_valid[x_valid >= self.boundaries.max()] = self.ncolors

        # Set middle colours
        if self.set_middle:
            in_mid = np.where((x_valid >= self.mid[0]) &
                              (x_valid < self.mid[-1]))
            col_valid[in_mid] = self.ncolors / 2

        col[not_masked] = col_valid

        return col
