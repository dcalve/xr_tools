import operator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.colors import ListedColormap

# Use scientific notation for 1e-3 <= n <= 1e3 in *default* tick labels
matplotlib.rcParams.update({'font.size': 10, 'axes.formatter.limits': (-3, 3)})


def figure(da, figsize=(8, 5), **kwargs):
    """"""

    plt.figure(figsize=figsize)
    LinePlot(da, **kwargs).panel()


class LinePlot(object):
    """"""

    def __init__(self, da,
                 title=None, xlabel=None, ylabel=None, xlims=None, clims=None, xname=None,
                 labels=None, color=None, linestyle=('solid',), linewidth=(1.5,), style_cycle=None,
                 subplot=111, **kwargs):
        """"""

        # Data args
        self.da = da

        # Axis args
        self.title = title
        self.xlims = xlims
        self.ylims = clims
        self.x, self.y, self.z = self.axes(xname)
        self.xlabel = xlabel if xlabel else self.x
        self.ylabel = ylabel if ylabel else self.y

        # Line plot args
        self.labels = labels if labels else self.da[self.z].values if self.z else None
        self.set_line_styles(color, linestyle, linewidth, style_cycle)

        # Other plotting args
        self.pm = (self.ylims[1] > 0 > self.ylims[0]) if self.ylims else (self.da.max() > 0 > self.da.min())
        self.subplot = subplot
        self.kwargs = kwargs

    def panel(self):
        """"""

        # Plotting
        self.plot()

        # Adjust axes before adding further annotation
        plt.gcf().tight_layout(rect=(.05, .05, .95, .95))

        # Legend
        if self.labels is not None:
            self.legend()

        # Further annotation
        self.annotate()

    def plot(self):
        """"""

        self.da.plot.line(x=self.x, hue=self.z)

    def annotate(self):
        """"""

        # Title
        if self.title:
            plt.title(self.title, fontsize=12)

        # Set limits
        if self.xlims:
            plt.gca().set_xlim(*self.xlims)
        if self.ylims:
            plt.gca().set_ylim(*self.ylims)

        # Set labels
        if self.xlabel:
            plt.gca().set_xlabel(self.xlabel)
        if self.ylabel:
            plt.gca().set_ylabel(self.ylabel)

        # Gridlines
        plt.gca().grid(True, which='both', axis='x')

        # Zero line
        if self.pm:
            plt.gca().axhline(0, linestyle='dotted', color='grey', linewidth=1.5)

    def axes(self, xname):
        """"""

        # Default x coordinate is the last dimension in the DataArray
        x = xname if xname else self.da.dims[-1]

        # y coordinate is always the data
        y = self.da.name

        # z coordinate is the other dimension not represented by the x coordinate.
        # If the x coordinate is multidimensional (occurs when the coordinate diverges over n_file) 
        # then the z coordinate is taken to be the smallest dimension represented by the x coordinate.
        x_dim = self.da[x].dims[self.da[x].shape.index(max(self.da[x].shape))]
        z = self.da[{x_dim: 0}].dims[0] if self.da.ndim == 2 else None

        return x, y, z

    # TODO: Make independent function
    def legend(self):
        """"""

        leg = plt.legend(self.labels, fontsize=8, loc=0, fancybox=True)
        leg.get_frame().set_alpha(0.5)

    def set_line_styles(self, color, linestyle, linewidth, style_cycle):
        """"""

        cycle = matplotlib.rcParams['axes.prop_cycle'].by_key()
        cycle_len = self.da[self.z].size

        if color:
            if len(color) == 1:
                cycle['color'] = self.read_cmap(color[0])(np.linspace(0, 1, 9))
            else:
                cycle['color'] = color
        if linestyle:
            cycle['linestyle'] = linestyle
        if linewidth:
            cycle['linewidth'] = linewidth

        cycle = {k: (v * cycle_len)[:cycle_len] for k, v in cycle.items()}

        if style_cycle:
            prop_cycle = [cycler(**{i: cycle.pop(i)}) for i in style_cycle[::-1] if i in cycle]
            prop_cycle = reduce(operator.mul, prop_cycle, 1)[:cycle_len]
            prop_cycle += cycler(**cycle)
        else:
            prop_cycle = cycler(**cycle)

        plt.rc('axes', prop_cycle=prop_cycle)

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
