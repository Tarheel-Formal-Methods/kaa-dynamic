from kaa.opts.kodiak import KodiakProd
from kaa.opts.bernstein import BernsteinProd

from kaa.templates import *

"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    use_parallel = False

    'The optimiation procedure to use in the bundle transformation. Optimization procedures are located in kaa.opts'
    OptProd = KodiakProd

    'The default template loading/unloading strategy to use during reachable set computations'
    DefaultStrat = StaticStrat

    'Suppress Output?'
    SuppressOutput = False

    'Number of samples to be used for volume estimation'
    VolumeSamples = 10000

    'seed for random.seed'
    RandSeed = 897987178


class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    plot_font = 15

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    default_fig_path = "/home/edward/work/kaa-dynamic/figures"

    'Figure dimensions'
    fig_size = (30,20)
