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
    OptProd = BernsteinProd

    'The default template loading/unloading strategy to use during reachable set computations'
    DefaultStrat = StaticStrat

    'Suppress Output?'
    SuppressOutput = False

class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    plot_font = 15

    'Toggle to save the figures to disk'
    save_fig = True

    'Path to save figures'
    default_fig_path = "/Users/edwardkim/Work/kaa-optimize/figures"

    'Figure dimensions'
    fig_size = (30,20)
