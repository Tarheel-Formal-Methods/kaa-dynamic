from kaa.opts.kodiak import KodiakProd
from kaa.opts.bernstein import BernsteinProd

"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    Parallelize = True

    'The optimiation procedure to use in the bundle transformation. Optimization procedures are located in kaa.opts'
    OptProd = KodiakProd

    'Suppress Output?'
    SuppressOutput = False

    'Number of samples to be used for volume estimation'
    VolumeSamples = 10000

    'seed for random.seed'
    RandSeed = 897987178

    'Use random sampling scheme for volume'
    RandVol = False

    'Save the flowpipe when error appears during transformation'
    SaveStateonError = True

    'Path for data directory to save all xlsx files from experiments.'
    DataDir = "/Users/edwardkim/Work/kaa-optimize/data"

    'Flag to trigger enveloping box threshold checking'
    UseThreshold = False

class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    plot_font = 15

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    default_fig_path = "/Users/edwardkim/Work/kaa-optimize/figures"

    'Figure dimensions'
    fig_size = (40,20)
