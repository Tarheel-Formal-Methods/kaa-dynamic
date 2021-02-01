from kaa.opts.kodiak import KodiakProd
from kaa.opts.bernstein import BernsteinProd

"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    Parallelize = False

    'The optimiation procedure to use in the bundle transformation. Optimization procedures are located in kaa.opts'
    OptProd = BernsteinProd

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

    'Use support points to sample inital points of trajectories'
    UseSuppPoints = True

    'Use pre-generated directions for all experiments.'
    UsePreGenDirs = False


class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    plot_font = 15

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    default_fig_path = "/Users/edwardkim/Work/kaa-optimize/figures"

    'Figure dimensions'
    fig_size = (30,20)
