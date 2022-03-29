import os
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    Parallelize = False

    '''
    The optimiation procedure to use in the bundle transformation. Optimization procedures are located in kaa.opts
    Takes the strings "Kodiak" or "Bernstein"
    '''
    OptProd = "Bernstein"

    'Suppress Output?'
    SuppressOutput = False

    'Number of samples to be used for volume estimation'
    VolumeSamples = 10000

    'Seed for random.seed'
    RandSeed = 897987178

    'Save the flowpipe when error appears during transformation'
    SaveStateonError = True

    'Path for data directory to save all xlsx files from experiments.'
    DataDir = os.path.join(ROOT_DIR, "data")

    'Flag to trigger enveloping box threshold checking'
    UseThreshold = False

    'Run Normalization method if condition number becomes too large.'
    NormalizeLinDir = True

    'Evaluate and simplify the polynomial after performing functional composition.'
    EvalFinalPoly = False

class DebugSettings:
    TimerProfileLabels = set()

class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    PlotFont = 13

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    default_fig_path = os.path.join(ROOT_DIR, "figures")

    'Figure dimensions'
    fig_size = (60, 60)