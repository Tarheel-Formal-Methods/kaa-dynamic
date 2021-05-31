
"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    Parallelize = False

    'The optimiation procedure to use in the bundle transformation. Optimization procedures are located in kaa.opts'
    OptProd = "Kodiak"

    'Suppress Output?'
    SuppressOutput = False

    'Number of samples to be used for volume estimation'
    VolumeSamples = 10000

    'seed for random.seed'
    RandSeed = 897987178

    'Save the flowpipe when error appears during transformation'
    SaveStateonError = True

    'Path for data directory to save all xlsx files from experiments.'
    DataDir = "/home/edward/work/kaa-dynamic/data"

    'Flag to trigger enveloping box threshold checking'
    UseThreshold = False

    "Number of threads to instantiate when running parallel routines."
    ThreadCount = 16

    'Run Normalization method if condition number becomes too large.'
    NormalizeLinDir = False

    DelFlowpipe = False

class DebugSettings:
    TimerProfileLabels = set()

class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    PlotFont = 19

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    default_fig_path = "/home/edward/work/kaa-dynamic/figures"

    'Figure dimensions'
    fig_size = (60, 60)

    NumSteps = 2
