from models.chem import Chem_UnitBox
from kaa.experiment import *
from kaa.experi_init import *
from kaa.timer import Timer
from kaa.bundle import BundleTransMode

def test_arch_Chem():
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    num_steps_one = 40
    num_steps_two = 40
    num_steps_three = 40

    delta_one = 0.05
    delta_two = 0.0075
    delta_three = 0.00005

    model_one = Chem_UnitBox(100, 1000, delta=delta_one)
    model_two = Chem_UnitBox(1000, 100000, delta=delta_two)
    model_three = Chem_UnitBox(1000, 10000000, delta=delta_three)

    pca_window_size_one = 10
    lin_window_size_one = 0

    pca_strat_one = SlidingPCAStrat(model_one, lifespan=pca_window_size_one)
    lin_strat_one = SlidingLinStrat(model_one, lifespan=lin_window_size_one)

    experi_input_one = dict(model=model_one,
                        strat=MultiStrategy(pca_strat_one, lin_strat_one),
                        label=f"ROBE21 Case 1",
                        num_steps=num_steps_one,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        trans_mode=BundleTransMode.AFO)

    pca_window_size_two = 5
    lin_window_size_two = 5

    pca_strat_two = SlidingPCAStrat(model_two, lifespan=pca_window_size_two)
    lin_strat_two = SlidingLinStrat(model_two, lifespan=lin_window_size_two)

    experi_input_two = dict(model=model_two,
                        strat=MultiStrategy(pca_strat_two, lin_strat_two),
                        label=f"ROBE21 Case 2",
                        num_steps=num_steps_two,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        trans_mode=BundleTransMode.AFO)

    pca_window_size_three = 0
    lin_window_size_three = 0

    pca_strat_three = SlidingPCAStrat(model_three, lifespan=pca_window_size_three)
    lin_strat_three = SlidingLinStrat(model_three, lifespan=lin_window_size_three)

    experi_input_three = dict(model=model_three,
                        strat=MultiStrategy(pca_strat_three, lin_strat_three),
                        label=f"ROBE21 Case 3",
                        num_steps=num_steps_three,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        trans_mode=BundleTransMode.AFO)

    experi = ProjectionPlotExperiment(experi_input_three, plot_total_width=True)
    experi.execute(ylims=(0.999,1.001))

    Timer.generate_stats()


def plot_arch_robe21():
    import matplotlib.pyplot as plot

    min_case_one = [0.9999902,  0.9999902,  0.9999902,  0.99999019, 0.99999018, 0.99999017,
                    0.99999016, 0.99999013, 0.99999009, 0.99999002, 0.99998995, 0.99998992,
                    0.99998989, 0.99998982, 0.99998976, 0.99998973, 0.99998973, 0.99998964,
                    0.99998959, 0.99998957, 0.99998954, 0.99998949, 0.99998942, 0.9999894,
                    0.99998941, 0.9999894,  0.9999894,  0.99998941, 0.99998937, 0.99998936,
                    0.99998943, 0.99998939, 0.99998938, 0.99998936, 0.99998934, 0.99998933,
                    0.99998929, 0.99998924, 0.99998925, 0.99998923, 0.99998921]

    max_case_one = [1.0000002,  1.0000002,  1.0000002,  1.0000002,  1.00000021, 1.00000021,
     1.00000023, 1.00000028, 1.00000029, 1.00000032, 1.00000035, 1.00000039,
     1.00000047, 1.00000054, 1.00000068, 1.0000007,  1.00000068, 1.00000077,
     1.00000071, 1.00000069, 1.00000068, 1.00000069, 1.0000007,  1.0000007,
     1.00000073, 1.00000073, 1.00000074, 1.00000073, 1.00000076, 1.00000078,
     1.00000079, 1.00000081, 1.00000082, 1.00000084, 1.00000085, 1.00000087,
     1.00000086, 1.00000091, 1.00000102, 1.00000112, 1.00000121]

    min_case_two = [0.9999902,  0.9999902,  0.99999018, 0.99999014, 0.99999016, 0.99999013,
                    0.99999012, 0.99999007, 0.99998996, 0.99998962, 0.99998955, 0.99998943,
                    0.99998914, 0.9999889,  0.99998869, 0.99998837, 0.99998776, 0.99998755,
                    0.99998726, 0.99998729, 0.9999868,  0.99998674, 0.99998627, 0.99998635,
                    0.99998586, 0.9999859,  0.99998559, 0.99998559, 0.99998529, 0.9999853,
                    0.999985,   0.99998496, 0.99998469, 0.99998462, 0.99998436, 0.99998427,
                    0.99998401, 0.99998392, 0.99998364, 0.99998352, 0.99998324]

    max_case_two = [1.0000002,  1.0000002 , 1.00000022, 1.00000022, 1.00000024, 1.00000027,
    1.00000027, 1.00000035, 1.00000043, 1.00000063, 1.00000068, 1.00000091,
    1.00000117, 1.00000142, 1.00000168, 1.00000198, 1.0000026,  1.0000027,
    1.0000031,  1.00000307, 1.00000361, 1.00000358, 1.00000407, 1.00000417,
    1.00000438, 1.00000443, 1.0000047,  1.00000459, 1.00000489, 1.00000485,
    1.00000515, 1.00000513, 1.00000545, 1.00000545, 1.00000573, 1.00000575,
    1.00000607, 1.0000061,  1.00000639, 1.00000644 ,1.00000676]


    step_two = 0.01
    step_two = 0.005

    fig, ax = plot.subplots()

    ax.fill_between(np.arange(41), min_case_one, max_case_one, color="C0", label="Case 1", alpha=0.4)
    ax.fill_between(np.arange(41), min_case_two, max_case_two, color="C1", label="Case 2", alpha=0.4)

    ax.set_ylim(0.999,1.001)

    plot.show()
