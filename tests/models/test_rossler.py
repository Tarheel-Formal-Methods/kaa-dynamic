from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.rossler import Rossler, Rossler_UnitBox
from kaa.trajectory import Traj
from kaa.experi_init import *
from kaa.timer import Timer

def test_sapo_Rossler():

    model = Rossler()
    num_steps = 150

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoRossler",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1,2)

def test_sapo_vol_Rossler():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    model = Rossler(delta=0.5)
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoRossler",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input)
    harosc.execute(1)

def test_init_reach_vol_vs_ran_Rossler():
    num_steps = 150
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 10
    lin_window_size = 10

    inputs = []
    for inc in range(5):
        inc /= 100

        box = ((0,0.1 + inc), (4.8 - inc,5), (0,0.1))

        unit_model = Rossler_UnitBox(init_box=box)
        model = Rossler(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"Rossler SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)

        inputs.append(experi_input_one)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = InitReachVSRandomPlotExperiment(*inputs, num_ran_temps=pca_window_size+lin_window_size, num_trials=10)
    experi.execute()


def test_sliding_strat_comb_Rossler():
    model = Rossler_UnitBox()
    test_sliding_strat_comb(model, 150, 4000, use_supp=True, use_pregen=False)
    Timer.generate_stats()

def test_skewed_sliding_strat_comb_Rossler():
    unit_model = Rossler_UnitBox()
    model = Rossler()
    test_skewed_sliding_strat_comb(unit_model, 150, 4000, num_temps=5, incre=1, use_supp=True, use_pregen=False, use_sapo=model)
    Timer.generate_stats()

def test_sliding_pca_Rossler():
    model = Rossler_UnitBox()
    test_sliding_pca(model, 20, 150, -1, use_supp=True, use_pregen=False)
    Timer.generate_stats()

def test_sliding_lin_Rossler():
    model = Rossler_UnitBox()
    test_sliding_lin(model, 20, 150, -1, use_supp=True, use_pregen=False)
    Timer.generate_stats()

def gen_save_dirs_Rossler():
    model = Rossler_UnitBox()
    gen_save_dirs(model, 150)
