from models.cvdp import CoupledVDP_UnitBox, CoupledVDP
from kaa.experi_init import *
from kaa.timer import Timer
from kaa.bundle import BundleTransMode
from kaa.temp.diag_static_strat import RandomDiagStaticStrat
from kaa.templates import MultiStrategy

def test_sapo_CVDP():
    num_steps = 40

    model = CoupledVDP(delta=0.08)
    experi_input = dict(model=model,
                        strat=None,
                        label=f"Sapo's Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1)
    Timer.generate_stats

def test_arch_CVDP():
    use_supp = True
    use_pregen = False

    num_trajs = 0

    delta_one = 0.1
    delta_two = 0.1

    num_steps_one = int(6 / delta_one) - 1
    num_steps_two = int(6 / delta_two) - 1
    init_box_one = ((1.25, 1.55), (2.35, 2.45), (1.25, 1.55), (2.35, 2.45))
    init_box_two = ((1.55, 1.85), (2.35, 2.45), (1.55, 1.85), (2.35, 2.45))

    model_one = CoupledVDP_UnitBox(delta=delta_one, init_box=init_box_one , mu=1)
    model_two = CoupledVDP_UnitBox(delta=delta_two, init_box=init_box_two , mu=2)

    #pca_window_size = 10
    #lin_window_size = 0

    ran_static_strat_one = RandomDiagStaticStrat(model_one, 30)
    #pca_strat_one = SlidingPCAStrat(model_one, 10)

    experi_input_one = dict(model=model_one,
                            strat=ran_static_strat_one,
                            label=f"CVDP mu=1",
                            num_steps=num_steps_one,
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            trans_mode=BundleTransMode.AFO,
                            restrict_inter=(-10,10))

    ran_static_strat_two = RandomDiagStaticStrat(model_two, 30)
    #pca_strat_two = SlidingPCAStrat(model_two, 10)


    experi_input_two = dict(model=model_two,
                            strat=ran_static_strat_two,
                            label=f"CVDP mu=2",
                            num_steps=num_steps_two,
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            trans_mode=BundleTransMode.AFO,
                            restrict_inter=(-10,10))

    experi = PhasePlotExperiment(experi_input_one, experi_input_two, label="CVDP20")
    experi_data = experi.execute(0, 1)

    return experi_data

def test_sliding_phase_plot_CVDP():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    model = VanDerPol_UnitBox(delta=0.08)

    pca_window_size = 2
    lin_window_size = 0

    pca_strat = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0,1, plot_border_traj=False)
    Timer.generate_stats()


def test_ran_diag_static_CVDP():
    unit_model = CoupledVDP_UnitBox(delta=0.08)
    test_ran_diag_static(unit_model, 40, 5)

def test_skewed_sliding_strat_comb_CVDP():
    unit_model = CoupledVDP_UnitBox(delta=0.08)
    model = CoupledVDP(delta=0.08)
    test_skewed_sliding_strat_comb(unit_model, 50, 4000, num_temps=5, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

def test_ran_strat_CVDP():
  model = VanDerPol_UnitBox(delta=0.08)
  test_ran_strat(model, 70, 5000, use_supp=True, use_pregen=False)

def test_init_reach_ran_vol_CVDP():
    num_steps = 40
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 2
    lin_window_size = 18

    inputs = []
    for inc in range(5):
        inc /= 100

        box = ((1.25-inc,1.55), (2.25-inc,2.35),(1.25-inc,1.55), (2.25-inc,2.35))

        unit_model = CoupledVDP_UnitBox(delta=0.08, init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                              strat=MultiStrategy(pca_strat, lin_strat),
                              label=f"CVDP PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
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

    avg_ran_vals = [14502.5921434655, 17671.5240439977, 22296.9709601768, 32466.8538839619, 39276.8448263568]

    experi = InitReachVSRandomPlotExperiment(*inputs, num_ran_temps=pca_window_size+lin_window_size, num_trials=10, precalc_vals=avg_ran_vals)
    experi.execute()

def test_init_reach_vol_CVDP():
    num_steps = 40
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 1
    lin_window_size = 4

    inputs_one = []
    inputs_two = []
    for inc in range(1,6,1):
        inc /= 100

        box = ((1.25-inc,1.55), (2.25-inc,2.35),(1.25-inc,1.55), (2.25-inc,2.35))

        unit_model = CoupledVDP_UnitBox(delta=0.08, init_box=box)
        model = CoupledVDP(delta=0.08, init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                              strat=MultiStrategy(pca_strat, lin_strat),
                              label=f"CVDP PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
                              supp_mode = use_supp,
                              pregen_mode = use_pregen,
                              num_trajs=num_trajs,
                              num_steps=num_steps)

        experi_input_two = dict(model=model,
                              strat=None,
                              label=f"SapoCVDP",
                              supp_mode = use_supp,
                              pregen_mode = use_pregen,
                              num_trajs=num_trajs,
                              num_steps=num_steps)


        inputs_one.append(experi_input_one)
        inputs_two.append(experi_input_two)


    inputs = inputs_one + inputs_two

    if use_supp:
      file_identifier = "(SUPP)"
    elif use_pregen:
      file_identifier = f"(PREGEN: {num_trajs})"
    else:
      file_identifier = "(RAND)"

    experi = InitReachPlotExperiment(*inputs, log_scale=True)
    experi.execute()
