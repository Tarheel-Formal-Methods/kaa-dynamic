from models.cvdp import CoupledVDP_UnitBox, CoupledVDP
from kaa.experiment import *
from kaa.experi_init import *
from kaa.timer import Timer

def test_sapo_CVDP():
    num_steps = 40

    model = CoupledVDP(delta=0.08)

    experi_input = dict(model=model,
                        strat=None,
                        label=f"Sapo's Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1)
    Timer.generate_stats()


def test_skewed_sliding_strat_comb_CVDP():
    unit_model = CoupledVDP_UnitBox(delta=0.08)
    model = CoupledVDP(delta=0.08)
    test_skewed_sliding_strat_comb(unit_model, 50, 4000, num_temps=5, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

def test_ran_strat_CVDP():
  model = VanDerPol_UnitBox(delta=0.08)
  test_ran_strat(model, 70, 5000, use_supp=True, use_pregen=False)

def test_init_reach_vol_CVDP():
    num_steps = 35
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 1
    lin_window_size = 4

    inputs_one = []
    inputs_two = []
    for inc in range(10):
        inc /= 100

        box = ((1.25-inc,1.55), (2.25-inc,2.35),(1.25-inc,1.55), (2.25-inc,2.35))

        unit_model = CoupledVDP_UnitBox(delta=0.08, init_box=box)
        model = CoupledVDP(delta=0.08, init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                              strat=MultiStrategy(pca_strat, lin_strat),
                              label=f"VDP SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                              supp_mode = use_supp,
                              pregen_mode = use_pregen,
                              num_trajs=num_trajs,
                              num_steps=num_steps)

        experi_input_two = dict(model=model,
                              strat=None,
                              label=f"SapoVDP",
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

    experi = InitReachPlotExperiment(*inputs)
    experi.execute()
