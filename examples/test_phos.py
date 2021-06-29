from models.phos import Phosphorelay, Phosphorelay_UnitBox
from kaa.templates import MultiStrategy
from kaa.experi_init import *

from kaa.timer import Timer

def test_sapo_Phos():
    num_steps = 100
    model = Phosphorelay()

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoPhos",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0, 1, 2)

def test_sapo_vol_Phos():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 100

    model = Phosphorelay()
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoVDP",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input)
    harosc.execute(1)

def test_skewed_sliding_strat_comb_Phos():
    unit_model = Phosphorelay_UnitBox()
    model = Phosphorelay()
    test_skewed_sliding_strat_comb(unit_model, 100, 5000, num_temps=3, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_pca_Phos():
    model = Phosphorelay_UnitBox()
    test_sliding_pca(model, 20, 200, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_Phos():
    model = Phosphorelay_UnitBox()
    test_sliding_lin(model, 20, 200, 5000, use_supp=True, use_pregen=False)
