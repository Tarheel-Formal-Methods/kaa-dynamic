from kaa.experiment import *
from kaa.experi_init import *
from kaa.timer import Timer
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox

def test_sapo_LV():
    num_steps = 100

    model = LotkaVolterra()
    experi_input = dict(model=model,
                        strat=None,
                        label=f"Sapo's Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1)
    Timer.generate_stats()

def test_skewed_sliding_strat_comb_LV():
    unit_model = LotkaVolterra_UnitBox()
    model = LotkaVolterra()
    test_skewed_sliding_strat_comb(unit_model, 100, 5000, num_temps=4, incre=1, use_supp=True, use_pregen=False, use_sapo=model)
