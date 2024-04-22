import pytest
from pyomo.environ import value
from SappyCSTR.core import EASap  # Make sure to import your function

def test_optimize_cstr():
    initial_conditions = {
        'C_A_in': 1.0,  # Inlet Concentration of A (M)
        'C_B_in': 1.0,  # Inlet Concentration of B (M)
        'v0': 0.5,  # volumetric flow rate (L/s)
        'V': 5,  # Reactor volume (L)
    }
    rate_constant_303K = 0.18  # Example rate constant at 303K for ESTERIFICATION
    activation_energy = 48.4*1000  # Example activation energy in J/mol for ESTERIFICATION
    compound_names = {
        'reactants': {
            'A': {'formula': 'C4H8O2', 'stoichiometry': 1},  # Ethyl Acetate
            'B': {'formula': 'NaOH', 'stoichiometry': 1},     # Sodium Hydroxide
        },
        'products': {
            'C': {'formula': 'C2H5OH', 'stoichiometry': 1},  # Ethanol
            'D': {'formula': 'NaC2H3O2', 'stoichiometry': 1},   # Sodium Acetate
        }
    }

    model = EASap(initial_conditions, rate_constant_303K, activation_energy, compound_names)

    assert value(model.C_A_in) == pytest.approx(3.1170817473657824, 0.01)
    assert value(model.C_B_in) == pytest.approx(1.5222873831305916, 0.01)
    assert value(model.V) == pytest.approx(0.6234563738047896, 0.01)
    assert value(model.v0) == pytest.approx(0.6814717895958581, 0.01)
    assert value(model.T_in) == pytest.approx(332.0706816612416, 0.01)
