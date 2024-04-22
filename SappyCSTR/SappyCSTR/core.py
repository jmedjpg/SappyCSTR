

from pyomo.environ import ConcreteModel, Var, Constraint, Objective, SolverFactory, Param, NonNegativeReals, value, exp
from thermo import Chemical

def EASap(initial_conditions, rate_constant_303K, activation_energy, compound_names):

    model = ConcreteModel()
    def reaction_enthalpy(compound_names):
        delta_H_rxn_react = 0.0
        delta_H_rxn_prod = 0.0
        for compound in compound_names['reactants'].values():
            chemr = Chemical(compound['formula'])
            delta_H_rxn_react += compound['stoichiometry'] * chemr.Hf
        for compound in compound_names['products'].values():
            chemp = Chemical(compound['formula'])
            delta_H_rxn_prod += compound['stoichiometry'] * chemp.Hf
        return (delta_H_rxn_prod - delta_H_rxn_react)

    def specific_heat(compound_formula):
        compound = Chemical(compound_formula)
        compound.calculate(model.T())
        return compound.Cp

    R = 8.314
    model.v0 = Var(initialize=initial_conditions['v0'], within=NonNegativeReals, bounds=(0.01, 1))
    model.C_A_in = Var(within=NonNegativeReals, initialize=initial_conditions['C_A_in'], bounds=(0, 5))
    model.C_B_in = Var(within=NonNegativeReals, initialize=initial_conditions['C_B_in'], bounds=(0, 5))
    model.C_A = Var(within=NonNegativeReals, initialize=model.C_A_in, bounds=(0, 5))
    model.C_B = Var(within=NonNegativeReals, initialize=model.C_B_in, bounds=(0, 5))
    model.V = Var(initialize=initial_conditions['V'], within=NonNegativeReals, bounds=(0.1, 10))
    model.T_in = Var(initialize=303, within=NonNegativeReals, bounds=(300, 1000))
    model.T = Var(within=NonNegativeReals, initialize=model.T_in, bounds=(300, 1000))
    model.T_ref = Param(initialize=303)
    model.C_C_in = Param(initialize=0)
    model.C_D_in = Param(initialize=0)
    model.C_C = Var(within=NonNegativeReals, initialize=0)
    model.C_D = Var(within=NonNegativeReals, initialize=0)
    model.Q_dot = Param(initialize=0)
    model.W_dot = Param(initialize=0)
    model.k0 = Param(initialize=rate_constant_303K)
    model.Ea = Param(initialize=activation_energy)
    model.k = model.k0 * exp(-model.Ea / (R * model.T))
    model.reaction_rate = Var()
    model.reaction_rate_constraint = Constraint(expr=model.reaction_rate == model.k * model.C_A * model.C_B)

    def energy_balance_rule(model):
        Cp_A = specific_heat(compound_names['reactants']['A']['formula'])
        Cp_B = specific_heat(compound_names['reactants']['B']['formula'])
        Cp_C = specific_heat(compound_names['products']['C']['formula'])
        Cp_D = specific_heat(compound_names['products']['D']['formula'])
        heat_capacity_flow_inlet = (model.C_A_in * Cp_A + model.C_B_in * Cp_B) * (model.T - model.T_in)
        heat_capacity_flow_outlet = (model.C_C * Cp_C + model.C_D * Cp_D) * (model.T - model.T_ref)
        reaction_enthalpy_term = reaction_enthalpy(compound_names) * model.reaction_rate * model.V
        return (model.Q_dot - model.W_dot - heat_capacity_flow_inlet + heat_capacity_flow_outlet + reaction_enthalpy_term) == 0

    model.energy_balance = Constraint(rule=energy_balance_rule)

    def mole_balance_ruleA(model):
        return model.C_A_in * model.v0 - model.C_A * model.v0 - model.reaction_rate * model.V == 0
    model.mole_balance_A = Constraint(rule=mole_balance_ruleA)

    def mole_balance_ruleB(model):
        return model.C_B_in * model.v0 - model.C_B * model.v0 - model.reaction_rate * model.V == 0
    model.mole_balance_B = Constraint(rule=mole_balance_ruleB)

    def mole_balance_ruleC(model):
        return model.C_C_in * model.v0 - model.C_C * model.v0 + model.reaction_rate * model.V == 0
    model.mole_balance_C = Constraint(rule=mole_balance_ruleC)

    def mole_balance_ruleD(model):
        return model.C_D_in * model.v0 - model.C_D * model.v0 + model.reaction_rate * model.V == 0
    model.mole_balance_D = Constraint(rule=mole_balance_ruleD)

    def objective_rule(model):
        return (model.C_A_in - model.C_A) / model.C_A_in
    model.objective = Objective(rule=objective_rule, sense=1)

    solver = SolverFactory('ipopt')
    results = solver.solve(model)

    if (results.solver.status == 'ok') and (results.solver.termination_condition == 'optimal'):
        print('Found optimal solution.')
        print('C_A_in:', value(model.C_A_in), 'M')
        print('C_B_in:', value(model.C_B_in), 'M')
        print('V:', value(model.V), 'L')
        print('v0:', value(model.v0), 'L/s')
        print('T_in:', value(model.T_in), 'K')
    else: 
        print('No optimal solution found.')

    return model
