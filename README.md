
# SappyCSTR

The `SappyCSTR` package is designed to enhance the conversion of reactants in a Continuous Stirred Tank Reactor (CSTR) engaged in a saponification reaction. Utilizing the Python libraries `Pyomo` for mathematical optimization and `Thermo` for thermodynamic property calculations, this package focuses on optimizing key operational parameters such as inlet temperature, inlet concentration of reactants, reactor volume, and inlet volumetric flow rate. By providing a structured and efficient approach to adjust these parameters, `SappyCSTR` aims to achieve optimal reactant conversion, thereby enhancing the efficiency and output of chemical processes within CSTR systems.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

Run the following line in a shell. For example, in a Bash shell:

'pip install SappyCSTR'

## Usage 

The primary usage of `SappyCSTR` is to optimize reactor conditions for saponification reactions occuring in a CSTR in adiabatic, isobaric conditions.  

First, import the function like so:

from SappyCSTR.core import optimize_CSTR

The inputs necessary to use the function `optimize_CSTR` are a dictionary of initial conditions, a rate constant at a reference temperature (scalar), activation energy for the reaction (scaler), and a dictionary defining the compound names, their role (reactant or product), and their stoichiometry. 

Here is the structure for initial_conditions:  

initial_conditions = {
    'C_A_in': **value**  # Inlet Concentration of A (M)
    'C_B_in': **value**,  # Inlet Concentration of B (M)
     'v0': **value**,# volumwtric flow rate (L/s)
     'V': **value**, # Reactor volume (L)
}

And here is the structure for compound_names:  

compound_names = {
    'reactants': {
        'A': {'formula': **'condensed structural formula'**, 'stoichiometry': **value**},  
        'B': {'formula': **'condensed structural formula'**, 'stoichiometry': **value**},   
    },
    'products': {
        'C': {'formula': **'condensed structural formula'**, 'stoichiometry': **value**},  # Ethanol
        'D': {'formula': **'condensed structural formula'**, 'stoichiometry': **value**},   # Sodium Acetate
    }
}
