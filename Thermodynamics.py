import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import linprog, minimize, root_scalar
from dataclasses import asdict, dataclass, field, replace
import warnings

# Optionale Bibliotheken laden
try:
    import CoolProp.CoolProp as CP
    COOLPROP_AVAILABLE = True
except ImportError:
    COOLPROP_AVAILABLE = False
    warnings.warn("CoolProp ist nicht installiert. Verwende ideale Gasannahmen für kalte Gase.")

try:
    import cantera as ctpip 
    CANTERA_AVAILABLE = True
except ImportError:
    CANTERA_AVAILABLE = False
    warnings.warn("Cantera ist nicht installiert. Verwende vereinfachte Stöchiometrie für Verbrennung.")

# Konstanten
R_UNIV = 8314.4626  # J/(kmol*K)
P_REF = 101325.0    # [Pa]

# NASA-7-Polynome fuer Produktgase (GRI/NASA-Format, 200-6000 K).
# Die Form entspricht der in NASA CEA verwendeten Glenn-Polynomfamilie.
NASA7 = {
    "CO2": {
        "M": 44.0095,
        "T_mid": 1000.0,
        "low": [2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699548e-13, -4.83719697e4, 9.90105222],
        "high": [4.63659493, 2.74131985e-3, -9.95828531e-7, 1.60373011e-10, -9.16103468e-15, -4.90249341e4, -1.93489550],
    },
    "H2O": {
        "M": 18.01528,
        "T_mid": 1000.0,
        "low": [4.19864056, -2.03643410e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -3.02937267e4, -0.849032208],
        "high": [3.03399249, 2.17691804e-3, -1.64072518e-7, -9.70419870e-11, 1.68200992e-14, -3.00042971e4, 4.96677010],
    },
    "CO": {
        "M": 28.0101,
        "T_mid": 1000.0,
        "low": [3.57953347, -6.10353680e-4, 1.01681433e-6, 9.07005884e-10, -9.04424499e-13, -1.43440860e4, 3.50840928],
        "high": [2.71518561, 2.06252743e-3, -9.98825771e-7, 2.30053008e-10, -2.03647716e-14, -1.41518724e4, 7.81868772],
    },
    "H2": {
        "M": 2.01588,
        "T_mid": 1000.0,
        "low": [2.34433112, 7.98052075e-3, -1.94781510e-5, 2.01572094e-8, -7.37611761e-12, -917.935173, 0.683010238],
        "high": [3.33727920, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10, 2.00255376e-14, -950.158922, -3.20502331],
    },
    "O2": {
        "M": 31.9988,
        "T_mid": 1000.0,
        "low": [3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12, -1.06394356e3, 3.65767573],
        "high": [3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14, -1.08845772e3, 5.45323129],
    },
    "N2": {
        "M": 28.0134,
        "T_mid": 1000.0,
        "low": [3.53100528, -1.23660987e-4, -5.02999433e-7, 2.43530612e-9, -1.40881235e-12, -1.04697628e3, 2.96747468],
        "high": [2.95257626, 1.39690040e-3, -4.92631603e-7, 7.86010195e-11, -4.60755204e-15, -923.948688, 5.87188762],
    },
}

EQUILIBRIUM_SPECIES = ["CO2", "H2O", "CO", "H2", "O2", "N2"]
ELEMENT_ORDER = ["C", "H", "O", "N"]
SPECIES_ELEMENTS = {
    "CO2": {"C": 1, "O": 2},
    "H2O": {"H": 2, "O": 1},
    "CO": {"C": 1, "O": 1},
    "H2": {"H": 2},
    "O2": {"O": 2},
    "N2": {"N": 2},
}

ETHANOL_HF_LIQ = -277.69e6  # [J/kmol] bei 298.15 K
ETHANOL_CP_LIQ = 112400.0   # [J/(kmol*K)] einfache Sensible-Enthalpie-Naeherung
CH4_HF = -74.873e6          # [J/kmol] bei 298.15 K
CH4_CP = 35700.0            # [J/(kmol*K)] einfache Sensible-Enthalpie-Naeherung
BUTANE_M = 58.1222          # [kg/kmol]
BUTANE_HF_LIQ = -148.0e6    # [J/kmol] bei 298.15 K, n-Butan fluessig
BUTANE_CP_LIQ = 134000.0    # [J/(kmol*K)] einfache Sensible-Enthalpie-Naeherung
BUTANE_CP_GAS = 98490.0     # [J/(kmol*K)] bei 298.15 K
BUTANE_T_C = 425.125        # [K]
BUTANE_P_C = 3.796e6        # [Pa]
N2O_HF = 82.05e6            # [J/kmol] bei 298.15 K
N2O_CP = 38000.0            # [J/(kmol*K)] einfache Sensible-Enthalpie-Naeherung
T_REF = 298.15

FUEL_DATA = {
    "C2H5OH": {"formula": "C2H5OH", "M": 46.06844, "elements": {"C": 2.0, "H": 6.0, "O": 1.0, "N": 0.0}, "hf": ETHANOL_HF_LIQ, "cp": ETHANOL_CP_LIQ},
    "CH4": {"formula": "CH4", "M": 16.0425, "elements": {"C": 1.0, "H": 4.0, "O": 0.0, "N": 0.0}, "hf": CH4_HF, "cp": CH4_CP},
    "C4H10": {"formula": "C4H10", "M": BUTANE_M, "elements": {"C": 4.0, "H": 10.0, "O": 0.0, "N": 0.0}, "hf": BUTANE_HF_LIQ, "cp": BUTANE_CP_LIQ},
    "H2": {"formula": "H2", "M": 2.01588, "elements": {"C": 0.0, "H": 2.0, "O": 0.0, "N": 0.0}, "hf": 0.0, "cp": None},
}

@dataclass
class NozzleInput:
    """Eingabeparameter für die Düsenberechnung."""
    # Allgemein
    p_amb: float = 101325.0       # [Pa]
    T_amb: float = 293.15         # [K]
    mdot: float = 0.1             # [kg/s]
    p_c: float = 2e6              # [Pa]
    T_in: float = 293.15          # [K]
    eta_nozzle: float = 0.95      # [-] (Wirkungsgrad)
    C_d: float = 0.98             # [-] (Entladungskoeffizient)
    
    # Betriebsmodus: "cold_gas" oder "combustion"
    mode: str = "combustion"
    
    # Für cold_gas
    gas_medium: str = "O2"        # z.B. "O2", "Ethanol"
    
    # Für combustion
    fuel: str = "C2H5OH"          # Ethanol
    oxidizer: str = "Air"         # "Air", "O2" oder "N2O"
    phi: float = 1.0              # Äquivalenzverhältnis
    combustion_model: str = "equilibrium"  # "equilibrium" oder "complete"
    fuel_T: float = 293.15        # [K]
    oxidizer_T: float = 293.15    # [K]
    
    # Brennraumdimensionierung (Methode A oder B)
    chamber_method: str = "A"     # "A": Contraction Ratio, "B": L*
    contraction_ratio: float = 4.0
    L_star: float = 1.2           # [m]
    alpha_conv: float = 20.0      # [Grad]
    alpha_div: float = 12.0       # [Grad]

    # Duesenmodus: "design" berechnet A_t aus mdot, "given_geometry" nutzt Vorgaben
    nozzle_mode: str = "design"
    D_t_given: float = 0.004      # [m]
    D_e_given: float = 0.008      # [m]
    D_c_given: float = 0.012      # [m]
    L_c_given: float = 0.02       # [m]

    # Kammerquelle: "manual" nutzt Eingaben oben, "injector" koppelt Tanks/Injektor/Duese
    chamber_source: str = "manual"

    # Injektor/Tanks (vereinfachtes 0D-Modell)
    fuel_tank_p: float = 3.0e6
    fuel_tank_T: float = 293.15
    fuel_tank_V: float = 0.005
    oxidizer_tank_p: float = 5.0e6
    oxidizer_tank_T: float = 293.15
    oxidizer_tank_V: float = 0.005
    fuel_injector_d: float = 0.001
    fuel_injector_L: float = 0.002
    fuel_injector_count: float = 1.0
    fuel_injector_Cd: float = 0.77
    oxidizer_injector_d: float = 0.001
    oxidizer_injector_L: float = 0.002
    oxidizer_injector_count: float = 5.0
    oxidizer_injector_Cd: float = 0.77
    injector_pressure_margin: float = 1000.0

@dataclass
class GasState:
    """Thermodynamischer Zustand des Gases im Brennraum (Stagnationszustand)."""
    p: float
    T: float
    rho: float
    R_spec: float
    cp: float
    gamma: float
    M_w: float  # Molare Masse [kg/kmol]
    Z: float = 1.0  # Kompressibilitätsfaktor
    composition: dict = field(default_factory=dict)
    partial_pressures: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    
    @property
    def a(self) -> float:
        """Schallgeschwindigkeit [m/s]"""
        return np.sqrt(self.gamma * self.R_spec * self.T * self.Z)

@dataclass
class NozzleGeometry:
    """Geometrische Abmessungen der Düse."""
    D_c: float = 0.0
    D_t: float = 0.0
    D_e: float = 0.0
    A_c: float = 0.0
    A_t: float = 0.0
    A_e: float = 0.0
    L_c: float = 0.0
    L_conv: float = 0.0
    L_div: float = 0.0
    alpha_conv: float = 35.0  # [Grad]
    alpha_div: float = 12.0   # [Grad]

    @property
    def L_total(self) -> float:
        return self.L_conv + self.L_div

def _nasa_coeffs(species: str, T: float):
    data = NASA7[species]
    return data["low"] if T < data["T_mid"] else data["high"]

def nasa_cp_molar(species: str, T: float) -> float:
    a = _nasa_coeffs(species, T)
    cp_r = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4
    return R_UNIV * cp_r

def nasa_h_molar(species: str, T: float) -> float:
    a = _nasa_coeffs(species, T)
    h_rt = a[0] + a[1]*T/2 + a[2]*T**2/3 + a[3]*T**3/4 + a[4]*T**4/5 + a[5]/T
    return R_UNIV * T * h_rt

def nasa_s_molar(species: str, T: float) -> float:
    a = _nasa_coeffs(species, T)
    s_r = a[0]*np.log(T) + a[1]*T + a[2]*T**2/2 + a[3]*T**3/3 + a[4]*T**4/4 + a[6]
    return R_UNIV * s_r

def nasa_g_rt(species: str, T: float) -> float:
    return nasa_h_molar(species, T) / (R_UNIV * T) - nasa_s_molar(species, T) / R_UNIV

def normalize_oxidizer(name: str) -> str:
    key = name.strip().lower().replace(" ", "")
    aliases = {
        "air": "Air",
        "luft": "Air",
        "o2": "O2",
        "oxygen": "O2",
        "sauerstoff": "O2",
        "n2o": "N2O",
        "lachgas": "N2O",
        "nitrousoxide": "N2O",
        "distickstoffmonoxid": "N2O",
        "n2": "N2",
        "nitrogen": "N2",
        "stickstoff": "N2",
    }
    return aliases.get(key, name)

def normalize_fuel(name: str) -> str:
    key = name.strip().lower().replace(" ", "")
    aliases = {
        "c2h5oh": "C2H5OH",
        "ethanol": "C2H5OH",
        "ethylalcohol": "C2H5OH",
        "ch4": "CH4",
        "methane": "CH4",
        "methan": "CH4",
        "c4h10": "C4H10",
        "butane": "C4H10",
        "butan": "C4H10",
        "n-butane": "C4H10",
        "n-butan": "C4H10",
        "nbutane": "C4H10",
        "nbutan": "C4H10",
        "h2": "H2",
        "hydrogen": "H2",
        "wasserstoff": "H2",
    }
    return aliases.get(key, name)

def coolprop_fluid_name(fluid: str) -> str:
    fuel_name = normalize_fuel(fluid)
    oxidizer_name = normalize_oxidizer(fluid)
    names = {
        "C2H5OH": "Ethanol",
        "Ethanol": "Ethanol",
        "CH4": "Methane",
        "Methane": "Methane",
        "C4H10": "n-Butane",
        "Butane": "n-Butane",
        "H2": "Hydrogen",
        "Hydrogen": "Hydrogen",
        "N2O": "NitrousOxide",
        "O2": "Oxygen",
        "N2": "Nitrogen",
        "Air": "Air",
    }
    return names.get(fuel_name, names.get(oxidizer_name, fluid))

def is_gas_like_phase(phase: str) -> bool:
    return str(phase).lower() in {"gas", "supercritical", "supercritical_gas"}

def fuel_stoich_o2(fuel: str) -> float:
    data = FUEL_DATA[fuel]
    elements = data["elements"]
    return elements["C"] + elements["H"] / 4.0 - elements["O"] / 2.0

def fuel_enthalpy(fuel: str, T: float) -> float:
    if fuel == "H2":
        return nasa_h_molar("H2", T)
    data = FUEL_DATA[fuel]
    return data["hf"] + data["cp"] * (T - T_REF)

def thermochemistry_phi_bounds(inputs: NozzleInput):
    fuel = normalize_fuel(inputs.fuel)
    if fuel not in FUEL_DATA:
        return 0.2, 5.0
    data = FUEL_DATA[fuel]
    elements = data["elements"]
    stoich_o2 = fuel_stoich_o2(fuel)
    phi_max = 5.0
    carbon_oxygen_deficit = elements["C"] - elements["O"]
    if carbon_oxygen_deficit > 0 and stoich_o2 > 0:
        # Ohne kondensierten Kohlenstoff muss mindestens genug Sauerstoff fuer CO vorhanden sein.
        phi_max = min(phi_max, 0.98 * (2.0 * stoich_o2 / carbon_oxygen_deficit))
    return 0.2, max(0.2, phi_max)

def reactant_inventory(inputs: NozzleInput):
    fuel = normalize_fuel(inputs.fuel)
    if fuel not in FUEL_DATA:
        supported = ", ".join(FUEL_DATA.keys())
        raise ValueError(f"Die Gleichgewichtsrechnung unterstützt aktuell {supported} als Brennstoff.")

    fuel_data = FUEL_DATA[fuel]
    oxidizer = normalize_oxidizer(inputs.oxidizer)
    phi = max(float(inputs.phi), 1e-6)
    elements = dict(fuel_data["elements"])

    stoich_o2 = fuel_stoich_o2(fuel)
    if stoich_o2 <= 0 and oxidizer != "N2":
        raise ValueError(f"Für {fuel} konnte kein positiver stöchiometrischer Sauerstoffbedarf bestimmt werden.")

    if oxidizer == "Air":
        n_o2 = stoich_o2 / phi
        n_n2 = 3.76 * n_o2
        elements["O"] += 2.0 * n_o2
        elements["N"] += 2.0 * n_n2
        reactants = {fuel: 1.0, "O2": n_o2, "N2": n_n2}
        of_stoich = (stoich_o2 * (NASA7["O2"]["M"] + 3.76 * NASA7["N2"]["M"])) / fuel_data["M"]
        equation = f"{fuel_data['formula']} + {stoich_o2:.4g} O2 + {3.76 * stoich_o2:.4g} N2 -> stöchiometrische Produkte"
    elif oxidizer == "O2":
        n_o2 = stoich_o2 / phi
        elements["O"] += 2.0 * n_o2
        reactants = {fuel: 1.0, "O2": n_o2}
        of_stoich = (stoich_o2 * NASA7["O2"]["M"]) / fuel_data["M"]
        equation = f"{fuel_data['formula']} + {stoich_o2:.4g} O2 -> stöchiometrische Produkte"
    elif oxidizer == "N2O":
        n_n2o = 2.0 * stoich_o2 / phi
        elements["O"] += n_n2o
        elements["N"] += 2.0 * n_n2o
        reactants = {fuel: 1.0, "N2O": n_n2o}
        of_stoich = (2.0 * stoich_o2 * 44.0128) / fuel_data["M"]
        equation = f"{fuel_data['formula']} + {2.0 * stoich_o2:.4g} N2O -> stöchiometrische Produkte + N2"
    elif oxidizer == "N2":
        warnings.warn("Stickstoff (N2) liefert keinen Sauerstoff. Es findet keine Oxidation statt.")
        n_n2 = 1.0 / phi
        elements["N"] += 2.0 * n_n2
        reactants = {fuel: 1.0, "N2": n_n2}
        of_stoich = np.inf
        equation = f"{fuel_data['formula']} + N2 -> keine Oxidation ohne Sauerstoff"
    else:
        raise ValueError(f"Unbekannter Oxidator: {inputs.oxidizer}")

    return elements, reactants, {
        "fuel": fuel,
        "oxidizer": oxidizer,
        "stoich_of": float(of_stoich),
        "actual_of": float(of_stoich / phi) if np.isfinite(of_stoich) else np.inf,
        "equation_stoich": equation,
    }

def reactant_enthalpy(inputs: NozzleInput, reactants: dict) -> float:
    fuel = normalize_fuel(inputs.fuel)
    H = fuel_enthalpy(fuel, inputs.fuel_T)
    for species, moles in reactants.items():
        if species == fuel:
            continue
        if species == "N2O":
            H += moles * (N2O_HF + N2O_CP * (inputs.oxidizer_T - T_REF))
        else:
            H += moles * nasa_h_molar(species, inputs.oxidizer_T)
    return H

def element_matrix(species_list):
    matrix = []
    for element in ELEMENT_ORDER:
        row = []
        for species in species_list:
            row.append(SPECIES_ELEMENTS[species].get(element, 0.0))
        matrix.append(row)
    return np.array(matrix, dtype=float)

def complete_combustion_products(elements: dict):
    C = elements["C"]
    H = elements["H"]
    O = elements["O"]
    N = elements["N"]

    n_h2o = min(H / 2.0, max(O, 0.0))
    O_left = O - n_h2o
    n_co2 = min(C, max(O_left / 2.0, 0.0))
    O_left -= 2.0 * n_co2
    C_left = C - n_co2
    n_co = min(C_left, max(O_left, 0.0))
    O_left -= n_co
    C_left -= n_co
    n_h2 = max((H - 2.0 * n_h2o) / 2.0, 0.0)
    n_o2 = max(O_left / 2.0, 0.0)
    n_n2 = max(N / 2.0, 0.0)

    products = {"CO2": n_co2, "H2O": n_h2o, "CO": n_co, "H2": n_h2, "O2": n_o2, "N2": n_n2}
    if C_left > 1e-8:
        warnings.warn("Sehr fette Mischung: elementarer Kohlenstoff wird in diesem Modell nicht abgebildet.")
    return products

def equilibrium_moles(T: float, p: float, elements: dict, species_list=EQUILIBRIUM_SPECIES):
    A = element_matrix(species_list)
    b = np.array([elements[element] for element in ELEMENT_ORDER], dtype=float)
    feasible = linprog(np.ones(len(species_list)), A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not feasible.success:
        raise ValueError("Elementbilanz konnte für die Gleichgewichtsrechnung nicht erfüllt werden.")

    n0 = np.clip(feasible.x, 0.0, None)
    g0 = np.array([nasa_g_rt(species, T) for species in species_list])
    pressure_ratio = max(p / P_REF, 1e-12)

    def objective(n):
        n_safe = np.clip(n, 1e-30, None)
        n_total = np.sum(n_safe)
        y = n_safe / n_total
        return float(np.sum(n_safe * (g0 + np.log(np.clip(y * pressure_ratio, 1e-300, None)))))

    constraints = [{"type": "eq", "fun": lambda n, row=row, target=target: float(np.dot(row, n) - target)}
                   for row, target in zip(A, b)]
    result = minimize(
        objective,
        n0,
        method="SLSQP",
        bounds=[(0.0, None)] * len(species_list),
        constraints=constraints,
        options={"maxiter": 1000, "ftol": 1e-10, "disp": False},
    )

    if not result.success:
        residual = float(np.max(np.abs(A @ np.clip(result.x, 0.0, None) - b))) if result.x is not None else np.inf
        if residual > 1e-5 or result.x is None or not np.all(np.isfinite(result.x)):
            warnings.warn(
                "Gleichgewichts-Minimierung nicht belastbar konvergiert "
                f"({result.message}, Element-Restfehler={residual:.2e}). Nutze beste Näherung."
            )

    return {species: max(float(amount), 0.0) for species, amount in zip(species_list, result.x)}

def product_enthalpy(products: dict, T: float) -> float:
    return sum(moles * nasa_h_molar(species, T) for species, moles in products.items())

def product_properties(products: dict, T: float, p: float):
    total_moles = sum(products.values())
    if total_moles <= 0:
        raise ValueError("Keine Produktmole zur Auswertung vorhanden.")

    mole_fractions = {species: moles / total_moles for species, moles in products.items() if moles / total_moles > 1e-9}
    M_w_mix = sum(mole_fractions[species] * NASA7[species]["M"] for species in mole_fractions)
    cp_molar = sum(mole_fractions[species] * nasa_cp_molar(species, T) for species in mole_fractions)
    cp_mass = cp_molar / M_w_mix
    R_spec = R_UNIV / M_w_mix
    gamma = cp_mass / (cp_mass - R_spec)
    rho = p / (R_spec * T)
    partial_pressures = {species: y * p for species, y in mole_fractions.items()}
    return GasState(p, T, rho, R_spec, cp_mass, gamma, M_w_mix, 1.0, mole_fractions, partial_pressures)

def solve_enthalpy_temperature(balance, inputs: NozzleInput, metadata: dict, H_reactants: float):
    T_min = max(220.0, min(inputs.fuel_T, inputs.oxidizer_T, T_REF) - 80.0)
    T_max = 6500.0
    if T_min >= T_max:
        raise ValueError("Eintrittstemperaturen liegen ausserhalb des unterstuetzten Enthalpie-Suchbereichs.")
    anchor_points = np.array([300.0, 500.0, 800.0, 1200.0, 1800.0, 2500.0, 3200.0, 4200.0, 5200.0, 6000.0, T_max])
    grid = np.unique(np.concatenate([
        np.linspace(T_min, T_max, 36),
        np.clip(anchor_points, T_min, T_max),
    ]))

    values = []
    for T in grid:
        value = float(balance(T))
        values.append(value if np.isfinite(value) else np.nan)

    bracket = None
    clean = [(float(T), float(value)) for T, value in zip(grid, values) if np.isfinite(value)]
    for (T_low, f_low), (T_high, f_high) in zip(clean[:-1], clean[1:]):
        if abs(f_low) < 1e-6:
            return T_low, f_low, 0.0
        if f_low * f_high < 0:
            bracket = (T_low, T_high)
            break

    if bracket:
        T_c = root_scalar(balance, bracket=bracket, method="brentq").root
        residual = float(balance(T_c))
        return float(T_c), residual, abs(residual) / max(abs(H_reactants), 1.0)

    if not clean:
        raise ValueError("Enthalpiebilanz konnte in keinem Temperaturpunkt ausgewertet werden.")

    T_best, residual = min(clean, key=lambda item: abs(item[1]))
    scale = max(abs(H_reactants), max(abs(value) for _, value in clean), 1.0)
    relative_residual = abs(residual) / scale

    if relative_residual > 0.02:
        warnings.warn(
            "Keine physikalisch plausible Enthalpie-Nullstelle gefunden; nutze beste Temperatur-Naeherung. "
            f"Eingaben: Oxidator={metadata.get('oxidizer')}, phi={inputs.phi:.4g}, "
            f"p_c={inputs.p_c / 1e5:.3g} bar, Modell={inputs.combustion_model}, "
            f"Suchbereich={T_min:.0f}-{T_max:.0f} K, T_best={T_best:.0f} K, "
            f"relativer Restfehler={relative_residual:.2%}."
        )

    return float(T_best), float(residual), float(relative_residual)

def solve_combustion_equilibrium(inputs: NozzleInput):
    elements, reactants, metadata = reactant_inventory(inputs)
    H_reactants = reactant_enthalpy(inputs, reactants)
    p = inputs.p_c
    use_equilibrium = inputs.combustion_model.lower() == "equilibrium"
    cache = {}

    def products_at(T):
        key = round(float(T), 4)
        if key not in cache:
            if use_equilibrium:
                cache[key] = equilibrium_moles(T, p, elements)
            else:
                cache[key] = complete_combustion_products(elements)
        return cache[key]

    def balance(T):
        return product_enthalpy(products_at(T), T) - H_reactants

    T_c, enthalpy_residual, enthalpy_relative_residual = solve_enthalpy_temperature(balance, inputs, metadata, H_reactants)

    products = products_at(T_c)
    state = product_properties(products, T_c, p)
    state.metadata = {
        **metadata,
        "reactants_mol_per_mol_fuel": reactants,
        "products_mol_per_mol_fuel": products,
        "combustion_model": "Gibbs equilibrium (NASA-7 subset)" if use_equilibrium else "complete combustion",
        "product_species": EQUILIBRIUM_SPECIES,
        "enthalpy_residual_J_per_kmol_fuel": float(enthalpy_residual),
        "enthalpy_relative_residual": float(enthalpy_relative_residual),
    }
    return state

def orifice_area(diameter: float, count: float) -> float:
    return max(float(count), 0.0) * np.pi / 4 * max(float(diameter), 0.0)**2

def n2o_vapor_pressure(T: float) -> float:
    # Grobe Clausius-Clapeyron-Naehung, gefittet auf ca. 50.8 bar bei 293 K und kritischen Punkt.
    T_c = 309.57
    p_c = 7.245e6
    p_293 = 5.08e6
    if T >= T_c:
        return p_c
    B = np.log(p_c / p_293) / (1 / 293.15 - 1 / T_c)
    A = np.log(p_293) + B / 293.15
    return float(np.exp(A - B / T))

def butane_vapor_pressure(T: float) -> float:
    if T >= BUTANE_T_C:
        return BUTANE_P_C

    if T >= 272.66:
        A, B, C = 4.35576, 1175.581, -2.071
    else:
        A, B, C = 3.85002, 909.65, -36.146
    return float(1e5 * 10 ** (A - B / (T + C)))

def fluid_properties(fluid: str, p: float, T: float):
    fuel_name = normalize_fuel(fluid)
    name = normalize_oxidizer(fluid)
    warnings_list = []

    if COOLPROP_AVAILABLE:
        try:
            cp_name = coolprop_fluid_name(fluid)
            rho = CP.PropsSI("D", "P", p, "T", T, cp_name)
            cp = CP.PropsSI("C", "P", p, "T", T, cp_name)
            cv = CP.PropsSI("CVMASS", "P", p, "T", T, cp_name)
            M = CP.PropsSI("M", "P", p, "T", T, cp_name) * 1000
            phase = CP.PhaseSI("P", p, "T", T, cp_name)
            return {"rho": rho, "cp": cp, "gamma": cp / cv, "R_spec": R_UNIV / M, "M_w": M, "phase": phase, "warnings": warnings_list}
        except Exception as exc:
            warnings_list.append(f"CoolProp konnte {fluid} nicht auswerten ({exc}); nutze Fallback-Fluidmodell.")

    if fuel_name == "C2H5OH":
        rho = max(650.0, 789.0 - 0.85 * (T - 293.15))
        return {"rho": rho, "cp": 2440.0, "gamma": None, "R_spec": None, "M_w": 46.06844, "phase": "liquid", "warnings": warnings_list}

    if fuel_name == "C4H10":
        p_vap = butane_vapor_pressure(T)
        R_spec = R_UNIV / BUTANE_M
        gamma = BUTANE_CP_GAS / (BUTANE_CP_GAS - R_UNIV)
        if T > 0.98 * BUTANE_T_C:
            warnings_list.append("Butan liegt nahe am kritischen Punkt; Dichte und Injektorstroemung sind nur grob angenaehert.")
        if p > 1.05 * p_vap and T < BUTANE_T_C:
            rho = max(350.0, min(650.0, 584.0 - 1.1 * (T - 293.15)))
            phase = "liquid"
        elif p < 0.95 * p_vap:
            rho = p / (R_spec * T)
            phase = "gas"
        else:
            rho_liq = max(350.0, min(650.0, 584.0 - 1.1 * (T - 293.15)))
            rho_gas = p / (R_spec * T)
            rho = 0.7 * rho_liq + 0.3 * rho_gas
            phase = "two_phase"
            warnings_list.append("Butan ist nahe der Siedelinie; Injektor-Massenstrom ist als homogenes Zweiphasen-Fallback modelliert.")
        cp_mass = (BUTANE_CP_GAS if phase == "gas" else BUTANE_CP_LIQ) / BUTANE_M
        return {"rho": rho, "cp": cp_mass, "gamma": gamma, "R_spec": R_spec, "M_w": BUTANE_M, "phase": phase, "p_vap": p_vap, "warnings": warnings_list}

    if name == "N2O":
        T_c = 309.57
        p_vap = n2o_vapor_pressure(T)
        if T > 0.98 * T_c:
            warnings_list.append("N2O liegt nahe am kritischen Punkt; Dichte und Zweiphasenstroemung sind nur grob angenaehert.")
        if p > 1.05 * p_vap and T < T_c:
            rho = max(450.0, min(1000.0, 745.0 - 17.8 * (T - 293.15)))
            phase = "liquid"
        elif p < 0.95 * p_vap:
            R_spec = R_UNIV / 44.0128
            rho = p / (R_spec * T)
            phase = "gas"
        else:
            rho_liq = max(450.0, min(1000.0, 745.0 - 17.8 * (T - 293.15)))
            R_spec = R_UNIV / 44.0128
            rho_gas = p / (R_spec * T)
            rho = 0.65 * rho_liq + 0.35 * rho_gas
            phase = "two_phase"
            warnings_list.append("N2O ist nahe der Siedelinie; Injektor-Massenstrom ist als homogenes Zweiphasen-Fallback modelliert.")
        return {"rho": rho, "cp": 880.0, "gamma": 1.28, "R_spec": R_UNIV / 44.0128, "M_w": 44.0128, "phase": phase, "p_vap": p_vap, "warnings": warnings_list}

    gas_data = {
        "O2": (31.9988, 1.4),
        "Air": (28.97, 1.4),
        "N2": (28.0134, 1.4),
        "CH4": (16.0425, 1.31),
        "C4H10": (BUTANE_M, BUTANE_CP_GAS / (BUTANE_CP_GAS - R_UNIV)),
        "H2": (2.01588, 1.405),
    }
    M, gamma = gas_data.get(fuel_name, gas_data.get(name, (28.97, 1.4)))
    R_spec = R_UNIV / M
    return {"rho": p / (R_spec * T), "cp": gamma * R_spec / (gamma - 1), "gamma": gamma, "R_spec": R_spec, "M_w": M, "phase": "gas", "warnings": warnings_list}

def liquid_orifice_stream(label: str, rho: float, p_tank: float, T_tank: float, p_chamber: float, area: float, Cd: float):
    dp = max(p_tank - p_chamber, 0.0)
    mdot = Cd * area * np.sqrt(2 * rho * dp) if area > 0 and dp > 0 else 0.0
    velocity = mdot / (rho * area) if rho > 0 and area > 0 else 0.0
    return {
        "label": label,
        "mdot": float(mdot),
        "velocity": float(velocity),
        "T_exit": float(T_tank),
        "p_exit": float(p_chamber),
        "rho": float(rho),
        "phase": "liquid",
        "choked": False,
        "mach": None,
        "dp": float(dp),
    }

def gas_orifice_stream(label: str, props: dict, p_tank: float, T_tank: float, p_chamber: float, area: float, Cd: float):
    gamma = props["gamma"]
    R = props["R_spec"]
    if area <= 0 or p_tank <= p_chamber or gamma is None or R is None:
        return {"label": label, "mdot": 0.0, "velocity": 0.0, "T_exit": T_tank, "p_exit": p_chamber, "rho": props["rho"], "phase": props["phase"], "choked": False, "mach": 0.0, "dp": max(p_tank - p_chamber, 0.0)}

    critical_ratio = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    pressure_ratio = p_chamber / p_tank
    if pressure_ratio <= critical_ratio:
        mach = 1.0
        T_exit = T_tank * 2 / (gamma + 1)
        p_exit = p_tank * critical_ratio
        mdot = Cd * area * p_tank / np.sqrt(T_tank) * np.sqrt(gamma / R) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
        choked = True
    else:
        mach = np.sqrt(max(0.0, (2 / (gamma - 1)) * ((p_tank / p_chamber) ** ((gamma - 1) / gamma) - 1)))
        T_exit = T_tank / (1 + (gamma - 1) / 2 * mach**2)
        p_exit = p_chamber
        mdot = Cd * area * p_tank / np.sqrt(T_tank) * np.sqrt(gamma / R) * mach * (1 + (gamma - 1) / 2 * mach**2) ** (-(gamma + 1) / (2 * (gamma - 1)))
        choked = False
    velocity = mach * np.sqrt(gamma * R * T_exit)
    rho_exit = p_exit / (R * T_exit)
    return {
        "label": label,
        "mdot": float(mdot),
        "velocity": float(velocity),
        "T_exit": float(T_exit),
        "p_exit": float(p_exit),
        "rho": float(rho_exit),
        "phase": props["phase"],
        "choked": bool(choked),
        "mach": float(mach),
        "dp": float(max(p_tank - p_chamber, 0.0)),
    }

def injector_streams_for_pressure(inputs: NozzleInput, p_chamber: float):
    warnings_list = []
    fuel_area = orifice_area(inputs.fuel_injector_d, inputs.fuel_injector_count)
    ox_area = orifice_area(inputs.oxidizer_injector_d, inputs.oxidizer_injector_count)
    fuel_props = fluid_properties(inputs.fuel, inputs.fuel_tank_p, inputs.fuel_tank_T)
    ox_props = fluid_properties(inputs.oxidizer, inputs.oxidizer_tank_p, inputs.oxidizer_tank_T)
    warnings_list.extend(fuel_props.pop("warnings", []))
    warnings_list.extend(ox_props.pop("warnings", []))

    if is_gas_like_phase(fuel_props["phase"]):
        fuel = gas_orifice_stream("fuel", fuel_props, inputs.fuel_tank_p, inputs.fuel_tank_T, p_chamber, fuel_area, inputs.fuel_injector_Cd)
    else:
        fuel = liquid_orifice_stream("fuel", fuel_props["rho"], inputs.fuel_tank_p, inputs.fuel_tank_T, p_chamber, fuel_area, inputs.fuel_injector_Cd)
        fuel["phase"] = fuel_props["phase"]
    fuel.update({
        "diameter": float(inputs.fuel_injector_d),
        "length": float(inputs.fuel_injector_L),
        "count": float(inputs.fuel_injector_count),
        "Cd": float(inputs.fuel_injector_Cd),
        "area_total": float(fuel_area),
        "L_over_D": float(inputs.fuel_injector_L / inputs.fuel_injector_d) if inputs.fuel_injector_d > 0 else np.nan,
    })

    ox_name = normalize_oxidizer(inputs.oxidizer)
    if is_gas_like_phase(ox_props["phase"]):
        oxidizer = gas_orifice_stream("oxidizer", ox_props, inputs.oxidizer_tank_p, inputs.oxidizer_tank_T, p_chamber, ox_area, inputs.oxidizer_injector_Cd)
    else:
        oxidizer = liquid_orifice_stream("oxidizer", ox_props["rho"], inputs.oxidizer_tank_p, inputs.oxidizer_tank_T, p_chamber, ox_area, inputs.oxidizer_injector_Cd)
        oxidizer["phase"] = ox_props["phase"]
        if ox_name == "N2O" and "p_vap" in ox_props and p_chamber < ox_props["p_vap"]:
            warnings_list.append("N2O kann nach dem Injektor flashen, weil der Kammerdruck unter dem Sättigungsdruck liegt.")
    oxidizer.update({
        "diameter": float(inputs.oxidizer_injector_d),
        "length": float(inputs.oxidizer_injector_L),
        "count": float(inputs.oxidizer_injector_count),
        "Cd": float(inputs.oxidizer_injector_Cd),
        "area_total": float(ox_area),
        "L_over_D": float(inputs.oxidizer_injector_L / inputs.oxidizer_injector_d) if inputs.oxidizer_injector_d > 0 else np.nan,
    })

    return fuel, oxidizer, fuel_props, ox_props, warnings_list

def nozzle_capacity_for_pressure(inputs: NozzleInput, p_chamber: float, mdot_guess: float, phi_guess: float, fuel_T: float, oxidizer_T: float):
    local_inputs = replace(inputs, p_c=p_chamber, mdot=max(mdot_guess, 1e-9), phi=max(phi_guess, 1e-6), fuel_T=fuel_T, oxidizer_T=oxidizer_T, nozzle_mode="given_geometry")
    try:
        state = ThermoModel(local_inputs).get_state()
    except Exception as exc:
        if local_inputs.mode == "combustion" and local_inputs.combustion_model.lower() == "equilibrium":
            warnings.warn(
                "Gleichgewichtsmodell konnte fuer diese Injektor-Mischung nicht ausgewertet werden "
                f"({exc}); nutze vollstaendige Verbrennung als Fallback fuer die Massenbilanz."
            )
            local_inputs = replace(local_inputs, combustion_model="complete")
            state = ThermoModel(local_inputs).get_state()
        else:
            raise
    solver = NozzleSolver(local_inputs, state)
    A_t = np.pi / 4 * local_inputs.D_t_given**2
    return solver.compute_mdot_for_throat(A_t), state

def solve_injector_coupled(inputs: NozzleInput):
    if inputs.nozzle_mode != "given_geometry":
        inputs = replace(inputs, nozzle_mode="given_geometry")

    stoich = reactant_inventory(inputs)[2]["stoich_of"] if inputs.mode == "combustion" else 1.0
    phi_min, phi_max = thermochemistry_phi_bounds(inputs)
    upper = min(inputs.fuel_tank_p, inputs.oxidizer_tank_p) - max(inputs.injector_pressure_margin, 1.0)
    lower = max(inputs.p_amb * 1.01, 1000.0)
    if upper <= lower:
        raise ValueError("Tankdruecke muessen ueber Umgebungs-/Kammerdruck liegen, damit der Injektor stroemen kann.")

    cached = {}

    def evaluate(p_chamber):
        if p_chamber in cached:
            return cached[p_chamber]
        fuel, ox, fuel_props, ox_props, warnings_list = injector_streams_for_pressure(inputs, p_chamber)
        mdot_total = fuel["mdot"] + ox["mdot"]
        actual_of = ox["mdot"] / max(fuel["mdot"], 1e-12)
        phi_raw = stoich / actual_of if actual_of > 0 and np.isfinite(stoich) else inputs.phi
        if not np.isfinite(phi_raw) or phi_raw <= 0:
            phi_raw = inputs.phi
            warnings_list.append("Injektor-phi konnte nicht direkt bestimmt werden; nutze den bisherigen phi-Wert als Naeherung.")
        thermo_phi = float(np.clip(phi_raw, phi_min, phi_max))
        if abs(thermo_phi - phi_raw) > max(1e-6, 1e-3 * abs(phi_raw)):
            warnings_list.append(
                f"Auslegungswarnung: Der Injektor liefert phi={phi_raw:.3g} und damit eine Mischung ausserhalb "
                f"des reduzierten Thermochemie-Bereichs ({phi_min:.2g}-{phi_max:.2g}). "
                f"Die Brennraum-Thermochemie wird mit phi={thermo_phi:.3g} angenaehert. "
                "Fuer eine realistischere Rechnung Oxidatorbohrung/Tankdruck erhoehen oder Kraftstoffbohrung verkleinern."
            )
        mdot_out, state = nozzle_capacity_for_pressure(inputs, p_chamber, mdot_total, thermo_phi, fuel["T_exit"], ox["T_exit"])
        value = {
            "fuel": fuel,
            "oxidizer": ox,
            "mdot_total": mdot_total,
            "actual_of": actual_of,
            "phi": float(phi_raw),
            "thermochemistry_phi": thermo_phi,
            "mdot_out": mdot_out,
            "state": state,
            "warnings": warnings_list,
            "fuel_props": fuel_props,
            "oxidizer_props": ox_props,
        }
        cached[p_chamber] = value
        return value

    def residual(p_chamber):
        value = evaluate(float(p_chamber))
        return value["mdot_total"] - value["mdot_out"]

    grid = np.linspace(lower, upper, 28)
    values = [residual(p) for p in grid]
    bracket = None
    for p_low, p_high, f_low, f_high in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if f_low == 0 or f_low * f_high < 0:
            bracket = (p_low, p_high)
            break

    warnings_list = []
    if bracket:
        p_solution = root_scalar(residual, bracket=bracket, method="brentq").root
    else:
        idx = int(np.argmin(np.abs(values)))
        p_solution = float(grid[idx])
        warnings_list.append("Keine exakte Injektor-Duese-Massenbilanz gefunden; nutze beste Naeherung im Druckbereich.")

    final = evaluate(p_solution)
    warnings_list.extend(final["warnings"])
    updated_inputs = replace(
        inputs,
        chamber_source="injector",
        nozzle_mode="given_geometry",
        p_c=float(p_solution),
        mdot=float(final["mdot_total"]),
        phi=float(final["thermochemistry_phi"]),
        fuel_T=float(final["fuel"]["T_exit"]),
        oxidizer_T=float(final["oxidizer"]["T_exit"]),
    )

    injector_result = {
        "p_chamber": float(p_solution),
        "mdot_total": float(final["mdot_total"]),
        "mdot_nozzle_capacity": float(final["mdot_out"]),
        "actual_of": float(final["actual_of"]),
        "phi": float(final["phi"]),
        "thermochemistry_phi": float(final["thermochemistry_phi"]),
        "fuel": final["fuel"],
        "oxidizer": final["oxidizer"],
        "fuel_tank_mass_est": float(final["fuel_props"]["rho"] * inputs.fuel_tank_V),
        "oxidizer_tank_mass_est": float(final["oxidizer_props"]["rho"] * inputs.oxidizer_tank_V),
        "fuel_burn_time_est": float(final["fuel_props"]["rho"] * inputs.fuel_tank_V / max(final["fuel"]["mdot"], 1e-12)),
        "oxidizer_burn_time_est": float(final["oxidizer_props"]["rho"] * inputs.oxidizer_tank_V / max(final["oxidizer"]["mdot"], 1e-12)),
        "warnings": list(dict.fromkeys(warnings_list)),
    }
    return updated_inputs, injector_result

class ThermoModel:
    """Bestimmt die thermischen Zustände und Gasstoffeigenschaften."""
    
    def __init__(self, inputs: NozzleInput):
        self.inputs = inputs

    def get_state(self) -> GasState:
        if self.inputs.mode == "cold_gas":
            return self._compute_cold_gas_state()
        elif self.inputs.mode == "combustion":
            return self._compute_combustion_state()
        else:
            raise ValueError(f"Unbekannter Modus: {self.inputs.mode}")

    def _compute_cold_gas_state(self) -> GasState:
        p = self.inputs.p_c
        T = self.inputs.T_in
        gas = self.inputs.gas_medium
        fuel_name = normalize_fuel(gas)
        medium_name = fuel_name if fuel_name in FUEL_DATA else normalize_oxidizer(gas)
        
        if COOLPROP_AVAILABLE:
            try:
                cp_gas = coolprop_fluid_name(gas)
                if fuel_name in {"C2H5OH", "C4H10"}:
                    phase = CP.PhaseSI('P', p, 'T', T, cp_gas)
                    if not is_gas_like_phase(phase):
                        warnings.warn(f"Achtung: {gas} liegt bei {p/1e5} bar und {T} K nicht gasförmig vor! "
                                      f"(Phase: {phase}). Zweiphasenströmung wird nicht modelliert.")
                
                rho = CP.PropsSI('D', 'P', p, 'T', T, cp_gas)
                cp = CP.PropsSI('C', 'P', p, 'T', T, cp_gas)
                cv = CP.PropsSI('CVMASS', 'P', p, 'T', T, cp_gas)
                M_w = CP.PropsSI('M', 'P', p, 'T', T, cp_gas) * 1000  # kg/kmol
                Z = CP.PropsSI('Z', 'P', p, 'T', T, cp_gas)
                gamma = cp / cv
                R_spec = R_UNIV / M_w
                return GasState(p, T, rho, R_spec, cp, gamma, M_w, Z, {medium_name: 1.0})
            except Exception as e:
                warnings.warn(f"CoolProp Fehler ({e}). Nutze idealisiertes Fallback.")
        
        # Fallback (Ideal Gas)
        fallback = {
            "O2": (32.0, 1.4),
            "Air": (28.97, 1.4),
            "N2": (28.0134, 1.4),
            "N2O": (44.0128, 1.28),
            "Ethanol": (46.07, 1.13),
            "C2H5OH": (46.07, 1.13),
            "C4H10": (BUTANE_M, BUTANE_CP_GAS / (BUTANE_CP_GAS - R_UNIV)),
        }
        M_w, gamma = fallback.get(medium_name, fallback.get(gas, (32.0, 1.4)))
        R_spec = R_UNIV / M_w
        rho = p / (R_spec * T)
        cp = gamma * R_spec / (gamma - 1)
        return GasState(p, T, rho, R_spec, cp, gamma, M_w, 1.0, {medium_name: 1.0})

    def _compute_combustion_state(self) -> GasState:
        """Bestimmt adiabate Flammentemperatur und Stoffwerte."""
        if CANTERA_AVAILABLE:
            try:
                pass
            except Exception:
                pass

        return solve_combustion_equilibrium(self.inputs)

class NozzleSolver:
    """Löst die 1D isentropen Strömungsgleichungen."""
    
    def __init__(self, inputs: NozzleInput, state: GasState):
        self.inputs = inputs
        self.state = state

    def check_choking(self) -> bool:
        """Prüft, ob die Strömung den kritischen Zustand erreicht."""
        gamma = self.state.gamma
        p_star = self.state.p * (2 / (gamma + 1)) ** (gamma / (gamma - 1))
        if self.inputs.p_amb > p_star:
            warnings.warn(f"WARNUNG: p_amb ({self.inputs.p_amb} Pa) > p_star ({p_star:.1f} Pa). "
                          "Strömung ist im Hals nicht kritisch (unchoked)!")
            return False
        return True

    def compute_throat(self) -> float:
        """Berechnet den engsten Querschnitt A_t."""
        return self.inputs.mdot / self.choked_mass_flux()

    def choked_mass_flux(self) -> float:
        """Massenstrom pro Halsfläche für kritische Strömung."""
        gamma = self.state.gamma
        R = self.state.R_spec
        T0 = self.state.T
        p0 = self.state.p
        Cd = self.inputs.C_d
        
        term1 = Cd * p0 / np.sqrt(T0) * np.sqrt(gamma / R)
        term2 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
        return term1 * term2

    def compute_mdot_for_throat(self, A_t: float) -> float:
        return A_t * self.choked_mass_flux()

    def solve_exit_mach(self) -> float:
        """Löst iterativ nach der Austrittsmachzahl bei p_e = p_amb (Idealexpansion)."""
        gamma = self.state.gamma
        p_ratio = self.state.p / self.inputs.p_amb  # p0 / p_e
        # p0/p = (1 + (gamma-1)/2 M^2)^(gamma/(gamma-1))
        # Umgestellt nach M:
        M_e = np.sqrt((2 / (gamma - 1)) * (p_ratio ** ((gamma - 1) / gamma) - 1))
        return M_e

    def pressure_from_mach(self, M: float) -> float:
        gamma = self.state.gamma
        return self.state.p / (1 + (gamma - 1) / 2 * M**2) ** (gamma / (gamma - 1))

    def temperature_from_mach(self, M: float) -> float:
        gamma = self.state.gamma
        return self.state.T / (1 + (gamma - 1) / 2 * M**2)

    def velocity_from_mach(self, M: float) -> float:
        T = self.temperature_from_mach(M)
        return M * np.sqrt(self.state.gamma * self.state.R_spec * T)

    def area_ratio(self, M: float) -> float:
        """Berechnet das Flächenverhältnis A/A* für eine gegebene Machzahl."""
        gamma = self.state.gamma
        if M < 1e-5: return np.inf
        term1 = 2 / (gamma + 1)
        term2 = 1 + (gamma - 1) / 2 * M**2
        return (1 / M) * (term1 * term2) ** ((gamma + 1) / (2 * (gamma - 1)))

    def mach_from_area_ratio(self, area_ratio: float, supersonic: bool = True) -> float:
        """Findet die Machzahl für ein gegebenes Flächenverhältnis."""
        def obj_func(M):
            return self.area_ratio(M) - area_ratio
        
        try:
            if supersonic:
                sol = root_scalar(obj_func, bracket=[1.0001, 20.0], method='brentq')
            else:
                sol = root_scalar(obj_func, bracket=[1e-4, 0.9999], method='brentq')
            return sol.root
        except ValueError:
            return 1.0

class GeometryBuilder:
    """Erzeugt die Düsenkontur."""
    
    def __init__(self, inputs: NozzleInput, geo: NozzleGeometry):
        self.inputs = inputs
        self.geo = geo

    def build(self):
        self.geo.alpha_conv = self.inputs.alpha_conv
        self.geo.alpha_div = self.inputs.alpha_div

        if self.inputs.nozzle_mode == "given_geometry":
            self.geo.D_t = self.inputs.D_t_given
            self.geo.D_e = self.inputs.D_e_given
            self.geo.D_c = self.inputs.D_c_given
            self.geo.L_c = self.inputs.L_c_given
            self.geo.A_t = np.pi / 4 * self.geo.D_t**2
            self.geo.A_e = np.pi / 4 * self.geo.D_e**2
            self.geo.A_c = np.pi / 4 * self.geo.D_c**2
            self.geo.L_conv = (self.geo.D_c - self.geo.D_t) / (2 * np.tan(np.radians(self.geo.alpha_conv)))
            self.geo.L_div = (self.geo.D_e - self.geo.D_t) / (2 * np.tan(np.radians(self.geo.alpha_div)))
            return

        # Durchmesser
        self.geo.D_t = 2 * np.sqrt(self.geo.A_t / np.pi)
        self.geo.D_e = 2 * np.sqrt(self.geo.A_e / np.pi)
        
        # Brennraum
        if self.inputs.chamber_method == "A":
            self.geo.A_c = self.geo.A_t * self.inputs.contraction_ratio
            self.geo.D_c = 2 * np.sqrt(self.geo.A_c / np.pi)
            self.geo.L_c = self.geo.D_c  # Annahme: L_c = D_c für Methode A
        else: # Methode B
            V_c = self.inputs.L_star * self.geo.A_t
            # Annahme D_c für L_star Methode
            self.geo.D_c = self.geo.D_t * 3  # Heuristik
            self.geo.A_c = np.pi/4 * self.geo.D_c**2
            self.geo.L_c = V_c / self.geo.A_c
            
        # Längen (Konisch)
        self.geo.L_conv = (self.geo.D_c - self.geo.D_t) / (2 * np.tan(np.radians(self.geo.alpha_conv)))
        self.geo.L_div = (self.geo.D_e - self.geo.D_t) / (2 * np.tan(np.radians(self.geo.alpha_div)))

    def generate_contour(self, num_points=100) -> pd.DataFrame:
        """Erzeugt x-r-Koordinaten."""
        x_c = np.linspace(-self.geo.L_c, 0, num_points//3)
        r_c = np.full_like(x_c, self.geo.D_c / 2)
        
        x_conv = np.linspace(0, self.geo.L_conv, num_points//3)
        r_conv = np.linspace(self.geo.D_c / 2, self.geo.D_t / 2, num_points//3)
        
        x_div = np.linspace(self.geo.L_conv, self.geo.L_total, num_points//3)
        r_div = np.linspace(self.geo.D_t / 2, self.geo.D_e / 2, num_points//3)
        
        x = np.concatenate([x_c, x_conv, x_div])
        r = np.concatenate([r_c, r_conv, r_div])
        
        return pd.DataFrame({"x [m]": x, "r [m]": r})

def print_disclaimer():
    print("="*60)
    print("WICHTIGER HINWEIS: THEORETISCHE VORABSCHÄTZUNG")
    print("="*60)
    print("- Die Geometrie ist eine 1D-Vorabschätzung (Isentrope Strömung).")
    print("- Brennraumdurchmesser/Länge sind ohne Zusatzannahmen nicht eindeutig.")
    print("- Realgas- und Gleichgewichtsrechnungen (Dissoziation) müssen z.B. mit NASA CEA validiert werden.")
    print("- Die Verbrennungsrechnung bildet komplexe Kinetik, Reibung, Wärmeübergang und Grenzschichten nicht ab.")
    print("- Nicht für fertigungstaugliche oder sicherheitskritische Auslegungen verwenden!")
    print("="*60)

def _save_contour_plot(contour_df, geo, inputs, output_base, show_plot=False):
    plt.figure(figsize=(10, 4))
    plt.plot(contour_df['x [m]'], contour_df['r [m]'], 'b-', linewidth=2)
    plt.plot(contour_df['x [m]'], -contour_df['r [m]'], 'b-', linewidth=2)
    plt.axvline(0, color='gray', linestyle='--')
    plt.axvline(geo.L_conv, color='r', linestyle='--', label="Düsenhals")
    plt.fill_between(contour_df['x [m]'], contour_df['r [m]'], -contour_df['r [m]'], color='lightblue', alpha=0.3)
    plt.title(f"Rotationssymmetrische Düsenkontur ({inputs.mode})")
    plt.xlabel("Länge x [m]")
    plt.ylabel("Radius r [m]")
    plt.legend()
    plt.axis('equal')
    plt.grid(True)
    plt.savefig(f"{output_base}.png", dpi=200, bbox_inches="tight")
    if show_plot and "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close()

def compute_flowfield(contour_df: pd.DataFrame, geo: NozzleGeometry, solver: NozzleSolver):
    rows = []
    A_t = max(geo.A_t, 1e-18)

    for _, point in contour_df.iterrows():
        x = float(point["x [m]"])
        r = float(point["r [m]"])
        area = max(np.pi * r**2, A_t)
        ratio = max(area / A_t, 1.0)

        if x < geo.L_conv:
            M = solver.mach_from_area_ratio(ratio, supersonic=False)
        elif abs(ratio - 1.0) < 1e-6:
            M = 1.0
        else:
            M = solver.mach_from_area_ratio(ratio, supersonic=True)

        T = solver.temperature_from_mach(M)
        p = solver.pressure_from_mach(M)
        v_ideal = solver.velocity_from_mach(M)
        v = np.sqrt(solver.inputs.eta_nozzle) * v_ideal if x >= geo.L_conv else v_ideal
        rows.append({
            "x [m]": x,
            "r [m]": r,
            "A [m2]": float(area),
            "M [-]": float(M),
            "p [Pa]": float(p),
            "T [K]": float(T),
            "v [m/s]": float(v),
        })

    return pd.DataFrame(rows)

def calculate_nozzle(config_dict=None, num_points=100, write_files=False, make_plot=False, show_plot=False):
    """Berechnet Zustand, Performance, Geometrie und Kontur ohne Konsolenausgabe."""
    if config_dict:
        inputs = NozzleInput(**config_dict)
    else:
        inputs = NozzleInput()

    injector_result = None
    if inputs.chamber_source == "injector":
        inputs, injector_result = solve_injector_coupled(inputs)

    # 1. Thermodynamik
    thermo = ThermoModel(inputs)
    state = thermo.get_state()
    
    # 2. Strömungsmechanik
    solver = NozzleSolver(inputs, state)
    is_choked = solver.check_choking()

    # Zustand am Austritt
    if inputs.nozzle_mode == "given_geometry":
        geo = NozzleGeometry()
        builder = GeometryBuilder(inputs, geo)
        builder.build()
        A_t = geo.A_t
        A_e = geo.A_e
        A_ratio = A_e / A_t if A_t > 0 else np.nan
        if A_ratio < 1.0:
            warnings.warn("A_e/A_t ist kleiner als 1. Die Austrittsmachzahl wird auf den Halszustand begrenzt.")
            M_e = 1.0
        else:
            M_e = solver.mach_from_area_ratio(A_ratio, supersonic=is_choked)
        mdot_effective = solver.compute_mdot_for_throat(A_t)
        p_e = solver.pressure_from_mach(M_e)
    else:
        A_t = solver.compute_throat()
        M_e = solver.solve_exit_mach()
        A_ratio = solver.area_ratio(M_e)
        A_e = A_t * A_ratio
        p_e = inputs.p_amb
        mdot_effective = inputs.mdot

        # 3. Geometrie
        geo = NozzleGeometry(A_t=A_t, A_e=A_e)
        builder = GeometryBuilder(inputs, geo)
        builder.build()

    T_e = solver.temperature_from_mach(M_e)
    v_e_ideal = solver.velocity_from_mach(M_e)
    v_e = np.sqrt(inputs.eta_nozzle) * v_e_ideal
    pressure_thrust = (p_e - inputs.p_amb) * A_e
    thrust = mdot_effective * v_e + pressure_thrust
    
    # 4. Kontur
    contour_df = builder.generate_contour()
    if num_points != 100:
        contour_df = builder.generate_contour(num_points)
    flowfield_df = compute_flowfield(contour_df, geo, solver)

    result = {
        "inputs": asdict(inputs),
        "state": {
            "p": float(state.p),
            "T": float(state.T),
            "rho": float(state.rho),
            "R_spec": float(state.R_spec),
            "cp": float(state.cp),
            "gamma": float(state.gamma),
            "M_w": float(state.M_w),
            "Z": float(state.Z),
            "a": float(state.a),
            "composition": {key: float(value) for key, value in state.composition.items()},
            "partial_pressures": {key: float(value) for key, value in state.partial_pressures.items()},
            "metadata": state.metadata,
        },
        "performance": {
            "is_choked": bool(is_choked),
            "mdot": float(mdot_effective),
            "M_e": float(M_e),
            "A_ratio": float(A_ratio),
            "p_e": float(p_e),
            "T_e": float(T_e),
            "v_e_ideal": float(v_e_ideal),
            "v_e": float(v_e),
            "pressure_thrust": float(pressure_thrust),
            "thrust": float(thrust),
        },
        "geometry": {
            "D_c": float(geo.D_c),
            "D_t": float(geo.D_t),
            "D_e": float(geo.D_e),
            "A_c": float(geo.A_c),
            "A_t": float(geo.A_t),
            "A_e": float(geo.A_e),
            "L_c": float(geo.L_c),
            "L_conv": float(geo.L_conv),
            "L_div": float(geo.L_div),
            "L_total": float(geo.L_total),
            "alpha_conv": float(geo.alpha_conv),
            "alpha_div": float(geo.alpha_div),
        },
        "contour": contour_df.to_dict(orient="records"),
        "flowfield": flowfield_df.to_dict(orient="records"),
        "injector": injector_result,
        "files": {},
    }

    if write_files:
        output_base = f"nozzle_contour_{inputs.mode}"
        contour_df.to_csv(f"{output_base}.csv", index=False)
        result["files"]["csv"] = f"{output_base}.csv"

        if make_plot:
            _save_contour_plot(contour_df, geo, inputs, output_base, show_plot=show_plot)
            result["files"]["png"] = f"{output_base}.png"

    return result

def _print_result(result):
    inputs = result["inputs"]
    state = result["state"]
    performance = result["performance"]
    geo = result["geometry"]

    print("\n--- Thermodynamischer Zustand (Brennkammer) ---")
    print(f"p_c       : {state['p'] / 1e5:.2f} bar")
    print(f"T_c       : {state['T']:.2f} K")
    print(f"R_mix     : {state['R_spec']:.2f} J/(kg*K)")
    print(f"Molar Mass: {state['M_w']:.2f} kg/kmol")
    print(f"Gamma (γ) : {state['gamma']:.4f}")
    print(f"cp        : {state['cp']:.2f} J/(kg*K)")
    print(f"Z         : {state['Z']:.4f}")
    print("Zusammensetzung:", ", ".join([f"{k}: {v:.2%}" for k, v in state["composition"].items()]))

    print("\n--- Düsen-Performance ---")
    print(f"Choked?   : {'Ja' if performance['is_choked'] else 'Nein'}")
    print(f"mdot      : {performance['mdot']:.4f} kg/s")
    print(f"Mach Exit : {performance['M_e']:.2f}")
    print(f"eta_nozzle: {inputs['eta_nozzle']:.3f}")
    print(f"p_exit    : {performance['p_e']:.1f} Pa")
    print(f"v_exit    : {performance['v_e']:.2f} m/s")
    print(f"T_exit    : {performance['T_e']:.2f} K")
    print(f"Druckschub: {performance['pressure_thrust']:.2f} N")
    print(f"Schub (F) : {performance['thrust']:.2f} N")

    print("\n--- Geometrie ---")
    print(f"D_c (Kammer) : {geo['D_c']*1000:.1f} mm")
    print(f"D_t (Hals)   : {geo['D_t']*1000:.1f} mm")
    print(f"D_e (Austritt): {geo['D_e']*1000:.1f} mm")
    print(f"L_c (Kammer) : {geo['L_c']*1000:.1f} mm")
    print(f"L_conv       : {geo['L_conv']*1000:.1f} mm (Winkel: {geo['alpha_conv']}°)")
    print(f"L_div        : {geo['L_div']*1000:.1f} mm (Winkel: {geo['alpha_div']}°)")
    print(f"L_total_Düse : {geo['L_total']*1000:.1f} mm")

def main(config_dict=None):
    if config_dict:
        inputs = NozzleInput(**config_dict)
    else:
        inputs = NozzleInput()

    print(f"\nStarte Auslegung Modus: {inputs.mode.upper()}")

    result = calculate_nozzle(asdict(inputs), write_files=True, make_plot=True)

    # 4. Ausgabe
    _print_result(result)
    
    print_disclaimer()

# ==========================================
# Beispielaufrufe
# ==========================================
if __name__ == "__main__":
    # Beispiel 1: Stöchiometrische Verbrennung Ethanol-Luft
    print("\n" + "#"*40 + "\nBEISPIEL 1: ETHANOL + LUFT\n" + "#"*40)
    config_combustion = {
        "p_amb": 101325, "T_amb": 293.15,
        "mdot": 0.5, "p_c": 15e5, "mode": "combustion",
        "fuel": "C2H5OH", "oxidizer": "Air", "phi": 1.0,
        "chamber_method": "A", "contraction_ratio": 5.0
    }
    main(config_combustion)

    # Beispiel 2: Kaltes Gas (Sauerstoff)
    print("\n" + "#"*40 + "\nBEISPIEL 2: KALTGAS SAUERSTOFF\n" + "#"*40)
    config_cold = {
        "p_amb": 101325, "T_amb": 293.15,
        "mdot": 0.1, "p_c": 30e5, "T_in": 300, "mode": "cold_gas",
        "gas_medium": "O2",
        "chamber_method": "A", "contraction_ratio": 3.0
    }
    main(config_cold)
