import yaml
import sympy


def load_constants(path):
    """Parses a yaml file and builds a dictionary that maps SymPy symbols to
    floats."""
    with open(path, 'r') as f:
        constant_values_dict = yaml.load(f)

    for k, v in constant_values_dict.items():
        constant_values_dict[sympy.Symbol(k)] = v
        del constant_values_dict[k]

    return constant_values_dict


def map_values_to_autolev_symbols(constants):

    d = {}
    d['TrunkMass'] = constants['ma']
    d['TrunkInertia'] = constants['ia']
    d['TrunkCMy'] = constants['ya']
    d['ThighMass'] = constants['mb']
    d['ThighInertia'] = constants['ib']
    d['ThighCMy'] = constants['yb']
    d['ThighLen'] = constants['lb']
    d['ShankMass'] = constants['mc']
    d['ShankInertia'] = constants['ic']
    d['ShankCMy'] = constants['yc']
    d['ShankLen'] = constants['lc']
    d['FootMass'] = constants['md']
    d['FootInertia'] = constants['id']
    d['FootCMx'] = constants['xd']
    d['FootCMy'] = constants['yd']
    d['ContactY'] = constants['fyd']
    d['ContactHeelX'] = constants['hxd']
    d['ContactToeX'] = constants['txd']
    d['ContactStiff'] = constants['kc']
    d['ContactDamp'] = constants['cc']
    d['ContactV0'] = constants['vs']
    d['ContactFric'] = constants['mu']

    return d

# Set initial coniditions: 18 states

def controller(x, t):
    """Outputs the joint torques and hand of god as a function of state."""

    frequency = 1.0 # frequency of gait [hz]
    period = 1.0 / frequency

    time_in_gait_cycle = t / period

    # subtract the whole number and the remainder is the percent gait cycle

    # figure out phi from t and the gait frequency

    # interpolate the m_stra and K data based on phi

    return m_star[phi] + np.dot(K[phi], x)
