#!/usr/bin/env python3

import argparse, copy, itertools, os, pickle, shutil, subprocess, sys, time, xmltodict
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.linalg import cholesky
from scipy.optimize import fsolve, linprog
from scipy.stats import beta, chi2, gaussian_kde, kstest, multivariate_normal, shapiro
from scipy.interpolate import griddata
from sklearn.cluster import MeanShift, estimate_bandwidth
from time import sleep
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.text import TextPath
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Affine2D
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from itertools import product, combinations, chain
from sympy import factorint



# constants
k_constant = 8.314462618e-3



# for nice plots
converter = {
                'energy': {
                    'name': r'$\Delta H_{vap}$',
                    'unit': r'kJ mol$^{-1}$'
                },
                'volume': {
                    'name': r'$V_{M}$',
                    'unit': r'Å$^{3}$'
                },
                'rXX': {
                    'name': None,
                    'unit': r'Å'
                },
                'gXX': {
                    'name': None,
                    'unit': r'-'
                },
                'hbonddist': {
                    'name': r'$\langle r_{HB} \rangle$',
                    'unit': r'Å'
                },
                'hbondangle': {
                    'name': r'$\langle \theta_{HB} \rangle$',
                    'unit': r'deg'
                },
                'dielectric': {
                    'name': r'$\epsilon_{0}$',
                    'unit': r'-'
                },
                'diffusion': {
                    'name': r'$D$',
                    'unit': r'cm$^{2}$ s$^{-1}$'
                },
                'isobaricheatcapacity': {
                    'name': r'$C_{p}$',
                    'unit': r'J mol$^{-1}$ K$^{-1}$'
                },
                'thermalexpansioncoefficient': {
                    'name': r'$\alpha_{p}$',
                    'unit': r'K$^{-1}$'
                },
                'isothermalcompressibility': {
                    'name': r'$\kappa_{T}$',
                    'unit': r'm$^{3}$ J$^{-1}$'
                },
                'epsilon': {
                    'name': r'$\epsilon$',
                    'unit': r'kJ mol$^{-1}$'
                },
                'sigma': {
                    'name': r'$\sigma$',
                    'unit': r'nm'
                },
                'charge': {
                    'name': r'$q$',
                    'unit': r'e'
                },
                'vsite': {
                    'name': r'scaling factor',
                    'unit': r'-'
                }
            }



# useful dictionary for lookup
lookup = {
           'water':                       'water',
           'hydrogen-fluoride':           'fluorane',
           'hydrogen-chloride':           'chlorane',
           'hydrogen-bromide':            'bromane',
           'hydrogen-iodide':             'iodane',
           'energy':                      'dhvap',
           'volume':                      'volume',
           'rXX':                         'RDF-rXX',
           'gXX':                         'RDF-gXX',
           'hbonddist':                   'HB-distance',
           'hbondangle':                  'HB-angle',
           'dielectric':                  'dielectric',
           'diffusion':                   'diffusion',
           'isobaricheatcapacity':        'isobaric-heat-capacity',
           'thermalexpansioncoefficient': 'thermal-expansion',
           'isothermalcompressibility':   'isothermal-compressibility'
         }



def parse_arguments():

    parser  = argparse.ArgumentParser()

    parser.add_argument('-a',   '--action',           help='"submit" jobs, "analyze" results, or "statistics" from extra simulations',   required=True,  type=str, choices=['submit', 'analyze', 'statistics'])
    parser.add_argument('-o',   '--output',           help='path to output directory in ACTdata/HydrogenHalides/Bayesian/output',        required=True,  type=str)
    parser.add_argument('-p',   '--platform',         help='current platform (e.g. cluster name) used to run jobs',                      required=True,  type=str, choices=['local', 'csb', 'dardel', 'davinci'])
    parser.add_argument('-mol', '--molecule',         help='molecule in MolProps',                                                       required=False, type=str, choices=['water', 'hydrogen-fluoride', 'hydrogen-chloride', 'hydrogen-bromide', 'hydrogen-iodide'])
    parser.add_argument('-pt',  '--particletype',     help='particletype in .xml',                                                       required=False, type=str)
    parser.add_argument('-ff',  '--forcefield',       help='path to reference ACT .xml file',                                            required=False, type=str)
    parser.add_argument('-add', '--additional',       help='path to additional ACT .xml file(s) of similar models',                      required=False, type=str, nargs='+')
    parser.add_argument('-eq',  '--equilibration',    help='path to equilibration .dat file',                                            required=False, type=str)
    parser.add_argument('-sim', '--simulation',       help='path to simulation .dat file',                                               required=False, type=str)
    parser.add_argument('-pst', '--post',             help='path to post-production .dat file (only needed for "analyze" without grid)', required=False, type=str,   default=None)
    parser.add_argument('-pdb', '--system',           help='path to system .pdb file',                                                   required=False, type=str)
    parser.add_argument('-d',   '--distributions',    help='path to prior distributions .csv file',                                      required=False, type=str)
    parser.add_argument('-c',   '--coverage',         help='sampling coverage along each dimension',                                     required=False, type=float, default=0.7)
    parser.add_argument('-m',   '--samples',          help='number of samples',                                                          required=False, type=int,   default=10)
    parser.add_argument('-n',   '--points',           help='number of points per parameter',                                             required=False, type=int,   default=[10], nargs='+')
    parser.add_argument('-z',   '--initial',          help='initialize all n runs using parameters of initial force field',              action='store_true')
    parser.add_argument('-g',   '--grid',             help='perform a grid search for synthetic likelihood',                             action='store_true')
    parser.add_argument('-e',   '--epsilon',          help='include Van der Waals epsilon as parameter',                                 action='store_true')
    parser.add_argument('-s',   '--sigma',            help='include Van der Waals sigma as parameter',                                   action='store_true')
    parser.add_argument('-q',   '--charge',           help='include particle charges as parameter',                                      action='store_true')
    parser.add_argument('-vs',  '--vsite',            help='include vsite as parameter',                                                 action='store_true')
    parser.add_argument('-dhv', '--energy',           help='include enthalpy of vaporization as observable',                             action='store_true')
    parser.add_argument('-vol', '--volume',           help='include molecular volume as observable',                                     action='store_true')
    parser.add_argument('-rXX', '--rXX',              help='include peak distance "r" to the first RDF peak as observable',              action='store_true')
    parser.add_argument('-gXX', '--gXX',              help='include peak height "g(r)" of the first RDF peak as observable',             action='store_true')
    parser.add_argument('-hbd', '--hbonddist',        help='include hydrogen bond distance as observable',                               action='store_true')
    parser.add_argument('-hba', '--hbondangle',       help='include hydrogen bond angle as observable',                                  action='store_true')
    parser.add_argument('-v',   '--verbose',	      help='print more intermediates',                                                   action='store_true')
    args = parser.parse_args()

    return args



def invert_and_combine(properties):

    units_dict = {}
    for property in properties:
        units = converter[property]['unit']
        if units != '-':
            unit_parts = units.split()
            for unit_part in unit_parts:
                parts = unit_part.split('$')
                if len(parts) == 1:
                    unit      = parts[0]
                    dimension = -1
                elif len(parts) == 3:
                    unit      = parts[0]
                    dimension = -int(parts[1][2:-1])
                if unit in units_dict:
                    units_dict[unit] += dimension
                else:
                    units_dict[unit] = dimension
    if len(units_dict) == 0:
        return '-'
    else:
        return ' '.join([f'{key}' if value == 1 else f'{key}$^{{{value}}}$' for key, value in dict(sorted(units_dict.items(), key=lambda item: (-item[1], item[0]))).items()])



def xml2dict(path):

    with open(path, 'r') as xml:
        data = xml.read()
        dictionary = xmltodict.parse(data)
        return dictionary



def get_num_jobs(platform):

    if platform == 'davinci':
        command = "squeue -u nordman"
    else:
        return 0
    with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8", shell=True) as process:
        out, _ = process.communicate()
    return len(out.splitlines())



def burn_in_geweke(chain, window_fraction_sliding=0.1, window_fraction_reference=0.5):

    m, n               = chain.shape
    m_window_sliding   = int(window_fraction_sliding   * m)
    m_window_reference = int(window_fraction_reference * m)

    # reference window is constant
    window_reference         = chain[-m_window_reference:]
    window_reference_mean    = np.mean(window_reference, axis=0)
    window_reference_var     = np.var(window_reference,  axis=0, ddof=1)
    window_reference_std_err = np.sqrt(window_reference_var / m_window_reference)

    for m_start in range(m-m_window_sliding+1):

        # sliding window will change
        window_sliding         = chain[m_start:m_start+m_window_sliding]
        window_sliding_mean    = np.mean(window_sliding, axis=0)
        window_sliding_var     = np.var(window_sliding,  axis=0, ddof=1)
        window_sliding_std_err = np.sqrt(window_sliding_var / m_window_sliding)

        z_scores = (window_sliding_mean - window_reference_mean) / np.sqrt(window_sliding_std_err**2 + window_reference_std_err**2)

        # convergence criterium
        if np.all(np.abs(z_scores) <= 2.0):
            return m_start

    # not converged
    return m



def get_label_spacing(n, preferred_number=10):

    if n <= preferred_number:
        return 1
    else:
        factorization = factorint(n-1)
        factors = []
        for prime, exponent in factorization.items():
            factors.extend([prime] * exponent)
        all_combos = set()
        for r in range(1, len(factors) + 1):
            for combo in chain.from_iterable(combinations(factors, r) for r in range(1, len(factors) + 1)):
                product = 1
                for num in combo:
                    product *= num
                all_combos.add(product)
        combos = np.array(sorted(all_combos))
        mindex = np.argmin(np.abs(combos-preferred_number))
        return (n-1) // combos[mindex]



def epsilon_sigma_curve(eigenvectors, mu, i, j):

    slope    = eigenvectors[i,-1] / eigenvectors[j,-1]
    epsilon  = mu[i]
    sigma    = mu[j]
    r        = sigma*(((12*epsilon+sigma*slope)/(6*epsilon+sigma*slope))**(1/6))
    c        = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    sigmas   = np.linspace(mu[j]*(1-0.1), mu[j]*(1+0.1), 1000)
    epsilons = (c*r**12)/(4*((sigmas**12)-(r**6)*(sigmas**6)))
    keeps    = (epsilons > 0.5*epsilon) & (epsilons < 1.5*epsilon)
    sigmas   = sigmas[keeps]
    epsilons = epsilons[keeps]

    return epsilons, sigmas, r, c



def calculate_mode(samples, grid_size=20, zoom_factor=5, zoom_grid_size=10):

    # Initial coarse grid
    kde = gaussian_kde(samples.T, bw_method='scott')
    grid_ranges = [np.linspace(samples[:, i].min(), samples[:, i].max(), grid_size) for i in range(samples.shape[1])]
    grids = np.meshgrid(*grid_ranges, indexing='ij')
    grid_coords = np.vstack([grid.ravel() for grid in grids])
    kde_values = kde(grid_coords).reshape(grids[0].shape)

    # Find coarse mode
    max_idx = np.argmax(kde_values)
    coarse_mode = [grid.ravel()[max_idx] for grid in grids]

    # Zoomed-in grid around coarse mode
    zoomed_ranges = [
        np.linspace(center - (r[1] - r[0])/zoom_factor, center + (r[1] - r[0])/zoom_factor, zoom_grid_size)
        for center, r in zip(coarse_mode, grid_ranges)
    ]
    zoom_grids = np.meshgrid(*zoomed_ranges, indexing='ij')
    zoom_coords = np.vstack([g.ravel() for g in zoom_grids])
    zoom_kde_values = kde(zoom_coords).reshape(zoom_grids[0].shape)

    # Final mode
    max_idx_zoom = np.argmax(zoom_kde_values)
    final_mode = [grid.ravel()[max_idx_zoom] for grid in zoom_grids]

    return np.array(final_mode).reshape((1, samples.shape[1]))



def calculate_mode_meanshift(samples, quantile=0.2, n_samples=5000):
    # Estimate bandwidth (controls KDE kernel width)
    bandwidth = estimate_bandwidth(samples, quantile=quantile, n_samples=n_samples)

    # Apply MeanShift
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(samples)

    # Get the cluster centers (modes)
    cluster_centers = ms.cluster_centers_

    # Count members in each cluster and return the most populated one (i.e., global mode)
    labels, counts = np.unique(ms.labels_, return_counts=True)
    mode_idx = np.argmax(counts)
    mode = cluster_centers[mode_idx]

    return mode.reshape((1, samples.shape[1]))



def from_pdf(x, y):

    return np.sum(x*y)/np.sum(y), None



def from_rdf(x, y, dr, molecular_volume, coordination_number=4.96):

    rho = 1/molecular_volume # molecular volume in nm^3
    y_counts = y*(x**2)
    full = coordination_number/(4*np.pi*rho*dr)
    cumsum = np.cumsum(y_counts)
    i_above_full = min(np.argmax(cumsum > full), len(y)-1)
    x_cut = copy.deepcopy(x[:i_above_full+1])
    y_cut = copy.deepcopy(y_counts[:i_above_full+1])
    y_cut[-1] = (full-cumsum[i_above_full-1])
    y_sum = np.sum(y_cut)
    if y_sum > 0:
        p = y_cut/y_sum
        mean = np.sum(x_cut*p)/np.sum(p)
        var = np.sum(((x_cut-mean)**2)*p)/np.sum(p)
        if var == 0:
            var = dr**2
        x_peak = 0.5*(mean + np.sqrt(mean**2 + 8*var))
        y_peak = (y_sum*dr/np.sqrt(2*np.pi*var))*np.exp(-0.5*((x_peak-mean)**2)/var)/(x_peak**2)
    else:
        return x[-1], 0

    return x_peak, y_peak


def xvg2max(filename, distribution=None, delta_r=0.0002, molecular_volume=0.030):

    # Read XVG data
    data = np.loadtxt(filename, comments=["#", "@"])
    x    = data[:, 0]
    y    = data[:, 1]

    if np.all(np.isnan(y)):
        return x[-1], 0.0

    if distribution == 'RDF':
        x_max, y_max = from_rdf(x, y, delta_r, molecular_volume)
    elif distribution == 'PDF':
        x_max, y_max = from_pdf(x, y)
    else:
        sys.exit('Unknown distribution "{distribution}" for XVG analysis. Exiting...')

    return x_max, y_max


def create_numbered_marker(number, circle_size=20, font_size=24):

    circle = Path.unit_circle()
    transform = Affine2D().scale(circle_size)
    circle = transform.transform_path(circle)
    text_path = TextPath((0, 0), str(number), size=font_size)
    bbox = text_path.get_extents()
    text_width  = bbox.width
    text_height = bbox.height
    x_offset = (text_width) *(text_width -2*circle_size)/(2*circle_size)
    y_offset = (text_height)*(text_height-2*circle_size)/(2*circle_size)
    text_path = text_path.transformed(Affine2D().translate(x_offset, y_offset))
    vertices = circle.vertices.tolist() + text_path.vertices.tolist()
    codes = circle.codes.tolist() + text_path.codes.tolist()
    return Path(vertices, codes)



class Overhead:


    def __init__(self, args):

        # find where file was submitted from
        submit_directory = os.environ.get('SLURM_SUBMIT_DIR')
        if submit_directory is None:
            submit_directory = os.path.dirname(os.path.abspath(__file__))

        # define path to output directory
        output = os.path.join(submit_directory, 'output')

        # create output directory if it does not exist
        os.makedirs(output, exist_ok=True)

        # output directory name
        args['name'] = args['output']

        # store full output path
        args['output'] = os.path.join(output, args['output'])

        # path to this files directory
        args['base'] = submit_directory

        if args['action'] == 'submit':

            # create output directory or warn user if it already exists
            if os.path.isdir(args['output']):
                # TODO! Prevent resetting the grid
                if args['grid']:
                    sys.exit('It is not possible to resume from a grid search. Exiting...')
                else:
                    args['resume'] = True
            else:
                os.makedirs(args['output'])
                if args['grid']:
                    args['samples'] = 0
                    args['verbose'] = True
                args['resume'] = False

            # can not initialize a grid consisting of a single point!
            if args['initial'] and args['grid']:
                sys.exit(f'Can not initialize a grid consisting of a single point (option "grid" can not be coupled with "initial"). Exiting...')

            # summarize parameters
            parameters = []
            for parameter in ['epsilon', 'sigma', 'charge', 'vsite']:
                if args[parameter]:
                    parameters += [parameter]
                del args[parameter]
            if len(parameters) > 0:
                args['parameters'] = parameters
            else:
                sys.exit('Please provide at least one parameter to study using the flags on the command line. Exiting...')

            # summarize observables
            observables = []
            for observable in ['energy', 'volume', 'rXX', 'gXX', 'hbonddist', 'hbondangle']:
                if args[observable]:
                    if (observable in ['hbonddist', 'hbondangle']) and (args['molecule'] not in ['water']):
                        sys.exit(f'Can not compute "{observable}" for "{args['molecule']}"')
                    else:
                        observables += [observable]
                        del args[observable]
            if len(observables) > 0:
                args['observables'] = observables
            else:
                sys.exit('Please provide at least one observable to study using the flags on the command line. Exiting...')

            # starting points; still okay to initialize from a grid if MH and not just SL
            if args['initial'] and len(args['points']) != 1:
                sys.exit(f'Can not initialize from initial parameters without exactly one set of points. Exiting...')
            elif len(args['points']) == 1:
                args['points'] = args['points']*len(args['parameters'])
            elif len(args['points']) != len(args['parameters']):
                sys.exit(f'Please provide an equal amount of grid point numbers (now: {len(args["points"])}) and parameters (now: {len(args["parameters"])}) or a single point number to use for every parameter. Exiting...')

        else:

            output = args['output']
            if os.path.isdir(output):
                # load args from earlier
                pkl = os.path.join(output, 'args.pkl')
                with open(pkl, 'rb') as p:
                   arguments = pickle.load(p)
                for arg in args:
                    # always keep these!
                    if (not isinstance(args[arg], bool)) and (args[arg] is not None) and (arg not in ['points']):
                        arguments[arg] = args[arg]
                    # re-new these!
                    elif arg in ['additional']:
                        arguments[arg] = args[arg]
                    elif arg == 'post':
                        arguments['post'] = args['post']
                args = arguments
                #TODO! Remove in future; old runs need 'initial' and 'coverage'
                if 'initial' not in args:
                    args['initial'] = False
                if 'coverage' not in args:
                    args['coverage'] = 0.70
                if ('rXX' not in args['observables']) and ('rdfpeakdist' in args['observables']):
                    args['observables'][args['observables'].index('rdfpeakdist')] = 'rXX'
                if ('gXX' not in args['observables']) and ('rdfpeakintensity' in args['observables']):
                    args['observables'][args['observables'].index('rdfpeakintensity')] = 'gXX'
            else:
                sys.exit(f"The output directory {args['output']} does not exist. Exiting...")

        # check for necessary files are provided and exist
        for name in ['forcefield', 'equilibration', 'simulation', 'system', 'distributions']:
            path = args[name]
            # check if path is None and whether file exists or not
            if path is None:
                sys.exit(f'Please provide a path to the {name} file using the flags on the command line. Exiting...')
            elif os.path.isfile(path):
                print(f'{name.capitalize()} file found at {path}')
                args[name] = os.path.abspath(path)
            else:
                sys.exit(f'No {name} file found at {path}. Exiting...')
        # also check whether the path to the post-production file is correct
        if (args['action'] == 'analyze') and (args['post'] is not None):
            path = args['post']
            if os.path.isfile(path):
                print(f'Post-production file found at {path}')
                args['post'] = os.path.abspath(path)
            else:
                sys.exit(f'No post-production file found at {path}. Exiting...')

        # to handle older runs
        if 'additional' not in args:
            args['additional'] = None

        # check for necessary files are provided and exist
        additional = []
        if args['additional'] is not None:
            for additional_forcefield in args['additional']:
                # check if path is None and whether file exists or not
                if os.path.isfile(additional_forcefield):
                    additional += [os.path.abspath(additional_forcefield)]
                else:
                    sys.exit(f'No addtional forcefield file found at {addtional_forcefield}. Exiting...')
        args['additional'] = additional

        # check for particletype
        data = xml2dict(args['forcefield'])
        if not args['particletype'] in [particle['@identifier'] for particle in data['alexandria_chemistry']['particletypes']['particletype']]:
            sys.exit(f'Did not find particletype {args["particletype"]} in .xml file. Exiting...')

        # append experimental data file
        args['experiment'] = os.path.join(submit_directory, 'expdata.csv')

        # save args for later
        pkl = os.path.join(args['output'], 'args.pkl')
        with open(pkl, 'wb') as p:
            pickle.dump(args, p, pickle.HIGHEST_PROTOCOL)

        self.settings = args

        # modify labels if necessary
        element = ''.join([i for i in self.settings['particletype'] if not i.isdigit()]).upper()
        converter['rXX']['name'] = fr'$\tilde{{r}}_{{{element}{element},1}}$'
        converter['gXX']['name'] = fr'$\tilde{{g}}_{{{element}{element},1}}$'


    def get_temperature(self):

        temperature = 0.0
        with open(self.settings['simulation'], 'r') as dat:
            for line in dat:
                line_no_comment = line.split('#')[0]
                try:
                    key, value  = [v.strip() for v in line_no_comment.split('=')]
                    if key == 'temperature_c':
                        temperature_raw = value
                        temperature_parts = temperature_raw.split('*')
                        if len(temperature_parts) > 1:
                            if temperature_parts[1] == 'kelvin':
                                temperature = float(temperature_parts[0])
                            else:
                                sys.exit('Conversion from {temperature_parts[1]} to Kelvin not implemented. Exiting...')
                        else:
                            temperature = float(temperature_parts[0])
                except:
                    continue

        return temperature


    def get_experiment(self):

        temperature = self.get_temperature()

        columns = [lookup[observable] for observable in self.settings['observables']]

        with open(self.settings['experiment'], 'r') as ef:
            lines   = ef.readlines()
            indices = {column: c for c, column in enumerate(lines[0].split(',')) if column in columns}
            for line in lines[1:]:
                data = line.split(',')
                if data[0] == lookup[self.settings['molecule']] and float(data[1]) == temperature:
                    return np.array([float(data[indices[column]]) for column in columns])
            sys.exit(f'Could not find experimental data for {lookup[self.settings["molecule"]]} ')


    # get prior parameters
    def read_distributions(self):

        beta_values = {}

        # iterate through rows in file
        with open(self.settings['distributions'], 'r') as df:

             # get beta values
             for line in df.readlines():

                 contents = line.split(',')
                 beta_values[contents[0]] = {'a': float(eval(contents[1])), 'b': float(eval(contents[2])), 'min': float(eval(contents[3])), 'max': float(eval(contents[4]))}

        # check whether all beta values are assigned otherwise a = b = 1
        for parameter in self.settings['parameters']:

            if parameter not in beta_values:

                beta_values[parameter] = {'a': 1.0, 'b': 1.0, 'min': 0.5, 'max': 1.5}

        self.beta_values = beta_values


    def read_parameters(self):

        models = {}

        for ff, forcefield in enumerate([self.settings['forcefield']]+self.settings['additional']):

            name = forcefield.split('/')[-1].split('.')[0]

            # read from .xml
            data = xml2dict(forcefield)

            # to store parameter values
            parameters = np.zeros((len(self.settings['parameters']),))

            if ff == 0:
                # store all particles
                self.particletypes = {}
                for particle in data['alexandria_chemistry']['particletypes']['particletype']:
                    self.particletypes[particle['@identifier']] = {}
                    for option in particle['option']:
                        self.particletypes[particle['@identifier']][option['@key']] = option['@value']
                if self.settings['particletype'] not in self.particletypes:
                    sys.exit(f'Did not find particletype {self.settings["particletype"]} in .xml file. Could only find {", ".join([particle for particle in self.particletypes.keys()])}. Exiting...')
                linestyle = '-'
                if len(self.settings['additional']) > 0:
                    marker = create_numbered_marker(ff+1)
                else:
                    marker = 'o'
            else:
                linestyle = (0, (ff, ff))
                marker    = create_numbered_marker(ff+1)

            # start searching!
            for p, parameter in enumerate(self.settings['parameters']):
                if parameter in ['epsilon', 'sigma']:
                    # epsilon and sigma can be found under the Van der Waals section
                    for interaction in data['alexandria_chemistry']['interaction']:
                        if interaction['@type'] == 'VANDERWAALS':
                            for parameterlist in interaction['parameterlist']:
                                if parameterlist['@identifier'] == self.particletypes[self.settings['particletype']]['vdwtype']:
                                    for vdwparameter in parameterlist['parameter']:
                                        if parameter == vdwparameter['@type']:
                                            value = float(vdwparameter['@value'])
                elif parameter == 'charge':
                    # charge is found under particle types, but we need to find reference "q" first
                    for particle in data['alexandria_chemistry']['particletypes']['particletype']:
                        if self.settings['particletype'] == particle['@identifier']:
                            for particle_parameter in particle['parameter']:
                                if particle_parameter['@type'] == 'charge':
                                    value = float(particle_parameter['@value'])
                    # now we need to know how the other charges scale w.r.t. to "q"
                    for particle in data['alexandria_chemistry']['particletypes']['particletype']:
                        for particle_parameter in particle['parameter']:
                            if particle_parameter['@type'] == 'charge':
                                self.particletypes[particle['@identifier']]['q_scale'] = float(particle_parameter['@value'])/value
                elif parameter == 'vsite':
                    sys.exit('Accessing vsite is yet to be implemented. Exiting...')
                else:
                    sys.exit(f'Do not know how to access parameter {parameter}. Exiting...')
                parameters[p] = value

            if ff == 0:
                self.parameters = parameters
                self.model      = name

            models[name] = {'parameters': parameters, 'linestyle': linestyle, 'marker': marker}

        self.models = models


    def set_parameters(self, parameters, molecule, base):

        # update ACT forcefield file
        for p, parameter in enumerate(self.settings['parameters']):
            if parameter in ['epsilon', 'sigma']:
                command = [
                    'alexandria', 'edit_ff', '-force',
                    '-ff',  'act.xml',
                    '-o',   'act.xml',
                    '-p',   parameter,
                    '-val', f'{parameters[p]}',
                    '-a',   self.settings["particletype"]
                ]
                with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                    stdout, stderr = process.communicate()
                    if process.returncode != 0:
                        print(f'ALEXANDRIA ERROR: {stderr}\n')
            elif parameter == 'charge':
                for particletype in self.particletypes:
                    command = [
                        'alexandria', 'edit_ff', '-force',
                        '-ff',  'act.xml',
                        '-o',   'act.xml',
                        '-p',   parameter,
                        '-val', f'{self.particletypes[particletype]["q_scale"]*parameters[p]}',
                        '-a',   particletype
                    ]
                    with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                        stdout, stderr = process.communicate()
                        if process.returncode != 0:
                            print(f'ALEXANDRIA ERROR: {stderr}\n')
            elif parameter == 'vsite':
                sys.exit('Modifying vsite is yet to be implemented. Exiting...')
            else:
                sys.exit(f'Do not know how to modify {parameter}. Exiting...')

        # (forcefully) remove old sim.dat to prevent too many files
        command = ['rm', '-f', 'sim.dat']
        with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
            _, _ = process.communicate()

        # generate OpenMM forcefield file
        command = [
            'alexandria', 'gentop',
            '-ff',	'act.xml',
            '-db',	molecule,
            '-charges', f'{base}/MolProps/mp2-aug-cc-pvtz.xml',
            '-openmm', 'openmm.xml'
        ]
        with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                print(f'ALEXANDRIA ERROR: {stderr}\n')


    def compute_sigma_points(self, mu, Sigma, alpha=0.9, kappa=0, beta=2):

        n = len(mu)  # Dimensionality
        lambda_ = alpha**2 * (n + kappa) - n  # Scaling factor

        # Compute weights
        Wm = np.zeros(2 * n + 1)
        Wc = np.zeros(2 * n + 1)

        Wm[0] = lambda_ / (n + lambda_)
        Wc[0] = lambda_ / (n + lambda_) + (1 - alpha**2 + beta)

        for i in range(1, 2 * n + 1):
            Wm[i] = 1 / (2 * (n + lambda_))
            Wc[i] = 1 / (2 * (n + lambda_))

        # Compute square root of (n + lambda) * Sigma using Cholesky decomposition
        try:
            sqrt_matrix = np.linalg.cholesky((n + lambda_) * Sigma)
        except np.linalg.LinAlgError:
            sqrt_matrix = np.linalg.svd((n + lambda_) * Sigma)[0]  # Use SVD fallback if Cholesky fails

        # Compute sigma points
        sigma_points = np.zeros((2 * n + 1, n))
        sigma_points[0] = mu  # First sigma point is the mean

        for i in range(n):
            sigma_points[i + 1]     = mu + sqrt_matrix[:, i]
            sigma_points[n + i + 1] = mu - sqrt_matrix[:, i]

        return sigma_points, Wm, Wc


    def submit(self):

        # enter directory
        os.chdir(self.settings['output'])

        # get arguments to pass on to SLURM script
        action        = self.settings['action']
        molecule      = self.settings['molecule']
        particletype  = self.settings['particletype']
        base          = self.settings["base"]
        experiment    = self.settings['experiment']
        forcefield    = self.settings['forcefield']
        equilibration = self.settings['equilibration']
        simulation    = self.settings['simulation']
        system        = self.settings['system']
        distributions = self.settings['distributions']
        coverage      = self.settings['coverage']
        samples       = self.settings['samples']
        parameters    = self.settings['parameters']
        observables   = self.settings['observables']

        # actual element
        element = ''.join([i for i in self.settings['particletype'] if not i.isdigit()]).upper()

        # get output name to give a proper name to SLURM job
        name = self.settings['name']

        if not self.settings['resume']:
            os.makedirs('SIM')
        os.chdir('SIM')

        # get number of points per parameter
        points = self.settings['points']
        # positions to start in parameter space
        if self.settings['initial']:
            coordinate_indices_sets = [(point,) for point in range(1, points[0]+1)]
        else:
            grid_ranges        = [range(1, point+1) for point in points]
            coordinate_indices_sets = product(*grid_ranges)

        for ci, coordinate_indices in enumerate(coordinate_indices_sets):

            if self.settings['initial']:
                directory = f'{str(ci+1).zfill(len(str(len(coordinate_indices_sets))))}'
            else:
                directory_parts = []
                for p, parameter in enumerate(parameters):
                    directory_parts += [f'{parameter[:3]}{str(coordinate_indices[p]).zfill(len(str(points[p])))}']
                directory = '-'.join(directory_parts)

            if not self.settings['resume']:
                os.makedirs(directory)
                os.chdir(directory)
                # output directories
                if molecule in ['water']:
                    os.system(f'printf "2\n1\n" | gmx pdb2gmx -f {system} -o structure.gro -p topology.top -ignh > /dev/null 2>&1')
                    os.system(f'echo -e "integrator = md\ndt = 0.001\nnsteps = 100000" > dummy.mdp')
                    os.system(f'gmx grompp -f dummy.mdp -c structure.gro -p topology.top -o topology.tpr > /dev/null 2>&1')
                    os.system(f'printf "a OW\nname 3 {element}\na HW1 | a HW2\nname 4 H\nq\n" | gmx make_ndx -f structure.gro -o index.ndx > /dev/null 2>&1')
                else:
                    os.system(f'printf "a {element}\nname 3 {element}\nq\n" | gmx make_ndx -f {self.settings["system"]} -o index.ndx > /dev/null 2>&1')
            else:
                os.chdir(directory)

            pf = self.settings['platform']
            if samples == 0:
                duartion = '24:00:00'
            else:
                duration = '24:00:00'
            if pf == 'csb':
                command = f'sbatch -t {duration} --gres=gpu:1 -p CLUSTER-AMD --nodelist=compute-0-[32-36] --nodes=1 -J "{name}-{directory}" '
                related_jobs = [job.split() for job in os.popen(f"squeue --format='%.10i %100j' --me | grep '{name}-{directory}'").read().split('\n')[:-1]]
                if len(related_jobs) > 0:
                    last_job = sorted(related_jobs, key=lambda x: int(x[0]))[-1][0]
                    command += f'--dependency=afterany:{last_job} '
            elif pf == 'dardel':
                command = f'sbatch -A naiss2023-5-531 --nodes=1 -t {duration} -p gpu -J "{name}-{directory}" '
                related_jobs = [job.split() for job in os.popen(f"squeue --format='%.10i %100j' --me | grep '{name}-{directory}'").read().split('\n')[:-1]]
                if len(related_jobs) > 0:
                    last_job = sorted(related_jobs, key=lambda x: int(x[0]))[-1][0]
                    command += f'--dependency=afterany:{last_job} '
            elif pf == 'davinci':
                command = f'sbatch -n 1 -t {duration} --gres=gpu:1 -p regular -C RTX2080Ti -J "{name}-{directory}" --exclude=a002 '
                related_jobs = [job.split() for job in os.popen(f"squeue --format='%.10i %100j' --me | grep '{name}-{directory}'").read().split('\n')[:-1]]
                if len(related_jobs) > 0:
                    last_job = sorted(related_jobs, key=lambda x: int(x[0]))[-1][0]
                    command += f'--dependency=afterany:{last_job} '
            else:
                command = ''
            command += (f'{base}/metropolis_hastings.py '
                      + f'-mol {molecule} '
                      + f'-pt {particletype} '
                      + f'-b {base} '
                      + f'-exp {experiment} '
                      + f'-ff {forcefield} '
                      + f'-eq {equilibration} '
                      + f'-sim {simulation} '
                      + f'-pdb {system} '
                      + f'-d {distributions} '
                      + f'-c {coverage} '
                      + f'-m {samples} '
                      + f'-n {" ".join([str(point) for point in points])} '
                      + f'-i {" ".join([str(coordinate_index) for coordinate_index in coordinate_indices])} '
                      + f'-p {" ".join(parameters)} '
                      + f'-o {" ".join(observables)}')
            for flag in ['initial', 'verbose', 'resume']:
                if self.settings[flag]:
                      command += f' --{flag}'
            num_lines = get_num_jobs(pf)
            while num_lines > 90:
                sleep(60.0)
                num_lines = get_num_jobs(pf)
            os.system(command)
            os.chdir('..')


    def analyze_grid(self):

        # enter directory
        os.chdir(self.settings['output'])

        # observables and parameters
        observables = self.settings['observables']
        parameters  = self.settings['parameters']

        # get number of observables and parameters
        m = len(observables)
        n = len(parameters)

        # to store observables, parameters and likelihoods
        observable_means  = np.empty((0,m))
        observable_stds   = np.empty((0,m))
        parameter_samples = np.empty((0,n))

        # get number of points per parameter
        points = self.settings['points']

        # read experiment
        experiment = self.get_experiment()

        # get distributions and points
        self.read_distributions()
        self.read_parameters()

        # calculate min and max values
        self.min = self.parameters * np.array([self.beta_values[parameter]['min'] if self.parameters[p] > 0 else self.beta_values[parameter]['max'] for p, parameter in enumerate(parameters)])
        self.max = self.parameters * np.array([self.beta_values[parameter]['max'] if self.parameters[p] > 0 else self.beta_values[parameter]['min'] for p, parameter in enumerate(parameters)])

        # find point closest to original model
        labels         = []
        original_point = []
        for i in range(n):
            labels         += [[(self.min[i]+beta.ppf((1+k)/(1+points[i]),self.beta_values[parameters[i]]['a'],self.beta_values[parameters[i]]['b'])*(self.max[i]-self.min[i])) for k in range(points[i])]]
            original_point += [1 + np.argmin(np.abs(labels[i] - self.parameters[i]))]
        original_point = tuple(original_point)

        # miscellaneous plots
        os.makedirs('PDF/MISC', exist_ok=True)
        os.makedirs('PNG/MISC', exist_ok=True)
        # high dimension plots
        indices             = range(m)
        combinations        = {}
        likelihoods         = {}
        shapiro_wilks       = {}
        for i in range(1, m+1):
            os.makedirs(f'PDF/{i}D', exist_ok=True)
            os.makedirs(f'PNG/{i}D', exist_ok=True)
            combinations[i] = {}
            for combination in itertools.combinations(indices, i):
                combined = "-".join([observables[c] for c in combination])
                combinations[i][combined]     = list(combination)
                likelihoods[combined]         = np.empty((0,))
                shapiro_wilks[combined]       = np.empty((0,), dtype=bool)
                if i == m:
                    full_combined = combined

        os.chdir('SIM')

        # positions to start in parameter space
        grid_ranges = [range(1, point+1) for point in points]
        grid_points = product(*grid_ranges)

        for gi, grid_point in enumerate(grid_points):

            if grid_point == original_point:
                ogi = gi

            directory_parts = []
            for p, parameter in enumerate(parameters):
                directory_parts += [f'{parameter[:3]}{str(grid_point[p]).zfill(len(str(points[p])))}']
            directory = '-'.join(directory_parts)

            os.chdir(directory)
            # output
            successful = True
            try:
                mu           = np.load('mu.npy')
                Sigma        = np.load('Sigma.npy')
                measurements = np.loadtxt('observables.xvg')[:,1:].T # slicing removes the time values
            except:
                successful = False
                print(f'WARNING: could not find data for run {directory}')
            for i in combinations:
                for combined in combinations[i]:
                    if successful:
                        combination = combinations[i][combined]
                        try:
                            likelihood = multivariate_normal.pdf(experiment[combination], mean=mu[combination], cov=Sigma[combination][:,combination], allow_singular=True)
                        except:
                            likelihood = 0.0
                        shapiro_wilk	   = np.all([shapiro(measurements[index])[1]            > 0.05 for index in combination])
                    else:
                        likelihood         = 0.0
                        shapiro_wilk       = False
                    likelihoods[combined]         = np.append(likelihoods[combined],         likelihood)
                    shapiro_wilks[combined]       = np.append(shapiro_wilks[combined],       shapiro_wilk)

            if successful:
                observable_mean = np.zeros((1,m))
                observable_std  = np.zeros((1,m))
                for o, observable in enumerate(observables):
                    observable_mean[0,o] = mu[o]
                    observable_std[0,o]  = np.sqrt(Sigma[o,o])
            else:
                observable_mean = np.array([[None for _ in range(m)]], dtype=float)
                observable_std  = np.array([[None for _ in range(m)]], dtype=float)
            parameter_sample = np.array([(self.min[g]+beta.ppf(gp/(1+points[g]),self.beta_values[parameters[g]]['a'],self.beta_values[parameters[g]]['b'])*(self.max[g]-self.min[g])) for g, gp in enumerate(grid_point)]).reshape((1,n))
            # append
            observable_means  = np.append(observable_means, observable_mean, axis=0)
            observable_stds   = np.append(observable_stds,  observable_std,  axis=0)
            parameter_samples = np.append(parameter_samples, parameter_sample, axis=0)
            os.chdir('..')
        os.chdir('..')

        # add points of interest
        pois = {}
        # add original model too
        original_forcefield = self.settings['forcefield'].split('/')[-1].split('.')[0]
        index               = np.ravel_multi_index([int(op-1) for op in original_point], points)
        label               = ", ".join([f'{converter[parameters[p]]["name"]}={parameter:.3f}' for p, parameter in enumerate(parameter_samples[index])])
        pois['original']    = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'{original_forcefield} ({label})', 'color': 'r', 'marker': 'o'}
        # and others
        for combined in likelihoods:
            index          = np.where(likelihoods[combined]==likelihoods[combined].max())[0][0]
            label          = ", ".join([f'{converter[parameters[p]]["name"]}={parameter:.3f}' for p, parameter in enumerate(parameter_samples[index])])
            pois[combined] = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'highest in this space ({label})', 'color': '#8B4513', 'marker': 'X'}
            # add global maxmimum
            if combined == full_combined:
                pois['global'] = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'highest globally ({label})', 'color': 'g', 'marker': 's'}
            # add the other statistics
            if combined == 'energy-volume':
                pois['therm'] = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'highest {converter["energy"]["name"]}+{converter["volume"]["name"]} ({label})', 'color': 'b', 'marker': '^'}
            elif combined == 'rXX-gXX':
                pois['rdf'] = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'highest {converter["rXX"]["name"]}+{converter["gXX"]["name"]} ({label})', 'color': '#FFA500', 'marker': 'v'}
            elif combined == 'hbonddist-hbondangle':
                pois['hbond'] = {'index': index, 'indices': np.unravel_index(index, points), 'label': f'highest {converter["hbonddist"]["name"]}+{converter["hbondangle"]["name"]} ({label})', 'color': 'purple', 'marker': 'D'}

        # for KDE plots
        levels        = 11
        levels_values = np.linspace(0, 1, levels)

        if np.all([key in pois for key in ['therm', 'rdf', 'hbond']]):

            # for generating plots, relative to plot width = height = 1
            L  = 0.45  # left margin
            WS = 0.05  # width spacing between plots
            R  = 0.05  # right margin, between rightmost plot and colorbar
            CB = 0.75  # horizontal space of colorbar
            B  = 0.30  # bottom margin
            HS = 0.05  # height spacing between plots
            T  = 0.15  # top margin

            # full size
            figwidth  = L + 2 + WS + R + CB
            figheight = B + 2 + HS + T

            for i in range(n):
                parameter_i = parameters[i]
                points_i    = parameter_samples[:,i]
                min_i       = self.min[i]+beta.ppf(        1/(1+points[i]), self.beta_values[parameters[i]]['a'], self.beta_values[parameters[i]]['b'])*(self.max[i]-self.min[i])
                max_i       = self.min[i]+beta.ppf(points[i]/(1+points[i]), self.beta_values[parameters[i]]['a'], self.beta_values[parameters[i]]['b'])*(self.max[i]-self.min[i])
                grid_i      = np.linspace(min_i, max_i, points[i])
                for j in range(i+1, n):
                    parameter_j = parameters[j]
                    points_j    = parameter_samples[:,j]
                    min_j       = self.min[j]+beta.ppf(        1/(1+points[j]), self.beta_values[parameters[j]]['a'], self.beta_values[parameters[j]]['b'])*(self.max[j]-self.min[j])
                    max_j       = self.min[j]+beta.ppf(points[j]/(1+points[j]), self.beta_values[parameters[j]]['a'], self.beta_values[parameters[j]]['b'])*(self.max[j]-self.min[j])
                    grid_j      = np.linspace(min_j, max_j, points[j])
                    filename = f"TWOBYTWO-{parameter_i}-{parameter_j}"
                    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(2*figwidth, 2*figheight), dpi=300)
                    for ax_ij in ax.flat:
                        ax_ij.set_xticks([])
                        ax_ij.set_yticks([])
                        ax_ij.set_xlim(min_i, max_i)
                        ax_ij.set_ylim(min_j, max_j)
                    fig.supxlabel(f'{converter[parameter_i]["name"]} ({converter[parameter_i]["unit"]})', fontsize=16)
                    ax[1,0].xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    ax[1,0].tick_params(axis='x', labelsize=12)
                    ax[1,0].xaxis.set_tick_params(which='major', direction='out')
                    ax[1,1].xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    ax[1,1].tick_params(axis='x', labelsize=12)
                    ax[1,1].xaxis.set_tick_params(which='major', direction='out')
                    fig.supylabel(f'{converter[parameter_j]["name"]} ({converter[parameter_j]["unit"]})', fontsize=16)
                    ax[0,0].yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    ax[0,0].tick_params(axis='y', labelsize=12)
                    ax[0,0].yaxis.set_tick_params(which='major', direction='out')
                    ax[1,0].yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    ax[1,0].tick_params(axis='y', labelsize=12)
                    ax[1,0].yaxis.set_tick_params(which='major', direction='out')
                    grid_ii, grid_jj = np.meshgrid(grid_i, grid_j)
                    for p, combined in enumerate(['energy-volume', 'rXX-gXX', 'hbonddist-hbondangle', full_combined]):
                        likelihood_axis_sum = np.sum(likelihoods[combined].reshape(tuple(points)), axis=tuple([k for k in range(n) if ((k != i) and (k != j))])).flatten()
                        slice_i = eval(f"points_i.reshape(tuple(points))[{','.join([':' if ((k == i) or (k == j)) else '0' for k in range(n)])}].flatten()")
                        slice_j = eval(f"points_j.reshape(tuple(points))[{','.join([':' if ((k == i) or (k == j)) else '0' for k in range(n)])}].flatten()")
                        grid_y = griddata((slice_i, slice_j), likelihood_axis_sum, (grid_ii, grid_jj), method='nearest')
                        ax[p//2,p%2].set_facecolor([0.2, 0.0, 0.4])
                        ax[p//2,p%2].contourf(grid_ii, grid_jj, grid_y, levels=levels, cmap='rainbow')
                        #ax[p//2,p%2].scatter(points_i[shapiro_wilks[combined]],   points_j[shapiro_wilks[combined]],   c='w',                       marker='+',                        s=3)
                        ax[p//2,p%2].scatter(points_i[pois['original']['index']], points_j[pois['original']['index']], c=pois['original']['color'], marker=pois['original']['marker'], s=30, edgecolors='w', linewidth=1)
                        ax[p//2,p%2].annotate(chr(65+p), xy=(0.95, 0.95), xycoords="axes fraction", fontsize=20, fontweight="bold", ha="right", va="top", color="k", path_effects=[pe.withStroke(linewidth=5, foreground="white")])
                    norm = mcolors.BoundaryNorm(boundaries=levels_values, ncolors=256)
                    sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
                    sm.set_array([])
                    fig.subplots_adjust(left=L/figwidth, wspace=WS/figwidth, right=(figwidth-(R+CB))/figwidth, bottom=B/figheight, hspace=HS/figheight, top=(figheight-T)/figheight)
                    cbar_ax = fig.add_axes([(figwidth-CB)/figwidth, B/figheight, CB/(3*figwidth), (figheight-(T+B))/figheight])
                    cbar = fig.colorbar(sm, cax=cbar_ax, ticks=levels_values)
                    cbar.ax.tick_params(labelsize=12)
                    cbar.set_label("posterior density", fontsize=16)
                    fig.align_labels()
                    fig.savefig(f'PDF/MISC/{filename}-supporting.pdf')
                    fig.savefig(f'PNG/MISC/{filename}-supporting.png')

        for i in range(n):

            parameter_i = parameters[i]
            points_i    = parameter_samples[:,i]
            min_i	= self.min[i]+beta.ppf(        1/(1+points[i]), self.beta_values[parameters[i]]['a'], self.beta_values[parameters[i]]['b'])*(self.max[i]-self.min[i])
            max_i	= self.min[i]+beta.ppf(points[i]/(1+points[i]), self.beta_values[parameters[i]]['a'], self.beta_values[parameters[i]]['b'])*(self.max[i]-self.min[i])
            grid_i	= np.linspace(min_i, max_i, 1000)

            for combined in likelihoods:

                unique_vals, indices = np.unique(points_i, return_inverse=True)
                summed_likelihoods   = np.bincount(indices, weights=likelihoods[combined])
                dx                   = np.min(unique_vals[1:]-unique_vals[:-1])

                filename = f'BAR-{parameter_i}-{combined}'
                fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                ax.bar(unique_vals, summed_likelihoods/np.sum(summed_likelihoods*dx), width=dx)
                ax.set_xlabel(f'{converter[parameter_i]["name"]} ({converter[parameter_i]["unit"]})', fontsize=11)
                ax.set_ylabel(f'Likelihood ({invert_and_combine([parameter_i])})', fontsize=11)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.tick_params(axis='x', labelsize=9)
                ax.tick_params(axis='y', labelsize=9)
                fig.tight_layout()
                fig.savefig(f'PDF/{len(combined.split("-"))}D/{filename}.pdf')
                fig.savefig(f'PNG/{len(combined.split("-"))}D/{filename}.png')

                filename = f'LOG-BAR-{parameter_i}-{combined}'
                fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                ax.bar(unique_vals, np.log(summed_likelihoods/np.sum(summed_likelihoods*dx)), width=dx)
                ax.set_xlabel(f'{converter[parameter_i]["name"]} ({converter[parameter_i]["unit"]})', fontsize=11)
                ax.set_ylabel(f'Log-likelihood ({invert_and_combine([parameter_i])})', fontsize=11)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.tick_params(axis='x', labelsize=9)
                ax.tick_params(axis='y', labelsize=9)
                fig.tight_layout()
                fig.savefig(f'PDF/{len(combined.split("-"))}D/{filename}.pdf')
                fig.savefig(f'PNG/{len(combined.split("-"))}D/{filename}.png')

            for j in range(i+1, n):

                parameter_j = parameters[j]
                points_j    = parameter_samples[:,j] # could retrieve more points than desired with three or more parameters!
                min_j       = self.min[j]+beta.ppf(        1/(1+points[j]), self.beta_values[parameters[j]]['a'], self.beta_values[parameters[j]]['b'])*(self.max[j]-self.min[j])
                max_j       = self.min[j]+beta.ppf(points[j]/(1+points[j]), self.beta_values[parameters[j]]['a'], self.beta_values[parameters[j]]['b'])*(self.max[j]-self.min[j])
                grid_j      = np.linspace(min_j, max_j, 1000)

                for combined in likelihoods:

                    combined_parameters     = [parameter_i, parameter_j]
                    combined_inverted_units = invert_and_combine(combined_parameters)

                    try:
                        # 2D KDE
                        filename = f'KDE-{parameter_i}-{parameter_j}-{combined}'
                        fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                        ax.set_facecolor([0.2, 0.0, 0.4])
                        kde_plot = sns.kdeplot(x=points_i, y=points_j, weights=likelihoods[combined], fill=True, cmap='rainbow', ax=ax, levels=levels)
                        norm = mcolors.BoundaryNorm(boundaries=levels_values, ncolors=256)
                        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
                        sm.set_array([])
                        cbar = fig.colorbar(sm, ax=ax, ticks=levels_values)
                        cbar.set_label(f"posterior density")
                        #ax.scatter(points_i[shapiro_wilks[combined]],       points_j[shapiro_wilks[combined]],       c='w', marker='+', s=int(15000/(points[i]*points[j])))
                        ax.scatter(points_i[pois['original']['index']], points_j[pois['original']['index']], c=pois['original']['color'], marker=pois['original']['marker'], s=100, edgecolors='w', linewidth=1, label=pois['original']['label'])
                        ax.set_xlabel(f'{converter[parameter_i]["name"]} ({converter[parameter_i]["unit"]})', fontsize=11)
                        ax.set_ylabel(f'{converter[parameter_j]["name"]} ({converter[parameter_j]["unit"]})', fontsize=11)
                        ax.set_xlim(min_i, max_i)
                        ax.set_ylim(min_j, max_j)
                        ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                        ax.tick_params(axis='x', labelsize=9)
                        ax.tick_params(axis='y', labelsize=9)
                        fig.tight_layout()
                        fig.savefig(f'PDF/{len(combined.split("-"))}D/{filename}.pdf')
                        fig.savefig(f'PNG/{len(combined.split("-"))}D/{filename}.png')
                    except Exception as error:
                        print(f"An exception occurred when generating the KDE plot: {type(error).__name__} - {error}")
                    try:
                        # grid
                        likelihood_axis_sum = np.sum(likelihoods[combined].reshape(tuple(points)), axis=tuple([k for k in range(n) if ((k != i) and (k != j))])).flatten()
                        slice_i = eval(f"points_i.reshape(tuple(points))[{','.join([':' if ((k == i) or (k == j)) else '0' for k in range(n)])}].flatten()")
                        slice_j = eval(f"points_j.reshape(tuple(points))[{','.join([':' if ((k == i) or (k == j)) else '0' for k in range(n)])}].flatten()")
                        grid_ii, grid_jj = np.meshgrid(grid_i, grid_j)
                        grid_y = griddata((slice_i, slice_j), likelihood_axis_sum, (grid_ii, grid_jj), method='cubic')
                        # 2D contour
                        filename = f'CONTOUR-{parameter_i}-{parameter_j}-{combined}'
                        fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                        ax.set_facecolor([0.2, 0.0, 0.4])
                        norm = mcolors.BoundaryNorm(boundaries=levels_values, ncolors=256)
                        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
                        sm.set_array([])
                        cbar = fig.colorbar(sm, ax=ax, ticks=levels_values)
                        cbar.set_label(f"posterior density")
                        ax.contourf(grid_ii, grid_jj, grid_y, levels=levels, cmap='rainbow')
                        #ax.scatter(points_i[shapiro_wilks[combined]],       points_j[shapiro_wilks[combined]],       c='w', marker='+', s=int(15000/(points[i]*points[j])))
                        ax.scatter(points_i[pois['original']['index']], points_j[pois['original']['index']], c=pois['original']['color'], marker=pois['original']['marker'], s=100, edgecolors='w', linewidth=1, label=pois['original']['label'])
                        ax.set_xlabel(f'{converter[parameter_i]["name"]} ({converter[parameter_i]["unit"]})', fontsize=11)
                        ax.set_ylabel(f'{converter[parameter_j]["name"]} ({converter[parameter_j]["unit"]})', fontsize=11)
                        ax.set_xlim(min_i, max_i)
                        ax.set_ylim(min_j, max_j)
                        ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                        fig.tight_layout()
                        fig.savefig(f'PDF/{len(combined.split("-"))}D/{filename}.pdf')
                        fig.savefig(f'PNG/{len(combined.split("-"))}D/{filename}.png')
                    except Exception as error:
                        print(f"An exception occurred when generating the contour plot: {type(error).__name__} - {error}")

        # get optimal point and later generate 2D plots that intersects it
        reshaped = likelihoods[full_combined].reshape(tuple(points))
        optimum  = tuple(indices[0] for indices in np.where(reshaped==reshaped.max()))
        indices  = np.arange(len(likelihoods[full_combined])).reshape(tuple(points))
        # now, let's take a look at these 2D slices!
        for i in range(n):
            xlabels = [f'{value:.3f}' for value in labels[i]]
            for j in range(i+1, n):
                subset = indices[:]
                slices = ','.join([':' if ((k == i) or (k == j)) else str(optimum[k]) for k in range(n)])
                subset = eval(f'subset[{slices}]').flatten()
                observable_means_subset = observable_means[subset]
                observable_stds_subset  = observable_stds[subset]
                ylabels    = [f'{value:.3f}' for value in labels[j][::-1]]
                for o, observable in enumerate(observables):
                    custom_name = copy.deepcopy(converter[observable]["name"])
                    custom_name = "".join([character for character in custom_name if character != '$'])
                    sims = np.flip(observable_means_subset[:,o].reshape((points[i],points[j])).T, axis=0)
                    filename = f'OBS-1D-{observable}-{parameters[i]}-{parameters[j]}'
                    pdiffs = 100.0*((sims-experiment[o])/experiment[o])
                    fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                    ax.set_facecolor([0.0, 0.0, 0.0])
                    im = ax.imshow(pdiffs, cmap='seismic', vmin=-50, vmax=50)
                    ax.scatter(pois['original']['indices'][i], points[j]-pois['original']['indices'][j]-1, c=pois['original']['color'], marker=pois['original']['marker'], s=100, edgecolors='w', linewidth=1, label=pois['original']['label'])
                    ax.scatter(pois[observable]['indices'][i], points[j]-pois[observable]['indices'][j]-1, c=pois[observable]['color'], marker=pois[observable]['marker'], s=100, edgecolors='w', linewidth=1, label=pois[observable]['label'])
                    ax.scatter(pois['global']['indices'][i],   points[j]-pois['global']['indices'][j]-1,   c=pois['global']['color'],   marker=pois['global']['marker'],   s=100, edgecolors='w', linewidth=1, label=pois['global']['label'])
                    for key in ['therm', 'rdf', 'hbond']:
                        if key in pois:
                            ax.scatter(pois[key]['indices'][i], points[j]-pois[key]['indices'][j]-1, c=pois[key]['color'], marker=pois[key]['marker'], s=100, edgecolors='w', linewidth=1, label=pois[key]['label'])
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                    ax.set_title(f'$\\frac{{({custom_name})_{{calc}}-({custom_name})_{{exp}}}}{{({custom_name})_{{exp}}}}$', fontsize=14, pad=10)
                    ax.set_xlabel(f'{converter[parameters[i]]["name"]} ({converter[parameters[i]]["unit"]})', fontsize=11)
                    ax.set_ylabel(f'{converter[parameters[j]]["name"]} ({converter[parameters[j]]["unit"]})', fontsize=11)
                    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                    cbar_ticks = [-50, -25, 0, 25, 50]
                    cbar_labels = ['<-50%', '-25%', '0%', '25%', '>50%']
                    cbar.set_ticks(cbar_ticks)
                    cbar.set_ticklabels(cbar_labels)
                    fig.tight_layout()
                    fig.savefig(f'PDF/MISC/{filename}.pdf')
                    fig.savefig(f'PNG/MISC/{filename}.png')

        for i in range(m):

            observable_i  = self.settings['observables'][i]
            points_i      = observable_means[:,i]-experiment[i]
            min_i         = -0.5*experiment[i]
            max_i         =  0.5*experiment[i]
            custom_name_i = copy.deepcopy(converter[observable_i]["name"])
            custom_name_i = "".join([character for character in custom_name_i if character != '$'])

            for j in range(i+1, m):

                observable_j  = self.settings['observables'][j]
                points_j      = observable_means[:,j]-experiment[j]
                min_j         = -0.5*experiment[j]
                max_j         =  0.5*experiment[j]
                custom_name_j = copy.deepcopy(converter[observable_j]["name"])
                custom_name_j = "".join([character for character in custom_name_j if character != '$'])

                # shared key
                combined = f'{observable_i}-{observable_j}'

                filename = f'OBS-2D-{observable_i}-{observable_j}'
                fig, ax = plt.subplots(dpi=600, figsize=(6, 6))
                ax.scatter(points_i,          points_j,          c='k', marker='.', label='simulated', s=25)
                ax.scatter(points_i[pois['original']['index']], points_j[pois['original']['index']], c=pois['original']['color'], marker=pois['original']['marker'], s=100, edgecolors='w', linewidth=1, label=pois['original']['label'])
                ax.scatter(points_i[pois[combined]['index']],   points_j[pois[combined]['index']],   c=pois[combined]['color'],   marker=pois[combined]['marker'],   s=100, edgecolors='w', linewidth=1, label=pois[combined]['label'])
                ax.scatter(points_i[pois['global']['index']],   points_j[pois['global']['index']],   c=pois['global']['color'],   marker=pois['global']['marker'],   s=100, edgecolors='w', linewidth=1, label=pois['global']['label'])
                for key in ['therm', 'rdf', 'hbond']:
                    if key in pois:
                        ax.scatter(points_i[pois[key]['index']], points_j[pois[key]['index']], c=pois[key]['color'], marker=pois[key]['marker'], s=100, edgecolors='w', linewidth=1, label=pois[key]['label'])
                ax.hlines(0, min_i, max_i, colors='k', linestyles='--', alpha=0.5)
                ax.vlines(0, min_j, max_j, colors='k', linestyles='--', alpha=0.5)
                ax.set_xlim(min_i, max_i)
                ax.set_ylim(min_j, max_j)
                ax.set_title(f'Calculated minus experiment, {converter[observable_j]["name"]} vs {converter[observable_i]["name"]}', fontsize=14)
                ax.set_xlabel(f'$({custom_name_i})_{{calc}}-({custom_name_i})_{{exp}}$ ({converter[observable_i]["unit"]})', fontsize=11)
                ax.set_ylabel(f'$({custom_name_j})_{{calc}}-({custom_name_j})_{{exp}}$ ({converter[observable_j]["unit"]})', fontsize=11)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5]))
                ax.tick_params(axis='x', labelsize=9)
                ax.tick_params(axis='y', labelsize=9)
                ax.set_box_aspect(1)
                ax.legend()
                fig.tight_layout()
                fig.savefig(f'PDF/MISC/{filename}.pdf')
                fig.savefig(f'PNG/MISC/{filename}.png')


    def analyze_sampling(self, minimum_burn_in_period=None, maximum_burn_in_period=None, threshold=2):

        # enter directory
        os.chdir(self.settings['output'])

        # observables and parameters
        observables = self.settings['observables']
        parameters  = self.settings['parameters']

        # get number of observables and parameters
        m = len(observables)
        n = len(parameters)

        # to store measurements and samples
        measurements = np.empty((0,m))
        samples      = np.empty((0,n))

        # get distributions and points
        self.read_distributions()
        self.read_parameters()

        # read experiment
        experiment = self.get_experiment()

        # get number of points per parameter
        points = self.settings['points']
        # positions to start in parameter space
        if self.settings['initial']:
            coordinate_indices_sets = [(point,) for point in range(1, points[0]+1)]
        else:
            grid_ranges        = [range(1, point+1) for point in points]
            coordinate_indices_sets = product(*grid_ranges)

        os.chdir('SIM')
        print('\nRemoving burn-in:')
        for ci, coordinate_indices in enumerate(coordinate_indices_sets):

            if self.settings['initial']:
                directory = f'{str(ci+1).zfill(len(str(len(coordinate_indices_sets))))}'
            else:
                directory_parts = []
                for p, parameter in enumerate(parameters):
                    directory_parts += [f'{parameter[:3]}{str(coordinate_indices[p]).zfill(len(str(points[p])))}']
                directory = '-'.join(directory_parts)

            os.chdir(directory)
            blockmeasurementsfile = 'measurements.dat'
            subsamplesfile        = 'samples.dat'
            if os.path.isfile(subsamplesfile):
                try:
                    subsamples      = np.atleast_2d(np.loadtxt(subsamplesfile,      delimiter=',', skiprows=1))[:,1:] # skriprows removes header, slicing the iteration no.
                    blockmeasurements = np.atleast_2d(np.loadtxt(blockmeasurementsfile, delimiter=',', skiprows=1))[:,1:] # skriprows removes header, slicing the iteration no.
                    burn_in_period  = burn_in_geweke(subsamples)
                    reason          = 'Geweke'
                    if minimum_burn_in_period is not None:
                        if minimum_burn_in_period > burn_in_period:
                            burn_in_period = minimum_burn_in_period
                            reason         = 'minimum'
                    if maximum_burn_in_period is not None:
                        if maximum_burn_in_period < burn_in_period:
                            burn_in_period = maximum_burn_in_period
                            reason         = 'maximum'
                    print(f'* {directory}: removed {burn_in_period}/{subsamples.shape[0]} ({reason})')
                    # no point creating it now if we throw away all other data anyway...
                    if burn_in_period < subsamples.shape[0]:
                        samples      = np.append(samples,      subsamples[burn_in_period:],      axis=0)
                        measurements = np.append(measurements, blockmeasurements[burn_in_period:], axis=0)
                except:
                    print(f'WARNING! No samples found for run {directory}! (file exists but is empty)')
            else:
                print(f'WARNING! No samples found for run {directory}! (file does not exist)')
            os.chdir('..')
        os.chdir('..')

        os.makedirs('PDF', exist_ok=True)
        os.makedirs('PNG', exist_ok=True)

        # for KDE plots
        levels        = 11
        levels_values = np.linspace(0, 1, levels)

        # condifence levels
        confidence_levels = {
            '1': {
                'confidence': 0.642, 'linestyle': 'dashed', 'space':  5
            },
            '2': {
                'confidence': 0.954, 'linestyle': 'dotted', 'space': 10
            }
        }

        # for generating plots, relative to plot width = height = 1
        L  = 0.70  # left margin
        WS = 0.00  # width spacing between plots
        R  = 0.05  # right margin, between rightmost plot and colorbar
        CB = 0.85  # horizontal space of colorbar
        B  = 0.70  # bottom margin
        HS = 0.00  # height spacing between plots
        T  = 0.15  # top margin

        print('\n###################### PARAMETERS ######################')

        # before subtracting reference, we need actual mean to calculate sigma points
        mu    = np.mean(samples, axis=0).reshape((1,n))
        Sigma = np.atleast_2d(np.cov(samples, rowvar=False))

        if self.settings['post'] is not None:
            # from settings
            system   = self.settings['system']
            molecule = self.settings['molecule']
            element  = ''.join([i for i in self.settings['particletype'] if not i.isdigit()]).upper()
            # post-processing directories
            os.makedirs('POST', exist_ok=True)
            os.chdir('POST')
            # find sigma points from mu and Sigma
            sigma_points, Wm, Wc = self.compute_sigma_points(mu.reshape((n,)), Sigma, kappa=(3-n))
            # save weights for "statistics" step
            np.save("Wm.npy", Wm)
            np.save("Wc.npy", Wc)
            # calculate mode
            mode = calculate_mode_meanshift(samples)
            # append mode
            extra_points = np.append(sigma_points, mode, axis=0)
            # directory names
            directories = ['MEAN' + suffix for suffix in ([''] + ['+' + str(num) for num in range(1, n+1)] + ['-' + str(num) for num in range(1, n+1)])] + ['MODE']
            for d, directory in enumerate(directories):
                if os.path.isdir(directory):
                    shutil.rmtree(directory)
                os.makedirs(directory)
                os.chdir(directory)
                # fill directory
                with subprocess.Popen(['cp', self.settings["forcefield"], 'act.xml'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                    _, _ = process.communicate()
                # modify force field
                if d < (2*n+1):
                    description = 'sigma point'
                else:
                    description = 'mode'
                print(f'Point #{d}: {extra_points[d]} ({description})')
                self.set_parameters(extra_points[d], molecule, self.settings['base'])
                if molecule in ['water']:
                    os.system(f'printf "2\n1\n" | gmx pdb2gmx -f {system} -o structure.gro -p topology.top -ignh > /dev/null 2>&1')
                    os.system(f'{self.settings["base"]}/openmm2gromacs.py openmm.xml topology.top')
                    #os.system(f'echo -e "integrator = md\ndt = 0.001\nnsteps = 100000" > dummy.mdp')
                    #os.system(f'gmx grompp -f dummy.mdp -c structure.gro -p topology.top -o topology.tpr > /dev/null 2>&1')
                    os.system(f'printf "a OW\nname 3 {element}\na HW1 | a HW2\nname 4 H\nq\n" | gmx make_ndx -f structure.gro -o index.ndx > /dev/null 2>&1')
                else:
                    os.system(f'printf "a {element}\nname 3 {element}\nq\n" | gmx make_ndx -f {self.settings["system"]} -o index.ndx > /dev/null 2>&1')
                # launch job
                #pf = self.settings['platform']
                #if pf == 'csb':
                #    command = f'sbatch -t 24:00:00 --gres=gpu:1 -p CLUSTER-AMD --nodelist=compute-0-[32-36] --nodes=1 -J "{directory}" '
                #elif pf == 'dardel':
                #    command = f'sbatch -A naiss2023-5-531 --nodes=1 -t 24:00:00 -p gpu -J "{directory}" '
                #elif pf == 'davinci':
                #    command = f'sbatch -n 1 -t 24:00:00 --gres=gpu:1 -p regular -C RTX2080Ti -J "{directory}" '
                #else:
                #    command = ''
                #command += (f'{self.settings["base"]}/post-production.py '
                #          + f'-d   {self.settings["post"]} '
                #          + f'-s   {self.settings["system"]} '
                #          + f'-t   {self.get_temperature()} '
                #          + f'-mol {self.settings["molecule"]} '
                #          + f'-pt  {self.settings["particletype"]}')
                #os.system(command)
                os.chdir('..')
            os.chdir('..')

        # parameter errors
        mu            = np.mean(samples, axis=0)
        model_strings = [f'{value:.3f}' for value in self.parameters]
        print(f'Paramater coefficient of variation (CV), using {self.model} = ({", ".join(model_strings)}) as reference:')
        for p, parameter in enumerate(parameters):
            print(f'* {parameter}:')
            cvxx = np.sqrt(np.mean(((samples[:,p]-mu[p])/mu[p])**2))
            print(f'  + CV(X):   {cvxx*100:.1f} %')
            rbxx = (self.parameters[p]-mu[p])/mu[p]
            print(f'  + RB(X,x): {rbxx*100:.1f}')

        # subtract original parameter values
        delta_samples = samples - self.parameters.reshape((1,n))
        delta_mu      = np.mean(delta_samples, axis=0).reshape((1,n))
        Sigma         = np.atleast_2d(np.cov(samples, rowvar=False))

        # full size
        figwidth  = L + n + WS*(n-1) + R + CB
        figheight = B + n + HS*(n-1) + T

        # full KDE figure (6x6)
        filename_param = "parameters"
        fig_param, ax_param = plt.subplots(nrows=n, ncols=n, figsize=(2*figwidth, 2*figheight), dpi=int(600/n))
        for ax_ij in ax_param.flat:
            ax_ij.set_xticks([])
            ax_ij.set_yticks([])

        # iterate over parameters
        for i in range(n):
            parameter_i = parameters[i]
            print(f'* {parameter_i}')
            # get the correct statistics
            keeps_i      = np.zeros(delta_mu.shape, dtype=bool)
            keeps_i[0,i] = True
            mu_i         = delta_mu[keeps_i]
            Sigma_i	 = Sigma[keeps_i.T @ keeps_i]
            # use the diagonal subplot from ax_param as ax_ii
            ax_ii = ax_param[i, i]
            # plot data
            ax_ii.hist(delta_samples[:,i], bins=100, color='r', alpha=1.0)
            # plot data mean
            ax_ii.axvline(mu_i, c='k', linestyle='solid', linewidth=1)
            # plot data confidence intervals
            confidence_width_min =  np.inf
            confidence_width_max = -np.inf
            for confidence_level in confidence_levels:
                std_dev_1d  = np.sqrt(Sigma_i)
                chi2_val_1d = chi2.ppf(confidence_levels[confidence_level]['confidence'], df=1)
                cri_lower   = mu_i - np.sqrt(chi2_val_1d) * std_dev_1d
                cri_upper   = mu_i + np.sqrt(chi2_val_1d) * std_dev_1d
                ax_ii.axvline(cri_lower, c='k', linestyle=confidence_levels[confidence_level]['linestyle'], linewidth=1)
                ax_ii.axvline(cri_upper, c='k', linestyle=confidence_levels[confidence_level]['linestyle'], linewidth=1)
                confidence_width_min = np.min((confidence_width_min, cri_lower[0]))
                confidence_width_max = np.max((confidence_width_max, cri_upper[0]))
            # determine limits for ith observable
            if np.abs(confidence_width_min) > np.abs(confidence_width_max):
                limit_i = np.abs(confidence_width_min)
            else:
                limit_i = np.abs(confidence_width_max)
            ticks = [round(value, 4) for value in np.linspace(-4*limit_i/3, 4*limit_i/3, 5)]
            for j in range(n):
                if i != j:
                    ax_param[i,j].set_ylim(-2.00*limit_i, 2.00*limit_i)
                    for tick in ticks:
                        if j < i:
                            ax_param[i,j].axhline(tick, color='k', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                            ax_param[j,i].axvline(tick, color='w', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                        else:
                            ax_param[i,j].axhline(tick, color='w', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                            ax_param[j,i].axvline(tick, color='k', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                ax_param[j,i].set_xlim(-2.00*limit_i, 2.00*limit_i)
            # set axes describing observables and reapply ticks
            ax_param[-1,i].set_xlabel(rf'$\Delta {converter[parameter_i]["name"][1:]} ({converter[parameter_i]["unit"]})', fontsize=20)
            ax_param[-1,i].set_xticks(ticks)
            ax_param[-1,i].set_xticklabels([f"{tick:.3f}" for tick in ticks], rotation=45, fontsize=16, ha="right", va="top")
            ax_param[ i,0].set_ylabel(rf'$\Delta {converter[parameter_i]["name"][1:]} ({converter[parameter_i]["unit"]})', fontsize=20)
            # special case for the first plot
            if i == 0:
                ax_param[0,1].set_yticks(ticks)
                y_pixel_positions = ax_param[0,1].transData.transform(np.column_stack((np.zeros_like(ticks), ticks)))[:,1]
                new_ticks = ax_param[0,0].transData.inverted().transform(np.column_stack((np.zeros_like(y_pixel_positions), y_pixel_positions)))[:,1]
                # Set yticks on ax[0,0] with proper labels
                ax_param[0,0].set_yticks(new_ticks)
                ax_param[0,0].set_yticklabels([f"{tick:.3f}" for tick in ticks], fontsize=16)
                ax_param[0,1].set_yticks([])
            else:
                ax_param[i,0].set_yticks(ticks)
                ax_param[i,0].set_yticklabels([f"{tick:.3f}" for tick in ticks], fontsize=16)
            # off-diagonal plots
            for j in range(i+1, n):
                parameter_j  = parameters[j]
                print(f'* {parameter_i}+{parameter_j}')
                keeps_ij      = np.zeros(delta_mu.shape, dtype=bool)
                # get the correct statistics
                keeps_ij[0,i] = True
                keeps_ij[0,j] = True
                mu_ij         = delta_mu[keeps_ij].reshape((1,2))
                Sigma_ij      = Sigma[keeps_ij.T @ keeps_ij].reshape((2,2))
                theta         = np.linspace(0, 2 * np.pi, 100)
                eigenvalues_ij, eigenvectors_ij = np.linalg.eigh(Sigma_ij)
                # 2D KDE
                ax_ij = ax_param[i,j]
                ax_ij.set_facecolor([0.2, 0.0, 0.4])
                kde_plot = sns.kdeplot(x=delta_samples[:,j], y=delta_samples[:,i], fill=True, cmap='rainbow', ax=ax_ij, levels=levels)
                ax_ij.scatter(mu_ij[:,1], mu_ij[:,0], c='w', s=30, zorder=2)
                # 2D points
                ax_ji = ax_param[j,i]
                ax_ji.scatter(delta_samples[:, i], delta_samples[:, j], c='r', alpha=0.010)
                ax_ji.scatter(mu_ij[:,0], mu_ij[:,1], c='k', s=30, zorder=2)
                # ellipse
                for confidence_level in confidence_levels:
                    chi2_val_2d     = chi2.ppf(confidence_levels[confidence_level]['confidence'], df=2)  # Degrees of freedom is 2 for 2D
                    radii_ij        = np.sqrt(eigenvalues_ij * chi2_val_2d)
                    ellipse_i       = radii_ij[0] * np.cos(theta)
                    ellipse_j       = radii_ij[1] * np.sin(theta)
                    ellipse         = np.column_stack((ellipse_i, ellipse_j))
                    ellipse_rotated = ellipse @ eigenvectors_ij.T + mu_ij
                    ax_ij.plot(ellipse_rotated[:, 1], ellipse_rotated[:, 0], c='w', linewidth=1)
                    ax_ji.plot(ellipse_rotated[:, 0], ellipse_rotated[:, 1], c='k', linewidth=1)
        norm = mcolors.BoundaryNorm(boundaries=levels_values, ncolors=256)
        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
        sm.set_array([])
        fig_param.subplots_adjust(left=L/figwidth, wspace=WS/figwidth, right=(figwidth-(R+CB))/figwidth, bottom=B/figheight, hspace=HS/figheight, top=(figheight-T)/figheight)
        cbar_ax = fig_param.add_axes([(figwidth-CB)/figwidth, B/figheight, CB/(3*figwidth), (figheight-(T+B))/figheight])
        cbar = fig_param.colorbar(sm, cax=cbar_ax, ticks=levels_values)
        cbar.ax.tick_params(labelsize=20)
        cbar.set_label("posterior density", fontsize=20)
        # automatically align all labels
        fig_param.align_labels()
        # save figure
        fig_param.savefig(f'PNG/{filename_param}.png')
        fig_param.savefig(f'PDF/{filename_param}.pdf')

        print('\n###################### OBSERVABLES #####################')

        # observable errors
        mu          = np.mean(measurements, axis=0)
        exp_strings = [f'{value:.3f}' for value in experiment]
        print(f'Observable errors, with references = ({", ".join(exp_strings)}):')
        for o, observable in enumerate(observables):
            print(f'* {lookup[observable]}:')
            cvyy = np.sqrt(np.mean(((measurements[:,o]-mu[o])/mu[o])**2))
            print(f'  + CV(Y):   {cvyy*100:.1f} %')
            rbyy = (experiment[o]-mu[o])/mu[o]
            print(f'  + RB(Y,y): {rbyy*100:.1f} %')

        # subtract experiment
        delta_measurements = measurements - experiment.reshape((1,m))
        delta_mu           = np.mean(delta_measurements, axis=0).reshape((1,m))
        Sigma              = np.atleast_2d(np.cov(measurements, rowvar=False))

        # full size
        figwidth  = L + m + WS*(m-1) + R + CB
        figheight = B + m + HS*(m-1) + T

        # full KDE figure (6x6)
        filename_obs = "observables"
        fig_obs, ax_obs = plt.subplots(nrows=m, ncols=m, figsize=(2*figwidth, 2*figheight), dpi=int(600/m))
        for ax_ij in ax_obs.flat:
            ax_ij.set_xticks([])
            ax_ij.set_yticks([])

        # iterate over observables
        for i in range(m):
            observable_i = observables[i]
            print(f'* {observable_i}')
            # get the correct statistics
            keeps_i	 = np.zeros(delta_mu.shape, dtype=bool)
            keeps_i[0,i] = True
            mu_i         = delta_mu[keeps_i]
            Sigma_i	 = Sigma[keeps_i.T @ keeps_i]
            # use the diagonal subplot from ax_obs as ax_ii
            ax_ii = ax_obs[i, i]
            # plot data
            ax_ii.hist(delta_measurements[:,i], bins=100, color='r', alpha=1.0)
            # plot data mean
            ax_ii.axvline(mu_i, c='k', linestyle='solid', linewidth=1)
            # plot data confidence intervals
            confidence_width_min =  np.inf
            confidence_width_max = -np.inf
            for confidence_level in confidence_levels:
                std_dev_1d  = np.sqrt(Sigma_i)
                chi2_val_1d = chi2.ppf(confidence_levels[confidence_level]['confidence'], df=1)
                cri_lower   = mu_i - np.sqrt(chi2_val_1d) * std_dev_1d
                cri_upper   = mu_i + np.sqrt(chi2_val_1d) * std_dev_1d
                ax_ii.axvline(cri_lower, c='k', linestyle=confidence_levels[confidence_level]['linestyle'], linewidth=1)
                ax_ii.axvline(cri_upper, c='k', linestyle=confidence_levels[confidence_level]['linestyle'], linewidth=1)
                confidence_width_min = np.min((confidence_width_min, cri_lower[0]))
                confidence_width_max = np.max((confidence_width_max, cri_upper[0]))
            # determine limits for ith observable
            if np.abs(confidence_width_min) > np.abs(confidence_width_max):
                limit_i = np.abs(confidence_width_min)
            else:
                limit_i = np.abs(confidence_width_max)
            ticks = [round(value, 3) for value in np.linspace(-4*limit_i/3, 4*limit_i/3, 5)]
            for j in range(m):
                if i != j:
                    ax_obs[i,j].set_ylim(-2.00*limit_i, 2.00*limit_i)
                    for tick in ticks:
                        if j < i:
                            ax_obs[i,j].axhline(tick, color='k', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                            ax_obs[j,i].axvline(tick, color='w', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                        else:
                            ax_obs[i,j].axhline(tick, color='w', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                            ax_obs[j,i].axvline(tick, color='k', linewidth=0.8, linestyle=(0, (5, 10)), zorder=2)
                ax_obs[j,i].set_xlim(-2.00*limit_i, 2.00*limit_i)
            # set axes describing observables and reapply ticks
            ax_obs[-1,i].set_xlabel(rf'$\Delta {converter[observable_i]["name"][1:]} ({converter[observable_i]["unit"]})', fontsize=20)
            ax_obs[-1,i].set_xticks(ticks)
            ax_obs[-1,i].set_xticklabels([f"{tick:.2f}" for tick in ticks], rotation=45, fontsize=16, ha="right", va="top")
            ax_obs[ i,0].set_ylabel(rf'$\Delta {converter[observable_i]["name"][1:]} ({converter[observable_i]["unit"]})', fontsize=20)
            # special case for the first plot
            if i == 0:
                ax_obs[0,1].set_yticks(ticks)
                y_pixel_positions = ax_obs[0,1].transData.transform(np.column_stack((np.zeros_like(ticks), ticks)))[:,1]
                new_ticks = ax_obs[0,0].transData.inverted().transform(np.column_stack((np.zeros_like(y_pixel_positions), y_pixel_positions)))[:,1]
                # Set yticks on ax[0,0] with proper labels
                ax_obs[0,0].set_yticks(new_ticks)
                ax_obs[0,0].set_yticklabels([f"{tick:.2f}" for tick in ticks], fontsize=16)
                ax_obs[0,1].set_yticks([])
            else:
                ax_obs[i,0].set_yticks(ticks)
                ax_obs[i,0].set_yticklabels([f"{tick:.2f}" for tick in ticks], fontsize=16)
            # off-diagonal plots
            for j in range(i+1, m):
                observable_j  = observables[j]
                print(f'* {observable_i}+{observable_j}')
                keeps_ij      = np.zeros(delta_mu.shape, dtype=bool)
                # get the correct statistics
                keeps_ij[0,i] = True
                keeps_ij[0,j] = True
                mu_ij         = delta_mu[keeps_ij].reshape((1,2))
                Sigma_ij      = Sigma[keeps_ij.T @ keeps_ij].reshape((2,2))
                theta         = np.linspace(0, 2 * np.pi, 100)
                eigenvalues_ij, eigenvectors_ij = np.linalg.eigh(Sigma_ij)
                # 2D KDE
                ax_ij = ax_obs[i,j]
                ax_ij.set_facecolor([0.2, 0.0, 0.4])
                kde_plot = sns.kdeplot(x=delta_measurements[:,j], y=delta_measurements[:,i], fill=True, cmap='rainbow', ax=ax_ij, levels=levels)
                ax_ij.scatter(mu_ij[:,1], mu_ij[:,0], c='w', s=30, zorder=2)
                # 2D points
                ax_ji = ax_obs[j,i]
                ax_ji.scatter(delta_measurements[:, i], delta_measurements[:, j], c='r', alpha=0.010)
                ax_ji.scatter(mu_ij[:,0], mu_ij[:,1], c='k', s=30, zorder=2)
                # ellipse
                for confidence_level in confidence_levels:
                    chi2_val_2d     = chi2.ppf(confidence_levels[confidence_level]['confidence'], df=2)  # Degrees of freedom is 2 for 2D
                    radii_ij        = np.sqrt(eigenvalues_ij * chi2_val_2d)
                    ellipse_i       = radii_ij[0] * np.cos(theta)
                    ellipse_j       = radii_ij[1] * np.sin(theta)
                    ellipse         = np.column_stack((ellipse_i, ellipse_j))
                    ellipse_rotated = ellipse @ eigenvectors_ij.T + mu_ij
                    ax_ij.plot(ellipse_rotated[:, 1], ellipse_rotated[:, 0], c='w', linewidth=1)
                    ax_ji.plot(ellipse_rotated[:, 0], ellipse_rotated[:, 1], c='k', linewidth=1)
        norm = mcolors.BoundaryNorm(boundaries=levels_values, ncolors=256)
        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
        sm.set_array([])
        fig_obs.subplots_adjust(left=L/figwidth, wspace=WS/figwidth, right=(figwidth-(R+CB))/figwidth, bottom=B/figheight, hspace=HS/figheight, top=(figheight-T)/figheight)
        cbar_ax = fig_obs.add_axes([(figwidth-CB)/figwidth, B/figheight, CB/(3*figwidth), (figheight-(T+B))/figheight])
        cbar = fig_obs.colorbar(sm, cax=cbar_ax, ticks=levels_values)
        cbar.ax.tick_params(labelsize=20)
        cbar.set_label("posterior predictive density", fontsize=20)
        # automatically align all labels
        fig_obs.align_labels()
        # save figure
        fig_obs.savefig(f'PNG/{filename_obs}.png')
        fig_obs.savefig(f'PDF/{filename_obs}.pdf')

        print('\nAnalysis done.')


    def statistics(self, num_blocks=5, block_size=10000):

        # enter directory
        os.chdir(self.settings['output'])

        # and into POST
        os.chdir('POST')

        # update observables
        self.settings['observables'] = [
            'energy',    'volume',
            'rXX',       'gXX',
            'hbonddist', 'hbondangle',
            'dielectric',
            'diffusion',
            'isobaricheatcapacity',
            'thermalexpansioncoefficient',
            'isothermalcompressibility'
        ]

        # properties in log space
        log_observables = [] #'diffusion']

        # observables and parameters
        observables = self.settings['observables']
        parameters  = self.settings['parameters']

        # get number of observables and parameters
        m = len(observables)
        n = len(parameters)

        # number of molecules
        with open(self.settings['system'], 'r') as pdb:
            for line in pdb.readlines():
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    number_of_molecules = int(line[22:26].strip())

        # read experiment
        experiment = self.get_experiment()

        # get temperature
        temperature = self.get_temperature()

        # to store data later
        measurements      = np.zeros((2*n+1, m))
        blockmeasurements = np.zeros((2*n+1, num_blocks, m))

        # directories with extra data
        directories_all   = ['MEAN' + suffix for suffix in ([''] + ['+' + str(num) for num in range(1, n+1)] + ['-' + str(num) for num in range(1, n+1)])] + ['MODE']
        directories_found = [o for o in os.listdir('.') if os.path.isdir(o)]
        directories       = [directory for directory in directories_all if directory in directories_found]

        for d, directory in enumerate(directories):
            os.chdir(directory)
            unfinished = True
            while unfinished:
                command = [
                    "gmx", "energy",
                    "-f",  "md.edr",
                    "-fluct_props",
                    "-nmol", str(number_of_molecules)
                ]
                with subprocess.Popen(command, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                    stdout, stderr = process.communicate(input="Potential\nTemperature\nPressure\nVolume\npV\nEnthalpy\n0\n")
                    if process.returncode != 0:
                        print(f'GROMACS ERROR: {stderr}\n')
                        time.sleep(0.1)
                    else:
                        unfinished = False
            # get data from GROMACS output
            measurement = {}
            for line in stdout.split('\n'):
                contents = line.split()
                if line.startswith('Potential'):
                    measurement['energy'] = k_constant * temperature - float(contents[1])
                elif line.startswith('Volume'):
                    if len(contents) == 6:
                        measurement['volume'] = 1000 * float(contents[1]) / number_of_molecules
                elif line.startswith('Coefficient of Thermal Expansion'):
                    measurement['thermalexpansioncoefficient'] = float(contents[6])
                elif line.startswith('Isothermal Compressibility'):
                    measurement['isothermalcompressibility'] = float(contents[4])
                elif line.startswith('Heat capacity at constant pressure'):
                    measurement['isobaricheatcapacity'] = float(contents[7])
            max_distance, max_intensity = xvg2max('rdf.xvg', distribution='RDF', delta_r=0.0002, molecular_volume=measurement['volume']/1000)
            measurement['rXX'] = 10 * max_distance # convert to Å
            measurement['gXX'] = max_intensity
            max_distance, _ = xvg2max('dist.xvg', distribution='PDF')
            measurement['hbonddist'] = 10 * max_distance # convert to Å
            max_angle, _ = xvg2max('angle.xvg', distribution='PDF')
            measurement['hbondangle'] = max_angle
            measurement['dielectric'] = np.loadtxt('epsilon.xvg', comments=["#", "@"])[-1,1]
            with open('msd.xvg', 'r') as msd:
                for line in msd.readlines():
                    if line.startswith("@ s0 legend"):
                        measurement['diffusion'] = float(line.split()[6]) * 1e-5
                        break
            if directory != 'MODE':
                measurements[d] = np.array([measurement[observable] for observable in observables])
            else:
                modeblockmeasurements = np.zeros((num_blocks, m))
            for b in range(num_blocks):
                unfinished = True
                while unfinished:
                    command = [
                        "gmx", "energy",
                        "-f",  "md.edr",
                        "-fluct_props",
                        "-b",    str(b*block_size),
                        "-e",    str((b+1)*block_size),
                        "-nmol", str(number_of_molecules)
                    ]
                    with subprocess.Popen(command, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                        stdout, stderr = process.communicate(input="Potential\nTemperature\nPressure\nVolume\npV\nEnthalpy\n0\n")
                        if process.returncode != 0:
                            print(f'GROMACS ERROR: {stderr}\n')
                            time.sleep(0.1)
                        else:
                            unfinished = False
                # get data from GROMACS output
                blockmeasurement = {}
                for line in stdout.split('\n'):
                    contents = line.split()
                    if line.startswith('Potential'):
                        blockmeasurement['energy'] = k_constant * temperature - float(contents[1])
                    elif line.startswith('Volume'):
                        if len(contents) == 6:
                            blockmeasurement['volume'] = 1000 * float(contents[1]) / number_of_molecules
                    elif line.startswith('Coefficient of Thermal Expansion'):
                        blockmeasurement['thermalexpansioncoefficient'] = float(contents[6])
                    elif line.startswith('Isothermal Compressibility'):
                        blockmeasurement['isothermalcompressibility'] = float(contents[4])
                    elif line.startswith('Heat capacity at constant pressure'):
                        blockmeasurement['isobaricheatcapacity'] = float(contents[7])
                max_distance, max_intensity = xvg2max(f'rdf_block{b}.xvg', distribution='RDF', delta_r=0.0002, molecular_volume=blockmeasurement['volume']/1000)
                blockmeasurement['rXX'] = 10 * max_distance # convert to Å
                blockmeasurement['gXX'] = max_intensity
                max_distance, _ = xvg2max(f'dist_block{b}.xvg', distribution='PDF')
                blockmeasurement['hbonddist'] = 10 * max_distance # convert to Å
                max_angle, _ = xvg2max(f'angle_block{b}.xvg', distribution='PDF')
                blockmeasurement['hbondangle'] = max_angle
                blockmeasurement['dielectric'] = np.loadtxt(f'epsilon_block{b}.xvg', comments=["#", "@"])[-1,1]
                with open(f'msd_block{b}.xvg', 'r') as msd:
                    for line in msd.readlines():
                        if line.startswith("@ s0 legend"):
                            blockmeasurement['diffusion'] = float(line.split()[6]) * 1e-5
                            break
                if directory != 'MODE':
                    blockmeasurements[d,b] = np.array([blockmeasurement[observable] for observable in observables])
                else:
                    modeblockmeasurements[b] = np.array([blockmeasurement[observable] for observable in observables])
            os.chdir('..')

        print('Block averages:')
        for d, directory in enumerate(directories):
            print(f'* {directory}:')
            for o, observable in enumerate(observables):
                if directory == 'MEAN':
                    with open('blocks.csv', 'w') as bf:
                         bf.write(','.join(observables) + '\n')
                         for b in range(num_blocks):
                             bf.write(','.join([str(value) for value in blockmeasurements[d,b].flatten()]) + '\n')
                    blockmeasurements[d,b]
                if directory != 'MODE':
                    print(f'  + {observable}: {np.mean(blockmeasurements[d,:,o])}')
                else:
                    print(f'  + {observable}: {np.mean(modeblockmeasurements[:,o])}')

        # have all directories been found to perform unscented transform?
        if np.all([directory in directories_found for directory in directories_all]):

            # number of sigma points is the number of directories minus one (= the mode)
            num_sigma_points = len(directories)-1

            Wm = np.load('Wm.npy').reshape(1, num_sigma_points)
            Wc = np.load('Wc.npy').reshape(1, num_sigma_points)

            print(f'Mean weights:     {Wm}')
            print(f'Variance weights: {Wc}')

            print(np.tensordot(Wm, blockmeasurements, axes=(1, 0)).shape)

            # blockmeasurements.shape = (num_sigma_points, num_blocks, m), where m is the number of observables
            measurements   = np.mean(blockmeasurements, axis=1)
            mean           = (Wm @ measurements).reshape(1, -1)
            variance       = (Wc @ (measurements - mean)**2).flatten()
            standard_error = (np.std(np.tensordot(Wm, blockmeasurements, axes=(1, 0)), ddof=1, axis=1) / np.sqrt(num_blocks)).flatten()

            # flatten mean to be safe
            mean = mean.flatten()

            # convert some of the errors from log space back into linear space
            processed_mean               = np.array([np.exp(mean[o]+0.5*variance[o])                                                       if (observable in log_observables) else mean[o]              for o, observable in enumerate(observables)])
            processed_standard_deviation = np.array([0.5*(np.exp(mean[o] + np.sqrt(variance[o])) - np.exp(mean[o] + np.sqrt(variance[o]))) if (observable in log_observables) else np.sqrt(variance[o]) for o, observable in enumerate(observables)])
            processed_standard_error     = np.array([0.5*(np.exp(mean[o] + standard_error[o])    - np.exp(mean[o] + standard_error[o]))    if (observable in log_observables) else standard_error[o]    for o, observable in enumerate(observables)])
            #experiment = np.array([np.log(experiment[o]) if (observable in log_observables) else experiment[o] for o, observable in enumerate(observables)])

            print('MEANS:')
            print(mean)
            print('VARIANCE:')
            print(variance)
            print('STANDARD ERRORS:')
            print(standard_error)

            for o, observable in enumerate(observables):
                cv  = processed_standard_deviation[o]   / processed_mean[o]
                rb  = (experiment[o]-processed_mean[o]) / processed_mean[o]
                rse = processed_standard_error[o]       / processed_mean[o]
                print(f'* {lookup[observable]}:')
                print(f'  + y:            {experiment[o]}')
                print(f'  + E[Y]:         {processed_mean[o]}')
                print(f'  + sqrt(Var(Y)): {processed_standard_deviation[o]}')
                print(f'  + SE(Y):        {processed_standard_error[o]}')
                print(f'  + CV(Y):        {cv*100:.1f} %')
                print(f'  + RB(Y,y):      {rb*100:.1f} %')
                print(f'  + RSE(Y):       {rse*100:.1f} %')

            print('\nExtra statistics from unscented transform sigma point calculations in GROMACS done.')


    def run(self):

        action = self.settings['action']

        if action == 'submit':
            self.submit()
        elif action == 'analyze':
            if self.settings['grid']:
                self.analyze_grid()
            else:
                self.analyze_sampling(maximum_burn_in_period=None)
        elif action == 'statistics':
            self.statistics()
        else:
            sys.exit(f"Unknown action {action}. Exiting...")



def main():

    # parse arguments
    args = vars(parse_arguments())

    # intiialize the overhead
    OH = Overhead(args)

    # start!
    OH.run()



if __name__ == "__main__":

    main()
