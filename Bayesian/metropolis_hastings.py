#!/usr/bin/env python3

from act_openmm import ActOpenMMSim
from scipy.integrate import quad
from scipy.optimize import curve_fit, fmin, fsolve, OptimizeWarning
from scipy.stats import beta, ks_2samp, mannwhitneyu, multivariate_normal, norm, ttest_ind
from scipy.special import erf
import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt
import argparse, copy, os, pickle, subprocess, sys, time, warnings, xmltodict



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
           'water':             'water',
           'hydrogen-fluoride': 'fluorane',
           'hydrogen-chloride': 'chlorane',
           'hydrogen-bromide':  'bromane',
           'hydrogen-iodide':   'iodane',
           'energy':            'dhvap',
           'volume':            'volume',
           'rXX':               'RDF-rXX',
           'gXX':               'RDF-gXX',
           'hbonddist':         'HB-distance',
           'hbondangle':        'HB-angle'
         }



covariance_type = {
                    'energy':     'single',
                    'volume':     'single',
                    'rXX':        'pair',
                    'gXX':        'pair',
                    'hbonddist':  'pair',
                    'hbondangle': 'pair'
                  }



def parse_arguments():

    parser  = argparse.ArgumentParser()

    parser.add_argument('-mol', '--molecule',      help='molecule in MolProps',                                          required=True, type=str)
    parser.add_argument('-pt',  '--particletype',  help='particletype in .xml',                                          required=True, type=str)
    parser.add_argument('-b',   '--base',          help='path to base directory',                                        required=True, type=str)
    parser.add_argument('-exp', '--experiment',    help='path to experimental reference data .csv file',                 required=True, type=str)
    parser.add_argument('-ff',  '--forcefield',    help='path to original ACT .xml file',                                required=True, type=str)
    parser.add_argument('-eq',  '--equilibration', help='path to equilibration .dat file',                               required=True, type=str)
    parser.add_argument('-sim', '--simulation',    help='path to simulation .dat file',                                  required=True, type=str)
    parser.add_argument('-pdb', '--system',        help='path to system .pdb file',                                      required=True, type=str)
    parser.add_argument('-d',   '--distributions', help='path to prior distributions .csv file',                         required=True, type=str)
    parser.add_argument('-c',   '--coverage',	   help='sampling coverage along each dimension',                        required=True, type=float)
    parser.add_argument('-m',   '--samples',       help='number of samples to sample',                                   required=True, type=int)
    parser.add_argument('-n',   '--points',        help='number of points per parameter',                                required=True, type=int, nargs='+')
    parser.add_argument('-i',   '--indices',       help='indices of the parameters on the grid',                         required=True, type=int, nargs='+')
    parser.add_argument('-z',   '--initial',       help='initialize all n runs using parameters of initial force field', action='store_true')
    parser.add_argument('-p',   '--parameters',    help='parameters to include',                                         required=True, type=str, nargs='+')
    parser.add_argument('-o',   '--observables',   help='observables to include',                                        required=True, type=str, nargs='+')
    parser.add_argument('-r',   '--resume',        help='whether this job resumes a previous job',                       action='store_true')
    parser.add_argument('-v',   '--verbose',       help='print more intermediates',                                      action='store_true')
    args = parser.parse_args()

    return args



def xml2dict(path):

    with open(path, 'r') as xml:
        data = xml.read()
        dictionary = xmltodict.parse(data)
        return dictionary



def integrand(x, sigma, min, max, value):

    if min == 0.0:
        return norm.cdf(np.log(value*max/x)/sigma)/(value*(max))
    else:
        return (norm.cdf(np.log(value*max/x)/sigma)-norm.cdf(np.log(value*min/x)/sigma))/(value*(max-min))



def initial_sigma(sigma, min, max, value, area):

    probability = quad(integrand, value*min, value*max, args=(sigma, min, max, value))[0]

    return (probability - area)



class Model:


    # initialize an ACT OpenMM model which can be modified
    def __init__(self, args, output=None):

        # store args as settings
        self.settings = args

        # to print during simulation
        self.output = output

        # print more
        self.verbose = self.settings['verbose']

        # read hyperparameters of simulation from datfile
        self.read_datfile()

        # get temperature
        temperature_raw = self.hyperparameters['temperature_c']
        temperature_parts = temperature_raw.split('*')
        if len(temperature_parts) > 1:
            if temperature_parts[1] == 'kelvin':
                self.temperature = float(temperature_parts[0])
            else:
                sys.exit('Conversion from {temperature_parts[1]} to Kelvin not implemented. Exiting...')
        else:
            self.temperature = float(temperature_parts[0])

        # get timestep in picoseconds
        timestep       = self.hyperparameters['dt']
        timestep_parts = timestep.split('*')
        if len(timestep_parts) > 1:
            if timestep_parts[1] == 'picoseconds':
                self.dt = float(timestep_parts[0])
            else:
                sys.exit('Conversion from {timestep_parts[1]} to picoseconds not implemented. Exiting...')
        else:
            self.dt = float(timestep_parts[0])

        # get number of molecules
        with open(self.settings['system'], 'r') as pdb:
            for line in pdb.readlines():
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    number_of_molecules = int(line[22:26].strip())
            self.number_of_molecules = number_of_molecules
            if self.verbose and self.output is not None:
                self.output.write(f'There are a total of {self.number_of_molecules} molecules in the box.\n\n')

        # read from ACT forcefield field
        self.read_parameters()

        # read from expdata.csv
        self.read_experiment()


    def read_datfile(self):

        self.hyperparameters = {}
        with open(self.settings['simulation'], 'r') as dat:
            for line in dat:
                line_no_comment = line.split('#')[0]
                try:
                    key, value                = [v.strip() for v in line_no_comment.split('=')]
                    self.hyperparameters[key] = value
                except:
                    continue


    def read_parameters(self):

        # read from .xml
        data = xml2dict(self.settings["forcefield"])

        # to store parameter values
        parameters = np.zeros((len(self.settings['parameters']),))

        # store all particles
        self.particletypes = {}
        for particle in data['alexandria_chemistry']['particletypes']['particletype']:
            self.particletypes[particle['@identifier']] = {}
            for option in particle['option']:
                self.particletypes[particle['@identifier']][option['@key']] = option['@value']
        if self.settings['particletype'] not in self.particletypes:
            sys.exit(f'Did not find particletype {self.settings["particletype"]} in .xml file. Could only find {", ".join([particle for particle in self.particletypes.keys()])}. Exiting...')

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

        self.parameters = parameters


    def read_experiment(self):

        columns = [lookup[observable] for observable in self.settings['observables']]

        with open(self.settings['experiment'], 'r') as ef:
            lines   = ef.readlines()
            indices = {column: c for c, column in enumerate(lines[0].split(',')) if column in columns}
            for line in lines[1:]:
                data = line.split(',')
                if data[0] == lookup[self.settings['molecule']] and float(data[1]) == self.temperature:
                    self.experiment = np.array([float(data[indices[column]]) for column in columns])
                    return True
            sys.exit(f'Could not find experimental data for {lookup[self.settings["molecule"]]} ')


    def get_parameters(self):

        return self.parameters


    # update ACT forcefield
    def set_parameters(self, parameters, molecule, base):

        if self.output is not None:
            self.output.write(f'Setting parameters to {parameters}\n')

        # add these to current parameters
        self.parameters = parameters

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
                        if self.verbose and self.output is not None:
                            self.output.write(f'ALEXANDRIA ERROR: {stderr}\n')
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
                            if self.verbose and self.output is not None:
                                self.output.write(f'ALEXANDRIA ERROR: {stderr}\n')
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
            '-ff',      'act.xml',
            '-db',      molecule,
            '-charges', f'{base}/MolProps/mp2-aug-cc-pvtz.xml',
            '-openmm', 'openmm.xml'
        ]
        with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                if self.verbose and self.output is not None:
                    self.output.write(f'ALEXANDRIA ERROR: {stderr}\n')



    def xvg2obs(self, iteration=0):

        # observables, for which we expect to find XVG files
        xvg_observables = ['energy', 'volume', 'rXX', 'gXX', 'hbonddist', 'hbondangle']

        # store for plotting (if verbose)
        popts            = {}
        running_averages = {}

        # check outputs of sampled observables
        for o, observable in enumerate(self.settings['observables']):

           if observable in xvg_observables:

               data = np.loadtxt(f'out_{observable}.xvg', dtype=float, skiprows=1).T

               # convert interaction energy to enthalpy of vaporization
               if observable == 'energy':
                   data[1] = - data[1] / self.number_of_molecules + k_constant * self.temperature
               elif observable == 'volume':
                   data[1] = 1000 * data[1] / self.number_of_molecules # convert to Å^3
               elif observable in ['rXX', 'hbonddist']:
                   data[1] = 10 * data[1] # convert to Å

               times  = data[0]
               values = data[1]

               try:
                   observables[o+1] = values
               except:
                   observables    = np.zeros((len(self.settings['observables'])+1, len(values)))
                   observables[0] = times
                   observables[1] = values

           else:

               sys.exit(f'Not sure how to collect data for {observable}. Exiting...')

        # save all computed observables
        self.observables = observables


    def simulate(self, equilibrate=False, iteration=0):

        # from settings
        if equilibrate:
            datfile = self.settings['equilibration']
        else:
            datfile = self.settings['simulation']
        system = self.settings['system']

        # output files
        system_chk             = 'checkpoint.chk'
        system_csv             = 'energies.csv'
        system_txt             = 'output.txt'
        system_pdb             = 'trajectory.pdb'
        system_xtc             = 'trajectory.xtc'
        system_energy_xvg      = f'out_energy.xvg'
        system_temperature_xvg = f'out_temperature.xvg'
        system_volume_xvg      = f'out_volume.xvg'
        system_density_xvg     = f'out_density.xvg'
        system_fnl             = 'final.pdb'

        # remove old .xtc file to prevent crashes
        command = ['rm', '-f', system_xtc]
        with subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
            _, _ = process.communicate()

        # create system
        system_opt = ActOpenMMSim(system, datfile,
            actfile     = None,
            xmlfile     = 'openmm.xml',
            chkfile     = system_chk,
            enefile     = system_csv,
            txtfile     = system_txt,
            pdbtraj     = system_pdb,
            xtctraj     = system_xtc
        )
        system_opt.set_monomer_energy(0) # TODO! Get minimum energy from .xml file
        system_opt.run()
        if not equilibrate:
            system_opt.log_to_xvg(system_energy_xvg,      [ 'Potential Energy (kJ/mole)'])
            system_opt.log_to_xvg(system_temperature_xvg, [ 'Temperature (K)' ])
            system_opt.log_to_xvg(system_volume_xvg,      [ "Box Volume (nm^3)" ])
            system_opt.log_to_xvg(system_density_xvg,     [ 'Density (g/mL)' ])
        system_opt.write_coordinates(system_fnl)

        # remove object to be on the safe side
        del system_opt


    def from_rdf(self, x, y, dr, molecular_volume, coordination_number=4.96):

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


    def from_pdf(self, x, y):

        return np.sum(x*y)/np.sum(y), None


    def xvg2max(self, filename, distribution=None, delta_r=0.0002, molecular_volume=0.030):

        # Read XVG data
        data = np.loadtxt(filename, comments=["#", "@"])
        x    = data[:, 0]
        y    = data[:, 1]

        if np.all(np.isnan(y)):
            return x[-1], 0.0

        if distribution == 'RDF':
            x_max, y_max = self.from_rdf(x, y, delta_r, molecular_volume)
        elif distribution == 'PDF':
            x_max, y_max = self.from_pdf(x, y)
        else:
            sys.exit('Unknown distribution "{distribution}" for XVG analysis. Exiting...')

        return x_max, y_max


    def xtc2xvgs(self, iteration=0, scratch='/scratch/burst/nordman'):

        # simulation parameters
        savetime    = self.dt*int(self.hyperparameters['saveXtc'])
        saved_steps = int(self.hyperparameters['steps']) // int(self.hyperparameters['saveXtc'])

        # prefix to keep track of scratch files TODO! add directory name to prefix
        prefix = f"{scratch}/{'-'.join([str(index) for index in self.settings['indices']])}"

        # path to XTC file
        xtcfile = 'trajectory.xtc'

        # run GROMACS and get observable data
        observables = {}
        for observable in self.settings['observables']:
            if observable in ['rXX', 'gXX', 'hbonddist', 'hbondangle']:
                observables[observable] = np.zeros((saved_steps,))
        if ('rXX' in observables) or ('gXX' in observables):
            molecular_volumes = np.loadtxt(f'out_volume.xvg', comments=["#", "@"])[:,1]/self.number_of_molecules
        if not os.path.isfile(xtcfile):
            print(f"Could not find the XTC file for iteration {iteration}...")
        elif len(observables) > 0:
            unfinished = True
            while unfinished:
                command = [
                    "gmx",   "trjconv",
                    "-f",     xtcfile,
                    "-s",    "topology.tpr",
                    "-o",   f"{prefix}-frame.g96",
                    "-sep"
                ]
                with subprocess.Popen(command, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                    stdout, stderr = process.communicate(input=f"0\n")
                    if process.returncode != 0:
                        if self.verbose and self.output is not None:
                            self.output.write(f'GROMACS ERROR: {stderr}\n')
                        time.sleep(0.1)
                    else:
                        unfinished = False
            element = ''.join([i for i in self.settings['particletype'] if not i.isdigit()]).upper()
            for step in range(0, saved_steps):
                frame = f"{prefix}-frame{step}.g96"
                if ('rXX' in observables) or ('gXX' in observables):
                    rdffile = f"{prefix}-rdf.xvg"
                    delta_r = 0.0002 # default = 0.002
                    unfinished = True
                    while unfinished:
                        command = [
                            "gmx",    "rdf",
                            "-f",      frame,
                            "-n",     "index.ndx",
                            "-o",      rdffile,
                            "-bin",    str(delta_r),
                            "-rmax",  "0.8",
                        ]
                        with subprocess.Popen(command, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                            stdout, stderr = process.communicate(input=f"{element}\n{element}\n")
                            if process.returncode != 0:
                                if self.verbose and self.output is not None:
                                    self.output.write(f'GROMACS ERROR: {stderr}\n')
                                time.sleep(0.1)
                            else:
                                unfinished = False
                    max_distance, max_intensity = self.xvg2max(rdffile, distribution='RDF', delta_r=delta_r, molecular_volume=molecular_volumes[step])
                    if 'rXX' in observables:
                        observables['rXX'][step] = max_distance
                    if 'gXX' in observables:
                        observables['gXX'][step] = max_intensity
                    with subprocess.Popen(["rm", "-f", rdffile], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                        _, _ = process.communicate()
                if (self.settings["molecule"] in ["water"]) and (('hbonddist' in observables) or ('hbondangle' in observables)):
                    numfile   = f"{prefix}-num.xvg"
                    distfile  = f"{prefix}-dist.xvg"
                    anglefile = f"{prefix}-angle.xvg"
                    unfinished = True
                    while unfinished:
                        command = [
                            "gmx",    "hbond-legacy", "-nomerge", "-noda",
                            "-f",      frame,
                            "-s",     "topology.tpr",
                            "-num",    numfile,
                            "-dist",   distfile,
                            "-ang",    anglefile,
                            "-abin",  "0.1",    # default = 1.0
                            "-rbin",  "0.0005", # default = 0.005
                            "-a",     "50",
                            "-r",     "0.26",
                        ]
                        with subprocess.Popen(command, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                            stdout, stderr = process.communicate(input=f"Water\nWater\n")
                            if process.returncode != 0:
                                if self.verbose and self.output is not None:
                                    self.output.write(f'GROMACS ERROR: {stderr}\n')
                                time.sleep(0.1)
                            else:
                                unfinished = False
                    if 'hbonddist' in observables:
                        max_distance, _ = self.xvg2max(distfile, distribution='PDF')
                        observables['hbonddist'][step] = max_distance
                    if 'hbondangle' in observables:
                        max_angle, _ = self.xvg2max(anglefile, distribution='PDF')
                        observables['hbondangle'][step] = max_angle
                    with subprocess.Popen(["rm", "-f", numfile, distfile, anglefile], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                        _, _ = process.communicate()
                with subprocess.Popen(["rm", "-f", frame], text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                    _, _ = process.communicate()
            # save observable data XVG files
            xtarget = "Time (ps)"
            for observable in observables:
                xvgfile = f'out_{observable}.xvg'
                with open(xvgfile, "w") as outfile:
                    outfile.write("@ xaxis label \"%s\"\n" % xtarget)
                    for step in range(1, saved_steps+1):
                        outfile.write(f"{str(np.around(step*savetime, decimals=1)).rjust(10)}  {str(np.around(observables[observable][step-1], decimals=5)).rjust(10)}\n")


    # get statistics from simulations
    def statistics(self, iteration, coordination_number=4.96):

        if self.settings['samples'] == 0:
            np.savetxt('observables.xvg', self.observables.T, fmt='%10.4f', delimiter=' ')

        mu    = np.mean(self.observables[1:], axis=1)
        Sigma = np.atleast_2d(np.cov(self.observables[1:]))*self.number_of_molecules # all observables should be scaled by number of molecules

        for o1 in range(len(self.settings['observables'])):
            if covariance_type[self.settings['observables'][o1]] == 'single':
                o1_scaling = 1/(coordination_number+1)
            else:
                o1_scaling = 0.5*coordination_number/(2*coordination_number-1)
            Sigma[o1,o1] *= o1_scaling
            for o2 in range(o1+1, len(self.settings['observables'])):
                if covariance_type[self.settings['observables'][o2]] == 'single':
                    o2_scaling = 1/(coordination_number+1)
                else:
                    o2_scaling = 0.5*coordination_number/(2*coordination_number-1)
                Sigma[o1,o2] *= np.sqrt(o1_scaling*o2_scaling)
                Sigma[o2,o1] *= np.sqrt(o1_scaling*o2_scaling)

        return mu, Sigma


    # check whether energy and volume has converged
    def check_convergence(self, current_step=0):

        convergence_observables = ['energy', 'volume']

        for observable in convergence_observables:
            self.convergence_data[observable]['current'] = np.loadtxt(f'out_{observable}.xvg', dtype=float, skiprows=1)[:,1]

        if current_step > 0: # earlier run -> compare to current run

            p_values = []
            for observable in convergence_observables:
                p_values += [
                    mannwhitneyu(self.convergence_data[observable]['earlier'], self.convergence_data[observable]['current']                 )[1],
                    ks_2samp(    self.convergence_data[observable]['earlier'], self.convergence_data[observable]['current']                 )[1],
                    ttest_ind(   self.convergence_data[observable]['earlier'], self.convergence_data[observable]['current'], equal_var=False)[1]
                ]

            if self.verbose:
               self.output.write(f"  + step {current_step}: {', '.join([str(p_value) for p_value in p_values])}\n")

            if np.all(np.array(p_values) > 0.05): # earlier run IS similar to current run -> converged
                return True
            else: # earlier run IS NOT similar to current run -> not converged and save current run
                for observable in convergence_observables:
                    self.convergence_data[observable]['earlier'] = self.convergence_data[observable]['current']
                return False

        else: # no earlier run -> not converged and save current run
            for observable in convergence_observables:
                self.convergence_data[observable]['earlier'] = self.convergence_data[observable]['current']
            return False


    # take a step
    def step(self, iteration=0, max_steps=10, maximum_retry_count=3):

        # setup
        converged             = False
        current_step          = 0
        self.convergence_data = {
            'energy': {'earlier': None, 'current': None},
            'volume': {'earlier': None, 'current': None}
        }

        if self.verbose:
            self.output.write(f"* iteration #{iteration} checking for convergence; Mann-Whitney U, Kolmogorov-Smirnov, and t-Test p-values for energy and volume respectively:\n")

        # simulate until results are converged
        retry_count = 0
        while not converged:

            if current_step > max_steps:
                if self.output is not None:
                    self.output.write(f"Observables did not converged within {max_steps} steps. Forced rejection.\n")
                return None, None

            try:
                # simulate
                self.simulate(iteration=iteration)
            except Exception as error:
                if str(error).startswith('Error initializing CUDA'):
                    if retry_count < maximum_retry_count:
                        retry_count += 1
                        time.sleep(10)
                        continue
                    else:
                        self.output.write(f"Simulating failed after maximum={maximum_retry_count} retries: {type(error).__name__} - {error}.\n")
                        sys.exit()
                elif self.output is not None:
                    self.output.write(f"An exception occurred when simulating: {type(error).__name__} - {error}. Forced rejection.\n")
                    return None, None
                else:
                    return None, None

            converged     = self.check_convergence(current_step=current_step)
            current_step += 1

        # extra observables
        self.xtc2xvgs(iteration=iteration)
        # get observables data and check for convergence
        self.xvg2obs(iteration=iteration)

        # get statistics from sample
        mu, Sigma = self.statistics(iteration)

        if self.verbose and self.output is not None:
            np.set_printoptions(edgeitems=np.inf, linewidth=np.inf, precision=3, suppress=False, formatter={'float': lambda x: f"{x: .3e}"})
            self.output.write('* Averages of sample:\n')
            self.output.write(f'  {mu.T}\n')
            self.output.write('* Covariance matrix of sample:\n')
            for i in range(Sigma.shape[0]):
                self.output.write(f'  {Sigma[i]}\n')
            np.set_printoptions(edgeitems=3, linewidth=75, precision=None, suppress=True, formatter=None)

        # store most recent statistics
        self.mu    = mu
        self.Sigma = Sigma

        return mu, Sigma



class MetropolisHastings:


    # initialize Metropolis-Hastings class
    def __init__(self, args, output='mh.out'):

        # store args as settings
        self.settings = args

        # to print during simulation
        if not os.path.isfile(output):
            self.output = open(output, 'w')
            self.output.write('Welcome to Metropolis-Hastings!\n\n')
        else:
            self.output = open(output, 'a')
            self.output.write('\n############################################\n')
            self.output.write('Welcome back to Metropolis-Hastings!\n\n')

        # print more
        self.verbose = self.settings['verbose']

        # more likely to succeed using a known initial model rather than grid point
        if self.settings['initial']:
            equilibration_tries = 10
        else:
            equilibration_tries = 3

        # copy forcefield file if it has not been copied yet
        if not os.path.isfile('act.xml'):
            with subprocess.Popen(['cp', self.settings["forcefield"], 'act.xml'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                _, _ = process.communicate()

        # set up OpenMM model
        self.model = Model(args, output=self.output)

        # get beta values
        self.read_distributions()

        # set boundaries
        self.set_boundaries()

        # keep track of progress
        self.accepted  = 0
        self.iteration = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.std = np.array([fsolve(initial_sigma, 0.1, args=(self.beta_values[parameter]['min'], self.beta_values[parameter]['max'], self.model.get_parameters()[p], self.settings['coverage']))[0] for p, parameter in enumerate(self.settings['parameters'])])

        # to store information
        self.parameters     = {'accepted': np.nan, 'proposed': np.nan}
        self.log_likelihood = {'accepted': np.nan, 'proposed': np.nan}
        self.log_prior_sum  = {'accepted': np.nan, 'proposed': np.nan}
        self.measurements   = {'accepted': np.nan, 'proposed': np.nan}

        if not self.settings['resume']:

            if self.settings['samples'] != 0:
                # save samples
                with open('samples.dat', 'w') as sdat:
                    sdat.write("iteration," + ",".join([parameter for parameter in self.settings['parameters']]) + '\n')
                # save measurements
                with open('measurements.dat', 'w') as odat:
                    odat.write("iteration," + ",".join([observable for observable in self.settings['observables']]) + '\n')

            # equlibration
            self.output.write(f'Preparing equilibration.\n\n')

            # calculate deterministic parameters corresponding to grid point value
            self.model.set_parameters(self.deterministic(), self.settings['molecule'], self.settings['base'])

            # equilibrate with deterministic parameters and retry if it crashes
            successful        = False
            equilibration_try = 1
            while not successful and equilibration_try <= equilibration_tries:
                try:
                    self.model.simulate(equilibrate=True)
                    self.output.write(f'\nSuccessful equilibration after {equilibration_try} deterministic tries.\n\n')
                    successful = True
                except Exception as error:
                    self.output.write(f"An exception occurred when equilibrating (deterministic): {type(error).__name__} - {error}\n")
                    equilibration_try += 1

            # randomize parameters if deterministic parameters fail and equilibrate again
            while not successful and (equilibration_try <= (equilibration_tries+10)) and (self.settings['samples'] != 0):
                try:
                    self.model.set_parameters(self.stochastic(), self.settings['molecule'], self.settings['base'])
                    self.model.simulate(equilibrate=True)
                    self.output.write(f'\nSuccessful equilibration after {equilibration_try-equilibration_tries} stochastic tries.\n\n')
                    successful = True
                except Exception as error:
                    self.output.write(f"An exception occurred when equilibrating (stochastic): {type(error).__name__} - {error}\n")
                    equilibration_try += 1

            # if all fails...
            if not successful:
                sys.exit('Could not equilibrate. Exiting...')

            # simulate
            self.output.write('Running first simulation for initial statistics...\n')
            self.mu, self.Sigma = self.model.step(iteration=self.iteration)

            # save setting
            save_setting = 1

        else:

            # set last accepted parameter values
            self.model.set_parameters(np.load("parameters.npy"), self.settings['molecule'], self.settings['base'])

            # keep track of iteration
            if os.path.isfile('samples.dat'):
                self.std      = np.load("std.npy")
                self.accepted = np.load("accepted.npy")
                with open('samples.dat', 'r') as sdat:
                    self.iteration = len(sdat.readlines())-1

            # get statistics from where we left off
            self.mu    = np.load("mu.npy")
            self.Sigma = np.load("Sigma.npy")

            # and set the same statistics for the model
            self.model.mu    = self.mu
            self.model.Sigma = self.Sigma

            # save setting
            save_setting = 0

        # set up Metropolis-Hastings
        self.parameters['proposed'] = self.model.get_parameters()
        self.update_log_likelihood()
        self.update_log_prior_sum()
        self.measurements['proposed'] = self.mu

        # store proposed as accepted
        self.accept(save_setting=save_setting)

        # confirm that setup is OK
        self.output.write('Setup completed.\n\n')

        # modify labels if necessary
        element = ''.join([i for i in self.settings['particletype'] if not i.isdigit()]).upper()
        converter['rXX']['name'] = fr'$\tilde{{r}}_{{{element}{element},1}}$'
        converter['gXX']['name'] = fr'$\tilde{{g}}_{{{element}{element},1}}$'


    # upon Metropolis-Hastings deletion
    def __del__(self):

        if self.output is not None:
            self.output.write('Thank you for running Metropolis-Hastings.\n')
            self.output.close()


    # deterministic initial parameters from prior distributions and indices
    def deterministic(self):

        if self.settings['initial']:
            return self.model.get_parameters()
        else:
            return np.array([self.min[p]+beta.ppf(self.settings['indices'][p]/(self.settings['points'][p]+1), self.beta_values[parameter]['a'], self.beta_values[parameter]['b']) * (self.max[p]-self.min[p]) for p, parameter in enumerate(self.settings['parameters'])])


    # stochastic initial parameters from prior distributions
    def stochastic(self):

        return np.array([self.min[p]+beta.ppf(np.random.rand(), self.beta_values[parameter]['a'], self.beta_values[parameter]['b'])*(self.max[p]-self.min[p]) for p, parameter in enumerate(self.settings['parameters'])])


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


    # set boundaries or read from disk
    def set_boundaries(self):

        parameters = self.model.get_parameters()

        if os.path.isfile('min.pkl'):
            with open('min.pkl', 'rb') as p:
                self.min = pickle.load(p)
        else:
            min = parameters * np.array([self.beta_values[parameter]['min'] if parameters[p] > 0 else self.beta_values[parameter]['max'] for p, parameter in enumerate(self.settings['parameters'])])
            with open('min.pkl', 'wb') as p:
                pickle.dump(min, p, pickle.HIGHEST_PROTOCOL)
            self.min = min

        if os.path.isfile('max.pkl'):
            with open('max.pkl', 'rb') as p:
                self.max = pickle.load(p)
        else:
            max = parameters * np.array([self.beta_values[parameter]['max'] if parameters[p] > 0 else self.beta_values[parameter]['min'] for p, parameter in enumerate(self.settings['parameters'])])
            with open('max.pkl', 'wb') as p:
                pickle.dump(max, p, pickle.HIGHEST_PROTOCOL)
            self.max = max


    # log-likelihood function using the normal distribution for observables
    def update_log_likelihood(self):

        # calculate log-likelihood
        if (self.mu is None) or (self.Sigma is None):
            log_likelihood = -np.inf
        else:
            log_likelihood = multivariate_normal.logpdf(self.model.experiment, mean=self.mu, cov=self.Sigma, allow_singular=True)

        self.log_likelihood['proposed'] = log_likelihood

        if self.verbose and self.output is not None:
            self.output.write('* log likelihood:\n')
            self.output.write(f'  {self.log_likelihood["proposed"]}\n')


    # scaled log-prior distribution for parameters
    def update_log_prior_sum(self):

        parameters = self.parameters['proposed']

        # if parameters are out of bounds (NB! log(0) = -inf)
        if np.any(parameters <= self.min) or np.any(parameters >= self.max):

            # sum of constants and -infs is inf
            self.log_prior_sum['proposed'] = -np.inf

        # else calculate like normally
        else:

            # scale x to [0, 1]
            parameters_scaled = (parameters - self.min) / (self.max - self.min)

            # iterate over parameters
            log_prior_sum = 0
            for p, parameter in enumerate(self.settings['parameters']):

                # calculate log_priors
                log_prior_sum += beta.logpdf(parameters_scaled[p], self.beta_values[parameter]['a'], self.beta_values[parameter]['b'])

            self.log_prior_sum['proposed'] = log_prior_sum

        if self.verbose and self.output is not None:
            self.output.write('* log prior sum:\n')
            self.output.write(f'  {self.log_prior_sum["proposed"]}\n')


    # proposal function on a log-scale as all parameters are positive
    def propose(self):

        # generate random variates
        random_variates = norm.rvs(scale=self.std)

        # calculate log proposal sum
        self.log_proposal_sum = np.sum(random_variates)

        # get how to factor old parameters
        proposal_ratios = np.exp(random_variates)

        # multiply old parameters to get new factors
        self.parameters['proposed'] = self.parameters['accepted'] * proposal_ratios

        # update forcefield file
        self.model.set_parameters(self.parameters['proposed'], self.settings['molecule'], self.settings['base'])

        if self.verbose and self.output is not None:
            self.output.write('* log proposal sum:\n')
            self.output.write(f'  {self.log_proposal_sum}\n')


    # log acceptance with contributions from log-likelihoods and log-priors
    def acceptance(self):

        log_proposed_sum = self.log_likelihood['proposed'] + self.log_prior_sum['proposed']
        log_accepted_sum = self.log_likelihood['accepted'] + self.log_prior_sum['accepted']

        log_acceptance = log_proposed_sum - log_accepted_sum + self.log_proposal_sum

        if log_acceptance > 0:
            return 1
        else:
            return np.exp(log_acceptance)


    # accept a parameter change
    def accept(self, save_setting=2):

        # accept the change in parameters and assign proposed as new accepted
        self.parameters['accepted']     = self.parameters['proposed']
        self.log_likelihood['accepted'] = self.log_likelihood['proposed']
        self.log_prior_sum['accepted']  = self.log_prior_sum['proposed']
        self.measurements['accepted']   = self.measurements['proposed']

        if save_setting > 0:

            # save to resume later
            np.save("parameters.npy", self.parameters['proposed'])
            np.save("std.npy",        self.std)
            np.save("mu.npy",         self.mu)
            np.save("Sigma.npy",      self.Sigma)

            if save_setting > 1:

                # keep track of acceptance rate
                self.accepted += 1
                np.save("accepted.npy", self.accepted)

                # confirm accepted
                if self.output is not None:
                    self.output.write('* accepted\n\n')


    # reject a parameter change
    def reject(self):

        # reset to last checkpoint
        with subprocess.Popen(['cp', 'backup.chk', 'checkpoint.chk'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
            _, _ = process.communicate()

        # confirm rejected
        if self.output is not None:
            self.output.write('* rejected\n\n')


    # sample from simulations
    def sample(self, update_frequency=10):

        for i in range(self.iteration, self.iteration+self.settings['samples']):

            # current iteration
            self.iteration += 1

            # back up checkpoint
            with subprocess.Popen(['cp', 'checkpoint.chk', 'backup.chk'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as process:
                _, _ = process.communicate()

            # propose new parameters
            self.propose()

            # check if any proposals are out-of-bounds
            if np.any(self.parameters['proposed'] <= self.min) or np.any(self.parameters['proposed'] >= self.max):

                alpha = 0

            else:

                # simulate
                self.mu, self.Sigma = self.model.step(iteration=self.iteration)

                # store measurements
                self.measurements['proposed'] = self.mu

                # update log-likelihood of proposed parameters
                self.update_log_likelihood()

                # calculate log-priors of new parameters (the accepted is already calculated)
                self.update_log_prior_sum()

                # get log-acceptance ratio
                alpha = self.acceptance()

            # accept or reject
            if np.random.rand() < alpha:
                self.accept()
            else:
                self.reject()

            # save samples
            with open('samples.dat', 'a') as sdat:
                sdat.write(f"{self.iteration}," + ",".join([str(value) for value in self.parameters['accepted'].flatten()]) + '\n')

            # save measurements
            with open('measurements.dat', 'a') as odat:
                odat.write(f"{self.iteration}," + ",".join([str(value) for value in self.measurements['accepted'].flatten()]) + '\n')

        if self.settings['samples'] > 0:

            # print some statistics
            if self.output is not None:
                self.output.write(f'Standard deviation: {self.std}\n')
                with open('samples.dat', 'r') as sdat:
                    no_samples = len(sdat.readlines())-1 # subtract one for header
                    self.output.write(f'Total acceptance ratio:  {self.accepted / no_samples}\n')



def main():

    # parse args
    args = vars(parse_arguments())

    # initialize MetropolisHastings
    MH = MetropolisHastings(args)

    # start sampling!
    MH.sample()



if __name__ == "__main__":

    main()
