# This script defines classes to implement the calvados model
# All credit to Tesei, Schulze, Crehuet and Lindorff-Larsen
# Read Tesei et al, PNAS (2021). DOI:10.1073/pnas.2111696118
# All errors due to David De Sancho's intervention

import os, sys
import numpy as np
import pandas as pd

from Bio import SeqIO

import openmm
from openmm import unit
from openmm import app
from openmm import XmlSerializer

import mdtraj as md

eps_factor = 0.2

def read_das_sequences(file="das_seqs.dat"):
    '''
    Reads protein sequences from Das et al, PNAS (2013) 

    Parameters
    ----------
    file :str
        File with sequences as list

    Raises
    ------
    IOError            

    '''
    try:
        das_proteins = np.loadtxt(file, skiprows=3, delimiter=',', \
                   dtype={'names': ('seq', 'kappa', 'scd'), \
                  'formats':('S50', 'f', 'f')})
    except IOError as e:
        print (e)
        sys.exit()

    proteins = {}
    for i, prot in enumerate(das_proteins):
        name, sequence = "das%i"%i, str(prot['seq'], 'utf-8')
        proteins[name] = {}
        proteins[name]['name'] = name
        proteins[name]['fasta'] = list([x for x in sequence])
        proteins[name]['eps_factor'] = 0.2
        proteins[name]['pH'] = 7
        proteins[name]['ionic'] = 0.15
    proteins = proteins
    proteins_df = pd.DataFrame.from_dict(proteins, orient='index')
    return proteins_df

def read_fasta_sequences(file_fasta):
    '''
    Reads protein sequences in fasta format

    Parameters
    ----------
    file_fasta :str
        File with sequences in fasta format

    Raises
    ------
    IOError            

    '''
    try:
        fasta_sequences = SeqIO.parse(open(file_fasta), 'fasta')
    except IOError as e:
        print (e)
        pass

    proteins = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        proteins[name] = {}
        proteins[name]['name'] = name
        proteins[name]['fasta'] = list(fasta)
        proteins[name]['eps_factor'] = 0.2
        proteins[name]['pH'] = 7
        proteins[name]['ionic'] = 0.15
    proteins = proteins
    proteins_df = pd.DataFrame.from_dict(proteins, orient='index')
    return proteins_df

def genParamsLJ(df, prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['X','MW'] += 2
    r.loc['Z','MW'] += 16
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    types = list(np.unique(fasta))
    MWs = [r.loc[a,'MW'] for a in types]
    lj_eps = eps_factor*4.184
    return lj_eps, fasta, types, MWs

def genParamsDH(df, prot, temp):
    pH = 7
    ionic = 0.15
    kT = 8.3145*temp*1e-3
    fasta = prot.fasta.copy()
    r = df.copy()
    # Set the charge on HIS based on the pH of the protein solution
    r.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    r.loc['X','q'] = r.loc[prot.fasta[0],'q'] + 1.
    r.loc['Z','q'] = r.loc[prot.fasta[-1],'q'] - 1.
    # Calculate the prefactor for the Yukawa potential
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    yukawa_eps = [r.loc[a].q*np.sqrt(lB*kT) for a in fasta]
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    return yukawa_eps, yukawa_kappa


class CalvadosModel(object):
    '''
    A class for building models with Calvados and running in OpenMM

    '''
    def __init__(self, file_res="residues.csv", temp=300.):
        '''
        Parameters
        ----------
        file_res : str
            File with information about residues
        temp : float
            Temperature for generating parameters.
            Default temperature is 300 K

        '''
        self.read_residues(file_res)
        self.temp = temp 

    def read_residues(self, file_res):
        '''
        Reads residues from file

        Parameters
        ----------
        file_res : str
            File with information about residues

        Raises
        ------
        IOError            

        '''
        try:
            residues = pd.read_csv(file_res).set_index('three', \
                    drop=False)
        except IOError as e:
            print (e)
            return
        self.residues = residues.set_index('one')

    def add_proteins(self, prot, file_fasta="sequences.fasta", \
            file_das='das_seqs.dat'):
        '''
        Adds proteins to object

        Parameters
        ----------
        prot : str, list
            The name of the proteins to simulate
        file_fasta :str
            File with sequences in fasta format

        '''
        #proteins_df = read_fasta_sequences(file_fasta)
        proteins_df = read_das_sequences(file_das)

        try:
            if isinstance(prot, str):
                self.prot = [proteins_df.loc[prot]]
            elif isinstance(prot, list):
                self.prot = [proteins_df.loc[x] for x in prot]
            else:
                print ("Unrecognized data type")
            self.set_params()
        except KeyError as e:
            print (e)
            print ("The protein requested was not found in the proteins database\n")
            print ("Available values are:", [x for x in proteins_df['name']]) 

    def set_params(self):
        '''
        Defines LJ and DH parameters

        '''
        residues = self.residues
        prot = self.prot
        temp = self.temp
        print (prot)
        self.paramsLJ = [genParamsLJ(residues, x) for x in prot]
        self.paramsDH = [genParamsDH(residues, x, temp) for x in prot]

class OMMsystem(object):
    '''
    A class for running simulation using the Calvados model using OpenMM

    '''
    def __init__(self, model, n_chains=1, name=None, box=None):
        '''
        Parameters
        ----------
        model : object
            Calvados models for running the simulations

        n_chains : int, list
            Number of chains. Either a single integer value for all proteins
            or a list of values with length equal to the number of proteins.
        name : str
            Root name for output
        box : int/list
            Length for cubic boxes and x,y,z values for slabs

        '''
        system = openmm.System()
        self.system = system

        self.model = model
        try:
            assert len(self.model.prot) == len(n_chains)
            self.n_chains = n_chains
        except AssertionError:
            self.n_chains = [n_chains for p in range(len(self.model.prot))]
        except TypeError:
            self.n_chains = [n_chains for p in range(len(self.model.prot))]
         
        if not name:
            self.name = '-'.join("%s_%i"%(p.name, self.n_chains[i]) for i,p in enumerate(model.prot))
        else:
            self.name = name
        print ("Setting root for filenames: %s"%self.name)

        self.set_box_vectors(box)

        self.set_config(box)
        self.set_forcefield()

    def set_box_vectors(self, box):
        '''
        Sets box vectors. Uses length of sequence as proxy for box size

        Parameters
        ----------
        box : int/list
            Length for cubic boxes and x,y,z values for slabs

        '''
        N = np.max([len(p.fasta) for p in self.model.prot])
        L = N*0.38+10

        a = unit.Quantity(np.zeros([3]), unit.nanometers)
        b = unit.Quantity(np.zeros([3]), unit.nanometers)
        c = unit.Quantity(np.zeros([3]), unit.nanometers)
        if not box:
            print (" Automatically setting cubic box length: ", L)
            a[0] = L * unit.nanometers
            b[1] = L * unit.nanometers
            c[2] = L * unit.nanometers
        else:
            print (" Manually setting box size: ", box)
            try:
                x,y,z = box
            except TypeError:
                x = box; y = box; z = box
            a[0] = x * unit.nanometers
            b[1] = y * unit.nanometers
            c[2] = z * unit.nanometers

        self.system.setDefaultPeriodicBoxVectors(a, b, c)

    def set_config(self, box):
        '''
        Defines initial configuration 

        Parameters
        ----------
        box : int/list
            Length for cubic boxes and x,y,z values for slabs

        '''
        if isinstance(box, int):
            self.build_cube(box)
        else:
            self.build_slab(box)

    def build_cube(self, box):
        '''
        Builds cubic configuration

        Parameters
        ----------
        box : list
            Length for cubic boxes and x,y,z values for slabs

        '''
        N = np.max([len(p.fasta) for p in self.model.prot])
        margin = 1
        L = box

        n_chains = self.n_chains 
        prot = self.model.prot

        top = md.Topology()
        pos = []

        for i,p in enumerate(prot):
            print (" Creating XY coordinates for %i chains of protein %s\n"%(n_chains[i], p.name))
            ni = 0 
            xyz = np.empty(0)
            while True:
                x,y,z = np.random.rand(3)*(L-margin)-(L-margin)/2
                x1 = x-L if x>0 else x+L
                y1 = y-L if y>0 else y+L
                z1 = z-L if z>0 else z+L
                try:
                    if np.all(np.linalg.norm(xyz-[x,y,z],axis=1)>.6):
                        if np.all(np.linalg.norm(xyz-[x1,y,z],axis=1)>.6):
                            if np.all(np.linalg.norm(xyz-[x,y1,z],axis=1)>.6):
                                if np.all(np.linalg.norm(xyz-[x,y,z1],axis=1)>.6):
                                    xyz = np.append(xyz,[x,y,z]).reshape((-1,3))
                                    ni += 1
                except ValueError:
                    xyz = np.append(xyz, np.random.rand(3)*(L-margin)-(L-margin)/2).reshape((-1,3))
                    ni += 1
                if ni >= n_chains[i]:
                    ni = 0
                    break

            for x, y, z in xyz[:n_chains[i]]:
                chain = top.add_chain()
                N = len(p.fasta)
                pos.append([[x,y,z +(i-N/2.)*.38] for i in range(N)])
                for j,resname in enumerate(p.fasta):
                    resname3 = self.model.residues.loc[resname, 'three']
                    residue = top.add_residue(resname3, chain, resSeq=j+1)
                    top.add_atom('CA', element=md.element.carbon, \
                            residue=residue, serial=j+1)
                for j in range(chain.n_atoms-1):
                    top.add_bond(chain.atom(j), chain.atom(j+1))

        self.top = top
        try:
            assert(os.path.isdir('data'))
        except AssertionError as e:
            print (e)
            os.makedirs('data')
        md.Trajectory(np.vstack(pos), top, 0, \
              [L, L, L], [90,90,90]).save_pdb('data/%s_top.pdb'%self.name)

    def build_slab(self, box):
        '''
        Builds slab configuration

        Parameters
        ----------
        box : list
            Length for cubic boxes and x,y,z values for slabs

        '''
        N = np.max([len(p.fasta) for p in self.model.prot])
        L = N*0.38+10
        margin = 1
        if box:
            L = np.min(box)

        n_chains = self.n_chains 
        prot = self.model.prot

        top = md.Topology()
        pos = []

        for i,p in enumerate(prot):
            print (" Creating XY coordinates for %i chains of protein %s\n"%(n_chains[i], p.name))
            ni = 0 
            xy = np.empty(0)
            while True:
                x,y = np.random.rand(2)*(L-margin)-(L-margin)/2
                x1 = x-L if x>0 else x+L
                y1 = y-L if y>0 else y+L
                try:
                    if np.all(np.linalg.norm(xy-[x,y],axis=1)>.6):
                        if np.all(np.linalg.norm(xy-[x1,y],axis=1)>.6):
                            if np.all(np.linalg.norm(xy-[x,y1],axis=1)>.6):
                                xy = np.append(xy,[x,y]).reshape((-1,2))
                                ni += 1

                except ValueError:
                    xy = np.append(xy,np.random.rand(2)*(L-margin)-(L-margin)/2).reshape((-1,2))
                    ni += 1

                if ni >= n_chains[i]:
                    ni = 0
                    break

            for x, y in xy[:n_chains[i]]:
                chain = top.add_chain()
                N = len(p.fasta)
                pos.append([[x,y,L/2+(i-N/2.)*.38] for i in range(N)])
                for j,resname in enumerate(p.fasta):
                    resname3 = self.model.residues.loc[resname, 'three']
                    residue = top.add_residue(resname3, chain, resSeq=j+1)
                    top.add_atom('CA', element=md.element.carbon, \
                            residue=residue, serial=j+1)
                for j in range(chain.n_atoms-1):
                    top.add_bond(chain.atom(j), chain.atom(j+1))

        self.top = top
        try:
            assert(os.path.isdir('data'))
        except AssertionError as e:
            print (e)
            os.makedirs('data')
        md.Trajectory(np.vstack(pos), top, 0, \
              [box[0], box[1], box[2]], [90,90,90]).save_pdb('data/%s_top.pdb'%self.name)

    def set_forcefield(self):
        '''
        Adds particles to system and creates interactions

        '''
        n_chains = self.n_chains
        prot = self.model.prot
        residues = self.model.residues
        paramsDH = self.model.paramsDH
        paramsLJ = self.model.paramsLJ

        for k,p in enumerate(prot):
            for _ in range(n_chains[k]):
                self.system.addParticle((residues.loc[p.fasta[0]].MW+2)*unit.amu)
                for a in p.fasta[1:-1]:
                    self.system.addParticle(residues.loc[a].MW*unit.amu) 
                self.system.addParticle((residues.loc[p.fasta[-1]].MW+16)*unit.amu)

        hb = openmm.openmm.HarmonicBondForce()
        energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6),4*eps*((s/r)^12-(s/r)^6)+eps*(1-l))'
        ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=0.5*(l1+l2)')
        yu = openmm.openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r - exp(-kappa*4)/4); q=q1*q2')
        
        yukawa_kappa = paramsDH[0][1]
        yu.addGlobalParameter('kappa', yukawa_kappa/unit.nanometer)
        yu.addPerParticleParameter('q')

        lj_eps = paramsLJ[0][0]
        ah.addGlobalParameter('eps', lj_eps*unit.kilojoules_per_mole)
        ah.addPerParticleParameter('s')
        ah.addPerParticleParameter('l')
        
        N_prev = 0
        for k, p in enumerate(prot):
            N = len(p.fasta)
            for j in range(n_chains[k]):
                begin = j*N + N_prev
                end = j*N+N + N_prev 
                yukawa_eps = paramsDH[k][0]
                for a,e in zip(p.fasta, yukawa_eps):
                    yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
                    ah.addParticle([residues.loc[a].sigmas*unit.nanometer, \
                            residues.loc[a].lambdas*unit.dimensionless])
                for i in range(begin,end-1):
                    hb.addBond(i, i+1, 0.38*unit.nanometer, \
                            8033.28*unit.kilojoules_per_mole/(unit.nanometer**2))
                    yu.addExclusion(i, i+1)
                    ah.addExclusion(i, i+1)
            N_prev = end
        
        yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
        ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
        hb.setUsesPeriodicBoundaryConditions(True)
        yu.setCutoffDistance(4*unit.nanometer)
        ah.setCutoffDistance(4*unit.nanometer)
        
        self.system.addForce(hb)
        self.system.addForce(yu)
        self.system.addForce(ah)
        
        serialized_system = XmlSerializer.serialize(self.system)
        outfile = open('system.xml','w')
        outfile.write(serialized_system)
        outfile.close()

class OMMrunner(object):
    def __init__(self, system, platform='CPU'):
        '''
        Parameters
        ----------
        system : object
            The OMMsystem object
        platform : str
            The platform for running simulations. Options are CPU (default) and CUDA

        '''
        self.system = system

        temp = system.model.temp
        name = system.name

        self.pdb = app.pdbfile.PDBFile('data/%s_top.pdb'%(name))

        integrator = openmm.openmm.LangevinIntegrator(temp*unit.kelvin,\
                0.01/unit.picosecond, 0.005*unit.picosecond)

        platform = openmm.Platform.getPlatformByName(platform)

        self.simulation = app.simulation.Simulation(self.pdb.topology, system.system, \
                integrator, platform)

        self.gen_file_names(name, temp)

        self.initialize()


    def gen_file_names(self, name, temp):
        '''
        Generates file names for runs

        Parameters
        ----------
        name: str
            Name for output files
        temp: float
            Temperature in K

        '''
        self.cpt = 'data/%s_T%gK.chk'%(name, temp)
        self.dcd = 'data/%s_T%gK.dcd'%(name, temp)
        self.log = 'data/%s_T%gK.log'%(name, temp)

    def initialize(self):
        '''
        Initializes the simulation either from checkpoints or 
        from initial coordinates and generates reporters

        '''
        if os.path.isfile(self.cpt):
            self.simulation.loadCheckpoint(self.cpt)
            self.simulation.reporters.append(app.dcdreporter.DCDReporter(self.dcd, \
                                    int(1e3), append=True))
        else:
            self.simulation.context.setPositions(self.pdb.positions)
            self.simulation.minimizeEnergy()
            self.simulation.reporters.append(app.dcdreporter.DCDReporter(self.dcd, \
                                    int(1e3)))
        self.simulation.reporters.append(app.statedatareporter.StateDataReporter(self.log, \
                int(1e3), potentialEnergy=True, temperature=True, step=True, \
                  speed=True, elapsedTime=True,separator='\t'))

    def run(self, time=0.1):
        '''
        Runs MD simulation

        Parameters
        ----------
        time : float
            Length of simulation run (in hours)

        '''
        self.simulation.runForClockTime(time*unit.hour, checkpointFile=self.cpt, \
                           checkpointInterval=1*unit.hour)
