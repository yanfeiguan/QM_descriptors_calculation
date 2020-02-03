#!/usr/bin/python
from __future__ import print_function, absolute_import

from multiprocessing import Process, Manager
#######################################################################
# Molecular conformer generator in progress
# genConf.py -isdf file_input.sdf -osdf file_output.sdf
# -n number_of_conformers (optional, if not specified is based
# on the nomber of rotable bonds) -rtpre rms_threshold_pre_opt(optional)
# -rtpost rms_threshold_post_opt(optional) -e energy_window (optional, Kcal/mol)
# -t number_of_threads (if not specify 1)
#######################################################################

## known issues / to-do list
## logging doesn't work properly
## need to print out rotatable bond atoms for reference
## comparison of printing vs. old full monte - more options
## comparison of printing vs. GoodVibes
## tabulation of results in addition to sdf output file
## pythonify where possible
## currently RMS but torsion_list would be better?

from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent import futures
import argparse, logging, os, sys, time, copy
import gzip
import pandas as pd

# PHYSICAL CONSTANTS
GAS_CONSTANT, PLANCK_CONSTANT, BOLTZMANN_CONSTANT, SPEED_OF_LIGHT, AVOGADRO_CONSTANT, AMU_to_KG, atmos = 8.3144621, 6.62606957e-34, 1.3806488e-23, 2.99792458e10, 6.0221415e23, 1.66053886E-27, 101.325
# UNIT CONVERSION
j_to_au = 4.184 * 627.509541 * 1000.0

# version number
__version__ = "1.0.1"

# Formatted output to command line and log file

# algorithm to generate nc conformations
def _genConf(s, nc, rms, efilter, rmspost, return_dict, name):
    m = Chem.MolFromSmiles(s)
    if not m:
        return
    try:
        AllChem.EmbedMolecule(m, AllChem.ETKDG())
        m = Chem.AddHs(m, addCoords=True)
    except:
        return

    nr = int(AllChem.CalcNumRotatableBonds(m))
    #m = Chem.AddHs(m)
    #Chem.AssignAtomChiralTagsFromStructure(m, replaceExistingTags=True)
    if not nc: nc = 3**nr

    if not rms: rms = -1
    ids=AllChem.EmbedMultipleConfs(m, numConfs=200, pruneRmsThresh=rms, randomSeed=1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)

    if len(ids)== 0:
        ids = m.AddConformer(m.GetConformer(), assignID=True)

    diz = []
    diz2 = []
    diz3 = []

    try:
        for id in ids:
            prop = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
            ff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=id)
            ff.Minimize()
            en = float(ff.CalcEnergy())
            econf = (en, id)
            diz.append(econf)
    except:
        return_dict['return'] = (None, None, None)
        return
    
    if efilter != "Y":
        n, diz2 = energy_filter(m, diz, efilter)
    else:
        n = m
        diz2 = diz

    if rmspost != None and n.GetNumConformers() > 1:
        o, diz3 = postrmsd(n, diz2, rmspost)
    else:
        o = n
        diz3 = diz2

    return_dict['return'] = (o, diz3, nr)


# wrap the genConf in process so that the genConf can be stopped
class genConf:
    def __init__(self, m, args):
        chembl_id, SMILES = m
        self.s = SMILES
        self.name = chembl_id
        self.nc = args.nconf
        self.rms = args.rmspre
        self.efilter = args.E_cutoff
        self.rmspost = args.rmspost
        self.timeout = args.timeout
        
    def __call__(self):
        self.return_dict = Manager().dict()
        self.process = Process(target=_genConf, args=(self.s, self.nc, self.rms, self.efilter, self.rmspost, self.return_dict, self.name))

        self.process.start()
        self.process.join(self.timeout)
        if 'return' in self.return_dict:
            return self.return_dict['return']
        else:
            self.terminate()
            return (None, None, None)

    def terminate(self):
        self.process.terminate()

# filter conformers based on relative energy
def energy_filter(m, diz, efilter):
    diz.sort()
    mini = float(diz[0][0])
    sup = mini + efilter
    n = Chem.Mol(m)
    n.RemoveAllConformers()
    n.AddConformer(m.GetConformer(int(diz[0][1])))
    nid = []
    ener = []
    nid.append(int(diz[0][1]))
    ener.append(float(diz[0][0])-mini)
    del diz[0]
    for x,y in diz:
        if x <= sup:
            n.AddConformer(m.GetConformer(int(y)))
            nid.append(int(y))
            ener.append(float(x-mini))
        else:
            break
    diz2 = list(zip(ener, nid))
    return n, diz2

# filter conformers based on geometric RMS
def postrmsd(n, diz2, rmspost):
    diz2.sort(key=lambda x: x[0])
    o = Chem.Mol(n)
    confidlist = [diz2[0][1]]
    enval = [diz2[0][0]]
    nh = Chem.RemoveHs(n)
    del diz2[0]
    for z,w in diz2:
        confid = int(w)
        p=0
        for conf2id in confidlist:
            rmsd = AllChem.GetBestRMS(nh, nh, prbId=confid, refId=conf2id)
            if rmsd < rmspost:
                p=p+1
                break
        if p == 0:
            confidlist.append(int(confid))
            enval.append(float(z))
    diz3 = list(zip(enval, confidlist))
    return o, diz3


# conformational search / handles parallel threads if more than one structure is defined
def csearch(supp, total, args, logger):
    if not os.path.isdir(args.MMFF_conf_folder):
        os.mkdir(args.MMFF_conf_folder)

    conf_sdfs = []
    with futures.ProcessPoolExecutor(max_workers=args.MMFF_threads) as executor:
        n_tasks = args.MMFF_threads if args.MMFF_threads < total else total
        tasks = [genConf(next(supp), args) for m in range(n_tasks)]
        running_pool = {task.name: executor.submit(task) for task in tasks}

        while True:
            if len(running_pool) == 0:
                break

            for mol_id in list(running_pool):
                future = running_pool[mol_id]
                if future.done():
                    mol, ids, nr = future.result(timeout=0)
                    if mol:
                        lowest_en, lowest_id = ids[0]
                        mol.SetProp('_Name', mol_id)
                        mol.SetProp('ConfId', str(lowest_id))
                        mol.SetProp('ConfEnergies', str(lowest_en) + ' kcal/mol')
                        writer = Chem.SDWriter(os.path.join(args.MMFF_conf_folder, '{}.sdf'.format(mol_id)))
                        writer.write(mol, confId=lowest_id)
                        conformers_found = len(ids)
                        logger.info('conformer searching for {} completed: '
                                    '{} conformers found, keep the lowest one'.format(mol_id, conformers_found))
                        conf_sdfs.append('{}.sdf'.format(mol_id))
                    else:
                        logger.info('conformer searching for {} failed.'.format(mol_id))
                        pass

                    # add new task
                    del(running_pool[mol_id])
                    
                    try:
                        task = genConf(next(supp), args)
                    except StopIteration:
                        # reach end of the supp
                        logger.info('MMFF conformer searching finished')
                        pass
                    else:
                        running_pool[task.name] = executor.submit(task)
    return conf_sdfs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Molecular conformer generator')
    parser.add_argument('-ismiles', required=True, 
                        help='sdf input file')
    parser.add_argument('-osdf', required=True, 
                        help='sdf output file')
    parser.add_argument('-nconf', type=int, required=False, 
                        help='number of conformers')
    parser.add_argument('-rmspre', type=float, required=False, 
                        help='rms threshold pre optimization')
    parser.add_argument('-rmspost', type=float, required=False, default=0.4, 
                        help='rms threshold post minimization')
    parser.add_argument('-cutoff', type=float, required=False, default=10.0, 
                        help='energy window')
    parser.add_argument('-printproperty', action='store_true', default=True, 
                        help='Print molecule properties (energy and rotable bond number)')
    parser.add_argument('-threads', type=int, required=False, default=40, 
                        help='number of threads')
    parser.add_argument('-timeout', required=False, default=600,
                        help = 'time out to kill sub processors')
    args = parser.parse_args()

    # Check that the input structure exists and has the correct format   
    if os.path.exists(args.ismiles):
        inp = args.ismiles; out = args.osdf
        filename = os.path.splitext(inp)[0]

    # Define input and outputs
    df = pd.read_csv('100k.csv.gz')

    supp = (x for x in df[['chembl_id', 'SMILES']].values)    
    csearch(supp, args)


