import sys,os
sys.path.append('/home/yanfei/bin')

from g16_log import G16Log
from tqdm import tqdm
import pickle

target = sys.argv[1]
dirs = [os.path.join(target, x) for x in os.listdir(target) if os.path.isdir(x)]
logs = []
for dir in dirs:
    logs.extend([os.path.join(dir, x) for x in os.listdir(dir) if '.log' in x])

def ReadLog(log):
    log = G16Log(log)
    if log.termination:
        QMs = {}
    else:
        return

    if log.formal_charge == 0:
        jtype = 'neutral'
    elif log.formal_charge == 1:
        jtype = 'plus1'
    elif log.formal_charge == -1:
        jtype = 'minus1'

    QMs['type'] = jtype
    if jtype == 'neutral':
        QMs['Coords'] = log.Coords
        QMs['Atoms'] = log.AtomsType
        QMs['mulliken_charge'] = log.mulliken_charge
        QMs['mulliken_dipole_moment'] = log.mulliken_dipole_moment
        QMs['hirshfeld_charges'] = log.hirshfeld_charges
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles
        QMs['NPA_Charge'] = log.NPA_Charge
        QMs['electron_configuration'] = log.electron_configuration
        QMs['bond_index_matrix'] = log.bond_index_matrix
        QMs['lone_pairs'] = log.lone_pairs
        QMs['bond_lewis'] = log.bond_lewis
        QMs['bond_non_lewis'] = log.bond_non_lewis
        QMs['bond_lewis_contribution'] = log.bond_lewis_contribution
        QMs['bond_non_lewis_contribution'] = log.bond_non_lewis_contribution
        QMs['SCF'] = log.SCF
        QMs['NMR'] = log.NMR
        QMs['homo'] = log.homo
        QMs['lumo'] = log.lumo
    elif jtype == 'minus1':
        QMs['mulliken_charge'] = log.mulliken_charge
        QMs['mulliken_spin_density'] = log.mulliken_spin_density
        QMs['hirshfeld_charges'] = log.hirshfeld_charges
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles
        QMs['NPA_Charge'] = log.NPA_Charge
        QMs['SCF'] = log.SCF
    elif jtype == 'plus1':
        QMs['mulliken_charge'] = log.mulliken_charge
        QMs['mulliken_spin_density'] = log.mulliken_spin_density
        QMs['hirshfeld_charges'] = log.hirshfeld_charges
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles
        QMs['NPA_Charge'] = log.NPA_Charge
        QMs['SCF'] = log.SCF

    return QMs

QM_descriptors = [ReadLog(x) for x in tqdm(logs, total=len(logs))]

with open('{}.pickle'.format(target), 'wb') as out:
    pickle.dump(QM_descriptors, out)




