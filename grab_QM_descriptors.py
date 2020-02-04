from g16_log import G16Log
import numpy as np


def read_log(log):
    log = G16Log(log)
    if log.termination:
        QMs = {}
    else:
        return None

    if log.formal_charge == 0:
        jtype = 'neutral'
    elif log.formal_charge == 1:
        jtype = 'plus1'
    elif log.formal_charge == -1:
        jtype = 'minus1'

    QMs['type'] = jtype
    if jtype == 'neutral':
        QMs['Coords'] = log.Coords.astype(np.float) if log.Coords is not None else np.nan
        QMs['Atoms'] = log.AtomsType
        QMs['mulliken_charge'] = log.mulliken_charge.astype(np.float) if log.mulliken_charge is not None else np.nan
        QMs['mulliken_dipole_moment'] = log.mulliken_dipole_moment.astype(np.float) if log.mulliken_dipole_moment is not None else np.nan
        QMs['hirshfeld_charges'] = log.hirshfeld_charges.astype(np.float) if log.hirshfeld_charges is not None else np.nan
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density.astype(np.float) if log.hirshfeld_spin_density is not None else np.nan
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles.astype(np.float) if log.hirshfeld_dipoles is not None else np.nan
        QMs['NPA_Charge'] = log.NPA_Charge.astype(np.float) if log.NPA_Charge is not None else np.nan
        QMs['electron_configuration'] = log.electron_configuration.astype(np.float) if log.electron_configuration is not None else np.nan
        QMs['bond_index_matrix'] = log.bond_index_matrix.astype(np.float) if log.bond_index_matrix is not None else np.nan
        QMs['lone_pairs'] = log.lone_pairs.astype(np.float) if log.lone_pairs is not None else np.nan
        QMs['bond_lewis'] = log.bond_lewis.astype(np.float) if log.bond_lewis is not None else np.nan
        QMs['bond_non_lewis'] = log.bond_non_lewis.astype(np.float) if log.bond_non_lewis is not None else np.nan
        QMs['bond_lewis_contribution'] = log.bond_lewis_contribution.astype(np.float) if log.bond_lewis_contribution is not None else np.nan
        QMs['bond_non_lewis_contribution'] = log.bond_non_lewis_contribution.astype(np.float) if log.bond_non_lewis_contribution is not None else np.nan
        QMs['SCF'] = log.SCF if log.SCF is not None else np.nan
        QMs['NMR'] = log.NMR if log.NMR is not None else np.nan
        QMs['homo'] = log.homo if log.homo is not None else np.nan
        QMs['lumo'] = log.lumo if log.lumo is not None else np.nan
    elif jtype == 'minus1':
        QMs['mulliken_charge'] = log.mulliken_charge.astype(np.float) if log.mulliken_charge is not None else np.nan
        QMs['mulliken_spin_density'] = log.mulliken_spin_density.astype(np.float) if log.mulliken_spin_density is not None else np.nan
        QMs['hirshfeld_charges'] = log.hirshfeld_charges.astype(np.float) if log.hirshfeld_charges is not None else np.nan
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density.astype(np.float) if log.hirshfeld_spin_density is not None else np.nan
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles.astype(np.float) if log.hirshfeld_dipoles is not None else np.nan
        QMs['NPA_Charge'] = log.NPA_Charge.astype(np.float) if log.NPA_Charge is not None else np.nan
        QMs['SCF'] = log.SCF if log.SCF is not None else np.nan
    elif jtype == 'plus1':
        QMs['mulliken_charge'] = log.mulliken_charge.astype(np.float) if log.mulliken_charge is not None else np.nan
        QMs['mulliken_spin_density'] = log.mulliken_spin_density.astype(np.float) if log.mulliken_spin_density is not None else np.nan
        QMs['hirshfeld_charges'] = log.hirshfeld_charges.astype(np.float) if log.hirshfeld_charges is not None else np.nan
        QMs['hirshfeld_spin_density'] = log.hirshfeld_spin_density.astype(np.float) if log.hirshfeld_spin_density is not None else np.nan
        QMs['hirshfeld_dipoles'] = log.hirshfeld_dipoles.astype(np.float) if log.hirshfeld_dipoles is not None else np.nan
        QMs['NPA_Charge'] = log.NPA_Charge.astype(np.float) if log.NPA_Charge is not None else np.nan
        QMs['SCF'] = log.SCF if log.SCF is not None else np.nan

    return QMs




