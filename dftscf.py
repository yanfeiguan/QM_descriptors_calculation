from rdkit import Chem
import os
import shutil
import subprocess
from file_parser import mol2xyz, xyz2com
from grab_QM_descriptors import read_log


def dft_scf(folder, sdf, g16_path, level_of_theory, n_procs, logger):
    basename = os.path.basename(sdf)

    parent_folder = os.getcwd()
    os.chdir(folder)

    file_name = os.path.splitext(basename)[0]

    xyz = mol2xyz(Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0])

    pwd = os.getcwd()

    g16_command = os.path.join(g16_path, 'g16')
    QM_descriptors = {}
    for jobtype in ['neutral', 'plus1', 'minus1']:
        if not os.path.isdir(jobtype):
            os.mkdir(jobtype)

        if jobtype == 'neutral':
            charge = 0
            mult = 1
            head = '%chk={}.chk\n%nprocshared={}\n# b3lyp/def2svp nmr=GIAO scf=(maxcycle=512, xqc) ' \
                   'pop=(full,mbs,hirshfeld,nbo6read)\n'.format(file_name, n_procs)
        elif jobtype == 'plus1':
            charge = 1
            mult = 2
            head = '%chk={}.chk\n%nprocshared=20\n# b3lyp/def2svp scf=(maxcycle=512, xqc) ' \
                   'pop=(full,mbs,hirshfeld,nbo6read)\n'.format(file_name, n_procs)
        elif jobtype == 'minus1':
            charge = -1
            mult = 2
            head = '%chk={}.chk\n%nprocshared=20\n# b3lyp/def2svp scf=(maxcycle=512, xqc) ' \
                   'pop=(full,mbs,hirshfeld,nbo6read)\n'.format(file_name, n_procs)


        os.chdir(jobtype)
        comfile = file_name + '.gjf'
        xyz2com(xyz, head=head, comfile=comfile, charge=charge, mult=mult, footer='$NBO BNDIDX $END\n')

        logfile = file_name + '.log'
        outfile = file_name + '.out'
        with open(outfile, 'w') as out:
            #subprocess.run('{} < {} >> {}'.format(g16_command, comfile, logfile), shell=True, stdout=out, stderr=out)
            QM_descriptors[jobtype] = read_log(logfile)
        os.chdir(pwd)

    QM_descriptor_return = QM_descriptors['neutral']

    # charges and fukui indices
    for charge in ['mulliken_charge', 'hirshfeld_charges', 'NPA_Charge']:
        QM_descriptor_return['{}_plus1'.format(charge)] = QM_descriptors['plus1'][charge]
        QM_descriptor_return['{}_minus1'.format(charge)] = QM_descriptors['minus1'][charge]

        QM_descriptor_return['{}_fukui_elec'.format(charge)] = QM_descriptors['neutral'][charge] - \
                                                               QM_descriptors['minus1'][charge]
        QM_descriptor_return['{}_fukui_neu'.format(charge)] = QM_descriptors['plus1'][charge] - \
                                                               QM_descriptors['neutral'][charge]

    # spin density
    for spin in ['mulliken_spin_density', 'hirshfeld_spin_density']:
        QM_descriptor_return['{}_plus1'.format(spin)] = QM_descriptors['plus1'][spin]
        QM_descriptor_return['{}_minus1'.format(charge)] = QM_descriptors['minus1'][spin]

    # SCF
    QM_descriptor_return['SCF_plus1'] = QM_descriptors['plus1']['SCF']
    QM_descriptor_return['SCF_minus1'] = QM_descriptors['minus1']['SCF']

    os.remove(sdf)

    os.chdir(parent_folder)

    return QM_descriptor_return
