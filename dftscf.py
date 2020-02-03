from rdkit import Chem
import os
import shutil
import subprocess
from file_parser import mol2xyz, xyz2com


def dft_scf(folder, sdf, g16_path, level_of_theory, n_procs, logger):
    basename = os.path.basename(sdf)

    parent_folder = os.getcwd()
    os.chdir(folder)

    file_name = os.path.splitext(basename)[0]

    xyz = mol2xyz(Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0])

    pwd = os.getcwd()

    g16_command = os.path.join(g16_path, 'g16')
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

        comfile = os.path.join(jobtype, file_name + '.gjf')
        logfile = os.path.join(jobtype, file_name + '.log')

        xyz2com(xyz, head=head, comfile=comfile, charge=charge, mult=mult, footer='$NBO BNDIDX $END\n')
        os.chdir(jobtype)

        with open(logfile, 'w') as out:
            subprocess.call([g16_command, '<', comfile], stdout=out, stderr=out)

    os.remove(sdf)

    os.chdir(parent_folder)
