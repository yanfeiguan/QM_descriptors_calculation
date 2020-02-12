from rdkit import Chem
import os
import shutil
import subprocess

import numpy as np

from .g16_log import XtbLog


def xtb_optimization(folder, sdf, xtb_path, logger):
    basename = os.path.basename(sdf)

    pwd = os.getcwd()

    os.chdir(folder)

    try:
        file_name = os.path.splitext(basename)[0]

        xtb_command = os.path.join(xtb_path, 'xtb')
        with open('{}_xtb_opt.log'.format(file_name), 'w') as out:
            print(xtb_command, '{}.sdf'.format(file_name))
            subprocess.call([xtb_command, '{}.sdf'.format(file_name), '-opt'],
                            stdout=out, stderr=out)
            shutil.move('xtbopt.sdf', '{}_opt.sdf'.format(file_name))
            os.remove('{}.sdf'.format(file_name))

        with open(file_name + '_freq.log', 'w') as out:
            subprocess.call([xtb_command, '{}_opt.sdf'.format(file_name), '-ohess'], stdout=out,
                            stderr=out)

            os.remove('hessian')
            os.remove('vibspectrum')

        log = XtbLog('{}_freq.log'.format(file_name))
    finally:
        os.chdir(pwd)

    if log.termination:
        peaks = log.wavenum
        if np.min(peaks) < 0:
            raise RuntimeError('imaginary frequency found for {}'.format(file_name))
        else:
            return '{}_opt.sdf'.format(file_name)
    else:
        raise RuntimeError('xtb optimization did not finish for {}'.format(file_name))
