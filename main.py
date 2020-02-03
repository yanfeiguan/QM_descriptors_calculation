from argparse import ArgumentParser, Namespace
import os
import shutil

import pandas as pd

from utils import create_logger
from genConf import csearch
from xtb_optimization import xtb_optimization
from dftscf import dft_scf

XTB_PATH = '/Users/yanfei/bin/xtb_exe/bin'
G16_PATH = ''

parser = ArgumentParser()
parser.add_argument('--ismiles', type=str, required=False,
                    help='input smiles included in a .csv file')
parser.add_argument('--output', type=str, default='QM_descriptors.pickle',
                    help='output as a .pickle file')
# conformer searching
parser.add_argument('--MMFF_conf_folder', type=str, default='MMFF_conf',
                    help='folder for MMFF searched conformers')
parser.add_argument('--nconf', type=int, default=500,
                    help='number of MMFF conformers')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False,
                        help='rms threshold pre optimization')
parser.add_argument('--rmspost', type=float, required=False, default=0.4,
                    help='rms threshold post MMFF minimization')
parser.add_argument('--E_cutoff', type=float, required=False, default=10.0,
                    help='energy window for MMFF minimization')
parser.add_argument('--MMFF_threads', type=int, required=False, default=40,
                    help='number of process for the MMFF conformer searching')
parser.add_argument('--timeout', required=False, default=600,
                    help='time window for each MMFF conformer searching sub process')
# xtb optimization
parser.add_argument('--xtb_folder', type=str, default='XTB_opt',
                    help='folder for XTB optimization')

# DFT calculation
parser.add_argument('--DFT_folder', type=str, default='DFT',
                    help='folder for DFT calculation')
parser.add_argument('--DFT_theory', type=str, default='b3lyp/def2svp',
                    help='level of theory for the DFT calculation')

args = parser.parse_args()
args.ismiles = '100k.csv'

name = os.path.splitext(args.ismiles)[0]
logger = create_logger(name=name)

df = pd.read_csv(args.ismiles, index_col=0)

# conformer searching
logger.info('starting MMFF conformer searching')
supp = (x for x in df[['id', 'smiles']].values)
conf_sdfs = csearch(supp, len(df), args, logger)

# xtb optimization
logger.info('starting GFN2-XTB structure optimization for the lowest MMFF conformer')
if not os.path.isdir(args.xtb_folder):
    os.mkdir(args.xtb_folder)

opt_sdfs = []
for conf_sdf in conf_sdfs:
    try:
        shutil.copyfile(os.path.join(args.MMFF_conf_folder, conf_sdf),
                        os.path.join(args.xtb_folder, conf_sdf))
        opt_sdf = xtb_optimization(args.xtb_folder, conf_sdf, XTB_PATH, logger)
        opt_sdfs.append(opt_sdf)
    except Exception as e:
        logger.error(f'XTB optimization for {os.path.splitext(conf_sdf)[0]} failed: {e}')
    else:
        logger.info(f'XTB optimization for {os.path.splitext(conf_sdf)[0]} completed. '
                    f'Structure saved in {args.MMFF_conf_folder}.')

# G16 DFT calculation
if not os.path.isdir(args.DFT_folder):
    os.mkdir(args.DFT_folder)

for opt_sdf in opt_sdfs:
    try:
        shutil.copyfile(os.path.join(args.xtb_folder, opt_sdf),
                        os.path.join(args.DFT_folder, opt_sdf))
        dft_log = dft_scf(args.DFT_folder, opt_sdf, G16_PATH, args.DFT_theory, logger)

