from rdkit import Chem
import os,subprocess, sys, shutil
from tqdm import tqdm

target = sys.argv[1]

pwd = os.getcwd()

sdf = target + '.sdf'

mols = Chem.SDMolSupplier(sdf, removeHs=False)

mol_name = sdf.split('.')[0]
os.chdir(target)

for conf in tqdm(mols, total=len(mols)):
    conf_name = conf.GetProp('_Name') + '_' + conf.GetProp('ConfId')

    writer = Chem.SDWriter('{}.sdf'.format(conf_name))
    writer.write(conf)
    writer.close()

    try:
        with open(conf_name+'.log', 'w') as out:
            subprocess.call(['/data/yanfei/software/xtb_exe/bin/xtb', '{}.sdf'.format(conf_name), '-opt'], stdout=out, stderr=out)
    except:
        continue

    shutil.move('xtbopt.sdf', '{}_opt.sdf'.format(conf_name))

os.chdir(pwd)
