import os, sys,re
sys.path.append('/home/yanfei/bin')

from file_parser import xyz2com
from tqdm import tqdm

xyzs = [x for x in os.listdir('.') if 'xyz' in x]

types = ['neutral', 'plus1', 'minus1']

for jobtype in types:
    os.mkdir(jobtype)

    for xyz in tqdm(xyzs, total=len(xyzs)):
        with open(xyz, 'r') as handle:
            xyz_string = handle.read()

        name = xyz.split('.')[0]
        if jobtype == 'neutral':
            charge = 0
            mult = 1
            head = '%chk={}.chk\n%nprocshared=20\n# b3lyp/def2svp nmr=GIAO scf=(maxcycle=512, xqc) pop=(full,mbs,hirshfeld,nbo6read)\n'.format(name)
        elif jobtype == 'plus1':
            charge = 1
            mult = 2
            head = '%chk={}.chk\n%nprocshared=20\n# b3lyp/def2svp scf=(maxcycle=512, xqc) pop=(full,mbs,hirshfeld,nbo6read)\n'.format(name)
        elif jobtype == 'minus1':
            charge = -1
            mult = 2
            head = '%chk={}.chk\n%nprocshared=20\n# b3lyp/def2svp scf=(maxcycle=512, xqc) pop=(full,mbs,hirshfeld,nbo6read)\n'.format(name)


        comfile = os.path.join(jobtype, name + '.gjf')
    
        xyz2com(xyz_string, head=head,
                comfile=comfile, charge=charge, mult=mult, footer='$NBO BNDIDX $END\n')


