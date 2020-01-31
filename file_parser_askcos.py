from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

def mol2xyz(mol, comment=None):
    c = mol.GetConformers()[0]
    coords = c.GetPositions()
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]

    xyz = '{}\n{}\n'.format(len(atoms), comment)
    for a, c in zip(atoms, coords):
        xyz += '{0}     {1:14.9f}    {2:14.9f}    {3:14.9f}\n'.format(a, *c)

    return xyz

def xyz2mol(xyz, smiles):
    lines = xyz.splitlines()
    N_atoms = int(lines[0])
    comments = lines[1]

    if N_atoms != len(lines[2:]):
        raise ValueError('Number of atoms does not match')

    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol = Chem.AddHs(mol, addCoords=True)
    try:
        conf = mol.GetConformers()[0]
    except:
        AllChem.EmbedMultipleConfs(mol, numConfs=1, pruneRmsThresh=0.5, randomSeed=1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        try:
            conf = mol.GetConformers()[0]
        except:
            return None, None

    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for i,coord in enumerate(lines[2:]):
        coord = coord.split()
        
        if atoms[i] != coord[0]:
            raise ValueError('Atom does not match')

        conf.SetAtomPosition(i, np.array(coord[1:]).astype('float'))

    mol.SetProp('comments', comments)
    return mol, comments    

def xyz2com(xyz, head, footer, comfile, charge=0, mult=1):
    title = xyz.splitlines()[1]
    coords = [x+'\n' for x in xyz.splitlines()[2:]]

    with open(comfile, 'w') as com:
        com.write(head)
        com.write('\n')
        com.write(title+'\n')
        com.write('\n')
        com.write('{} {}\n'.format(charge, mult))
        com.writelines(coords)
        com.write('\n')
        com.write(footer)
        com.write('\n\n\n')

def NBO2csv(out, acsv):
    #TODO bcsv
    with open(out, 'r') as outhandle:
        txt = outhandle.readlines()

    txt = [x.strip() for x in txt]
    df,txt = _GetNPACharge(txt)

    return df    
             

def _GetNPACharge(txt):
    columns = 'Atom No    Charge        Core      Valence    Rydberg      Total'
    start_id = txt.index(columns)+2
    end_id = start_id + txt[start_id:].index('====================================================================')

    NPACharge = txt[start_id:end_id] 

    NPACharge = [x.split() for x in NPACharge]

    df = pd.DataFrame(NPACharge, columns=columns.split())
    
    return df, txt[end_id:]











