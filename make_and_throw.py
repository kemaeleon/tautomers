"""import libraries"""
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops 
from rdkit.Chem import rdchem  
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
import argparse
from PIL import Image, ImageFont, ImageDraw, ImageEnhance

def _initialise_neutralisation_reactions():
    """from rdkit documentation"""
    patts = (
        # Imidazoles
        ('[n+;H]', 'n'),
        # Amines
        ('[N+;!H0]', 'N'), ('[N-;H1]', 'N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]', 'O'),
        # Ketones
        ('c=[O+;H1]', 'c=O'), ('C=[O+;H1]', 'C=O'),
        # Thiols
        ('[S-;X1]', 'S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]', 'N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]', 'N'),
        # Tetrazoles
        ('[n-]', '[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]', 'S'),
        # Amides
        ('[$([N-]C=O)]', 'N'),
        )
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]


def has_cn_cation(mol):
    patts = ['[c+]', '[C+]', '[c-]', '[C-]', '[n+;H0]', '[N+;H0]']
    #['O+;H2]'
    for p in patts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(p)):
            return True
    return False


_reactions = None


def neutralise_charges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions = _initialise_neutralisation_reactions()
        reactions = _reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return Chem.MolToSmiles(mol, True), True
    else:
        return smiles, False


def save_image_and_label(image, label, out_name):
    img = image
    cropped = img.crop((0, -34, img.width, img.height))
    draw = ImageDraw.Draw(cropped)
    draw.rectangle((0, -34, img.width, 34), fill="white")
    draw.text((100, 20), label, font=ImageFont.truetype("arial.ttf", 14),
              fill='#000000')
    canvas = Image.new('RGBA', cropped.size, (255, 255, 255, 255))
    canvas.paste(cropped, mask=cropped)
    canvas.save(out_name)



def smi_to_taut(smi, name, make_more=True):
    m = Chem.MolFromSmiles(smi)
    suppl = rdchem.ResonanceMolSupplier(m,
                                        rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION |
                                        rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS |
                                        rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS |
                                        rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS)
    tautomers = []
    count = 0
    for g in suppl:
        try:
            s = Chem.MolToSmiles(g)
            (s, neutralised) = neutralise_charges(s)
            z = Chem.MolFromSmiles(s)
            z1 = Chem.RemoveHs(z)
            Chem.Kekulize(z1)
            s = Chem.MolToSmiles(z1,isomericSmiles = True)
            g = Chem.MolFromSmiles(s)
            if not has_cn_cation(z1):
                try:
                    if make_more:
                        smi_to_taut(s, name+"_g2_"+str(count), make_more=False)
                        count += 1
                except:                
                    print("resonance structure generation from res struct failed")
                tautomers.append(g)
        except:
            print(str("{} could not be processed.").format(s))
    w = Chem.SDWriter(str(name) + ".sdf")
    fh = open(str(name) + ".smi", 'w')
    for t in tautomers:
        th = Chem.AddHs(t)
        AllChem.EmbedMolecule(th, AllChem.ETKDG())
        w.write(th)
        fh.write(AllChem.MolToSmiles(th) + "\n")
    mcs = rdFMCS.FindMCS(tautomers, atomCompare=rdFMCS.AtomCompare.CompareElements)
    mcsp = Chem.MolFromSmarts(mcs.smartsString)
    AllChem.EmbedMolecule(mcsp, useRandomCoords=True)
    #print(Chem.MolToMolBlock(mcsp))
    #print(mcsp)
    #for t in tautomers:
    #    AllChem.GenerateDepictionMatching2DStructure(t, mcsp)
    #print(mcs.smartsString)
    img = Draw.MolsToGridImage(tautomers, molsPerRow=3, legends=None)
    save_image_and_label(img, smi, str(name) + ".png")
    return tautomers

def main():
    """example: python --smiles C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O --name vitC"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles', dest='smiles', type=str,
                        default='Nc1[n-]c([O-])c2[nH+]cnc2[nH+]1',
                        help='smiles input molecule')
    parser.add_argument('--name', dest='name', type=str,
                        default='output', help='what to call the output files')
    args = parser.parse_args()
    smi_to_taut(args.smiles, args.name)


if __name__ == '__main__':
    main()
