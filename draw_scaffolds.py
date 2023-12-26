import os.path
import re
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from cairosvg import svg2png
from rdkit.Chem import rdDepictor


def draw(smiles, scale, scaffold_id, cluster_id) -> str:
    group_id = f"<g id = '{scaffold_id}' class = 'scaffold' cluster = '{cluster_id}'>"
    mol = Chem.MolFromSmiles(smiles)
    color = (1, 1, 1, 1)

    atoms = []
    for a in mol.GetAtoms():
        atoms.append(a.GetIdx())

    bonds = []
    for bond in mol.GetBonds():
        aid1 = atoms[bond.GetBeginAtomIdx()]
        aid2 = atoms[bond.GetEndAtomIdx()]
        bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

    atom_cols = {}
    for i, at in enumerate(atoms):
        atom_cols[at] = color

    ri = mol.GetRingInfo()
    largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
    if largest_ring_size >= 12:
        rdDepictor.SetPreferCoordGen(True)
    else:
        rdDepictor.SetPreferCoordGen(False)

    AllChem.Compute2DCoords(mol)

    sizer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    sizer.drawOptions().scalingFactor = scale * 20  # see https://github.com/rdkit/rdkit/discussions/6562, default is 20
    sizer.drawOptions().padding = 0  # only structure, no empty space
    sizer.drawOptions().bondLineWidth = 1.5 + scale * 0.03
    sizer.drawOptions().maxFontSize = - 1
    sizer.drawOptions().minFontSize = - 1
    sizer.drawOptions().baseFontSize = 0.5 + scale * 0.00000001
    sizer.drawOptions().fillHighlights = True
    sizer.drawOptions().setHighlightColour(color)
    sizer.drawOptions().highlightBondWidthMultiplier = 20
    size = sizer.GetMolSize(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.drawOptions().scalingFactor = scale
    drawer.drawOptions().padding = 0
    drawer.drawOptions().bondLineWidth = 1.5 + scale * 0.03
    drawer.drawOptions().maxFontSize = - 1
    drawer.drawOptions().minFontSize = - 1
    drawer.drawOptions().baseFontSize = 0.5 + scale * 0.00000001
    drawer.drawOptions().fillHighlights = True
    drawer.drawOptions().setHighlightColour(color)
    drawer.drawOptions().highlightBondWidthMultiplier = 20
    drawer.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors=atom_cols)
    drawer.FinishDrawing()
    drawing_text = drawer.GetDrawingText()
    drawing_text = re.sub(r'<rect.*</rect>', group_id, drawing_text)
    drawing_text = re.sub(r'</svg>', r'</g>\n</svg>', drawing_text)
    return drawing_text


def main(png=False):
    # insert path to the working directory here
    working_directory =
    common_path = working_directory + "\\scaffolds_svg"
    # insert path to the scaffold table here
    table = 
    with open(table) as file:
        lines = file.readlines()
    table_dict = {k: v for v, k in enumerate(lines[0].strip().split('\t'))}
    if not os.path.exists(common_path + f'\\png') and png:
        os.mkdir(common_path + f'\\png')
    for line in lines[1:]:
        parsed_line = line.split('\t')
        scaffold_id = parsed_line[table_dict['Scaffold ID']]
        scaffold_smiles = parsed_line[table_dict['Scaffold SMILES']]
        cluster_id = parsed_line[table_dict['Cluster ID']]
        scale = float(parsed_line[table_dict['Score']])

        if not os.path.exists(common_path + f'\\{cluster_id}'):
            os.mkdir(common_path + f'\\{cluster_id}')

        with open(common_path +
                  f'\\{cluster_id}' +
                  f'\\{cluster_id}_'
                  f'{int(round(scale, 3) * 1000)}_'
                  f'{scaffold_id}.svg', 'w') as file:
            svg = draw(scaffold_smiles, scale, scaffold_id, cluster_id)
            file.write(svg)
        if png:
            svg2png(svg, write_to=common_path + f'\\png\\{cluster_id}_'
                                                f'{int(round(scale, 3) * 1000)}_'
                                                f'{scaffold_id}.png')


if __name__ == '__main__':
    main()
