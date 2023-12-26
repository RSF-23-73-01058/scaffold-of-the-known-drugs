from rdkit.ML.Cluster.Butina import ClusterData
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import AllChem
import numpy as np


def get_data(lines, generic=True) -> dict:
    table_dict = {k: v for v, k in enumerate(lines[0].strip().split('\t'))}
    data = {}
    for line in lines[1:]:
        parsed_line = line.split('\t')
        scaffold_smiles = parsed_line[table_dict[f'{"G" if generic else "Non-g"}eneric Scaffold SMILES']].strip()
        scaffold_id = parsed_line[table_dict[f'{"G" if generic else "Non-g"}eneric Scaffold ID']].strip()

        if scaffold_smiles not in data:
            data[scaffold_smiles] = {}
            data[scaffold_smiles]['scaffold_id'] = scaffold_id
            data[scaffold_smiles]['frequency'] = 1
        else:
            data[scaffold_smiles]['frequency'] += 1
    return data


def cluster(data, threshold=0.7, distance='tanimoto') -> dict:
    def euclidian_distance(pi, pj):
        dv = np.array(pi) - np.array(pj)
        dv = sum(dv * dv)
        return np.sqrt(dv)

    def tanimoto_distance(pi, pj):
        dv = np.array(pi) + np.array(pj)
        numerator = sum(1 for n in dv if n == 2)
        denominator = sum(1 for n in dv if n > 0)
        return 1 - numerator / denominator

    scaffolds = list(data.keys())

    points = []
    for scaffold in scaffolds:
        mol = MolFromSmiles(scaffold)
        points.append(list(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=1024)))
    if distance == 'tanimoto':
        clusters = ClusterData(points,
                               nPts=len(scaffolds),
                               distThresh=threshold,
                               distFunc=lambda pi, pj: tanimoto_distance(pi=pi, pj=pj))
    else:
        clusters = ClusterData(points,
                               nPts=len(scaffolds),
                               distThresh=threshold,
                               distFunc=lambda pi, pj: euclidian_distance(pi=pi, pj=pj))

    for cluster_id in range(len(clusters)):
        for scaffold in clusters[cluster_id]:
            data[scaffolds[scaffold]]['cluster_id'] = cluster_id + 1
    return data


def scale(data, left=1, right=10) -> dict:
    frequencies = [data[scaffold]['frequency'] for scaffold in data.keys()]
    max_frequency = max(frequencies)
    min_frequency = min(frequencies)

    angular_coefficient = (right - left) / (max_frequency - min_frequency)
    bias = right - angular_coefficient * max_frequency

    for scaffold in data.keys():
        data[scaffold]['score'] = round(data[scaffold]['frequency'] * angular_coefficient + bias, ndigits=8)
    return data


def write_tsv(data, output):
    output.write('Scaffold ID\t'
                 'Scaffold SMILES\t'
                 'Cluster ID\t'
                 'Score\n')
    for scaffold in data.keys():
        output.write(f"{data[scaffold]['scaffold_id']}\t"
                     f"{scaffold}\t"
                     f"{data[scaffold]['cluster_id']}\t"
                     f"{data[scaffold]['score']}\n")


def main():
    # insert input path here
    path = 
    with open(path) as file:
        lines = file.readlines()
    data = get_data(lines, generic=True)
    data = cluster(data, distance='tanimoto', threshold=0.7)
    data = scale(data, left=1, right=10)

    # insert output path here
    output =
    with open(output, 'w') as file:
        write_tsv(data, file)


if __name__ == '__main__':
    main()
