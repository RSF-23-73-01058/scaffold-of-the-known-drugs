import json


def get_data(indications) -> list:
    indications_table_dict = {k: v for v, k in enumerate(indications[0].strip().split('\t'))}

    data = []
    for line in indications[1:]:
        parse_line = line.strip().split('\t')

        drug_id = parse_line[indications_table_dict['Parent Compound']]
        name = parse_line[indications_table_dict['Preferred Name']]
        smiles = parse_line[indications_table_dict['SMILES']]
        scaffold_id = parse_line[indications_table_dict['Generic Scaffold ID']]
        scaffold_smiles = parse_line[indications_table_dict['Generic Scaffold SMILES']]
        mesh_indications = parse_line[indications_table_dict['MESH Indications']]

        mesh_indications = ' | '.join(mesh_indications.split('|'))

        data.append([drug_id, name, smiles, scaffold_id, scaffold_smiles, mesh_indications])
    return data


def main():
    # insert path to the drug table here
    indications_path = 
    indications = open(indications_path).readlines()

    data = get_data(indications)
    # insert path to the output file here
    output =
    with open(output, 'w') as file:
        json.dump(data, file)


if __name__ == '__main__':
    main()
