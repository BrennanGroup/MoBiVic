"""
Script to functionalise the molecules in vehicle with substituents.

Written by Matthew Holland on 20 September 2024
"""
import json
import csv
from rdkit import Chem
from hcie.molecule import Molecule
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

with open('../Data/vehicle_by_regid.json', 'r') as json_file:
    vehicle_by_regid = json.load(json_file)


vehicle_smiles = [item['smiles'] for item in vehicle_by_regid.values()]
with open('monofunctionalised_smiles.csv', 'r') as csv_file:
    reader = csv.reader(csv_file)
    monofunctionalised_vehicle = list(reader)
    monofunctionalised_vehicle = [entry[0] for entry in monofunctionalised_vehicle]


carbon_substituents = [
    'C',
    'Cl',
    'OC',
    'F',
    'O',
    'N',
    'C(F)(F)(F)'  # CF3 needs to be written explicitly like this so that the carbon atom has the lowest atom idx.
]

nitrogen_substituents = [
    'C',
    'CC(F)(F)(F)'
]

# It is necessary to create RDKit mol objects for each substituent
carbon_substituent_mols = [Chem.MolFromSmiles(smiles) for smiles in carbon_substituents]
nitrogen_substituent_mols = [Chem.MolFromSmiles(smiles) for smiles in nitrogen_substituents]


def functionalise_heterocycle(smiles: str):
    """
       Generate all functionalised derivatives of a given heterocycle.

       For a given input SMILES string, this function identifies all functionalisable positions
       (exit vectors) on the molecule—defined as either C–H bonds on aromatic carbons or N–H bonds
       on pyrrole-like nitrogens. It then attaches each substituent from a predefined set to each
       exit vector, one at a time, generating a new molecule for each combination. Substituents are
       selected based on whether the exit vector is a carbon or nitrogen atom.

       Parameters:
           smiles (str): The SMILES string of the input heterocyclic scaffold.

       Returns:
           List[str]: A list of SMILES strings representing the functionalised derivatives.
       """
    substituted_smiles = []

    vehicle_mol = Molecule(smiles)

    for vector in vehicle_mol.exit_vectors:
        if vehicle_mol.mol.GetAtomWithIdx(vector[0]).GetSymbol() == 'N':
            for substituent in nitrogen_substituent_mols:
                substituted = replace_with_substituent(vehicle_mol.mol,
                                                              substituent,
                                                              vector[0])
                substituted_smiles.append(substituted)
        else:
            for substituent in carbon_substituent_mols:
                substituted = replace_with_substituent(vehicle_mol.mol,
                                                       substituent,
                                                       vector[0])
                substituted_smiles.append(substituted)

    return substituted_smiles


def replace_with_substituent(molecule: Chem.Mol,
                             substituent: Chem.Mol,
                             atom_id: int
                             ) -> str:
    """
    Adds a substituent to the molecule, bonding it to the atom at position given by atom_idx.
    The process is a bit convoluted - there may well be a better way of doing this, but it involves removing Hs from the
    parent molecule, then forming a combined molecule with both the parent and subsituent molecule.
    This is then converted to an editable molecule, and a bond created between the atom id and the substituent.

    Note that the first specified atom in the substituent SMILES string is the one that is bonded to the atom with id
    atom_id

    This new molecule is then sanitised and returned as a SMILES string.
    :param molecule: rdkit mol of parent molecule to functionalise
    :param substituent: The substituent to bond to the parent molecule
    :param atom_id: The atom id of the parent molecule to create the bond to the substituent from.
    :return: SMILES string of functionalised molecule.
    """
    mol = Chem.RemoveHs(molecule)

    if mol.GetAtomWithIdx(atom_id).GetSymbol() == 'N':
        mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)

    combined_mols = Chem.CombineMols(mol, substituent)

    editable_mol = Chem.RWMol(combined_mols)
    editable_mol.AddBond(atom_id, mol.GetNumAtoms(), Chem.BondType.SINGLE)

    substituted_mol = editable_mol.GetMol()
    Chem.SanitizeMol(substituted_mol)

    return Chem.MolToSmiles(substituted_mol)



def functionalise_parallel(smiles_list = vehicle_smiles):
    """
    Parallelise the functionalising of heterocycles for efficiency.
    :param smiles_list: list of all the library SMILES to functionalise
    :return: list of all the functionalised SMILES
    """
    functionalised_smiles = multiprocessing.Manager().list()

    with ProcessPoolExecutor() as executor:
        futures = [
        executor.submit(functionalise_heterocycle, smiles) for smiles in smiles_list
            ]

        for future in as_completed(futures):
            result = future.result()
            for smiles in result:
                functionalised_smiles.append(smiles)

    return list(functionalised_smiles)


def mono_main():
    """
    Monofunctionalise an unfunctionalised library
    :return: .csv file of monofunctionalised SMILES
    """
    functionalised_smiles = functionalise_parallel()
    with open('monofunctionalised_smiles.csv', 'w') as csv_file:
        for smiles in functionalised_smiles:
            csv_file.write(smiles+'\n')

def bi_main():
    """
    Bifunctionalise a monofunctionalised library
    :return: .csv file of bifunctionalised SMILES
    """
    bifunctionalised_smiles = functionalise_parallel(monofunctionalised_vehicle)
    print(len(bifunctionalised_smiles))
    with open('bifunctionalised_smiles.csv', 'w') as csv_file:
        for smiles in bifunctionalised_smiles:
            csv_file.write(smiles+'\n')
