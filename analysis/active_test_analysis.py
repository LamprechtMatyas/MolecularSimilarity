#!/usr/bin/env python3
""""
Takes active molecules and active molecules from test set and compare their fragments.
"""
import argparse
import json
import copy

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit

import inputoutput_utils
import analysis.ecfpfcfpanalysis_utils as utils


def _main():
    configuration = read_configuration()
    active_molecules = utils.read_molecules(configuration["active_file"])
    test_molecules = utils.read_molecules(configuration["test_file"])
    inputoutput_utils.create_parent_directory(configuration["output_directory"] + "/0")
    options = {
        "kekule": False,
        "isomeric": False
    }
    active_fragments = []
    for active_molecule in active_molecules:
        molecule_smiles = active_molecule.strip("\"")
        molecule = Chem.MolFromSmiles(molecule_smiles)
        mol_indexes = extract_neighbourhood_fragments(molecule, 6, options)
        active_fragments.append(mol_indexes)
    active_fragments1 = copy.deepcopy(active_fragments)
    test_fragments = []
    for test_molecule in test_molecules:
        molecule_smiles = test_molecule.strip("\"")
        molecule = Chem.MolFromSmiles(molecule_smiles)
        mol_indexes = extract_neighbourhood_fragments(molecule, 6, options)
        test_fragments.append(mol_indexes)
    test_fragments1 = copy.deepcopy(test_fragments)
    file1 = configuration["output_directory"] + "/difference.txt"
    file2 = configuration["output_directory"] + "/diff_index.txt"
    file3 = configuration["output_directory"] + "/no_same.txt"
    _prepare_files(file1)
    _prepare_files(file2)

    with open(configuration["nbit_file"], "r", encoding="utf-8") as nbit_stream:
        new_line = nbit_stream.read()
        line = json.loads(new_line)
        nbits = line["nbits"]
        for nbit in nbits:
            print(nbit)
            for mol_num, molecule in enumerate(active_fragments):
                for num, fragment in enumerate(molecule):
                    active_fragments[mol_num][num][0] = int(active_fragments1[mol_num][num][0] % nbit)
            for mol_num, molecule in enumerate(test_fragments):
                for num, fragment in enumerate(molecule):
                    test_fragments[mol_num][num][0] = test_fragments1[mol_num][num][0] % nbit

            sorted_active_fragments = []
            for molecule in active_fragments:
                sorted_active_fragments.append(_sort_indexes(molecule))
            sorted_test_fragments = []
            for molecule in test_fragments:
                sorted_test_fragments.append(_sort_indexes(molecule))
            equvilence_class1 = _equivalence_classes(sorted_active_fragments)
            equvilence_class2 = _equivalence_classes(sorted_test_fragments)
            same_equivalence_classes = utils.same_equivalence_class(equvilence_class1, equvilence_class2)
            difference1 = utils.difference_of_equivalence(equvilence_class1, same_equivalence_classes, 0)
            difference2 = utils.difference_of_equivalence(equvilence_class2, same_equivalence_classes, 1)
            with open(file1, "a", encoding="utf-8") as output_stream:
                output_stream.write("nbit:  " + str(nbit) + "\n")
                output_stream.write("difference1: \n")
                for item in difference1:
                    output_stream.write(str(item) + "\n")
                output_stream.write("difference2:\n")
                for item in difference2:
                    output_stream.write(str(item) + "\n")
                output_stream.write("-------------------------------------\n")
            same_content_and_diff_indexes = []
            for item in same_equivalence_classes:
                if item[0] != item[1]:
                    same_content_and_diff_indexes.append(item)
            with open(file2, "a", encoding="utf-8") as output_stream:
                output_stream.write("nbit:  " + str(nbit) + "\n")
                for item in same_content_and_diff_indexes:
                    output_stream.write(str(item) + "\n")
                output_stream.write("-------------------------------------\n")
    act_frag = []
    for molecule_items in active_fragments1:
        for index in molecule_items:
            act_frag.append(index)
    test_frag = []
    for molecule_items in test_fragments1:
        for index in molecule_items:
            test_frag.append(index)
    same_class_equivalence = utils.same_equivalence_class(act_frag, test_frag)
    print(same_class_equivalence)
    difference1 = utils.difference_of_equivalence(act_frag, same_class_equivalence, 0)
    difference2 = utils.difference_of_equivalence(test_frag, same_class_equivalence, 1)
    with open(file3, "w", encoding="utf-8") as output_stream:
        output_stream.write("actives: \n")
        for item in difference1:
            output_stream.write(str(item) + "\n")
        output_stream.write("inactives: \n")
        for item in difference2:
            output_stream.write(str(item) + "\n")

        output_stream.write(str(difference2) + "\n")


def read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="ecfp/fcfp analysis")
    parser.add_argument("-a", type=str, dest="active_file", required=True)
    parser.add_argument("-t", type=str, dest="test_file", required=True)
    parser.add_argument("-n", type=str, dest="nbit_file", help="json file with nbit configurations",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_directory", help="directory where to store"
                                                                      " output file",
                        required=True)
    return vars(parser.parse_args())


# Extract and return circular fragments.
def extract_neighbourhood_fragments(molecule: object, diameter: int, options: dict) -> list:
    output = []
    info = {}
    size = diameter // 2
    mol_index = []
    AllChem.GetMorganFingerprint(molecule, size, bitInfo=info)
    i = 0
    for element in info:
        for item in info[element]:
            # item = [rooted atom, radius]
            # assemble fragments into atom
            env = rdkit.Chem.FindAtomEnvironmentOfRadiusN(
                molecule, item[1], item[0])
            atoms = set()
            for bidx in env:
                atoms.add(molecule.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(molecule.GetBondWithIdx(bidx).GetEndAtomIdx())
            # check if we have some atoms
            if atoms == set():
                atoms.add(item[0])
            try:
                # kekuleSmiles - we may lost some information
                # about aromatic atoms, but if we do not kekulize
                # we can get invalid smiles
                smiles = rdkit.Chem.MolFragmentToSmiles(
                    molecule, atomsToUse=list(atoms), bondsToUse=env,
                    rootedAtAtom=item[0], kekuleSmiles=options["kekule"],
                    isomericSmiles=options["isomeric"])
            except Exception:
                print("Wrong")
            mol_index.append([element, smiles])
            i += 1
            #if i >= 3:
            #    return mol_index
    return mol_index


def _prepare_files(output_file: str):
    with open(output_file, "w", encoding="utf-8"):
        pass


def _equivalence_classes(mol_indexes: list) -> list:
    list_of_equivalence = []
    for molecule in mol_indexes:
        i = 0
        while i < len(molecule):
            index = [molecule[i][0], molecule[i][1]]
            while (i+1 < len(molecule)) and (molecule[i][0] == molecule[i+1][0]) :
                i +=1
                index.append(molecule[i][1])
            i += 1
            if index not in list_of_equivalence:
                list_of_equivalence.append(index)

    return list_of_equivalence


def _sort_indexes(mol_indexes: list) -> list:
    for i in range(len(mol_indexes) - 1):
        for j in range(len(mol_indexes) - i - 1):
            if mol_indexes[j+1][1] < mol_indexes[j][1]:
                tmp = mol_indexes[j+1]
                mol_indexes[j+1] = mol_indexes[j]
                mol_indexes[j] = tmp
    return mol_indexes


if __name__ == "__main__":
    _main()

