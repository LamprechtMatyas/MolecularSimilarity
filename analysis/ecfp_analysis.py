#!/usr/bin/env python3
""""
Takes 2 nbits as input arguments and input file with molecules that we want to compare.
We use ecfp fingerprints as bit vector then we print the differences of representation of
each molecule to the output file.
"""
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    molecules = _read_molecules(configuration["input_file"])
    inputoutput_utils.create_parent_directory(configuration["output_file"])
    options = {
        "kekule": False,
        "isomeric": False
    }
    with open(configuration["output_file"], "w", encoding="utf-8") as output_stream:
        for active_molecule in molecules:
            molecule_smiles = active_molecule.strip("\"")
            molecule = Chem.MolFromSmiles(molecule_smiles)
            fragments1 = extract_neighbourhood_fragments(molecule, 6, options,
                                                         int(configuration["first_nbit"]))
            fragments2 = extract_neighbourhood_fragments(molecule, 6, options,
                                                         int(configuration["second_nbit"]))
            equivalence_class1 = _equivalence_class(fragments1)
            equivalence_class2 = _equivalence_class(fragments2)
            same_equivalence = _same_equivalence_class(equivalence_class1, equivalence_class2)
            difference1 = _difference_of_equivalence(equivalence_class1, same_equivalence, 0)
            difference2 = _difference_of_equivalence(equivalence_class2, same_equivalence, 1)
            output_stream.write("molecule:  " + molecule_smiles + "\n")
            output_stream.write(str(difference1) + "\n")
            output_stream.write(str(difference2) + "\n")
            output_stream.write(50*"-" + "\n")


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="ecfp analysis")
    parser.add_argument("-f", type=str, dest="first_nbit", help="write number",
                        required=True)
    parser.add_argument("-s", type=str, dest="second_nbit", help="write number",
                        required=True)
    parser.add_argument("-i", type=str, dest="input_file", help="input smi file with molecules",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file", help="output txt file",
                        required=True)
    return vars(parser.parse_args())


def _read_molecules(file: str) -> list:
    molecules = []
    with open(file, mode="r", encoding="utf-8") as input_sream:
        for new_line in input_sream:
            atributes = new_line.split("\t")
            molecules.append(atributes[0])
    return molecules


# Extract and return circular fragments.
def extract_neighbourhood_fragments(molecule: object, diameter: int, options: dict, nbits) -> list:
    output = []
    info = {}
    size = diameter // 2
    AllChem.GetMorganFingerprintAsBitVect(molecule, size, nBits=nbits, bitInfo=info)
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
            output.append({
                "smiles": smiles,
                "index": element,
                "type": "ECFP",
                "size": size
            })
    return output


def _equivalence_class(output):
    last = -1
    num = -1
    list_of_equivalence = []
    for i in range(len(output)):
        if last == int(output[i]["index"]):
            list_of_equivalence[num].append(output[i]["smiles"])
        else:
            list_of_equivalence.append([output[i]["index"]])
            num += 1
            list_of_equivalence[num].append(output[i]["smiles"])
            last = output[i]["index"]
    return list_of_equivalence


def _same_equivalence_class(clekv1, clekv2):
    same = []
    num = -1
    for item1 in clekv1:
        for item2 in clekv2:
            if item1[1:] == item2[1:]:
                same.append([item1[0]])
                num += 1
                same[num].append(item2[0])
                same[num].append(item1[1:])
                clekv2.remove(item2)
                break
    return same


def _difference_of_equivalence(ekv1, intersection, position: int):
    difference = []
    for item in ekv1:
        founded = False
        for intersect in intersection:
            if item[0] == intersect[position]:
                founded = True
                break
        if founded is False:
            difference.append(item)
    return difference


if __name__ == "__main__":
    _main()

