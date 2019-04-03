#!/usr/bin/env python3
""""
Takes 1 nbit as input argument and input file with molecules. It takes all molecules from input
file together and computes equivalence classes. Each equivalence class in represented only once.
It writes all equivalence classes that have got more than 2 items.
"""
import argparse
import json

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit

import inputoutput_utils
import analysis.ecfpfcfpanalysis_utils as utils


def _main():
    configuration = _read_configuration()
    molecules = utils.read_molecules(configuration["input_file"])
    inputoutput_utils.create_parent_directory(configuration["output_file"])
    options = {
        "kekule": False,
        "isomeric": False
    }
    equivalence_classes = []

    for active_molecule in molecules:
        molecule_smiles = active_molecule.strip("\"")
        molecule = Chem.MolFromSmiles(molecule_smiles)
        fragments1 = extract_neighbourhood_fragments(molecule, 6, options,
                                                     int(configuration["first_nbit"]))
        equivalence_class1 = utils.equivalence_class(fragments1)
        if not equivalence_classes:
            equivalence_classes.extend(equivalence_class1)
        else:
            for new_item in equivalence_class1:
                if new_item not in equivalence_classes:
                    equivalence_classes.append(new_item)
    eq_classes = []
    for item in equivalence_classes:
        if len(item) > 2:
            eq_classes.append(item[1:])
    with open(configuration["output_file"], "w", encoding="utf-8") as output_stream:
        json.dump(eq_classes, output_stream)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="prints all ecfp equivalence classes ")
    parser.add_argument("-f", type=str, dest="first_nbit", help="write number",
                        required=True)
    parser.add_argument("-i", type=str, dest="input_file", help="input smi file with molecules",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file", help="output json file", required=True)
    return vars(parser.parse_args())


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


if __name__ == "__main__":
    _main()

