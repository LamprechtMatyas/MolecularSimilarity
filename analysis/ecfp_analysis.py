#!/usr/bin/env python3
""""
Takes 2 nbits as input arguments and input file with molecules that we want to compare.
We use ecfp fingerprints as bit vector then we print the differences of representation of
each molecule to the output file.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit

import inputoutput_utils
import analysis.ecfpfcfpanalysis_utils as utils


def _main():
    configuration = utils.read_configuration()
    molecules = utils.read_molecules(configuration["input_file"])
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
            equivalence_class1 = utils.equivalence_class(fragments1)
            equivalence_class2 = utils.equivalence_class(fragments2)
            same_equivalence = utils.same_equivalence_class(equivalence_class1, equivalence_class2)
            difference1 = utils.difference_of_equivalence(equivalence_class1, same_equivalence, 0)
            difference2 = utils.difference_of_equivalence(equivalence_class2, same_equivalence, 1)
            output_stream.write("molecule:  " + molecule_smiles + "\n")
            output_stream.write(str(difference1) + "\n")
            output_stream.write(str(difference2) + "\n")
            output_stream.write(50*"-" + "\n")


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

