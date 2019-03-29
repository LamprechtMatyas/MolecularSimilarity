#!/usr/bin/env python3
""""
Takes 2 nbits as input argument and input file with molecules that we want to compare.
We use topological torsion fingerprints that are hashed and then we print the comparison of
each molecule to the output file.
"""
import argparse

from rdkit import Chem
from rdkit.Chem.AtomPairs import Torsions

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    molecules = _read_molecules(configuration["input_file"])
    inputoutput_utils.create_parent_directory(configuration["output_file"])
    with open(configuration["output_file"], "w", encoding="utf-8") as output_stream:
        for active_molecule in molecules:
            molecule_smiles = active_molecule.strip("\"")
            molecule = Chem.MolFromSmiles(molecule_smiles)
            torsion1 = Torsions.GetHashedTopologicalTorsionFingerprint(molecule,
                                                                       nBits=int(configuration["first_nbit"]))
            d = torsion1.GetNonzeroElements()
            torsion2 = Torsions.GetHashedTopologicalTorsionFingerprint(molecule,
                                                                       nBits=int(configuration["second_nbit"]))
            d2 = torsion2.GetNonzeroElements()
            path_score1 = []
            for item in d:
                path_score1.append(Torsions.ExplainPathScore(item))
            path_score2 = []
            for item in d2:
                path_score2.append(Torsions.ExplainPathScore(item))

            output_stream.write(molecule_smiles + "\n\n")
            output_stream.write(
                configuration["first_nbit"] + 69*" " + configuration["second_nbit"] + "\n")
            for i in range(min(len(path_score1), len(path_score2))):
                output_stream.write(str(path_score1[i]) + "     " + str(path_score2[i]) + "\n")
            if len(path_score2) > len(path_score1):
                for i in range(len(path_score1), len(path_score2)):
                    output_stream.write(44*" " + str(path_score2[i]) + "\n")
            if len(path_score2) < len(path_score1):
                for i in range(len(path_score2), len(path_score1)):
                    output_stream.write(str(path_score1[i]) + "\n")
            output_stream.write(69*"-" + "\n\n")


def _read_molecules(file: str) -> list:
    molecules = []
    with open(file, mode="r", encoding="utf-8") as input_sream:
        for new_line in input_sream:
            atributes = new_line.split("\t")
            molecules.append(atributes[0])
    return molecules


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="tt analysis")
    parser.add_argument("-f", type=str, dest="first_nbit", help="write number",
                        required=True)
    parser.add_argument("-s", type=str, dest="second_nbit", help="write number",
                        required=True)
    parser.add_argument("-i", type=str, dest="input_file", help="input smi file with molecules",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file", help="output txt file", required=True)
    return vars(parser.parse_args())


if __name__ == "__main__":
    _main()

