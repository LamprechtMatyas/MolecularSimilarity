#!/usr/bin/env python3
""""
Extract fragments for molecules in given SDF/SMILES files and save then into JSON files.

Usage:
    python extract_fragments.py
        -i (comma separated 3 input files with inactive molecules, first with active molecules, second with inactive
            molecules, third with test molecules)
        -o {comma separated path to 3 output files, first with fragments from active molecules, second with fragments
            from inactive molecules, third with molecules from test molecules}
        -f {optional, comma separated list of fragment types to extract}
        -p {type of input files, "sdf", "smi". Default is "sdf"}
        --kekule {generated kekule form of SMILES for fragments}
        --isomeric {put stereochemistry information into fragments SMILES}
Fragments type:
    - tt.{SIZE}
    - ecfp.{SIZE}
    - fcfp.{SIZE}
    - ap
where {SIZE} should be replaced by required fragment size. {SIZE} for ecfp and fcfp are only even numbers.
Usage example:
    tt.3,ecfp.2
default value:
    ecfp.6
Kekule smiles form has no aromatic bonds. Use of --kekule option thus may
reduce the number of generated unique fragments.
This file can be also imported as a python script. In such case please
use the extract_fragments method.

"""

import argparse
import logging
import json
from typing import Dict, Any, List, Union, TextIO

import rdkit
import rdkit.Chem
from rdkit.Chem import AllChem
import rdkit.Chem.AtomPairs.Utils
from rdkit.Chem.AtomPairs import Pairs

import inputoutput_utils


def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(module)s - %(message)s",
        datefmt="%H:%M:%S")
    configuration = _read_configuration()

    # Prepare configuration for the extraction.
    extraction_options = {
        "kekule": configuration["kekule"],
        "isomeric": configuration["isomeric"],
        "fragments": configuration["fragments"]
    }
    input_files = configuration["input"]
    extract_fragments(input_files, configuration["input_type"],
                      configuration["output"], extraction_options)


# Get and return application settings.
def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Extract molecular fragments. "
                    "See file header for more details.")
    parser.add_argument("-i", type=str, dest="input",
                        help="3 comma separated input files with actives, inactives and test molecules",
                        required=True)
    parser.add_argument("-o", type=str, dest="output",
                        help="3 comma separated output files", required=True)
    parser.add_argument("-f", type=str, dest="fragments",
                        help="fragment type - ecfp, fcfp, ap, tt", required=False)
    parser.add_argument("-p", type=str, dest="input_type",
                        help="type of input file smi/sdf, default sdf", default="sdf")
    parser.add_argument("--kekule", dest="kekule", help="kekule option",
                        action="store_true", required=False)
    parser.add_argument("--isomeric", dest="isomeric", help="isomeric option",
                        action="store_true", required=False)

    configuration = vars(parser.parse_args())

    if configuration["fragments"] is None:
        configuration["fragments"] = "ecfp.6"
    parsed_types = []
    for item in configuration["fragments"].split(","):
        item_split = item.split(".")
        if item_split[0] != "ap":
            if not len(item_split) == 2:
                logging.error("Invalid fragment type: %s", item)
                logging.info("Expected format {TYPE}.{SIZE} or ap")
                exit(1)
            parsed_types.append({
                "name": item_split[0],
                "size": int(item_split[1])
            })
        else:
            parsed_types.append({
                "name": item_split[0],
            })
    input_files = []
    for file in configuration["input"].split(","):
        input_files.append(file)
    if len(input_files) != 3:
        logging.info("Wrong number of input files")
        exit(1)
    configuration["input"] = input_files
    output_files = []
    for file in configuration["output"].split(","):
        output_files.append(file)
    if len(output_files) != 3:
        logging.info("Wrong number of output files")
        exit(1)
    configuration["output"] = output_files
    configuration["fragments"] = parsed_types
    configuration["input_type"] = configuration["input_type"].lower()
    return configuration


# Extract fragments from molecules and write them to output JSON files.
# The extraction_options["fragments"] must be a list with objects describing
# fragments to extract, see _read_configuration for more details.
def extract_fragments(input_files: list, input_type: str, output_files: list,
                      extraction_options: dict):
    # The write_molecule_json need some static info.

    # Count some statistics.
    total_fragments = 0
    for file in output_files:
        inputoutput_utils.create_parent_directory(file)

    for file_num, path in enumerate(input_files):
        holder = {"first": True}
        with open(output_files[file_num], "w", encoding="utf-8") as output_stream:
            for molecule in _LOAD_FUNCTIONS[input_type](path):
                item = {
                    "name": molecule.GetProp("_Name"),
                    "smiles": rdkit.Chem.MolToSmiles(molecule),
                    "fragments": _extract_fragments_from_molecule(
                        molecule, extraction_options["fragments"],
                        extraction_options)
                }
                total_fragments += len(item["fragments"])
                # Append to output.
                _append_object_to_jsonlines(output_stream, item, holder)
    logging.info("\tfragments total: %d", total_fragments)


# Generate molecules from SDF file.
def _load_sdf(path: str) -> iter:
    logging.info("Loading (SDF): %s" % path)
    for molecule in rdkit.Chem.SDMolSupplier(path):
        if molecule is None:
            logging.error("Invalid molecule detected.")
            continue
        yield molecule


# Generate molecules from SMI file.
def _load_smi(path: str) -> iter:
    logging.info("Loading (SMI): %s" % path)
    with open(path, "r", encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            line_parts = line.split("\t")
            molecule = rdkit.Chem.MolFromSmiles(line_parts[0])
            if molecule is None:
                logging.error("Invalid molecule detected.")
                continue
            # Molecules created from SMILES does not have any name
            molecule.SetProp("_Name", line_parts[1])
            yield molecule


_LOAD_FUNCTIONS = {
    "sdf": _load_sdf,
    "smi": _load_smi
}


# Return fragments for given molecule.
def _extract_fragments_from_molecule(molecule: object, types: list, options: dict) -> list:
    output = []
    for item in types:
        if item["name"] == "tt":
            output.extend(extract_path_fragments(
                molecule, item["size"], options))
        elif item["name"] == "ecfp":
            _control_size(item["size"], item["name"])
            output.extend(extract_neighbourhood_fragments(
                molecule, item["size"], options, True))
        elif item["name"] == "fcfp":
            _control_size(item["size"], item["name"])
            output.extend(extract_neighbourhood_fragments(
                molecule, item["size"], options, False))
        elif item["name"] == "ap":
            output.extend(extract_atompair_fragments(
                molecule))
    return output


def extract_path_fragments(molecule: object, size: int, options: dict) -> list:
    output = []
    pattern = rdkit.Chem.MolFromSmarts("*" + ("~*" * (size - 1)))
    for atoms in molecule.GetSubstructMatches(pattern):
        smiles = rdkit.Chem.MolFragmentToSmiles(
            molecule, atomsToUse=list(atoms),
            kekuleSmiles=options["kekule"],
            isomericSmiles=options["isomeric"])
        output.append({
            "smiles": smiles,
            "index": _score_path(molecule, atoms, size),
            "type": "TT",
            "size": size
        })
    return output


def _score_path(molecule: object, path: list, size: int) -> int:
    codes = [None] * size
    for i in range(size):
        if i == 0 or i == (size - 1):
            sub = 1
        else:
            sub = 2
        codes[i] = _get_atom_code(molecule.GetAtomWithIdx(path[i]), sub)

    # We scan the vector for both sides, we want to make sure that
    # the begging is less or equal to the end.

    # "canonize" the code vector:
    beg = 0
    end = len(codes) - 1
    while beg < end:
        if codes[beg] > codes[end]:
            codes.reverse()
            break
        elif codes[beg] == codes[end]:
            beg += 1
            end -= 1
        else:
            break

    # Just add all together.
    accum = 0
    for i in range(size):
        accum |= (codes[i]) << (ATOM_CODE["bits"]["total"] * i)
    return accum


ATOM_CODE = {
    "bits": {
        "type": 4,
        "pi": 2,
        "branch": 4,
        "total": 10
    }
}


def _get_atom_code(atom: object, branch_subtract: int) -> int:
    # Constants;
    num_type_bits = ATOM_CODE["bits"]["type"]
    num_pi_bits = ATOM_CODE["bits"]["pi"]
    num_branch_bits = ATOM_CODE["bits"]["branch"]

    max_num_branches = (1 << num_branch_bits) - 1
    max_num_pi = (1 << num_pi_bits) - 1

    atom_number_types = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 43,
                         0]
    # Number of non-hydrogen? neighbor
    if atom.GetDegree() > branch_subtract:
        num_branches = atom.GetDegree() - branch_subtract
    else:
        num_branches = 0
    code = num_branches % max_num_branches
    # Number of bonding pi-electrons.
    n_pi = rdkit.Chem.AtomPairs.Utils.NumPiElectrons(atom) % max_num_pi
    code |= n_pi << num_branch_bits

    # If atom.getAtomicNum() is in atomNumberTypes then return
    # exact match. Otherwise return smallest bigger value.
    type_idx = 0
    n_types = 1 << num_type_bits
    while type_idx < n_types:
        if atom_number_types[type_idx] == atom.GetAtomicNum():
            break
        elif atom_number_types[type_idx] > atom.GetAtomicNum():
            type_idx = n_types
            break
        else:
            type_idx += 1

    # Make sure we do not point outside the array.
    if type_idx == n_types:
        type_idx -= 1

    # Atom type.
    code |= type_idx << (num_branch_bits + num_pi_bits)

    return code


def _control_size(size: int, name: str):
    if int(size) % 2 == 1:
        print("Incorrect input, size in ", name, " must be even!")
        exit(1)


# Extract and return circular fragments.
def extract_neighbourhood_fragments(molecule: object, diameter: int, options: dict, ecfp: bool) -> list:
    output = []
    info = {}
    size = diameter // 2
    if ecfp:
        AllChem.GetMorganFingerprint(molecule, radius=size, bitInfo=info)
    else:
        AllChem.GetMorganFingerprint(molecule, radius=size, bitInfo=info, useFeatures=True)
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
                logging.exception("Invalid fragment detected.")
                logging.info("Molecule: %s", molecule.GetProp("_Name"))
                logging.info("Atoms: %s", ",".join([str(x) for x in atoms]))
            if ecfp:
                output.append({
                    "smiles": smiles,
                    "index": element,
                    "type": "ECFP",
                    "size": size
                })
            else:
                output.append({
                    "smiles": smiles,
                    "index": element,
                    "type": "FCFP",
                    "size": size
                })
    return output


def extract_atompair_fragments(molecule: object) -> list:
    output = []
    pairFps = Pairs.GetAtomPairFingerprint(molecule)
    d = pairFps.GetNonzeroElements()
    for pair in d:
        atom1 = rdkit.Chem.AtomFromSmarts(Pairs.ExplainPairScore(pair)[0][0])
        atom2 = rdkit.Chem.AtomFromSmarts(Pairs.ExplainPairScore(pair)[2][0])
        smiles = (Pairs.ExplainPairScore(pair)[0][0] + Pairs.ExplainPairScore(pair)[2][0])
        atom1_type = atom1.GetAtomicNum()
        atom2_type = atom2.GetAtomicNum()
        atom1_num_pi_bonds = Pairs.ExplainPairScore(pair)[0][2]
        atom2_num_pi_bonds = Pairs.ExplainPairScore(pair)[2][2]
        atom1_num_neigh = Pairs.ExplainPairScore(pair)[0][1]
        atom2_num_neigh = Pairs.ExplainPairScore(pair)[2][1]
        atom1_property_value = 64 * atom1_type + 16 * atom1_num_pi_bonds + atom1_num_neigh
        atom2_property_value = 64 * atom2_type + 16 * atom2_num_pi_bonds + atom2_num_neigh
        dist = Pairs.ExplainPairScore(pair)[1] + 1
        atom_pair_key = min(atom1_property_value, atom2_property_value) + 1024 * (
                max(atom1_property_value, atom2_property_value) + 1024 * dist
        )
        num = (d[pair])
        for i in range(num):
            output.append({
                "smiles": smiles,
                "index": atom_pair_key,
                "type": "AP",
                "size": dist
            })
    return output


# Write given molecule as a JSON into stream.
def _append_object_to_jsonlines(output_stream: TextIO, item: dict, holder: Dict[str, bool]):
    if holder["first"]:
        holder["first"] = False
    else:
        output_stream.write("\n")
    json.dump(item, output_stream)


if __name__ == "__main__":
    _main()


