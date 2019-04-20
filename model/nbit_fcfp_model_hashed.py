#!/usr/bin/env python3
""""
Model that uses AllChem.GetHashedMorganFingerprint() function for fcfp fingerprintns
which hashes the data and we can set the store space as an input parameter as nbit.
Then it uses Tanimoto similarity metrics.
    {"model_name": "nbit_fcfp_model_hashed", "fragments": "fcfp.num1", "nbits": num}
    where num in natural number, num1 is a diameter of ECFP
"""

import json

from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import Chem

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils
import rdkitmodel_utils


class NbitFcfpModelHashed(IModel):
    model_name = "nbit_fcfp_model_hashed"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_smiles = rdkitmodel_utils.select_molecules(active_fragments)
        model = {
            "configuration": model_configuration,
            "data": {
                "active": active_smiles
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        model_data = model_configuration["data"]
        diameter = int(model_configuration["configuration"]["fragments"][0]["size"]) 
        if int(diameter) % 2 == 1:
            print("Incorrect input, size must be even!")
            exit(1)
        radius = diameter // 2
        active_molecules_fcfp = []
        nbits = model_configuration["configuration"]["nbits"]
        for active_molecule in model_data["active"]:
            molecule_smiles = active_molecule.strip("\"")
            molecule = Chem.MolFromSmiles(molecule_smiles)
            fcfp_fingerprint = AllChem.GetHashedMorganFingerprint(molecule, radius,
                                                                  useFeatures=True, nBits=nbits)
            active_molecules_fcfp.append(fcfp_fingerprint)

        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for new_line in input_stream:
                    line = json.loads(new_line)
                    test_molecule_input = line["smiles"]
                    test_molecule_smiles = test_molecule_input.strip("\"")
                    test_molecule = Chem.MolFromSmiles(test_molecule_smiles)
                    test_mol_fingerprint = AllChem.GetHashedMorganFingerprint(test_molecule, radius,
                                                                              useFeatures=True,
                                                                              nBits=nbits)
                    max_sim = max([DataStructs.TanimotoSimilarity(test_mol_fingerprint, fingerprint)
                                   for fingerprint in active_molecules_fcfp])
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


register_model(NbitFcfpModelHashed.model_name, lambda: NbitFcfpModelHashed())
