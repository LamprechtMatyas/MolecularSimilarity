"""
compute_descriptors without --fragments
"""

import json
import rdkit

from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import Chem
import pickle

from model_interface import IModel
from model_factory import register_model
import my_library


class RdkitEcfpModel(IModel):
    model_name = "rdkit_ecfp"

    def name(self):
        return self.model_name

    def create_model(self, active_descriptors: str, inactive_descriptors: str, model_name: str):
        active_smiles = _select_molecules(active_descriptors)
        model = {
            "metadata": {
                "name": model_name
            },
            "data": {
                "active": active_smiles,
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        my_library.save_to_json_file(output_file, model)

    def score_model(self, model_data: list, fragments_file: str, descriptors_file: str, output_file: str):
        my_library.create_parent_directory(output_file)
        first_line = True
        with open(output_file, "w") as output_stream:
            with open(fragments_file, "r") as input_stream:
                active_molecules_ecfp = []
                radius = 0
                for new_line in input_stream:
                    line = json.loads(new_line)
                    if first_line:
                        radius = line["fragments"][0]["size"]
                        for molecule in model_data["active"]:
                            molecule = molecule.strip("\"")
                            molecule = Chem.MolFromSmiles(molecule)
                            ecfp_fingerprint = AllChem.GetMorganFingerprint(molecule, radius)
                            active_molecules_ecfp.append(ecfp_fingerprint)

                    test_molecule = line["smiles"]
                    test_molecule = test_molecule.strip("\"")
                    test_molecule = Chem.MolFromSmiles(test_molecule)
                    test_mol_fingerprint = AllChem.GetMorganFingerprint(test_molecule, radius)
                    max_sim = 0
                    for fingerprint in active_molecules_ecfp:
                        sim = DataStructs.DiceSimilarity(test_mol_fingerprint, fingerprint)
                        if sim > max_sim:
                            max_sim = sim
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _select_molecules(input_file: str):
    smiles = []
    with open(input_file, "r") as stream:
        next(stream)
        for line in stream:
            line_parts = line.split(",")
            smiles.append(line_parts[0])
    return smiles


register_model(RdkitEcfpModel.model_name, RdkitEcfpModel())
