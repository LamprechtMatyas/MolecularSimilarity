#!/usr/bin/env python3
""""
Model that uses built in function from rdkit library for computing
topological torsion fingerprints and computing similarity by using
TanimotoSimilarity from DataStructs library.
Each molecule from test set is compared with all known active molecules and the highest
similarity is that assigned to a molecule.
input model_configuration should look like this:
    {"model_name": "rdkit_tt_model", "fragments": "tt.num1"}
    num1 is a natural number
"""

import json

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem.AtomPairs import Torsions

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils
import rdkitmodel_utils


class RdkitTTModel(IModel):
    model_name = "rdkit_tt_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        return rdkitmodel_utils.create_model(active_fragments, model_configuration)

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        model_data = model_configuration["data"]
        radius = int(model_configuration["configuration"]["fragments"][0]["size"]) 
        active_molecules_tt = []
        for active_molecule in model_data["active"]:
            molecule_smiles = active_molecule.strip("\"")
            molecule = Chem.MolFromSmiles(molecule_smiles)
            tt_fingerprint = Torsions.GetTopologicalTorsionFingerprintAsIntVect(molecule, radius)
            active_molecules_tt.append(tt_fingerprint)
        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for new_line in input_stream:
                    line = json.loads(new_line)
                    test_molecule_input = line["smiles"]
                    test_molecule_smiles = test_molecule_input.strip("\"")
                    test_molecule = Chem.MolFromSmiles(test_molecule_smiles)
                    test_mol_fingerprint = Torsions.GetTopologicalTorsionFingerprintAsIntVect(test_molecule, radius)
                    max_sim = max([DataStructs.TanimotoSimilarity(test_mol_fingerprint, fingerprint)
                                  for fingerprint in active_molecules_tt])
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


register_model(RdkitTTModel.model_name, lambda: RdkitTTModel())


