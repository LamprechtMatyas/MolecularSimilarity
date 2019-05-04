#!/usr/bin/env python3
""""
Model that does same hashing as is in MorganFingerprints - AllChem.GetHashedMorganFingerprint(),
so it hashes all active molecules and then it hashed all test molecules and for each test molecule
it computes similarity to all active hashed molecules and it takes the highest one.
Configuration file should look like this:
    {"model_name": "my_nbit_hashed_ap_model", "fragments": "ecfp.6", "nbits": num}
    where num in natural number
"""

import json

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils
import rdkitmodel_utils


class MyNbitHashedApModel(IModel):
    model_name = "my_nbit_hashed_ap_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = []
        nbit = int(model_configuration["nbits"])
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                molecule_indexes = []
                for item in line["fragments"]:
                    molecule_indexes.append(int(item["index"]) % nbit)
                active_indexes.append(molecule_indexes)
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "nbits": model_configuration["nbits"]
            },
            "data": {
                "active": active_indexes
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        nbit = int(model_configuration["configuration"]["nbits"])

        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for new_line in input_stream:
                    line = json.loads(new_line)
                    test_indexes = []
                    for item in line["fragments"]:
                        test_indexes.append(int(item["index"]) % nbit)

                    max_sim = 0
                    for active_indexes in model_configuration["data"]["active"]:
                        intersection = _intersection_of_two_arrays(active_indexes, test_indexes)
                        sim = len(intersection) / (
                                    len(test_indexes) + len(active_indexes) - len(intersection))
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


def _intersection_of_two_arrays(arr1: list, arr2: list):
    and_arr = []
    for item in arr1:
        for i in range(len(arr2)):
            if item == arr2[i]:
                and_arr.append(item)
                arr2 = arr2[:i] + arr2[(i+1):]
                break
    return and_arr


register_model(MyNbitHashedApModel.model_name, lambda: MyNbitHashedApModel())
