#!/usr/bin/env python3
""""
Model that does that if there is only one member of pair in molecule then it adds to the molecule
second member. This is done to actives molecules, that are separated and then to test molecules.
input model_configuration should look like this:
    {"model_name": "add_paired_model", "fragments": "ecfp.6", "pair": [num1,  num2]}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class AddPairedModel(IModel):
    model_name = "add_paired_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = []
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                active_indexes.append(_add_pair(line["fragments"], int(model_configuration["pair"][0]),
                                      int(model_configuration["pair"][1])))

        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "pair": model_configuration["pair"]
            },
            "data": {
                "active": active_indexes
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: list, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for new_line in input_stream:
                    line = json.loads(new_line)
                    test_active_indexes = _add_pair(line["fragments"],
                                                    int(model_configuration["configuration"]["pair"][0]),
                                                    int(model_configuration["configuration"]["pair"][1]))
                    max_sim = max([_compute_sim(item, test_active_indexes) for item in
                                   model_configuration["data"]["active"]])
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _add_pair(fragments: dict, pair1: int, pair2: int) -> list:
    molecule_indexes = []
    for item in fragments:
        if item["index"] not in molecule_indexes:
            molecule_indexes.append(item["index"])
    if (pair1 not in molecule_indexes) and (pair2 in molecule_indexes):
        molecule_indexes.append(pair1)
    elif (pair1 in molecule_indexes) and (pair2 not in molecule_indexes):
        molecule_indexes.append(pair2)
    return molecule_indexes


def _compute_sim(active_fragments: list, test_fragments: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
                summary += 1
    sim = summary / (len(active_fragments) + len(test_fragments) - summary)
    return sim


register_model(AddPairedModel.model_name, lambda: AddPairedModel())


