#!/usr/bin/env python3
""""
Model that takes pair of equivalence class, keeps active molecules separated (deletes multiple
same indexes - keeps only first occurrence) and applies this equivalence class to each active
molecule - so e.g. if pair is [1,2] and the molecule has both indexes then this pair replace each
index from pair, so it deletes [1], [2] and adds [1,2]. Then this same is applied on each molecule
in test set and then each molecule from test set in compared with each molecule from active
molecules and the maximum similarity for each test molecule is taken.
input model_configuration should look like this:
    {"model_name": "active_pair_model", "pair": [num1,  num2]}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class ActivePairModel(IModel):
    model_name = "active_pair_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = []
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                molecule_indexes = []
                for item in line["fragments"]:
                    if item["index"] not in molecule_indexes:
                        molecule_indexes.append(item["index"])
                if (int(model_configuration["pair"][0]) in molecule_indexes) and\
                    (int(model_configuration["pair"][1]) in molecule_indexes):
                    molecule_indexes.remove(int(model_configuration["pair"][0]))
                    molecule_indexes.remove(int(model_configuration["pair"][1]))
                    molecule_indexes.append(model_configuration["pair"])
                active_indexes.append(molecule_indexes)

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
                    test_active_indexes = []
                    for fragment in line["fragments"]:
                        if fragment["index"] not in test_active_indexes:
                            test_active_indexes.append(fragment["index"])
                    if (int(model_configuration["configuration"]["pair"][
                                0]) in test_active_indexes) and (
                            int(model_configuration["configuration"]["pair"][
                                    1]) in test_active_indexes):
                        test_active_indexes.remove(
                            int(model_configuration["configuration"]["pair"][0]))
                        test_active_indexes.remove(
                            int(model_configuration["configuration"]["pair"][1]))
                        test_active_indexes.append(model_configuration["configuration"]["pair"])
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


def _compute_sim(active_fragments: list, test_fragments: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
                summary += 1
    sim = summary / (len(active_fragments) + len(test_fragments) - summary)
    return sim


register_model(ActivePairModel.model_name, lambda: ActivePairModel())


