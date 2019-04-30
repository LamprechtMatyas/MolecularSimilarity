#!/usr/bin/env python3
""""
Model that takes groups of equivalence classes, keeps active molecules separated.
On each active molecule applies equivalence class in that way, that e.g.
classes: [[a, b, c], [d, e]] and we have got a, b, a, ,e g, g, h, then we will
keep a, d, g, g, h. That means that we keep only 1 occurence of the group and
it is first value.
Then we do the same thing to all test molecules and we compute max-fusion
similarity.
input model_configuration could like e.g like this:
    {"model_name": "delete_index_group_model", "groups": [[num1,  num2, num3], [num4, num5, num6]]}}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class DeleteIndexGroupModel(IModel):
    model_name = "delete_index_group_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        molecules_indexes = []
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                molecule_indexes = []
                for fragment in line["fragments"]:
                    index = fragment["index"]
                    for group in model_configuration["groups"]:
                        if index in group:
                            index = group[0]
                            break
                    if index not in molecule_indexes:
                        molecule_indexes.append(index)
                molecules_indexes.append(molecule_indexes)        

        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "groups": model_configuration["groups"]
            },
            "data": {
                "active": molecules_indexes
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
                        index = fragment["index"]
                        for group in model_configuration["configuration"]["groups"]:
                            if index in group:
                                index = group[0]
                                break
                        if index not in test_active_indexes:
                            test_active_indexes.append(index)
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

register_model(DeleteIndexGroupModel.model_name, lambda: DeleteIndexGroupModel())