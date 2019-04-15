#!/usr/bin/env python3
""""
Model that takes all indexes together and keeps each index occurrence once. We make groups
on these indexes and then we also make groups (if we can) on test molecules.
input model_configuration can look e.g. like this:
    {"model_name": "ecfp_group_model", "groups": [[num1,  num2], [num3, num4, num5]}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class EcfpGroupModel(IModel):
    model_name = "ecfp_group_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = []
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                for item in line["fragments"]:
                    if item["index"] not in active_indexes:
                        active_indexes.append(item["index"])
        groups_indexes = _make_groups(active_indexes, model_configuration["groups"])

        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "groups": model_configuration["groups"]
            },
            "data": {
                "active": groups_indexes
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
                    group_indexes = _make_groups(test_active_indexes,
                                                 model_configuration["configuration"]["groups"])
                    summary = 0
                    for test_index in group_indexes:
                        if test_index in model_configuration["data"]["active"]:
                            summary += 1
                    sim = summary / (len(model_configuration["data"]["active"]) +
                                     len(test_active_indexes) - summary)
                    score = {
                        "name": line["name"],
                        "score": sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _make_groups(indexes: list, groups: list) -> list:
    for group in groups:
        founded = True
        for item in group:
            if item not in indexes:
                founded = False
                break
        if founded:
            for item in group:
                indexes.remove(item)
            indexes.append(group)
    return indexes


register_model(EcfpGroupModel.model_name, lambda: EcfpGroupModel())


