#!/usr/bin/env python3
""""
Model that takes groups of equivalence classes, keeps active molecules separated (deletes multiple
same indexes - keeps only first occurrence) and applies these equivalence classes to each active
molecule - so e.g. if group is [1,2,3] and the molecule has these three indexes then this group replace each
index from group, so it deletes [1], [2], [3] and adds [1,2,3]. Then this same is applied on each molecule
in test set and then each molecule from test set in compared with each molecule from active
molecules and the maximum similarity for each test molecule is taken.
input model_configuration could look e.g. like this:
    {"model_name": "active_group_model", "groups": [[num1,  num2, num3], [num4, num5, num6]]}
    {"model_name": "active_group_model", "groups": [[num1,  num2, num3, num4]]}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class ActiveGroupModel(IModel):
    model_name = "active_group_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = []
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                molecule_indexes = _make_group_indexes(line["fragments"], model_configuration["groups"])
                active_indexes.append(molecule_indexes)
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "groups": model_configuration["groups"]
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
                    test_active_indexes = _make_group_indexes(line["fragments"],
                                                              model_configuration["configuration"]["groups"])
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


def _make_group_indexes(fragments: dict, groups: list) -> list:
    molecule_indexes = []
    for item in fragments:
        if item["index"] not in molecule_indexes:
            molecule_indexes.append(item["index"])
    for group in groups:
        founded = True
        for item in group:
            if item not in molecule_indexes:
                founded = False
                break
        if founded:
            for item in group:
                molecule_indexes.remove(item)
            molecule_indexes.append(group)
    return molecule_indexes


def _compute_sim(active_fragments: list, test_fragments: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
                summary += 1
    sim = summary / (len(active_fragments) + len(test_fragments) - summary)
    return sim


register_model(ActiveGroupModel.model_name, lambda: ActiveGroupModel())


