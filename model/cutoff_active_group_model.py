#!/usr/bin/env python3
""""
Model that takes all indexes together from all active molecules. Indexes that appear in equal
or more molecules than is cutoff these ones rest. Then we make groups from them
as in active_group_model.py and we evaluate the test molecules in a same way as
in active_group_model.py.
input model_configuration could look like this:
    {"model_name": "cutoff_active_pair_model", "cutoff":num, "groups": [[num1,  num2], [num3, num4]]}
    num has to be number between 0 and 100, because it is percent
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.active_pair_model as pair_model


class CutoffActiveGroupModel(IModel):
    model_name = "cutoff_active_group_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        indexes = {}
        num_of_active_mol = 0
        cutoff = int(model_configuration["cutoff"])
        if (cutoff < 0) or (cutoff > 100):
            print("Wrong input. Cutoff has to be between 0 and 100!!")
            exit(1)
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                num_of_active_mol += 1
                line = json.loads(new_line)
                indexes1 = []
                for item in line["fragments"]:
                    if item["index"] not in indexes1:
                        indexes1.append(item["index"])
                for item in indexes1:
                    if item not in indexes:
                        indexes[item] = 1
                    else:
                        indexes[item] += 1

        cut = num_of_active_mol * int(model_configuration["cutoff"]) / 100
        cutoff_indexes = []
        for item in indexes:
            if indexes[item] >= cut:
                cutoff_indexes.append(item)
        group_indexes = _make_groups(cutoff_indexes, model_configuration["groups"])
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "groups": model_configuration["groups"]
            },
            "data": {
                "active": group_indexes
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
                    max_sim = _compute_sim(model_configuration["data"]["active"], group_indexes)
                    score = {
                        "name": line["name"],
                        "score": max_sim
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


def _compute_sim(active_fragments: list, test_fragments: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
            summary += 1
    sim = summary / (len(active_fragments) + len(test_fragments) - summary)
    return sim


register_model(CutoffActiveGroupModel.model_name, lambda: CutoffActiveGroupModel())


