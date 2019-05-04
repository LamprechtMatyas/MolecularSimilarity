#!/usr/bin/env python3
""""
Model that keeps active molecules separately and for each active molecule constructs the groups (if
it is possible). Then we compare each molecule from test set with each molecule in
model and we compute similarity in this way: if test molecule has at least one index from group
then we add to the similarity sum +(length of the group). This means that the maximum similarity
can be more than 1, so we have to divide similarity by number of indexes in groups,
because we want similarity between 0 and 1.
input model_configuration could look like this:
    {"model_name": "group_and_add_model", "fragments": "ecfp.6", "groups": [[num1,num2], [num3, num4, num5]]}
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.active_group_model as group_model


class GroupAndAddModel(IModel):
    model_name = "group_and_add_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        new_model = group_model.ActiveGroupModel()
        return new_model.create_model(active_fragments, inactive_fragments, active_descriptors,
                                      inactive_descriptors, model_configuration)

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
                    groups = model_configuration["configuration"]["groups"]
                    test_active_indexes = _make_group_indexes(line["fragments"], groups)
                    max_sim = max([_compute_sim(item, test_active_indexes, groups) for item in
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
    test_active_indexes = []
    for fragment in fragments:
        if fragment["index"] not in test_active_indexes:
            test_active_indexes.append(fragment["index"])
    for group in groups:
        founded = True
        for item in group:
            if item not in test_active_indexes:
                founded = False
                break
        if founded:
            for item in group:
                test_active_indexes.remove(item)
            test_active_indexes.append(group)
    return test_active_indexes


def _compute_sim(active_fragments: list, test_fragments: list, groups: list) -> list:
    groups_len = 0
    for item in groups:
        groups_len += len(item)
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
                summary += 1
    for group in groups:
        if group in active_fragments:
            founded = False
            for item in group:
                if item in test_fragments:
                    founded = True
                    break
            if founded:
                summary += len(group)
    sim = summary / (groups_len*(len(active_fragments) + len(test_fragments) - summary))
    return sim


register_model(GroupAndAddModel.model_name, lambda: GroupAndAddModel())


