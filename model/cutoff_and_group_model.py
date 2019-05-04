#!/usr/bin/env python3
""""
Model that for each active molecule deletes multiple occurence of one index, adds its indexes
to dictionary if they are new or increase the counting. If active molecule does not have
all members from group then the counting of all other members from group is also increased.
Then is done cutoff, which means that model only takes indexes that has got
more percentage threshold than is inputted cutoff. Then this group thing is also applied on test
molecules and then the similarity is computed.
input model_configuration could look like this:
    {"model_name": "cutoff_add_paired_model", "fragments": "ecfp.6", "cutoff": num, "groups": [[num1,  num2, num3]]}
    num is between 0 and 100 and it simulates percentage
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.add_paired_model as paired_model


class CutoffAddGroupModel(IModel):
    model_name = "cutoff_and_group_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        if "cutoff" not in model_configuration:
            print("We need cutoff value")
            exit(1)
        cutoff = int(model_configuration["cutoff"])
        if (cutoff < 0) or (cutoff > 100):
            print("Wrong input. Cutoff has to be between 0 and 100!!")
            exit(1)
        indexes = {}
        num_of_active_mol = 0
        with open(active_fragments, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                num_of_active_mol += 1
                line = json.loads(new_line)
                indexes1 = _add_indexes(line["fragments"], model_configuration["groups"])
                for item in indexes1:
                    if item in indexes:
                        indexes[item] += 1
                    else:
                        indexes[item] = 1

        cut = num_of_active_mol * int(model_configuration["cutoff"]) / 100
        cutoff_indexes = []
        for item in indexes:
            if indexes[item] >= cut:
                cutoff_indexes.append(item)
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "groups": model_configuration["groups"]
            },
            "data": {
                "active": cutoff_indexes
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
                    test_active_indexes = _add_indexes(line["fragments"],
                                                       model_configuration["configuration"]["groups"])
                    max_sim = _compute_sim(model_configuration["data"]["active"], test_active_indexes)
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _add_indexes(fragments: dict, groups: list) -> list:
    molecule_indexes = []
    for item in fragments:
        if item["index"] not in molecule_indexes:
            molecule_indexes.append(item["index"])
    for group in groups:
        founded_pos = False
        founded_neg = False
        for item in group:
            if (founded_pos is False) and (item in molecule_indexes):
                founded_pos = True
            elif (founded_neg is False) and (item not in molecule_indexes):
                founded_neg = True
            if founded_pos and founded_neg:
                break
        if founded_pos and founded_neg:
            for item in group:
                if item not in molecule_indexes:
                    molecule_indexes.append(item)
    return molecule_indexes


def _compute_sim(active_fragments: list, test_fragments: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
            summary += 1
    sim = summary / (len(active_fragments) + len(test_fragments) - summary)
    return sim


register_model(CutoffAddGroupModel.model_name, lambda: CutoffAddGroupModel())


