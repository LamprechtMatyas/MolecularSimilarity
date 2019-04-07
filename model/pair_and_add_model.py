#!/usr/bin/env python3
""""
Model that keeps active molecules separately and for each active molecule constructs the pair (if
it is possible). Then we compare each molecule from test set with molecule each molecule in
model and we compute similarity in this way: if test molecule has only 1 index from pair then
we add +2 to the similarity. This means that the maximum similarity can be 2, so we have to divide
similarity by 2, because we want similarity between 0 and 1.
input model_configuration should look like this:
    {"model_name": "pair_and_add_model", "pair": [num1,num2]}
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.active_pair_model as pair_model


class PairAndAddModel(IModel):
    model_name = "pair_and_add_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        new_model = pair_model.ActivePairModel()
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
                    test_active_indexes = []
                    for fragment in line["fragments"]:
                        if fragment["index"] not in test_active_indexes:
                            test_active_indexes.append(fragment["index"])
                    pair1 = int(model_configuration["configuration"]["pair"][0])
                    pair2 = int(model_configuration["configuration"]["pair"][1])
                    pair = model_configuration["configuration"]["pair"]
                    if (pair1 in test_active_indexes) and (pair2 in test_active_indexes):
                        test_active_indexes.remove(pair1)
                        test_active_indexes.remove(pair2)
                        test_active_indexes.append(pair)
                    max_sim = max([_compute_sim(item, test_active_indexes, pair) for item in
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


def _compute_sim(active_fragments: list, test_fragments: list, pair: list) -> list:
    summary = 0
    for item in test_fragments:
        if item in active_fragments:
                summary += 1
    if pair in active_fragments:
        if (pair[0] in test_fragments) and (pair[1] not in test_fragments):
            summary += 2
        elif (pair[0] not in test_fragments) and (pair[1] in test_fragments):
            summary += 2
    sim = summary / (2*(len(active_fragments) + len(test_fragments) - summary))
    return sim


register_model(PairAndAddModel.model_name, lambda: PairAndAddModel())


