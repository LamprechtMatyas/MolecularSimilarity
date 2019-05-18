#!/usr/bin/env python3
""""
Baseline model for active_pair_model.py add_paired_model.py, so after we can compare which
paires are good to use.
input model_configuration should look like this:
    {"model_name": "baseline_active_pair_model", "fragments": "ecfp.6"}
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class BaselineActivePairModel(IModel):
    model_name = "baseline_active_pair_model"

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
                active_indexes.append(molecule_indexes)

        model = {
            "configuration": {
                "model_name": model_configuration["model_name"]
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
        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for new_line in input_stream:
                    line = json.loads(new_line)
                    test_active_indexes = []
                    for fragment in line["fragments"]:
                        if fragment["index"] not in test_active_indexes:
                            test_active_indexes.append(fragment["index"])
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


register_model(BaselineActivePairModel.model_name, lambda: BaselineActivePairModel())


