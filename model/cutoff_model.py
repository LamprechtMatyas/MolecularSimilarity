#!/usr/bin/env python3
""""
Model that takes all active molecules together with their indexes and keeps indexes
that appear in equal or more molecules than is cutoff. Then the score of test molecules is done
in same way as in baseline_model.py.
input model_configuration should look like this:
    {"model_name": "cutoff_model", "fragments": "ecfp.6", "cutoff": num}
    num is number of percent
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.baseline_model as baseline


class CutOffModel(IModel):
    model_name = "cutoff_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
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
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"]
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

        baseline_model = baseline.BaselineModel()
        baseline_model.score_model(model_configuration, fragments_file, descriptors_file,
                                   output_file)


register_model(CutOffModel.model_name, lambda: CutOffModel())


