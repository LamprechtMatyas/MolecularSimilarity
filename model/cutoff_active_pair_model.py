#!/usr/bin/env python3
""""
Model that takes all indexes together from all active molecules. Indexes that appear in equal
or more molecules than is cutoff this one rest. Then we pair them as in active_pair_model.py
and we evaluate the test molecules in a same way as in active_pair_model.py
input model_configuration should look like this:
    {"model_name": "cutoff_active_pair_model", "cutoff":num, "pair": [num1,  num2]}
    num has to be number between 0 and 100, because it is percent
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.active_pair_model as pair_model


class CutoffActivePairModel(IModel):
    model_name = "cutoff_active_pair_model"

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
                for item in line["fragments"]:
                    if item["index"] in indexes:
                        indexes[item["index"]] += 1
                    else:
                        indexes[item["index"]] = 1
        cut = num_of_active_mol * int(model_configuration["cutoff"]) / 100
        cutoff_indexes = []
        for item in indexes:
            if indexes[item] >= cut:
                cutoff_indexes.append(item)
        if (int(model_configuration["pair"][0]) in cutoff_indexes) and \
                (int(model_configuration["pair"][1]) in cutoff_indexes):
            cutoff_indexes.remove(int(model_configuration["pair"][0]))
            cutoff_indexes.remove(int(model_configuration["pair"][1]))
            cutoff_indexes.append(model_configuration["pair"])
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "pair": model_configuration["pair"]
            },
            "data": {
                "active": [cutoff_indexes]
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: list, fragments_file: str,
                    descriptors_file: str, output_file: str):
        model = pair_model.ActivePairModel()
        model.score_model(model_configuration, fragments_file, descriptors_file, output_file)


register_model(CutoffActivePairModel.model_name, lambda: CutoffActivePairModel())


