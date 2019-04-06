#!/usr/bin/env python3
""""
Model that for each active molecule adds its indexes to dictionary if they are new or increase the
counting. If in active molecule is only one member of pair then the counting of second member of
pair is also increased. Then is done cutoff, which means that model only takes indexes that has got
more percentage threshold than is inputted cutoff. Then this pair thing is also applied on test
molecules and then the similarity is computed.
input model_configuration should look like this:
    {"model_name": "cutoff_add_paired_model", "cutoff":num, "pair": [num1,  num2]}
    num is between 0 and 100 and it simulates percentage
"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils

import model.add_paired_model as paired_model


class CutoffAddPairedModel(IModel):
    model_name = "cutoff_add_paired_model"

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
                for item in line["fragments"]:
                    if item["index"] in indexes:
                        indexes[item["index"]] += 1
                    else:
                        indexes[item["index"]] = 1
                index1 = model_configuration["pair"][0]
                index2 = model_configuration["pair"][1]
                if (index1 in indexes) and (index2 not in indexes):
                    if index2 in indexes:
                        indexes[index2] += 1
                    else:
                        indexes[index2] = 1
                if (index1 not in indexes) and (index2 in indexes):
                    if index1 in indexes:
                        indexes[index1] += 1
                    else:
                        indexes[index1] = 1

        cut = num_of_active_mol * int(model_configuration["cutoff"]) / 100
        cutoff_indexes = []
        for item in indexes:
            if indexes[item] >= cut:
                cutoff_indexes.append(item)
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
        model = paired_model.AddPairedModel()
        model.score_model(model_configuration, fragments_file, descriptors_file, output_file)


register_model(CutoffAddPairedModel.model_name, lambda: CutoffAddPairedModel())

