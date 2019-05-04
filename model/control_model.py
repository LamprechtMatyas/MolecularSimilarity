#!/usr/bin/env python3
""""
Control model that controls that we compute correctly ecfp, fcfp, ap, tt.
We use indexes of active molecules and we simulate Tanimoto computation.
input model_configuration should look like this:
    {"model_name": "control_model", "fragments": "ecfp.6"}
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class ControlModel(IModel):
    model_name = "control_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str, model_configuration: dict):
        indexes = []
        with open(active_fragments, "r", encoding="utf-8") as stream:
            for input_line in stream:
                index = []
                line = json.loads(input_line)
                for fragment in line["fragments"]:
                    index.append(fragment["index"])
                indexes.append(index)
        model = {
            "configuration": model_configuration,
            "data": {
                "active": indexes
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
                    indexes = []
                    for fragment in line["fragments"]:
                        indexes.append(fragment["index"])
                    max_sim = 0
                    for active_indexes in model_configuration["data"]["active"]:
                        intersection = _intersection_of_two_arrays(active_indexes, indexes)
                        sim = len(intersection) / (len(indexes) + len(active_indexes) - len(intersection))
                        if sim > max_sim:
                            max_sim = sim
                    score = {
                        "name": line["name"],
                        "score": max_sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _intersection_of_two_arrays(arr1: list, arr2: list):
    and_arr = []
    for item in arr1:
        for i in range(len(arr2)):
            if item == arr2[i]:
                and_arr.append(item)
                arr2 = arr2[:i] + arr2[(i+1):]
                break
    return and_arr


register_model(ControlModel.model_name, lambda: ControlModel())
