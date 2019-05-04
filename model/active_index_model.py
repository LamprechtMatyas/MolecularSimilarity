#!/usr/bin/env python3
""""
Model based on active indexes.
For each molecule in test set we compare its indexes with the active indexes and we compute a sum.
If we find index in set of active indexes then we add +num1 to the sum else we add -num2 to the sum.
input model_configuration should look like this:
    {"model_name": "active_index_model", "fragments": "ecfp.6", "active_parameter": num1, "inactive_parameter": num2}
    where num1 and num2 are numbers
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class PositiveIndexModel(IModel):
    model_name = "active_index_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        indexes = []
        with open(active_fragments, "r") as stream:
            for line in stream:
                input_line = json.loads(line)
                for item in input_line["fragments"]:
                    indexes.append(item["index"])
        model = {
            "configuration": model_configuration,
            "data": indexes
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        active_parameter = int(model_configuration["configuration"]["active_parameter"])
        inactive_parameter = int(model_configuration["configuration"]["inactive_parameter"])
        with open(output_file, "w") as streamo:
            first = True
            with open(fragments_file, "r") as stream:
                for line in stream:
                    molecule = json.loads(line)
                    name = molecule["name"]
                    suma = 0
                    for fragment in molecule["fragments"]:
                        founded = False
                        for index in model_configuration["data"]:
                            if int(index) == int(fragment["index"]):
                                founded = True
                                break
                        if founded:
                            suma += active_parameter
                        else:
                            suma -= inactive_parameter
                    sim = suma / len(molecule["fragments"])
                    score = {
                        "name": name,
                        "score": sim
                    }
                    if first:
                        first = False
                    else:
                        streamo.write("\n")
                    json.dump(score, streamo)


register_model(PositiveIndexModel.model_name, lambda: PositiveIndexModel())


