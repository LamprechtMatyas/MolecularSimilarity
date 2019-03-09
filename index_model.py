""""
Model based on active indexes.
For each molecule in test set we compare its indexes with the active indexes and we compute a sum.
If we find index in set of active indexes then we add +1 to the sum else we add -1 to the sum.
"""
from model_interface import IModel
from model_factory import register_model
import my_library

import json


class PositiveIndexModel(IModel):
    model_name = "index_model"

    def name(self):
        return self.model_name

    def create_model(self, active_descriptors: str, inactive_descriptors: str, model_name: str) -> dict:
        indexes = []
        with open(active_descriptors, "r") as stream:
            for line in stream:
                line_parts = line.split(",")
                indexes.append(line_parts[1])
        indexes.pop(0)
        model = {
            "metadata": {
                "name": model_name
            },
            "data": indexes
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        my_library.save_to_json_file(output_file, model)

    def score_model(self, data: list, fragments_file: str, descriptors_file: str, output_file: str):
        my_library.create_parent_directory(output_file)
        with open(output_file, "w") as streamo:
            first = True
            with open(fragments_file, "r") as stream:
                for line in stream:
                    molecule = json.loads(line)
                    name = molecule["name"]
                    sum = 0
                    for fragment in molecule["fragments"]:
                        founded = False
                        for index in data:
                            if int(index) == int(fragment["index"]):
                                founded = True
                                break
                        if founded:
                            sum += 1
                        else:
                            sum -= 1
                    sim = sum / len(molecule["fragments"])
                    score = {
                        "name": name,
                        "score": sim
                    }
                    if first:
                        first = False
                    else:
                        streamo.write("\n")
                    json.dump(score, streamo)


register_model(PositiveIndexModel.model_name, PositiveIndexModel())





