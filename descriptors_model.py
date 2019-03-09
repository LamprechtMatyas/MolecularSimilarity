""""
For each molecule in test set we compare its descriptors with the active descriptors and we compute a sum.
If we find descriptor in active descriptors we add to the sum +1 and if we do not find it we add -1
"""

from model_interface import IModel
from model_factory import register_model
import my_library

import json


class DescriptorsModel(IModel):
    model_name = "descriptors_model"

    def name(self):
        return self.model_name

    def create_model(self, active_descriptors: str, inactive_descriptors: str, model_name: str) -> dict:
        descriptors = []
        with open(active_descriptors, "r") as stream:
            for line in stream:
                line_parts = line.split(",")
                descriptors.append(line_parts[2:])
        descriptors.pop(0)
        model = {
            "metadata": {
                "name": model_name
            },
            "data": descriptors
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        my_library.save_to_json_file(output_file, model)

    def score_model(self, data: list, fragments_file: str, descriptors_file: str, output_file: str):
        name_num = self._read_molecules(fragments_file)
        my_library.create_parent_directory(output_file)
        with open(output_file, "w") as streamo:
            first = True
            first_write = True
            with open(descriptors_file, "r") as stream:
                counter = 0
                molecule_num = 0
                sum = 0
                for line in stream:
                    if first:
                        first = False
                        continue
                    line_parts = line.split(",")
                    parts = line_parts[2:]
                    founded = False
                    for descriptors in data:
                        if descriptors == parts:
                            founded = True
                            break

                    if founded:
                        sum += 1
                    else:
                        sum -= 1
                    counter += 1
                    if counter == name_num[molecule_num]["fragments"]:
                        score = {
                            "name": name_num[molecule_num]["molecule"],
                            "score": sum / counter
                        }
                        counter = 0
                        sum = 0
                        molecule_num += 1
                        if first_write:
                            first_write = False
                        else:
                            streamo.write("\n")
                        json.dump(score, streamo)

    def _read_molecules(self, input_file: str):
        name_num = []
        with open(input_file, "r") as stream:
            for line in stream:
                molecule = json.loads(line)
                name = molecule["name"]
                num_fragments = len(molecule["fragments"])
                mol_frag = {
                    "molecule": name,
                    "fragments": num_fragments
                }
                name_num.append(mol_frag)
        return name_num


register_model(DescriptorsModel.model_name, DescriptorsModel())


