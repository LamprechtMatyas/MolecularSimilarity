#!/usr/bin/env python3
""""
For each molecule in test set we compare its descriptors
with the active descriptors and we compute a sum.
If we find descriptor in active descriptors we add to the sum +num1
and if we do not find it we add -num2
input model_configuration should look like this:
    {"model_name": "descriptors_model", "active_parameter": num1,"inactive_parameter": num2}
    where num1 and num2 are numbers
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class DescriptorsModel(IModel):
    model_name = "descriptors_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: str) -> dict:
        descriptors = []
        with open(active_descriptors, "r") as stream:
            line = stream.readline()
            line_parts = line.split(",")
            if line_parts[1] != "index":
                print("Wrong input, we want fragments descriptors not molecules")
                exit(1)
            for line in stream:
                line_parts = line.split(",")
                descriptors.append(line_parts[2:])
        model = {
            "configuration": model_configuration,
            "data": descriptors
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        name_num = _read_molecules(fragments_file)
        inputoutput_utils.create_parent_directory(output_file)
        active_parameter = int(model_configuration["configuration"]["active_parameter"])
        inactive_parameter = int(model_configuration["configuration"]["inactive_parameter"])
        with open(output_file, "w") as streamo:
            first_write = True
            with open(descriptors_file, "r") as stream:
                next(stream)
                counter = 0
                molecule_num = 0
                sum = 0
                for line in stream:
                    line_parts = line.split(",")
                    parts = line_parts[2:]
                    founded = False
                    for descriptors in model_configuration["data"]:
                        if descriptors == parts:
                            founded = True
                            break

                    if founded:
                        sum += active_parameter
                    else:
                        sum -= inactive_parameter
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


def _read_molecules(input_file: str):
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


register_model(DescriptorsModel.model_name, lambda: DescriptorsModel())


