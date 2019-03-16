""""
Model based on active and inactive indexes.
For each molecule in test set we compare its indexes with the active indexes
and then with inactive indexes and we compute a sum.
If the index is only in active in indexes then we add +num1 to the sum, if the index is only
in inactive indexes the we add -num3, otherwise we add +num2 to the sum.
input model_configuration should look like this:
    {"model_name": "active_inactive_index_model", "active_parameter": num1,
    "neutral_parameter": num2,"inactive_parameter": num3}
    where num1, num2, num3 are numbers
"""
import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class PosNegIndexModel(IModel):
    model_name = "active_inactive_index_model"

    def name(self):
        return self.model_name

    def create_model(self,  active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        active_indexes = _select_indexes(active_descriptors)
        inactive_indexes = _select_indexes(inactive_descriptors)
        model = {
            "configuration": model_configuration,
            "data": {
                "active": active_indexes,
                "inactive": inactive_indexes
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        active_parameter = int(model_configuration["configuration"]["active_parameter"])
        neutral_parameter = int(model_configuration["configuration"]["neutral_parameter"])
        inactive_parameter = int(model_configuration["configuration"]["inactive_parameter"])
        with open(output_file, "w") as output_stream:
            first = True
            with open(fragments_file, "r") as stream:
                for line in stream:
                    molecule = json.loads(line)
                    name = molecule["name"]
                    suma = 0
                    for fragment in molecule["fragments"]:
                        is_in_actives = False
                        for index in model_configuration["data"]["active"]:
                            if int(index) == int(fragment["index"]):
                                is_in_actives = True
                                break
                        is_in_inactives = False
                        for index in model_configuration["data"]["inactive"]:
                            if int(index) == int(fragment["index"]):
                                is_in_inactives = True
                                break
                        if is_in_actives & is_in_inactives is False:
                            suma += active_parameter
                        elif is_in_inactives & is_in_actives is False:
                            suma -= inactive_parameter
                        else:
                            suma += neutral_parameter
                    score = {
                        "name": name,
                        "score": suma / len(molecule["fragments"])
                    }
                    if first:
                        first = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def _select_indexes(input_file: str):
    indexes = []
    with open(input_file, "r") as stream:
        for line in stream:
            line_parts = line.split(",")
            indexes.append(line_parts[1])
    indexes.pop(0)
    return indexes


register_model(PosNegIndexModel.model_name, lambda: PosNegIndexModel())
