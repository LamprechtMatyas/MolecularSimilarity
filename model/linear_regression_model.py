#!/usr/bin/env python3
""""
Model that uses standard machine learning method - linear regression
It can be used for whole molecules or their fragments, but it has to be consistent.
If we use fragments then the similarity is computed as average similarity of the descriptors.
input model_configuration should look like this:
    {"model_name": "linear_regression_model", "molecules": 0/1}
    where 0 means that we use fragments, 1 means we use molecules
"""
import json

from sklearn import linear_model

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class LinearRegressionModel(IModel):
    model_name = "linear_regression_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict):
        act_descriptors = extract_descriptors(active_descriptors, model_configuration)
        inact_descriptors = extract_descriptors(inactive_descriptors, model_configuration)
        model = {
            "configuration": model_configuration,
            "data": {
                "active": act_descriptors,
                "inactive": inact_descriptors
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        reg = linear_model.LinearRegression()
        # get activity list
        actives = [1 for i in range(len(model_configuration["data"]["active"]))]
        inactives = [0 for i in range(len(model_configuration["data"]["inactive"]))]
        activity = actives + inactives

        reg.fit(model_configuration["data"]["active"] + model_configuration["data"]["inactive"],
                activity)
        test_descriptors = extract_descriptors(descriptors_file, model_configuration["configuration"])
        molecule_file = int(model_configuration["configuration"]["molecules"])
        prediction = (reg.predict(test_descriptors))
        if molecule_file == 1:
            first_line = True
            with open(output_file, "w", encoding="utf-8") as output_stream:
                with open(fragments_file, "r", encoding="utf-8") as input_stream:
                    for num_line, new_line in enumerate(input_stream):
                        line = json.loads(new_line)
                        score = {
                            "name": line["name"],
                            "score": prediction[num_line]
                        }
                        if first_line:
                            first_line = False
                        else:
                            output_stream.write("\n")
                        json.dump(score, output_stream)

        else:
            num_of_fragment = [0]
            names_of_molecules = []
            with open(fragments_file, "r", encoding="utf-8") as fragments_stream:
                suma = 0
                for new_line in fragments_stream:
                    line = json.loads(new_line)
                    fragment_length = len(line["fragments"])
                    suma += fragment_length
                    num_of_fragment.append(suma)
                    names_of_molecules.append(line["name"])
            first_line = True
            with open(output_file, "w", encoding="utf-8") as output_stream:
                for i in range(len(num_of_fragment) - 1):
                    prediction_of_molecule = prediction[num_of_fragment[i]:num_of_fragment[i+1]]
                    sim = sum(prediction_of_molecule) / len(prediction_of_molecule)
                    score = {
                        "name": names_of_molecules[i],
                        "score": sim
                    }
                    if first_line:
                        first_line = False
                    else:
                        output_stream.write("\n")
                    json.dump(score, output_stream)


def extract_descriptors(input_file: str, model_configuration: dict) -> list:
    descriptors = []
    with open(input_file, "r", encoding="utf-8") as stream:
        line = stream.readline()
        line_parts = line.split(",")
        if (((line_parts[1] == "index") & (int(model_configuration["molecules"]) == 1)) |
                ((line_parts[1] != "index") & (int(model_configuration["molecules"]) == 0))):
            print("Wrong input")
            exit(1)
        for line in stream:
            line_parts = line.split(",")
            descriptors.append(list(map(float, line_parts[1:])))
    return descriptors


register_model(LinearRegressionModel.model_name, lambda: LinearRegressionModel())
