""""
Model that uses standard machine learning method - decision trees
input model_configuration should look like this:
    {"model_name": "decision_tree_model"}
"""
import json

from sklearn import tree

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils


class DecisionTreeClassifier(IModel):
    model_name = "decision_tree_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict):
        act_descriptors = _extract_descriptors(active_descriptors)
        inact_descriptors = _extract_descriptors(inactive_descriptors)
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
        decision_tree = tree.DecisionTreeClassifier()
        # ger activity list
        actives = [1 for i in range(len(model_configuration["data"]["active"]))]
        inactives = [0 for i in range(len(model_configuration["data"]["inactive"]))]
        activity = actives + inactives

        decision_tree.fit(model_configuration["data"]["active"] + model_configuration["data"]["inactive"],
                          activity)
        test_descriptors = _extract_descriptors(descriptors_file)
        prediction = decision_tree.predict(test_descriptors)
        molecule_names = _extract_names(fragments_file)
        #write output
        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            for i in range(len(prediction)):
                score = {
                    "name": molecule_names[i],
                    "score": int(prediction[i])
                }
                if first_line:
                    first_line = False
                else:
                    output_stream.write("\n")
                json.dump(score, output_stream)


def _extract_names(input_file: str) -> tuple:
    names = []
    with open(input_file, "r", encoding="utf-8") as stream:
        for new_line in stream:
            line = json.loads(new_line)
            names.append(line["name"])
    return names


def _extract_descriptors(input_file: str) -> list:
    descriptors = []
    with open(input_file, "r", encoding="utf-8") as stream:
        line = stream.readline()
        line_parts = line.split(",")
        if line_parts[1] == "index":
            print("Wrong input, expected molecules not fragments")
            exit(1)
        for line in stream:
            line_parts = line.split(",")
            descriptors.append(list(map(float, line_parts[1:])))
    return descriptors


register_model(DecisionTreeClassifier.model_name, lambda: DecisionTreeClassifier())
