#!/usr/bin/env python3
"""
Model that obtains equivalence classes in input configuration. Then for each active and test molecule
it changes their equivalent classes according to input configuration and then it computes
Tanimoto similarity and takes the biggest one.
input model_configuration should look like this:
    {"model_name": "equivalent_ecfp_model", "equivalence": [[class1], [class2], [class3], ...]}
    classes can look e.g. like this:
        ["C[NH+]", "C[NH+]"]
        ["c(c)(c)F", "c1cc(N)cc(O)c1"]
        ["Fc", "Fc", "Fc"]

"""

import json

from model_interface import IModel
from model_factory import register_model
import inputoutput_utils
import rdkitmodel_utils


class EquivalentClassModel(IModel):
    model_name = "equivalent_model"

    def name(self):
        return self.model_name

    def create_model(self, active_fragments: str, inactive_fragments: str,
                     active_descriptors: str, inactive_descriptors: str,
                     model_configuration: dict) -> dict:
        indexes_without_duplication = _modify_input_fragments(active_fragments)
        model = {
            "configuration": {
                "model_name": model_configuration["model_name"],
                "radius": rdkitmodel_utils.find_radius(active_fragments),
                "equivalence": model_configuration["equivalence"]
            },
            "data": {
                "active": indexes_without_duplication
            }
        }
        return model

    def save_to_json_file(self, output_file: str, model: dict):
        inputoutput_utils.save_to_json_file(output_file, model)

    def score_model(self, model_configuration: dict, fragments_file: str,
                    descriptors_file: str, output_file: str):
        inputoutput_utils.create_parent_directory(output_file)
        model_data = model_configuration["data"]
        celkem = []
        # this makes from e.g. ["CO","CO", "NH"] this: [["CO","CO"],["NH"]]
        for equivalence_class1 in model_configuration["configuration"]["equivalence"]:
            equivalence_class = sorted(equivalence_class1)
            fr = []
            i = 0
            while i < len(equivalence_class):
                fr1 = [equivalence_class[i]]
                while (i+1 < len(equivalence_class)) and (equivalence_class[i] == equivalence_class[i+1]):
                    i += 1
                    fr1.append(equivalence_class[i])
                i += 1
                fr.append(fr1)
            celkem.append(fr)
        # we want only classes where are different values: not [["CO"], ["CO"]],
        # but yes: [["CO"], ["C"]]
        heterogenous_equivalence_class = []
        for item in celkem:
            if len(item) > 1:
                heterogenous_equivalence_class.append(item)

        # reading active data from model and deleting the information about their indexes
        active_data = []
        for molecule_fragment in model_data["active"]:
            mol_data = []
            for item in molecule_fragment:
                mol_data.append(item[1:])
            active_data.append(mol_data)

        active_data = _equivalence_classes(active_data, heterogenous_equivalence_class)

        molecule_fragment = _modify_input_fragments(fragments_file)
        test_data = []
        for molecule_fragment1 in molecule_fragment:
            mol_data = []
            for item in molecule_fragment1:
                mol_data.append(item[1:])
            test_data.append(mol_data)

        first_line = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            with open(fragments_file, "r", encoding="utf-8") as input_stream:
                for num, new_line in enumerate(input_stream):
                    line = json.loads(new_line)
                    max_sim = 0
                    test_molecule_fragments = test_data[num]
                    for active_indexes in active_data:
                        intersection = _intersection_of_two_arrays(active_indexes, test_molecule_fragments)
                        sim = len(intersection) / (len(test_molecule_fragments) +
                                                   len(active_indexes) - len(intersection))
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


def _modify_input_fragments(file: str) -> list:
    indexes = []
    with open(file, "r", encoding="utf-8") as stream:
        for input_line in stream:
            line = json.loads(input_line)
            index1 = []
            for fragment in line["fragments"]:
                index = [fragment["index"], fragment["smiles"]]
                index1.append(index)
            indexes.append(index1)
    indexes_without_duplication = []
    for molecule_fragments in indexes:
        actual = 0
        solo_indexes = [molecule_fragments[0]]
        for i in range(1, len(molecule_fragments)):
            if int(molecule_fragments[i][0]) == int(solo_indexes[actual][0]):
                solo_indexes[actual].append(molecule_fragments[i][1])
            else:
                solo_indexes.append(molecule_fragments[i])
                actual += 1
        indexes_without_duplication.append(solo_indexes)
    return indexes_without_duplication


def _equivalence_classes(active_data: list, heterogenous_equivalence_class: list) -> list:
    for active_molecule_fragments in active_data:
        for class_equivalence in heterogenous_equivalence_class:
            founded_all = True
            indexes = []
            for item in class_equivalence:
                founded = False
                for num, molecule_fragment in enumerate(active_molecule_fragments):
                    if molecule_fragment[0:] == item:
                        founded = True
                        indexes.append(num)
                        break
                if founded is False:
                    founded_all = False
            if founded_all:
                i = len(indexes)
                while i > 0:
                    index = indexes[i - 1]
                    del active_molecule_fragments[index]
                    i -= 1
                active_molecule_fragments.append(class_equivalence)
    return active_data


register_model(EquivalentClassModel.model_name, lambda: EquivalentClassModel())

