#!/usr/bin/env python3
""""
Runs whole process: from extract_fragments.py to compute_evaluation for 61 nbits values.
"""
import json
import argparse

import inputoutput_utils
import extract_fragments
import compute_evaluation
import add_activity
import model_factory
import make_configuration_input


def _main():
    configuration = _read_configuration()
    input_files = ["data/actives.smi", "data/inactives.smi", "data/test.smi"]
    make_configuration_input.make_configuration_input(configuration["model"],
                                                      configuration["output_directory"])
    for i in range(61):
        file_name = configuration["output_directory"] + "/configuration" + str(i) + ".json"
        with open(file_name, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                model_configuration = json.loads(new_line)
                def_directory = configuration["output_directory"]+ "/" + configuration["model"]\
                                + "/" + configuration["model"] + "0"
                directory = configuration["output_directory"] + "/" + configuration["model"]\
                            + "/" + configuration["model"] + str(i)
                if i == 0:
                    extraction_options = {
                        "kekule": None,
                        "isomeric": None,
                        "fragments": configuration["fragments"]
                    }
                    fragments_output_files = [directory + "/fragmentsa.json",
                                              directory + "/fragmentsi.json",
                                              directory + "/fragmentst.json"]
                    for file in fragments_output_files:
                        inputoutput_utils.create_parent_directory(file)
                    extract_fragments.extract_fragments(input_files, "smi",
                                                        fragments_output_files, extraction_options)

                # run create_model and score_molecules
                new_model = model_factory.create_model(model_configuration["model_name"])
                model = new_model.create_model(def_directory + "/fragmentsa.json",
                                               def_directory + "/fragmentsi.json", "", "",
                                               model_configuration)
                new_model.score_model(model, def_directory + "/fragmentst.json", "",
                                      directory + "/score.json")

                # run add_activity
                activity = add_activity.read_activity("data/test_activity.json")
                add_activity.add_activity_and_write_to_json(directory + "/score.json", activity,
                                                            directory + "/activity.json")

                #  run compute_evaluation
                score_act = compute_evaluation.read_file_with_score_and_activity(
                    directory + "/activity.json")
                activity = compute_evaluation.sort_activity(score_act)
                compute_evaluation.evaluation(activity, directory + "/output.json")


def _read_configuration():
    parser = argparse.ArgumentParser(description="Running whole process for more nbit configurations"
                                                 "See file header for more details.")
    parser.add_argument("-m", type=str, dest="model",
                        help="name of model", required=True)
    parser.add_argument("-f", type=str, dest="fragments", required=False)
    parser.add_argument("-o", type=str, dest="output_directory", help="output directory",
                        required=True)

    configuration = vars(parser.parse_args())
    parsed_types = []
    for item in configuration["fragments"].split(","):
        item_split = item.split(".")
        if item_split[0] != "ap":
            if not len(item_split) == 2:
                exit(1)
            parsed_types.append({
                "name": item_split[0],
                "size": int(item_split[1])
            })
        else:
            parsed_types.append({
                "name": item_split[0],
            })
    configuration["fragments"] = parsed_types
    return configuration


if __name__ == "__main__":
    _main()

