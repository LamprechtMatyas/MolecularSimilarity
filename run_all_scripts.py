#!/usr/bin/env python3
""""
runs all scripts

Usage:
    python run_all_scripts.py
        -a (optional input file with active molecules)
        -i (optional input file with inactive molecules)
        -t (optional input file with test molecules)
        -c (optional input file with activity of test molecules)
        -m (input configuration json file)
        -o (output file)
        -f (optional, comma seperated list of fragments)
        -p (optional, type of input molecules files, default sdf)
        -d (optional directory where to store intermediate results)
        --kekule {generated kekule form of SMILES for fragments}
        --isomeric {put stereochemistry information into fragments SMILES}
Fragments type:
    - tt.{SIZE}
    - ecfp.{SIZE}
    - fcfp.{SIZE}
    - ap

"""
import argparse
import logging
import tempfile
import json

import extract_fragments
import inputoutput_utils
import compute_descriptors
import model_factory
import add_activity
import compute_evaluation


def _main():
    # run extract_fragments
    configuration = _read_configuration()
    extraction_options = {
        "kekule": configuration["kekule"],
        "isomeric": configuration["isomeric"],
        "fragments": configuration["fragments"]
    }
    input_files = [configuration["input_actives"], configuration["input_inactives"],
                   configuration["test"]]
    directory = configuration["directory"]
    fragments_output_files = [directory+"/fragmentsa.json", directory+"/fragmentsi.json",
                              directory+"/fragmentst.json"]
    for file in fragments_output_files:
        inputoutput_utils.create_parent_directory(file)
    extract_fragments.extract_fragments(input_files, configuration["input_type"],
                                        fragments_output_files, extraction_options)

    # run extract_descriptors
    with open(configuration["model_configuration"], "r", encoding="utf-8") as input_stream:
        model_configuration = json.load(input_stream)
    descriptors_output_files = [directory+"/descriptorsa.csv", directory+"/descriptorsi.csv",
                                directory+"/descriptorst.csv"]
    for file in descriptors_output_files:
        inputoutput_utils.create_parent_directory(file)
    if (model_configuration["model_name"] == "decision_tree_model") |\
        ((model_configuration["model_name"] == "linear_regression_model") and (int(model_configuration["molecules"]) == 1)):
        compute_descriptors.compute_descriptors(fragments_output_files, descriptors_output_files,
                                                False)
    else:
        compute_descriptors.compute_descriptors(fragments_output_files, descriptors_output_files, True)

    # run create_model and score_molecules
    new_model = model_factory.create_model(model_configuration["model_name"])
    model = new_model.create_model(directory+"/fragmentsa.json", directory+"/fragmentsi.json",
                                   directory+"/descriptorsa.csv", directory+"/descriptorsi.csv",
                                   model_configuration)
    new_model.score_model(model, directory+"/fragmentst.json",
                          directory+"/descriptorst.csv", directory+"/score.json")

    # run add_activity
    activity = add_activity.read_activity(configuration["activity"])
    add_activity.add_activity_and_write_to_json(directory + "/score.json", activity,
                                                directory + "/activity.json")

    #  run compute_evaluation
    score_act = compute_evaluation.read_file_with_score_and_activity(directory + "/activity.json")
    activity = compute_evaluation.sort_activity(score_act)
    compute_evaluation.evaluation(activity, configuration["output"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="run all scripts "
                    "See file header for more details.")
    parser.add_argument("-a", type=str, dest="input_actives",
                        help="input file with active molecules", required=False,
                        default="data/actives.smi")
    parser.add_argument("-i", type=str, dest="input_inactives",
                        help="input file with inactive molecules", required=False,
                        default="data/inactives.smi")
    parser.add_argument("-t", type=str, dest="test",
                        help="input file with test molecules", required=False,
                        default="data/test.smi")
    parser.add_argument("-c", type=str, dest="activity",
                        help="input file with activity of test molecules", required=False,
                        default="data/test_activity.json")
    parser.add_argument("-m", type=str, dest="model_configuration",
                        help="input json file with model configuration", required=True)
    parser.add_argument("-o", type=str, dest="output",
                        help="output json file", required=True)
    parser.add_argument("-f", type=str, dest="fragments",
                        help="optional, comma separated list of fragments", required=False)
    parser.add_argument("-p", type=str, dest="input_type",
                        help="type of input files with molecules smi/sdf, default is smi",
                        default="smi")
    parser.add_argument("-d", dest="directory",
                        help="directory where to store intermediate results", required=False)
    parser.add_argument("--kekule", dest="kekule", help="kekule option",
                        action="store_true", required=False)
    parser.add_argument("--isomeric", dest="isomeric", help="isomeric option",
                        action="store_true", required=False)

    configuration = vars(parser.parse_args())
    if configuration["fragments"] is None:
        configuration["fragments"] = "tt.3"
    if configuration["directory"] is None:
        configuration["directory"] = tempfile.gettempdir()
    parsed_types = []
    for item in configuration["fragments"].split(","):
        item_split = item.split(".")
        if item_split[0] != "ap":
            if not len(item_split) == 2:
                logging.error("Invalid fragment type: %s", item)
                logging.info("Expected format {TYPE}.{SIZE} or ap")
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
    configuration["input_type"] = configuration["input_type"].lower()
    return configuration


if __name__ == "__main__":
    _main()
