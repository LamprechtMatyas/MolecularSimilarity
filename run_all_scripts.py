""""
runs all scripts

Usage:
    python run_all_scripts.py
        -a (input file with active molecules)
        -i (input file with inactive molecules)
        -t (input file with test molecules)
        -c (input file with activity of test molecules)
        -o (output file)
        -m (type of model)
        -f (optional, comma seperated list of fragments)
        -p (optional, type of input molecules files, default sdf)
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

import extract_fragments
import my_library
import compute_descriptors
import model_factory
import add_activity
import compute_evaluation


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="run all scripts "
                    "See file header for more details.")
    parser.add_argument("-a", type=str, dest="inputa", required=True)
    parser.add_argument("-i", type=str, dest="inputi", required=True)
    parser.add_argument("-t", type=str, dest="test", required=True)
    parser.add_argument("-c", type=str, dest="activity", required=True)
    parser.add_argument("-m", type=str, dest="model", required=True)
    parser.add_argument("-o", type=str, dest="output", required=True)
    parser.add_argument("-f", type=str, dest="fragments", required=False)
    parser.add_argument("-p", type=str, dest="input_type", default="sdf")
    parser.add_argument("--kekule", dest="kekule",
                        action="store_true", required=False)
    parser.add_argument("--isomeric", dest="isomeric",
                        action="store_true", required=False)

    configuration = vars(parser.parse_args())
    if "fragments" not in configuration or configuration["fragments"] is None:
        configuration["fragments"] = "tt.3"
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


def _main():
    # run extract_fragments
    configuration = _read_configuration()
    extraction_options = {
        "kekule": configuration["kekule"],
        "isomeric": configuration["isomeric"],
        "fragments": configuration["fragments"]
    }
    input_files = [configuration["inputa"], configuration["inputi"], configuration["test"]]
    fragments_output_files = ["data1/tmp/fragmentsa.json", "data1/tmp/fragmentsi.json", "data1/tmp/fragmentst.json"]
    for file in fragments_output_files:
        my_library.create_parent_directory(file)
    extract_fragments.extract_fragments(input_files, configuration["input_type"],
                                        fragments_output_files, extraction_options)

    # run extract_descriptors
    descriptors_output_files = ["data1/tmp/descriptorsa.csv", "data1/tmp/descriptorsi.csv", "data1/tmp/descriptorst.csv"]
    for file in descriptors_output_files:
        my_library.create_parent_directory(file)
    compute_descriptors.compute_descriptors(fragments_output_files, descriptors_output_files, True)

    # run create_model and score_molecules
    new_model = model_factory.create_model(configuration["model"])
    model = new_model.create_model("data1/tmp/descriptorsa.csv", "data1/tmp/descriptorsi.csv", configuration["model"])
    new_model.score_model(model["data"], "data1/tmp/fragmentst.json", "data1/tmp/descriptorst.csv", "data1/tmp/score.json")

    # run add_activity
    activity = add_activity.read_activity(configuration["activity"])
    add_activity.write_to_json("data1/tmp/score.json", activity, "data1/tmp/activity.json")

    #  run compute_evaluation
    score_act = compute_evaluation.read_file("data1/tmp/activity.json")
    activity = compute_evaluation.sort_activity(score_act)
    compute_evaluation.evaluation(activity, configuration["output"])


if __name__ == "__main__":
    _main()
