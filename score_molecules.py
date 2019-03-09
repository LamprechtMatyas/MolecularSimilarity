""""
Molecules scoring

Usage: python score_molecules.py
            -m (input json file with model)
            -f (input json file with fragments from the test set)
            -d (input csv file with descriptors from the test set)
            -o (output json file)
"""
import model_factory

import argparse
import json


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Molecule scoring "
                                             "See file header for more details.")
    parser.add_argument("-m", type=str, dest="model",
                        help="input model", required=True)
    parser.add_argument("-f", type=str, dest="fragments",
                        help="input test fragments", required=True)
    parser.add_argument("-d", type=str, dest="descriptors",
                        help="input test descriptors", required=True)
    parser.add_argument("-o", type=str, dest="output",
                        help="output file", required=True)
    return vars(parser.parse_args())


configuration = _read_configuration()
with open(configuration["model"]) as input_stream:
    model_content = json.load(input_stream)

model = model_factory.create_model(model_content["metadata"]["name"])
model.score_model(model_content["data"], configuration["fragments"], configuration["descriptors"],
                  configuration["output"])





