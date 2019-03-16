#!/usr/bin/env python3
""""
Creation of a specific model

Usage: python create_model
        -c (configuration json file of a model)
        -a (input json file with active fragments)
        -f (input json file with inactive fragments)
        -d (input csv file with descriptors from active molecules)
        -i (input csv file with descriptors from inactive molecules)
        -o (output json file with model)
"""
import argparse
import json

import model_factory
import inputoutput_utils


def _main():
    configuration = _read_configuration()
    with open(configuration["configuration"], "r", encoding="utf-8") as input_stream:
        model_configuration = json.load(input_stream)
    model_name = model_configuration["model_name"]

    new_model = model_factory.create_model(model_name)
    model = new_model.create_model(configuration["active_fragments"], configuration["inactive_fragments"],
                                   configuration["active_descriptors"], configuration["inactive_descriptors"],
                                   model_configuration)

    inputoutput_utils.create_parent_directory(configuration["output"])
    new_model.save_to_json_file(configuration["output"], model)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(
        description="Model creation.")
    parser.add_argument("-c", type=str, dest="configuration",
                        help="configuration input json fle", required=True)
    parser.add_argument("-a", type=str, dest="active_fragments",
                        help="input json file with active fragments", required=True)
    parser.add_argument("-f", type=str, dest="inactive_fragments",
                        help="input json file with inactive fragments", required=True)
    parser.add_argument("-d", type=str, dest="active_descriptors",
                        help="input csv file with active descriptors", required=True)
    parser.add_argument("-i", type=str, dest="inactive_descriptors",
                        help="input csv file with inactive descriptors", required=True)
    parser.add_argument("-o", type=str, dest="output",
                        help="output file", required=True)
    return vars(parser.parse_args())


if __name__ == "__main__":
    _main()
