""""
Creation of a specific model

Usage: python create_model
        -m (name of a model)
        -a (input file with descriptors from active molecules)
        -i (input fie with descriptors from inactive molecules)
        -o (output json file with model)
"""
import argparse

import model_factory
import my_library


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(
        description="Model creation.")
    parser.add_argument("-m", type=str, dest="model",
                        help="model type", required=True)
    parser.add_argument("-a", type=str, dest="active",
                        help="input csv file with active fragments", required=True)
    parser.add_argument("-i", type=str, dest="inactive",
                        help="input csv file with inactive fragments", required=True)
    parser.add_argument("-o", type=str, dest="output",
                        help="output file", required=True)
    return vars(parser.parse_args())


configuration = _read_configuration()
model_name = configuration["model"]

new_model = model_factory.create_model(model_name)
model = new_model.create_model(configuration["active"], configuration["inactive"], model_name)

my_library.create_parent_directory(configuration["output"])
new_model.save_to_json_file(configuration["output"], model)


