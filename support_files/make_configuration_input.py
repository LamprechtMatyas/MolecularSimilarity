import json
import argparse

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    make_configuration_input(configuration["model"], configuration["directory"])


def make_configuration_input(model: str, directory: str):
    inputoutput_utils.create_parent_directory(directory + "/0")
    for num, i in enumerate(range(1024, 16384+256, 256)):
        configuration = {
            "model_name": model,
            "nbits": i
        }
        output_file = directory + "/configuration" + str(num) + ".json"
        with open(output_file, "w") as output_stream:
            json.dump(configuration, output_stream)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Creation of configuration files.")
    parser.add_argument("-m", type=str, dest="model",
                        help="model_name", required=True)
    parser.add_argument("-d", type=str, dest="directory",
                        help="output directory", required=True)
    return vars(parser.parse_args())


if __name__ == "__main__":
    _main()

