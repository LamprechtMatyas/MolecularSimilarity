""""
Library for more used functions
"""

import os
import json


def create_parent_directory(path: str):
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)


def save_to_json_file(output_file: str, model: dict):
    create_parent_directory(output_file)
    with open(output_file, "w") as stream:
        json.dump(model, stream)

