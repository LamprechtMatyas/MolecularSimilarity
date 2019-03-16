#!/usr/bin/env python3
""""
Library for more used functions
"""

import os
import json


def create_parent_directory(path: str):
    dir_name = os.path.dirname(path)
    if dir_name != "":
        os.makedirs(dir_name, exist_ok=True)


def save_to_json_file(output_file: str, model: dict):
    create_parent_directory(output_file)
    with open(output_file, "w", encoding="utf-8") as stream:
        json.dump(model, stream)

