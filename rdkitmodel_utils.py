#!/usr/bin/env python3

import json


def create_model(active_fragments: str, model_configuration: dict) -> dict:
    active_smiles = select_molecules(active_fragments)
    model = {
        "configuration": model_configuration,
        "data": {
            "active": active_smiles
        }
    }
    return model


def select_molecules(input_file: str) -> list:
    molecules = []
    with open(input_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            molecules.append(line["smiles"])
    return molecules


