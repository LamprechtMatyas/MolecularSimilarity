""""
evaluation of model, computing AUC, EF...

Usage:
    python compute_evaluation.py
        -i (json file with name, score, activity)
        -o (json file with AUC, EF, ...)
"""
import argparse
import json

from rdkit.ML.Scoring import Scoring

import my_library


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="model evaluation "
                                                 "See file header for more details.")
    parser.add_argument("-i", type=str, dest="input", required=True)
    parser.add_argument("-o", type=str, dest="output", required=True)
    return vars(parser.parse_args())


def read_file(input_file: str) -> list:
    score_activity = []
    with open(input_file, "r") as stream:
        for line_num, line in enumerate(stream):
            molecule = json.loads(line)
            score_activity.append([molecule["score"]])
            score_activity[line_num].append(molecule["activity"])
    return score_activity


def sort_activity(score_activity: list) -> list:
    score_activity = sorted(score_activity, reverse=True)
    return [
        [score_activity[i][1]]
        for i in range(len(score_activity))
    ]


def evaluation(activity_arr: list, output_file: str):
    my_library.create_parent_directory(output_file)
    auc = Scoring.CalcAUC(activity_arr, 0)
    ef1 = Scoring.CalcEnrichment(activity_arr, 0, [0.01])
    ef5 = Scoring.CalcEnrichment(activity_arr, 0, [0.05])
    rie = Scoring.CalcRIE(activity_arr, 0, 100)
    bedroc = Scoring.CalcBEDROC(activity_arr, 0, 100)
    output = {
        "AUC": auc,
        "EF1": ef1[0],
        "EF5": ef5[0],
        "RIE": rie,
        "BEDROC": bedroc
    }
    with open(output_file, "w") as stream:
        json.dump(output, stream)


def _main():
    configuration = _read_configuration()
    score_act = read_file(configuration["input"])
    activity = sort_activity(score_act)
    evaluation(activity, configuration["output"])


if __name__ == "__main__":
    _main()
