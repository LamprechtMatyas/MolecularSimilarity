#!/usr/bin/env python3
""""
It takes files from different models from script add_activity.py and prints graph that shows
how each model scores active molecules.
Usage:
    python active_score_graph.py
        -a (input comma separated activity files)
        -n (input nicknames for files that will be printed in graphs)
        -o (output png file)
"""
import matplotlib.pyplot as plt

import argparse
import json
import inputoutput_utils


def _main():
    configuration = _read_configuration()
    _print_graph(configuration["input_activity"], configuration["nicknames"],
                 configuration["output_file"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Active molecule scoring by different models. "
                                                 "See file header for more details.")
    parser.add_argument("-a", type=str, dest="input_activity",
                        help="input comma separated activity files", required=False)
    parser.add_argument("-n", type=str, dest="nicknames",
                        help="input nicknames for files that will be printed in graphs",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file",
                        help="output png file", required=True)
    configuration = vars(parser.parse_args())

    input_files = []
    for file in configuration["input_activity"].split(","):
        input_files.append(file)
    configuration["input_activity"] = input_files

    input_nicknames = []
    for file in configuration["nicknames"].split(","):
        input_nicknames.append(file)
    configuration["nicknames"] = input_nicknames

    if len(configuration["input_activity"]) != len(configuration["nicknames"]):
        print("Wrong input. Number of input files must be equal to the number of nicknames.")
        exit(1)

    return configuration


def _print_graph(activity_files: list, nicknames: list, file_name: str):
    inputoutput_utils.create_parent_directory(file_name)
    for file_num, file in enumerate(activity_files):
        with open(file, "r", encoding="utf-8") as activity_file:
            names = []
            scores = []
            for new_line in activity_file:
                line = json.loads(new_line)
                if int(line["activity"]) == 1:
                    names.append(line["name"])
                    scores.append(line["score"])
            plt.plot(names, scores, marker="o")
    plt.legend(nicknames)
    plt.xlabel("molecule names")
    plt.ylabel("molecule scores")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(file_name)


if __name__ == "__main__":
    _main()
