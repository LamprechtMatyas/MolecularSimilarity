#!/usr/bin/env python3
""""
Script that prints graphs. For each key in output of compute_evaluation.py or run_all_scripts.py
it prints one graph.
Usage:
    python print_graph.py
        -e (comma separated input files(each file is an output of script: compute_evaluation.py or
            run_all_scrips.py))
        -n (comma separated nicknames - for each input file there must be its nickname, that will
            identify the input in a graph)
        -d (directory where to store graph)
"""

import argparse
import json

import print_graph


def _main():
    configuration = _read_configuration()
    keys = []
    with open(configuration["input_evaluation"][0], "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            keys = line.keys()
    for item in keys:
        print_graph.print_graph(configuration["input_evaluation"], configuration["directory"],
                                configuration["nicknames"], item)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="script that prints graphs from evaluation files"
                                                 "See file header for more details.")
    parser.add_argument("-e", type=str, dest="input_evaluation",
                        help="input comma separated files with evaluation", required=True)
    parser.add_argument("-n", type=str, dest="nicknames",
                        help="input nicknames for files that will be printed in some graphs",
                        required=True)
    parser.add_argument("-d", dest="directory",
                        help="directory where to store graphs results", required=False)

    configuration = vars(parser.parse_args())

    input_files = []
    for file in configuration["input_evaluation"].split(","):
        input_files.append(file)
    configuration["input_evaluation"] = input_files

    input_nicknames = []
    for file in configuration["nicknames"].split(","):
        input_nicknames.append(file)
    configuration["nicknames"] = input_nicknames

    if configuration["directory"] is None:
        configuration["directory"] = ""

    if len(configuration["input_evaluation"]) != len(configuration["nicknames"]):
        print("Wrong input, the number of nicknames must be equal to the number of input files")
        exit(1)
    return configuration


if __name__ == "__main__":
    _main()
