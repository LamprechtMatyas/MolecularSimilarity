#!/usr/bin/env python3
""""
Script that prints graph.
Usage:
    python print_graph.py
        -e (comma separated input files(each file is an output of script: compute_evaluation.py or
            run_all_scrips.py))
        -n (comma separated nicknames - for each input file there must be its nickname, that will
            identify the input in a graph)
        -t (type of graph that will be printed,
            it must be one of keys in input files e.g. AUC or EF1)
        -d (directory where to store graph)
"""
import argparse
import json

import matplotlib.pyplot as plt

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    print_graph(configuration["input_evaluation"], configuration["directory"],
                configuration["nicknames"], configuration["type"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="script that prints graph")
    parser.add_argument("-e", type=str, dest="input_evaluation",
                        help="input comma separated files with evaluation", required=True)
    parser.add_argument("-n", type=str, dest="nicknames",
                        help="input nicknames for files that will be printed in graph",
                        required=True)
    parser.add_argument("-t", type=str, dest="type",
                        help="type of graph", required=True)
    parser.add_argument("-d", dest="directory",
                        help="directory where to store graph", required=True)

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


def print_graph(activity_files: list, directory: str, nicknames: list, input_type: str):
    input_values = []
    for file in activity_files:
        with open(file, "r", encoding="utf-8") as activity_file:
            for new_line in activity_file:
                line = json.loads(new_line)
                input_values.append(line[input_type.upper()])
    plt.plot(nicknames, input_values, marker="o")
    if input_type.upper() == "EF1":
        plt.ylabel("EF 1%")
    elif input_type.upper() == "EF5":
        plt.ylabel("EF 5%")
    else:
        plt.ylabel(input_type.upper())

    if directory != "":
        file_name = directory + "/" + input_type.upper() + ".png"
    else:
        file_name = input_type.upper() + ".png"
    inputoutput_utils.create_parent_directory(file_name)
    plt.xticks(rotation=90, fontsize="x-small")
    plt.tight_layout()
    plt.savefig(file_name, dpi=150)
    plt.figure()


if __name__ == "__main__":
    _main()
