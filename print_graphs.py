#!/usr/bin/env python3
""""
For each input file prints 3 graphs:
    1) histogram for score distribution of all molecules
    2) histogram for score distribution of active molecules
    3) histogram for score distribution of inactive molecules

Usage:
    python print_graphs.py
        -a (input comma separated activity files - output from different model
            from script add_activity.py)
        -n (input comma separated nicknames for files that will be printed in graphs)
        -d (directory where to store output files)
"""
import matplotlib.pyplot as plt

import argparse
import json
import inputoutput_utils


def _main():
    configuration = _read_configuration()
    for i in range(len(configuration["input_activity"])):
        _print_histograms_with_scores(configuration["input_activity"][i],
                                      configuration["directory"], configuration["nicknames"][i])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Prints 3 histogram images"
                    "See file header for more details.")
    parser.add_argument("-a", type=str, dest="input_activity",
                        help="input comma separated activity files", required=False)
    parser.add_argument("-n", type=str, dest="nicknames",
                        help="input comma separated nicknames for files"
                             " that will be printed in graphs",
                        required=True)
    parser.add_argument("-d", type=str, dest="directory",
                        help="directory where to store output files", required=True)

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


def _print_histograms_with_scores(activity_file: str, directory: str, nick: str):
    score_actives = []
    score_inactives = []
    with open(activity_file, "r") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            if int(line["activity"]) == 1:
                score_actives.append(line["score"])
            else:
                score_inactives.append(line["score"])
    bins = [i/20 for i in range(20, )]
    plt.hist(score_actives+score_inactives, bins=bins)
    file_name = directory + "/" + nick + "_activeinactive.png"
    inputoutput_utils.create_parent_directory(file_name)
    plt.savefig(file_name)
    plt.figure()

    plt.hist(score_actives, bins=bins)
    file_name = directory + "/" + nick + "_actives.png"
    plt.savefig(file_name)
    plt.figure()

    plt.hist(score_inactives, bins=bins)
    file_name = directory + "/" + nick + "_inactives.png"
    plt.savefig(file_name)
    plt.figure()


if __name__ == "__main__":
    _main()
