#!/usr/bin/env python3
""""
Prints pie-chart graph that shows the ratio of number of actives and inactives molecules.
Usage:
    python active_inactive_graph.py
        -a (input file - output from add_activity.py)
        -o (output png file)
"""
import matplotlib.pyplot as plt

import argparse
import json
import os
import inputoutput_utils


def _main():
    configuration = _read_configuration()
    inputoutput_utils.create_parent_directory(configuration["output_file"])
    _print_graph(configuration["input_activity"], configuration["output_file"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="pie graph of actives and inactives molecules"
                    "See file header for more details.")
    parser.add_argument("-a", type=str, dest="input_activity",
                        help="input file", required=True)
    parser.add_argument("-o", type=str, dest="output_file",
                        help="output png file", required=True)

    configuration = vars(parser.parse_args())

    return configuration


def _print_graph(activity_file: str, output_file:str):
    actives = 0
    inactives = 0
    with open(activity_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            if int(line["activity"]) == 1:
                actives += 1
            else:
                inactives += 1
    act = ["actives", "inactives"]
    slices = [actives, inactives]
    colors = ["r", "b"]
    plt.pie(slices, labels=act, colors=colors)
    plt.legend()
    plt.title("number of actives vs inactives")
    plt.savefig(output_file)


if __name__ == "__main__":
    _main()
