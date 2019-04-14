#!/usr/bin/env python3
"""
For given evaluation files and baseline evaluation, it prints histogram for AUC
values or line for EF 1%, EF 5% for values that were higher than baseline value.
"""

import argparse
import json
from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    base_value = _baseline_val(configuration["baseline_output"], configuration["type"])
    values = _return_all_values(configuration["input_directory"], configuration["type"], base_value)
    _print_graph(base_value, values, configuration["type"], configuration["output_file"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="grpahs of pairs "
                                                 "See file header for more details.")
    parser.add_argument("-b", type=str, dest="baseline_output", required=True)
    parser.add_argument("-d", type=str, dest="input_directory",
                        help="directory to evaluation files", required=True)
    parser.add_argument("-t", type=str, dest="type", help="type of evaluation - AUC, EF1, EF5",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file",
                        help="output  png file", required=True)
    return vars(parser.parse_args())
    
    
def _baseline_val(file: str, type1: str) -> float:
    with open(file, "r", encoding="utf-8") as input_stream:
        new_line = input_stream.read()
        line = json.loads(new_line)
        print(line)
        print(line[type1])

        return float(line[type1])   
    
        
def _return_all_values(directory: str, type1: str, val: float) -> list:
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    values = []
    for file in files:
        with open(directory + "/" + file, "r", encoding="utf-8") as input_stream:
            new_line = input_stream.read()
            line = json.loads(new_line)
            if float(line[type1]) > val:
                values.append(line[type1])

    return values   
    
          
def _print_graph(val: float, values: list, type1: str, output: str):
    inputoutput_utils.create_parent_directory(output)
    textstr = "baseline " + type1 + ": " + "%.6f" %(val)
    if type1 == "AUC":
        min_val = min(values)
        max_val = max(values)
        diff = max_val - min_val
        diff_step = diff/10
        steps = []
        for i in range(10):
            steps.append(i*diff_step + min_val)
        steps.append(max_val)
        arr = plt.hist(values, bins=steps, color="blue")
        plt.xticks(steps, rotation=90)
        for i in range(10):
            plt.text(arr[1][i], arr[0][i], str(int(arr[0][i])), horizontalalignment="left")
        props = dict(boxstyle="round")
        plt.text(steps[7], int(arr[0][0]), textstr)
        plt.tight_layout()
        plt.savefig(output)
    else:
        num = {}
        values1 = sorted(values)
        for item in values1:
            if item in num:
                num[item] += 1
            else:
                num[item] = 1
        item = []
        values = []
        for l in num:
            item.append(l)
            values.append(num[l])
        xval = [val-2] + [val] + (item)
        plt.plot(item, values, marker="o")
        plt.xticks(xval)

        plt.axvline(val, color="red")
        plt.savefig(output)

    print(min(values))
    print(max(values))


if __name__ == "__main__":
    _main()
