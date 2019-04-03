#!/usr/bin/env python3
""""
Takes files from run_all_nbit_configurations.py and prints graphs with AUC,, EF,...
"""
import argparse

import print_evaluation_graphs


def _main():
    configuration = _read_configuration()
    files_names = []
    files_nicknames = []
    num = 0
    for i in range(1024, 16384+256, 257):
        file_name = configuration["input_directory"] + str(num) + "/output.json"
        files_names.append(file_name)
        files_nickname = configuration["nickname"] + str(i)
        files_nicknames.append(files_nickname)
        if i == 1024:
            num += 1
            file_name = configuration["input_directory"] + str(num) + "/output.json"
            files_names.append(file_name)
            files_nickname = configuration["nickname"] + str(1025)
            files_nicknames.append(files_nickname)
        num += 1
    configuration = {
        "input_evaluation": files_names,
        "nicknames": files_nicknames,
        "directory": configuration["output_directory"]
    }
    print_evaluation_graphs.print_graphs(configuration)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="Prints graphs for more nbits output files ")
    parser.add_argument("-i", type=str, dest="input_directory", help="directory to input files",
                        required=True)
    parser.add_argument("-n", type=str, dest="nickname", help="name that will be printed in graph",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_directory", help="directory where to store"
                                                                      " results",
                        required=True)
    return vars(parser.parse_args())


if __name__ == "__main__":
    _main()
