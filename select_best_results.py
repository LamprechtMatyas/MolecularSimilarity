#!/usr/bin/env python3
"""
From output file from pair_analysis.py or groups_analysis.py, it prints to
the output file only the given number of best results (based on AUC, EF 1% or
EF 5% which is given in input) 
"""
import argparse
import json
import inputoutput_utils


def _main():
    configuration = _read_configuration()
    inputoutput_utils.create_parent_directory(configuration["output"])
    values = [[], [], [], [], []]
    if configuration["type"] == "AUC":
        str1 = "EF1"
        str2 = "EF5"
    elif configuration["type"] == "EF1":
        str1 = "AUC"
        str2 = "EF5"
    elif configuration["type"] == "EF5":
        str1 = "AUC"
        str2 = "EF1"
    else:
        print("Wrong type!")
        print("It has to be: AUC. EF1 or EF5")
        exit(1)
    with open(configuration["input_file"], "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            values[0].append(line["groups"])
            values[1].append(line[configuration["type"]])
            values[2].append(line[str1])
            values[3].append(line[str2])
    if len(values[0]) < int(configuration["best"]):
        print("The input file does not have that much good results!")
        print("Number of results in input file: " + str(len(values[0])))
        print("Number you wanted to select: " + configuration["best"])
        print("Can not be like that!")
        exit(1)
    for i in range(len(values[1]) - 1):
        for j in range(len(values[1]) - i - 1):
            if values[1][j] < values[1][j+1]:
                tmp = values[1][j]
                values[1][j] = values[1][j+1]
                values[1][j+1] = tmp
                tmp = values[0][j]
                values[0][j] = values[0][j+1]
                values[0][j+1] = tmp
                tmp = values[2][j]
                values[2][j] = values[2][j + 1]
                values[2][j + 1] = tmp
                tmp = values[3][j]
                values[3][j] = values[3][j + 1]
                values[3][j + 1] = tmp
    with open(configuration["output"], "w", encoding="utf-8") as output_stream:
        for i in range(int(configuration["best"])):
            model = {
                "groups": values[0][i],
                configuration["type"]: values[1][i],
                str1: values[2][i],
                str2: values[3][i]
            }
            json.dump(model, output_stream)
            output_stream.write("\n")


def _read_configuration():
    parser = argparse.ArgumentParser(description="best results selecting "
                                                 "See file header for more details.")
    parser.add_argument("-f", type=str, dest="input_file",
                        help="file from which you want to select best", required=True)
    parser.add_argument("-b", type=str, dest="best", help="number of best results", required=True)
    parser.add_argument("-t", type=str, dest="type", help="type - AUC, EF1 or EF5", required=True)
    parser.add_argument("-o", type=str, dest="output", help="output_file", required=True)

    configuration = vars(parser.parse_args())
    return configuration


if __name__ == "__main__":
    _main()


