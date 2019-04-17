#!/usr/bin/env python3
""""
Division of groups to files based on AUC, EF 1%, EF 5%. We have some baseline model and
if our pair has better results(means > ) e.g. in AUC and EF 5% then we add this pair to file aucef5.
So outputfiles are:
auc, aucef1, aucaf5, aucef1ef5, ef1, ef1ef5, ef5, greater(takes all that were in
some aspect higher)
"""
import argparse
import json
from os import listdir
from os.path import isfile, join

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    baseline_results = _mean_results(configuration["baseline_input"])
    inputoutput_utils.create_parent_directory(configuration["output_directory"] + "/0")
    _prepare_files(configuration["output_directory"])
    files = [f for f in listdir(configuration["configuration"]) if
                 isfile(join(configuration["configuration"], f))]
    groups = []
    for file in files:
        with open(configuration["configuration"] + "/" + file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                groups.append(line["groups"])

    evaluation_files = [f for f in listdir(configuration["evaluations"]) if
                        isfile(join(configuration["evaluations"], f))]
    num = 0
    auc = baseline_results[0]
    ef1 = baseline_results[1]
    ef5 = baseline_results[2]
    for file in evaluation_files:
        with open(configuration["evaluations"] + "/" + file, "r", encoding="utf-8") as input_stream:
            new_line = input_stream.read()
            line = json.loads(new_line)
            output = {
                "groups": groups[num],
                "AUC": line["AUC"],
                "EF1": line["EF1"],
                "EF5": line["EF5"]
            }
            if (line["AUC"] > auc) and (line["EF1"] > ef1) and (line["EF5"] > ef5):
                with open(configuration["output_directory"] + "/aucef1ef5.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if (line["AUC"] > auc) and (line["EF1"] > ef1):
                with open(configuration["output_directory"] + "/aucef1.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if (line["AUC"] > auc) and (line["EF5"] > ef5):
                with open(configuration["output_directory"] + "/aucef5.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if (line["EF1"] > ef1) and (line["EF5"] > ef5):
                with open(configuration["output_directory"] + "/ef1ef5.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if line["AUC"] > auc:
                with open(configuration["output_directory"] + "/auc.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if line["EF5"] > ef5:
                with open(configuration["output_directory"] + "/ef5.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if line["EF1"] > ef1:
                with open(configuration["output_directory"] + "/ef1.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            if (line["AUC"] > auc) or (line["EF1"] > ef1) or (line["EF5"] > ef5):
                with open(configuration["output_directory"] + "/greater.json",
                          "a", encoding="utf-8") as output_stream:
                    json.dump(output, output_stream)
                    output_stream.write("\n")
            num += 1
    with open(configuration["output_directory"] + "/baseline.json", "w", encoding="utf-8") as output_stream:
        output = {
            "AUC": auc,
            "EF1": ef1,
            "EF5": ef5
        }
        json.dump(output, output_stream)


def _read_configuration():
    parser = argparse.ArgumentParser(description="group analysis "
                                                 "See file header for more details.")
    parser.add_argument("-b", type=str, dest="baseline_input",
                        help="file from pair_analysis or select best results or group analysis", required=True)
    parser.add_argument("-c", type=str, dest="configuration", help="directory with configuration files from"
                                                                   " run_all_groups", required=True)
    parser.add_argument("-e", type=str, dest="evaluations", help="directory to evaluation files,"
                                                                 " there must be stored only evaluation files",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_directory", help="output directory", required=True)

    configuration = vars(parser.parse_args())
    return configuration

    
def _mean_results(file: str) -> list:
    results = [0, 0, 0]
    num = 0
    with open(file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            num += 1
            line = json.loads(new_line)
            results[0] += float(line["AUC"])
            results[1] += float(line["EF1"])
            results[2] += float(line["EF5"])
    results[0] = results[0] / num
    results[1] = results[1] / num
    results[2] = results[2] / num
    return results  
    
      
def _prepare_files(output_directory: str):
    with open(output_directory + "/aucef1ef5.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/aucef1.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/aucef5.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/ef1ef5.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/auc.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/ef1.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/ef5.json",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/greater.json", "w", encoding="utf-8"):
        pass
  

if __name__ == "__main__":
    _main()
