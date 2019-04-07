#!/usr/bin/env python3
""""
Takes output from pair_analysis as one of input files and divides pairs to files based on:
if they appear in more test active molecules or more test inactive molecules
or they appear in same count of active and inactive test molecules or
they do not appear in active neither inactive test molecules.
"""
import argparse
import json

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    activity_list = read_activity(configuration["activity_file"])
    all_actives = 0
    all_inactives = 0
    inputoutput_utils.create_parent_directory(configuration["output_directory"] + "/0")
    output_directory = configuration["output_directory"]
    _prepare_files(output_directory)
    pair_active_inactive = []
    i = 0
    with open(configuration["analysis_file"], "r", encoding="utf-8") as analysis_stream:
        for new_analysis in analysis_stream:
            analysis = json.loads(new_analysis)
            pair = analysis["pair"]
            pair_active_inactive.append([pair[0]])
            pair_active_inactive[i].append(pair[1])
            pair_active_inactive[i].append(0)
            pair_active_inactive[i].append(0)
            i += 1
    # print(pair_active_inactive)
    with open(configuration["test_fragments"], "r", encoding="utf-8") as fragments_stream:
        for num, new_fragments in enumerate(fragments_stream):
            molecule = json.loads(new_fragments)
            for i in range(len(pair_active_inactive)):
                first_founded = False
                second_founded = False
                for fragment in molecule["fragments"]:
                    if int(fragment["index"]) == int(pair_active_inactive[i][0]):
                        first_founded = True
                    elif int(fragment["index"]) == int(pair_active_inactive[i][1]):
                        second_founded = True
                if (first_founded is True) and (second_founded is True):
                    if int(activity_list[num]) == 0:
                        pair_active_inactive[i][3] += 1
                        all_inactives += 1
                    elif int(activity_list[num]) == 1:
                        pair_active_inactive[i][2] += 1
                        all_actives += 1
            # print(num)

    for i in range(len(pair_active_inactive)):
        if pair_active_inactive[i][2] > pair_active_inactive[i][3]:
            with open(output_directory + "/actives.txt",
                        "a", encoding="utf-8") as output_stream:
                output_stream.write("\n" + str([pair_active_inactive[i][0], pair_active_inactive[i][1]]) + "\n")
                output_stream.write("number of active:    " + str(pair_active_inactive[i][2]) + "\n")
                output_stream.write("number of inactive:  " + str(pair_active_inactive[i][3]) + "\n")
                output_stream.write(22*"-" + "\n")
        elif  pair_active_inactive[i][2] < pair_active_inactive[i][3]:
            with open(output_directory + "/inactives.txt",
                      "a", encoding="utf-8") as output_stream:
                output_stream.write("\n" + str([pair_active_inactive[i][0], pair_active_inactive[i][1]]) + "\n")
                output_stream.write("number of active:    " + str(pair_active_inactive[i][2]) + "\n")
                output_stream.write("number of inactive:  " + str(pair_active_inactive[i][3]) + "\n")
                output_stream.write(22*"-" + "\n")
        elif (pair_active_inactive[i][2] == 0) and (pair_active_inactive[i][3] == 0):
            with open(output_directory + "/zeros.txt",
                      "a", encoding="utf-8") as output_stream:
                output_stream.write("\n" + str([pair_active_inactive[i][0], pair_active_inactive[i][1]]) + "\n")
                output_stream.write("number of active:    " + str(pair_active_inactive[i][2]) + "\n")
                output_stream.write("number of inactive:  " + str(pair_active_inactive[i][3]) + "\n")
                output_stream.write(22*"-" + "\n")
        else:
            with open(output_directory + "/same.txt",
                        "a", encoding="utf-8") as output_stream:
                output_stream.write("\n" + str([pair_active_inactive[i][0], pair_active_inactive[i][1]]) + "\n")
                output_stream.write("number of active:    " + str(pair_active_inactive[i][2]) + "\n")
                output_stream.write("number of inactive:  " + str(pair_active_inactive[i][3]) + "\n")
                output_stream.write(22*"-" + "\n")

    print("number of all actives:  " + str(all_actives))
    print("number of all inactives:  " + str(all_inactives))


def _prepare_files(output_directory: str):
    with open(output_directory + "/actives.txt",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/inactives.txt",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/zeros.txt",
              "w", encoding="utf-8"):
        pass
    with open(output_directory + "/same.txt",
              "w", encoding="utf-8"):
        pass


def read_activity(input_activity: str) -> list:
    with open(input_activity, "r", encoding="utf-8") as stream:
        for line in stream:
            activity = json.loads(line)
            return activity["activity"]


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="analysis of a lot of results")
    parser.add_argument("-a", type=str, dest="analysis_file", required=True)
    parser.add_argument("-c", type=str, dest="activity_file", help="test activity file",
                        required=True)
    parser.add_argument("-f", type=str, dest="test_fragments", required=True)
    parser.add_argument("-o", type=str, dest="output_directory",
                        help="output json lines file", required=True)
    return vars(parser.parse_args())


if __name__  == "__main__":
    _main()
