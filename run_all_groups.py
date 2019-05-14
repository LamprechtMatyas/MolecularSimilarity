#!/usr/bin/env python3
""""
Generates all possible groups from indexes and for all groups it computes the model, score, activity
and evaluation.
"""
import argparse
import json
import multiprocessing
import os
from os import listdir
from os.path import isfile, join

import inputoutput_utils
import model_factory
import add_activity
import compute_evaluation


def _main():
    cpu_counts = multiprocessing.cpu_count()
    configuration = _read_configuration()
    if int(configuration["num_cores"]) > cpu_counts:
        print("You have only " + str(cpu_counts) + " cores")
        exit(1)
    else:
        cpu_counts = int(configuration["num_cores"])   
    maximal_num = _make_configuration_files(configuration["groups"],
                                       configuration["output_directory"], configuration["model"],
                                       cpu_counts, configuration["cutoff"])

    for i in range(cpu_counts):
        process = multiprocessing.Process(target=_model_and_score_and_evaluate,
                                          args=(configuration["active_fragments"],
                                          configuration["test_fragments"], configuration["test_activity"],
                                          i, configuration["output_directory"], maximal_num))
        process.start()


def _read_configuration():
    parser = argparse.ArgumentParser(description="making groups from groups"
                                                 "See file header for more details.")
    parser.add_argument("-i", type=str, dest="active_fragments",
                        help="file with fragments from active molecules", required=True)
    parser.add_argument("-g", type=str, dest="groups",
                        help="file with fragments that we want to group", required=True)
    parser.add_argument("-t", type=str, dest="test_fragments",
                        help="file with fragments from test molecules", required=True)
    parser.add_argument("-a", type=str, dest="test_activity", 
                        help="activity of test molecules", required=True)
    parser.add_argument("-m", type=str, dest="model",
                        help="name of a model", required=True)
    parser.add_argument("-n", type=str, dest="num_cores",
                        help="number of cores that you want to use for program running", required=True)
    parser.add_argument("-c", type=str, dest="cutoff",
                        help="cutoff option", required=False)
    parser.add_argument("-o", type=str, dest="output_directory",
                        help="directory where to store results", required=True)

    configuration = vars(parser.parse_args())
    if configuration["cutoff"] is not None:
        configuration["cutoff"] = int(configuration["cutoff"])
        if (configuration["cutoff"] < 0) or (configuration["cutoff"] > 100):
            print("cutoff has to have value between 0 and 100!!")
            exit(1)
    else:
         configuration["cutoff"] = -1  
    return configuration


def _make_configuration_files(group_file: str, output_directory: str, model_name: str,
                              cpu_counts: int, cutoff_val: int) -> int:
    groups = []
    inputoutput_utils.create_parent_directory(output_directory + "/configurationfiles/0")
    inputoutput_utils.create_parent_directory(output_directory + "/evaluations/0")
    files = [f for f in listdir(output_directory + "/configurationfiles") if
                 isfile(join(output_directory + "/configurationfiles", f))]
    num_of_config = 0
    for file in files:
        with open(output_directory + "/configurationfiles/" + file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                num_of_config += 1
    config_files = []
    for file in files:
        first_part = file.split("_")[0]
        config_files.append(int(first_part[13:]))
    if num_of_config == 0:
        maximal_num = 0
    else:
        maximal_num = max(config_files)
    evaluation_files = [f for f in listdir(output_directory + "/evaluations") if
                        isfile(join(output_directory + "/evaluations", f))]
    if len(evaluation_files) != num_of_config:   
        num_of_max_num = 0
        for item in config_files:
            if item == maximal_num:
                num_of_max_num += 1
        if num_of_max_num != cpu_counts:
            print("Please run the program as before on " + str(num_of_max_num) + " cores")
            exit(1)
        else:
            return maximal_num
    
    with open(group_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            groups.append(line["groups"])

    num_of_groups = len(groups)
    group_list = []
    for i in range(len(groups)-1):
        for j in range(i+1, len(groups)):
            groups1 = groups[i].copy()
            groups2 = groups[j].copy()
            is_intersected = False
            while _control_intersection(groups1, groups2):
                fin = False
                is_intersected = True
                for group1 in groups1:
                    for item in group1:
                        for group2 in groups2:
                            if item in group2:
                                groups1.remove(group1)
                                groups2.remove(group2)
                                groups1.append(_new_group(group1, group2))
                                fin = True
                                break

                        if fin:
                            break
                    if fin:
                        break
            if is_intersected:
                groups11 = _one_group_intersection(groups1)
                new_group = groups11.copy()
                if groups2 != []:
                    new_group.extend(groups2)
                group_list.append(new_group)
            else:
                group_list.append(groups1 + groups2)

    for file in files:
        with open(output_directory + "/configurationfiles/" + file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                grupeto = line["groups"]
                for i in range(len(group_list)):
                    if grupeto == group_list[i]:
                        group_list.remove(grupeto)
                        break
    
    number = len(group_list) // cpu_counts
    ranges = []
    for i in range(cpu_counts):
        ranges.append(i * number)
    ranges.append(len(group_list))
    maximal_num += 1
    for i in range(cpu_counts):
        output_file = output_directory + "/configurationfiles/configuration" + str(maximal_num) + "_" + str(i) + ".json"
        first = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            for j in range(ranges[i], ranges[i+1]):
                if cutoff_val == -1:
                    model = {
                        "model_name": model_name,
                        "groups": group_list[j],
                        "evaluation": "evaluation" + str(maximal_num) + "_" + str(j) + ".json"
                    }
                else:
                    model = {
                        "model_name": model_name,
                        "cutoff": cutoff_val,
                        "groups": group_list[j],
                        "evaluation": "evaluation" + str(maximal_num) + "_" + str(j) + ".json"
                    }
                if first:
                    first = False
                else:
                    output_stream.write("\n")
                json.dump(model, output_stream)

    return maximal_num


def _control_intersection(groups1: list, groups2: list) -> bool:
    for group in groups1:
        for item in group:
            for group2 in groups2:
                if item in group2:
                    return True
    return False


def _make_intersection(groups1: list, groups2: list) -> list:
    new_groups = []
    for group in groups1:
        for item in group:
            founded = False
            for group2 in groups2:
                if item in group2:
                    founded = True

            if founded:
                new_groups.append(_new_group(group, group2))
            else:
                new_groups.append(group)
    return new_groups


def _new_group(group1: list, group2: list) -> list:
    new_group = group1.copy()
    for item in group2:
        if item not in new_group:
            new_group.append(item)
    return new_group


def _is_one_group_intersected(groups: list) -> bool:
    for i in range(len(groups) - 1):
        for j in range(i + 1, len(groups)):
            for item in groups[i]:
                if item in groups[j]:
                    return True
    return False


def _one_group_intersection(groups: list) -> list:
    while _is_one_group_intersected(groups):
        founded = False
        for i in range(len(groups) - 1):
            for j in range(i+1, len(groups)):
                for item in groups[i]:
                    if item in groups[j]:
                        groups.append(_new_group(groups[i], groups[j]))
                        groups.remove(groups[j])
                        groups.remove(groups[i])
                        founded = True
                        break
                if founded:
                    break
            if founded:
                break
    return groups
    

def _model_and_score_and_evaluate(active_fragments: str, test_fragments: str, test_activity: str,
                                  num: int, output_directory: str, maximal_num: int):
    inputoutput_utils.create_parent_directory(output_directory + "/scorefiles/0")
    inputoutput_utils.create_parent_directory(output_directory + "/activities/0")
    with open(output_directory + "/configurationfiles/configuration" + str(maximal_num) + "_" + str(num) + ".json", "r",
              encoding="utf-8") as input_file:
        for new_line in input_file:
            line = json.loads(new_line)
            if os.path.isfile(output_directory + "/evaluations/" + line["evaluation"]):
                continue 
            
            new_model = model_factory.create_model(line["model_name"])
            model = new_model.create_model(active_fragments, "", "", "",
                                           line)
            new_model.score_model(model, test_fragments, "",
                                  output_directory + "/scorefiles/score" + line["evaluation"])

            # run add_activity
            activity = add_activity.read_activity(test_activity)
            add_activity.add_activity_and_write_to_json(output_directory + "/scorefiles/score" + line["evaluation"],
                                                        activity,
                                                        output_directory + "/activities/activity" + line["evaluation"])

            # run compute_evaluation
            score_act = compute_evaluation.read_file_with_score_and_activity(output_directory + "/activities/activity"
                                                                            + line["evaluation"])
            activity = compute_evaluation.sort_activity(score_act)
            compute_evaluation.evaluation(activity, output_directory + "/evaluations/" + line["evaluation"])


if __name__ == "__main__":
    _main()
