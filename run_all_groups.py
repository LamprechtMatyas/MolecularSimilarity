#!/usr/bin/env python3
""""
Generates all possible groups from indexes and for all groups it computes the model, score, activity
and evaluation.
"""
import argparse
import json
import multiprocessing

import inputoutput_utils
import model_factory
import add_activity
import compute_evaluation


def _main():
    cpu_counts = multiprocessing.cpu_count()
    configuration = _read_configuration()
    ranges = _make_configuration_files(configuration["groups"],
                                    configuration["output_directory"], configuration["model"], cpu_counts)

    for i in range(cpu_counts):
        process = multiprocessing.Process(target=_model_and_score_and_evaluate,
                                          args=(configuration["active_fragments"],
                                          configuration["test_fragments"], configuration["test_activity"],
                                          i, configuration["output_directory"], ranges))
        process.start()


def _read_configuration():
    parser = argparse.ArgumentParser(description="model evaluation "
                                                 "See file header for more details.")
    parser.add_argument("-i", type=str, dest="active_fragments", required=True)
    parser.add_argument("-g", type=str, dest="groups",
                        help="file with fragments", required=True)
    parser.add_argument("-t", type=str, dest="test_fragments", required=True)
    parser.add_argument("-a", type=str, dest="test_activity", required=True)
    parser.add_argument("-m", type=str, dest="model", required=True)
    parser.add_argument("-o", type=str, dest="output_directory", required=True)

    configuration = vars(parser.parse_args())
    return configuration


def _make_configuration_files(group_file: str, output_directory: str, model_name: str, cpu_counts) -> list:
    groups = []
    inputoutput_utils.create_parent_directory(output_directory + "/configurationfiles/0")
    with open(group_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            groups.append(line["groups"])

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
                    new_group.append(groups2)
                group_list.append(new_group)
            else:
                group_list.append(groups1 + groups2)

    number = len(group_list) // cpu_counts
    ranges = []
    for i in range(cpu_counts):
        ranges.append(i * number)
    ranges.append(len(group_list))

    for i in range(cpu_counts):
        output_file = output_directory + "/configurationfiles/configuration" + str(i) + ".json"
        first = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            for j in range(ranges[i], ranges[i+1]):
                model = {
                    "model_name": model_name,
                    "groups": group_list[j]
                }
                if first:
                    first = False
                else:
                    output_stream.write("\n")
                json.dump(model, output_stream)

    return ranges


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
                                  num: int, output_directory: str, ranges: list):
    inputoutput_utils.create_parent_directory(output_directory + "/scorefiles/0")
    inputoutput_utils.create_parent_directory(output_directory + "/activities/0")
    inputoutput_utils.create_parent_directory(output_directory + "/evaluations/0")
    count = ranges[num]
    with open(output_directory + "/configurationfiles/configuration" + str(num) + ".json", "r",
              encoding="utf-8") as input_file:
        for new_line in input_file:
            line = json.loads(new_line)
            new_model = model_factory.create_model(line["model_name"])
            model = new_model.create_model(active_fragments, "", "", "",
                                           line)
            new_model.score_model(model, test_fragments, "",
                                  output_directory + "/scorefiles/score" + str(count) + ".json")

            # run add_activity
            activity = add_activity.read_activity(test_activity)
            add_activity.add_activity_and_write_to_json(output_directory + "/scorefiles/score" + str(count) + ".json",
                                                        activity,
                                                        output_directory + "/activities/activity" + str(count) + ".json")

            # run compute_evaluation
            score_act = compute_evaluation.read_file_with_score_and_activity(output_directory + "/activities/activity"
                                                                            + str(count) + ".json")
            activity = compute_evaluation.sort_activity(score_act)
            compute_evaluation.evaluation(activity, output_directory + "/evaluations/evaluation" + str(count) + ".json")
            count += 1


if __name__ == "__main__":
    _main()
