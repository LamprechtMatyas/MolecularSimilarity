#!/usr/bin/env python3
""""
Generates all possible pairs from indexes and for all pairs it computes the model, score, activity
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
    ranges = _make_configuration_files(configuration["active_fragments"],
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
    parser.add_argument("-i", type=str, dest="active_fragments",
                        help="file with fragments", required=True)
    parser.add_argument("-t", type=str, dest="test_fragments", required=True)
    parser.add_argument("-a", type=str, dest="test_activity", required=True)
    parser.add_argument("-m", type=str, dest="model", required=True)
    parser.add_argument("-o", type=str, dest="output_directory", required=True)

    configuration = vars(parser.parse_args())
    return configuration


def _make_configuration_files(input_file: str, output_directory: str, model_name: str, cpu_counts) -> list:
    active_indexes = []
    inputoutput_utils.create_parent_directory(output_directory + "/configurationfiles/0")
    with open(input_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            for item in line["fragments"]:
                if item["index"] not in active_indexes:
                    active_indexes.append(item["index"])

    pair_list = []
    for i in range(len(active_indexes)-1):
        for j in range(i+1, len(active_indexes)):
            pair_list.append([active_indexes[i], active_indexes[j]])

    number = len(pair_list) // cpu_counts
    ranges = []
    for i in range(cpu_counts):
        ranges.append(i * number)
    ranges.append(len(pair_list))

    for i in range(cpu_counts):
        output_file = output_directory + "/configurationfiles/configuration" + str(i) + ".json"
        first = True
        with open(output_file, "w", encoding="utf-8") as output_stream:
            for j in range(ranges[i], ranges[i+1]):
                model = {
                    "model_name": model_name,
                    "groups": [pair_list[j]]
                }
                if first:
                    first = False
                else:
                    output_stream.write("\n")
                json.dump(model, output_stream)

    return ranges


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
