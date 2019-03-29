#!/usr/bin/env python3
""""
Library for more used functions.
"""

import argparse


def read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="ecfp/fcfp analysis")
    parser.add_argument("-f", type=str, dest="first_nbit", help="write number",
                        required=True)
    parser.add_argument("-s", type=str, dest="second_nbit", help="write number",
                        required=True)
    parser.add_argument("-i", type=str, dest="input_file", help="input smi file with molecules",
                        required=True)
    parser.add_argument("-o", type=str, dest="output_file", help="output txt file",
                        required=True)
    return vars(parser.parse_args())


def read_molecules(file: str) -> list:
    molecules = []
    with open(file, mode="r", encoding="utf-8") as input_sream:
        for new_line in input_sream:
            atributes = new_line.split("\t")
            molecules.append(atributes[0])
    return molecules


def equivalence_class(output: list) -> list:
    last = -1
    num = -1
    list_of_equivalence = []
    for i in range(len(output)):
        if last == int(output[i]["index"]):
            list_of_equivalence[num].append(output[i]["smiles"])
        else:
            list_of_equivalence.append([output[i]["index"]])
            num += 1
            list_of_equivalence[num].append(output[i]["smiles"])
            last = output[i]["index"]
    return list_of_equivalence


def same_equivalence_class(clekv1: list, clekv2: list) -> list:
    same = []
    num = -1
    for item1 in clekv1:
        for item2 in clekv2:
            if item1[1:] == item2[1:]:
                same.append([item1[0]])
                num += 1
                same[num].append(item2[0])
                same[num].append(item1[1:])
                clekv2.remove(item2)
                break
    return same


def difference_of_equivalence(ekv1: list, intersection: list, position: int) -> list:
    difference = []
    for item in ekv1:
        founded = False
        for intersect in intersection:
            if item[0] == intersect[position]:
                founded = True
                break
        if founded is False:
            difference.append(item)
    return difference

