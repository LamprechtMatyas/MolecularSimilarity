#!/usr/bin/env python3
"""
We use hashing ecfp and we take pairs that improved performance and we compare
how many of pairs (that improved performance) are in hashed ecfp.
"""
import argparse
import json

from rdkit.Chem import AllChem
import rdkit


def _main():
    configuration = _read_configuration()
    make_all_things(configuration)


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="test hashed "
                                                 "See file header for more details.")
    parser.add_argument("-a", type=str, dest="input_file", help="smi file with active molecules",
                        required=True)
    parser.add_argument("-n", type=str, dest="nbits",
                        help="nbits to hash", required=True)
    parser.add_argument("-p", type=str, dest="input_pairs",
                        help="output file with pairs that improved the performance", required=True)
    configuration = vars((parser.parse_args()))
    return configuration
    
    
def make_all_things(configuration: dict):
    couples = []
    with open(configuration["input_file"], "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = new_line.strip()
            molecule = rdkit.Chem.MolFromSmiles(line)
            indexes = _extract_neighbourhood_fragments(molecule, 6, int(configuration["nbits"]))
            sorted_indexes = _sort_2d_list(indexes)
            groups = (_make_groups(sorted_indexes))
            couples.extend(_make_couples(groups[0]))
    num = 0
    with open(configuration["input_pairs"], "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            pair = line["pair"]
            if pair in couples:
                num += 1
    print(configuration["nbits"])
    print("Number of pairs included in input pairs: " + str(num))
    print("Number of all pairs after hash: " + str(len(couples)))


# Extract and return circular fragments.
def _extract_neighbourhood_fragments(molecule: object, diameter: int, nbit: int) -> list:
    info = {}
    size = diameter // 2
    AllChem.GetMorganFingerprint(molecule, size, bitInfo=info)
    indexes = [[], []]
    for num, element in enumerate(info):
        indexes[0].append(element)
        indexes[1].append(element % nbit)
    return indexes


def _sort_2d_list(indexes: list) -> list:
    for i in range(len(indexes[1]) - 1):
        for j in range(len(indexes[1]) - i - 1):
            if indexes[1][j] > indexes[1][j+1]:
                tmp = indexes[1][j]
                indexes[1][j] = indexes[1][j+1]
                indexes[1][j+1] = tmp
                tmp = indexes[0][j]
                indexes[0][j] = indexes[0][j + 1]
                indexes[0][j + 1] = tmp
    return indexes


def _make_groups(indexes: list) -> list:
    i = 0
    group_indexes = [[], []]

    while i < len(indexes[1]):
        new_list = [indexes[1][i]]
        new_list1 = [indexes[0][i]]
        item = indexes[1][i]
        i += 1
        while (i < len(indexes[1])) and (item == indexes[1][i]):
            new_list.append(indexes[1][i])
            new_list1.append(indexes[0][i])
            i += 1
        group_indexes[0].append(new_list1)
        group_indexes[1].append(new_list)
    return group_indexes


def _make_couples(indexes: list) -> list:
    couples = []
    for item in indexes:
        if len(item) > 1:
            for i in range(len(item) - 1):
                for j in range(i+1, len(item)):
                    couples.append([item[i], item[j]])
    return couples


if __name__ == "__main__":
    _main()
