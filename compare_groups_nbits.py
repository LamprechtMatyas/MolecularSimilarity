#!/usr/bin/env python3
"""
"""
import argparse
import json
import copy


def _main():
    configuration = _read_configuration()
    indexes = _read_indexes(configuration["input_file"])
    _make_all(indexes, configuration["nbits"], configuration["input_groups"], configuration["num"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="test hashed fingerprints")
    parser.add_argument("-f", type=str, dest="input_file", help="input file with"
                        "fragments from active molecules", required=True)
    parser.add_argument("-n", type=str, dest="nbits", help="nbits tha should be used to hash",
                        required=True)
    parser.add_argument("-m", type=str, dest="num", help="number of first groups",
                        required=True)
    parser.add_argument("-p", type=str, dest="input_groups", help="file with groups",
                        required=True)
    configuration = vars((parser.parse_args()))
    parsed = []
    for item in configuration["nbits"].split(","):
        parsed.append(int(item))
    configuration["nbits"] = parsed
    return configuration


def _read_indexes(input_file: str) -> list:
    indexes = []
    with open(input_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            indexes_molecule = []
            for fragment in line["fragments"]:
                if fragment["index"] not in indexes_molecule:
                    indexes_molecule.append(fragment["index"])
            indexes.append(indexes_molecule)
    return(indexes)    
            

def _make_all(indexes_copy: list, nbits: list, input_groups_file: str, num_of_groups: int):   
    num_of_iter = 0 
    num_of_groups = int(num_of_groups) * len(nbits)
    with open(input_groups_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            groups_file = line["groups"] 
            print("-------------------------------")
            print(groups_file)
            for nbit in nbits:
                indexes = copy.deepcopy(indexes_copy)
                for i in range(len(indexes)):
                    for j in range(len(indexes[i])):
                        indexes[i][j] = indexes[i][j] % nbit
                groups = _sort_and_make_groups(indexes, copy.deepcopy(indexes_copy))
                groups1 = []
                for group in groups:
                    if group not in groups1:
                        groups1.append(group)
                groups = groups1
                num = len(groups)
                count = 0
                num_groups_12 = 0
                   
                for group in groups_file: 
                    num_groups_12 += 1                   
                    for new_groups in groups:
                        founded_all = True
                        for item in group:
                            if item not in new_groups:
                                 founded_all = False
                        if founded_all:
                            count += 1
                            break 
                
                print(nbit)            
                print("Founded groups: " + str(count))
                num_of_iter += 1
                if int(num_of_iter) == int(num_of_groups):
                    print("-------------------------------")
                    exit()
      
                    
def _sort_and_make_groups(nbit_indexes: list, indexes: list) -> list:
    for num_mol, molecule_indexes in enumerate(nbit_indexes):
        for i in range(len(molecule_indexes) - 1):
            for j in range(len(molecule_indexes) - i - 1):
                if molecule_indexes[j] > molecule_indexes[j+1]:
                    tmp = molecule_indexes[j]
                    molecule_indexes[j] = molecule_indexes[j+1]
                    molecule_indexes[j+1] = tmp
                    tmp = indexes[num_mol][j]
                    indexes[num_mol][j] = indexes[num_mol][j+1]
                    indexes[num_mol][j+1] = tmp
    groups = []
    for num_mol, molecule in enumerate(nbit_indexes):
        for i in range(len(molecule)):
            group_index = [i]
            while (i+1 < len(molecule)) and (molecule[i] == molecule[i+1]):
                group_index.append(i+1)
                i += 1
            group = []
            if len(group_index) > 1:
                for index in group_index:
                    group.append(indexes[num_mol][index])
                groups.append(group)
    return groups
    
    
if __name__ == "__main__":
    _main()