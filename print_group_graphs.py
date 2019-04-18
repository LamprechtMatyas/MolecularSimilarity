#!/usr/bin/env python3
""""
Prints boxplots, histograms, tables for fragments.
"""

import argparse
import json
import operator

import matplotlib.pyplot as plt

import inputoutput_utils


def _main():
    configuration = _read_configuration()
    val = _get_baseline_value(configuration["baseline_val"])
    auc = _print_histogram(val, configuration["input_file"],
                           configuration["output_directory"])
    _fragments_distribution(configuration["input_file"], configuration["output_directory"])
    _print_boxplot(configuration["input_file"], configuration["output_directory"], auc)
    

def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="prints graphs and table for group and pairs"
                                                 " outputs")
    parser.add_argument("-b", type=str, dest="baseline_val", help="file with baseline value",
                        required=True)
    parser.add_argument("-i", type=str, dest="input_file", help="file from pair_analysis.py"
                                                                " or group_analysis.py", required=True)
    parser.add_argument("-o", type=str, dest="output_directory",
                        help="directory where to store output images", required=True)
    return vars(parser.parse_args())


def _get_baseline_value(file: str) -> float:
    with open(file, "r", encoding="utf-8") as input_stream:
        new_line = input_stream.read()
        line = json.loads(new_line)
        return line["AUC"]
                
                
def _print_histogram(baseline_val: int, input_file: str, output_directory: str) -> list:
    inputoutput_utils.create_parent_directory(output_directory + "/0")
    textstr = "baseline AUC: " + str(float(round(baseline_val, 5)))
    auc = []
    with open(input_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            auc.append(line["AUC"])
    min_val = min(auc)
    max_val = max(auc)
    diff = max_val - min_val
    diff_step = diff/10
    steps = []        
    for i in range(10):
        steps.append(i*diff_step + min_val)
    steps.append(max_val)
    arr = plt.hist(auc, bins=steps, color="blue")
    plt.xticks(steps, rotation=90)
    for i in range(10):
        plt.text(arr[1][i]+diff_step/2, arr[0][i], str(int(arr[0][i])), horizontalalignment="center")
    props = dict(boxstyle="round")
    plt.text(steps[7], int(arr[0][0]), textstr)
    plt.tight_layout()
    plt.savefig(output_directory + "/AUC.png", dpi=1000)
    plt.figure()
    return auc


def _fragments_distribution(input_file: str, output_directory: str):
    fragments = {}
    number_of_input = 0
    with open(input_file, "r", encoding="utf-8") as input_stream:
        for new_line in input_stream:
            line = json.loads(new_line)
            number_of_input += 1
            for group in line["groups"]:
                for item in group:
                    if item in fragments:
                        fragments[item] += 1
                    else:
                        fragments[item] = 1
    numbers = []
    for item in fragments:
        numbers.append(fragments[item])
    max_val = max(numbers)
    min_val = min(numbers)
    diff = max_val - min_val
    diff_step = diff/10
    steps = []        

    arr = plt.hist(numbers, color="blue")
    for i in range(10):
        plt.text(arr[1][i]+diff_step/2, arr[0][i], str(int(arr[0][i])), horizontalalignment="center")
    
    plt.xlabel("Number of fragments occurence")
    plt.ylabel("Number of fragments")
    plt.annotate("Number of groups: " + str(number_of_input), xy=(0.65, 0.95), xycoords='axes fraction')
    plt.tight_layout()
    plt.savefig(output_directory + "/fragments.png")
    plt.figure()
    

def _sort(fragments: dict) -> list:
    fragments_list = [[], []]
    for item in fragments:
        fragments_list[0].append(item)
        fragments_list[1].append(fragments[item])
    for i in range(len(fragments_list[1]) - 1):
        for j in range(len(fragments_list[1]) - i - 1):
            if fragments_list[1][j] < fragments_list[1][j+1]:
                tmp = fragments_list[1][j]
                fragments_list[1][j] = fragments_list[1][j+1]
                fragments_list[1][j+1] = tmp
    return fragments_list     


def _print_boxplot(input_file: str, output_directory: str, auc: list):
    indexes = {}
    index_auc = []
    num = 0
    with open(input_file, "r", encoding="utf-8") as input_stream:
         for new_line in input_stream:
             line = json.loads(new_line)
             for group in line["groups"]:
                 for item in group:
                     number = -1
                     if item in indexes:
                         number = indexes[item]
                     else:
                         indexes[item] = num
                         num += 1
                     if number == -1:
                         index_auc.append([item])
                         index_auc[len(index_auc)-1].append(line["AUC"])
                     else:
                         index_auc[number].append(line["AUC"])
    sorted_index_auc = _sort_by_length(index_auc)
    
    for item in sorted_index_auc:
        plt.figure()
        plt.boxplot(item[1:])
        plt.xticks([])
        plt.savefig(output_directory + "/" + str(item[0]) + ".png")
    table_things = []
    for i in range(15):
        table_things.append([sorted_index_auc[i][0], len(sorted_index_auc[i][1:]),
                             float(round(max(sorted_index_auc[i][1:]),5))])

    indexes_by_length = []
    for i in range(len(sorted_index_auc)):
        indexes_by_length.append(len(sorted_index_auc[i][1:]))
    nrows, ncols = 15, 3
    hcell, wcell = 0.23, 1.
    hpad, wpad = 0, 0 
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Index", "Number of groups", "Highest AUC")
    
    the_table = ax.table(cellText=table_things, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_groups_AUC.png", dpi=1000)   
    plt.figure()
    
    sorted_auc = sorted(auc, reverse=True)
    for i in range(15):
        for j in range(len(sorted_auc)):
            if table_things[i][2] == float(round(sorted_auc[j],5)):
                table_things[i][2] = j+1
                break 
    fig = plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Index", "Number of groups", "Ranking of AUC")
    
    the_table = ax.table(cellText=table_things, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_number_ranking.png", dpi=1000)
    plt.figure()
    
    tabul = []
    for i in range(15):
        tabul.append([table_things[i][0], table_things[i][1]])
    fig = plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Index", "Number of groups")
    
    the_table = ax.table(cellText=tabul, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_groups.png", dpi=1000)
    plt.figure()
    
    
    so_index_auc = _sort_by_max_auc(sorted_index_auc)
    table_things = []
    for i in range(15):
        table_things.append([so_index_auc[i][0], float(round(max(so_index_auc[i][1:]),5)), len(so_index_auc[i][1:])])
    
    fig = plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Index", "Highest AUC", "Number of groups")
    
    the_table = ax.table(cellText=table_things, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_auc_groups.png", dpi=1000)
    plt.figure()
    
    index_auc = []
    with open(input_file, "r", encoding="utf-8") as input_stream:
         for new_line in input_stream:
             line = json.loads(new_line)
             for group in line["groups"]:
                 for item in group:
                     index_auc.append([item])
                     index_auc[len(index_auc) - 1].append(line["AUC"])
    
    sorted_index_auc = _sort_by_auc(index_auc)
    table_st = []
    for i in range(15):
        table_st.append([sorted_index_auc[i][0], float(round(sorted_index_auc[i][1],5))])
    
    fig = plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Index", "AUC")
    
    the_table = ax.table(cellText=table_st, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_auc.png", dpi=1000)
    

def _sort_by_auc(auc: list) -> list:
    for i in range(len(auc) - 1):
        for j in range(len(auc) - i - 1):
            if auc[j][1] < auc[j+1][1]:
                tmp = auc[j]
                auc[j] = auc[j+1]
                auc[j+1] = tmp
    return auc
           
        
def _sort_by_max_auc(index_auc: list) -> list:
    for i in range(len(index_auc) - 1):
        for j in range(len(index_auc) - i - 1):
            if max(index_auc[j][1:]) < max(index_auc[j+1][1:]):
                tmp = index_auc[j]
                index_auc[j] = index_auc[j+1]
                index_auc[j+1] = tmp
    return index_auc
        
            
def _sort_by_length(index_auc: list) -> list:
    for i in range(len(index_auc) - 1):
        for j in range(len(index_auc) - i - 1):
            if len(index_auc[j]) < len(index_auc[j+1]):
                tmp = index_auc[j]
                index_auc[j] = index_auc[j+1]
                index_auc[j+1] = tmp
    return index_auc                     


if __name__ == "__main__":
    _main()

