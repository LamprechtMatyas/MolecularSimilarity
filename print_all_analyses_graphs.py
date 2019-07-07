#!/usr/bin/env python3
""""
Prints boxplots, histograms, tables for more evaluation files and groups.
"""

import argparse
import json

import matplotlib.pyplot as plt

import inputoutput_utils

def _main():
    configuration = _read_configuration()
    _print_histogram(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])
    _fragments_distribution(configuration["input_files"], configuration["output_directory"])
    _print_boxplots(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])
    _print_number_of_fragments(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])
    _print_highest_auc_each(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])
    _table_of_highest_auc(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])
    _highest_groups_auc(configuration["input_files"], configuration["nicknames"], configuration["output_directory"])


def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(description="prints graphs and tables more for"
                                                 " groups and pairs evaluations")
    parser.add_argument("-i", type=str, dest="input_files", help="comma separated input files",
                        required=True)
    parser.add_argument("-n", type=str, dest="nicknames", help="comma sepapated nicknames",
                         required=True)
    parser.add_argument("-o", type=str, dest="output_directory",help="output directory",
                         required=True)
    configuration = vars(parser.parse_args())
    configuration["input_files"] = configuration["input_files"].split(",")
    configuration["nicknames"] = configuration["nicknames"].split(",")
    if len(configuration["input_files"]) != len(configuration["nicknames"]):
        print("Number of input files must be equal to the number of nicknames!")
        exit(1)
    return configuration


def _print_histogram(input_files: list, nicknames: list, output_directory: str):
    all_auc = []
    inputoutput_utils.create_parent_directory(output_directory+"/0")
    for file in input_files:
        auc = []
        with open(file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                auc.append(line["AUC"])
        all_auc.append(auc)
    plt.hist(all_auc, stacked=True, density=False)
    plt.legend(nicknames, loc="upper right")
    plt.savefig(output_directory + "/AUC.png", dpi=1000)
    plt.figure()
    

def _fragments_distribution(input_files: list, output_directory: str):
    fragments = {}
    i = 0
    for file in input_files:
        i += 1
        with open(file, "r", encoding="utf-8") as input_file:
            for new_line in input_file:
                line = json.loads(new_line)
                for group in line["groups"]:
                    for item in group:
                        if item in fragments:
                            fragments[item] += 1
                        else:
                            fragments[item] = 1
    numbers = []
    for item in fragments:
        numbers.append(fragments[item])
    plt.hist(numbers)
    plt.savefig(output_directory + "/fragments.png", dpi=1000)
    plt.figure()


def _print_boxplots(input_files: list, nicknames: list, output_directory: str):
    indexes = {}
    index_auc = []
    num = 0
    for num_file,file in enumerate(input_files):
        with open(file, "r", encoding="utf-8") as input_stream:
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
                            for i in range(len(nicknames)):
                                index_auc[len(index_auc) - 1].append([])
                            index_auc[len(index_auc) - 1][num_file+1].append(line["AUC"])
                        else:
                            index_auc[number][num_file+1].append(line["AUC"])
 
    xaxis = []
    for i in range(len(nicknames)):
        xaxis.append(i+2)
    for index in index_auc:
        name = index[0]
        plt.boxplot(index[1:], labels=nicknames) 
        plt.ylabel("AUC")
        plt.tight_layout() 
        plt.savefig(output_directory + "/" +str(name), dpi=1000)
        plt.figure()
        plt.close("all")
        
        
def _print_number_of_fragments(input_files: list, nicknames: list, output_directory: list):
    indexes = {}
    index_num = []
    num = 0
    for num_file,file in enumerate(input_files):
        with open(file, "r", encoding="utf-8") as input_stream:
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
                            index_num.append([item])
                            for i in range(len(nicknames) + 1):
                                index_num[len(index_num) - 1].append(0)
                            index_num[len(index_num) - 1][num_file+1] = 1
                        else:
                            index_num[number][num_file+1] += 1
    all = 0
    for item in index_num:
        item[len(item) - 1] = sum(item[1:])
    sorted_index_num = _sort_by_max_num(index_num)
    nrows, ncols = 15, 3
    hcell, wcell = 0.23, 1.
    hpad, wpad = 0, 0 
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    col = ["Index"]
    
    for item in nicknames:
        col.append(item)
    col.append("Celkem")
    colLabels = tuple(col)   
    number_totally = sorted_index_num[:15]    
    the_table = ax.table(cellText=number_totally, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/index_number_totally.png", dpi=1000)
    plt.figure()
       

def _sort_by_max_num(index_num: list) -> list:
    vel = len(index_num[0])
    for i in range(len(index_num) - 1):
        for j in range(len(index_num) - i - 1):
            if index_num[j+1][vel-1] > index_num[j][vel-1]:
                tmp = index_num[j]
                index_num[j] = index_num[j+1]
                index_num[j+1] = tmp
    return index_num


def _print_highest_auc_each(input_files: list, nicknames: list, output_directory: list):
    highest = [0 for i in range(len(input_files))]
    for num_file,file in enumerate(input_files):
        with open(file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                if line["AUC"] > highest[num_file]:
                    highest[num_file] = line["AUC"]
    plt.plot(nicknames, highest, marker="o")
    plt.savefig(output_directory + "/highest_auc.png", dpi=1000)
    plt.figure()
    
 
def _table_of_highest_auc(input_files: list, nicknames: list, output_directory: list):
    auc_nick = []
    for num_file,file in enumerate(input_files):
        with open(file, "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
               line = json.loads(new_line)
               auc_nick.append([line["AUC"], nicknames[num_file]])
    sorted_auc_nick = _sort_by_auc(auc_nick)
    nrows, ncols = 15, 3
    hcell, wcell = 0.23, 1.
    hpad, wpad = 0, 0 
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels = ["AUC", "Nicknames"]
  
    number_totally = sorted_auc_nick[:15]    
    the_table = ax.table(cellText=number_totally, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/table_highest_auc.png", dpi=1000)
    plt.figure()
    
                    
def _sort_by_auc(auc_nick: list) -> list:
    for i in range(len(auc_nick) - 1):
        for j in range(len(auc_nick) - i - 1):
            if auc_nick[j][0] < auc_nick[j+1][0]:
                tmp = auc_nick[j]
                auc_nick[j] = auc_nick[j+1]
                auc_nick[j+1] = tmp
    return auc_nick


def _highest_groups_auc(input_files: list, nicknames: list, output_directory: list):
    indexes = {}
    index_group_auc = []
    num = 0
    for num_file,file in enumerate(input_files):
        with open(file, "r", encoding="utf-8") as input_stream:
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
                            index_group_auc.append([item])
                            for i in range(2):
                                index_group_auc[len(index_group_auc) - 1].append(0)
                            index_group_auc[len(index_group_auc) - 1][1] = 1
                            index_group_auc[len(index_group_auc) - 1][2] = line["AUC"]
                        else:
                            index_group_auc[number][1] += 1
                            if line["AUC"] > index_group_auc[number][2]:
                                index_group_auc[number][2] = line["AUC"]
    
    sorted_index_group_auc = _sort(index_group_auc, 1)
    nrows, ncols = 15, 3
    hcell, wcell = 0.23, 1.
    hpad, wpad = 0, 0 
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels = ["Index", "Number of groups", "Highest AUC"]
  
    number_totally = sorted_index_group_auc[:15]
    for item in number_totally:
         item[2] =  float(round(item[2],5))   
    the_table = ax.table(cellText=number_totally, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/table_groups_auc.png", dpi=1000)
    plt.figure()
    
    sorted_index_group_auc = _sort(index_group_auc, 2)
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels = ["Index", "Highest AUC", "Number of groups"]
  
    number_totally = sorted_index_group_auc[:15]
    for item in number_totally:
         tmp = float(round(item[2],5))  
         item[2] = item[1]
         item[1] = tmp  
    the_table = ax.table(cellText=number_totally, colLabels=colLabels, cellLoc="center", loc="center")
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)
    plt.savefig(output_directory + "/table_auc_groups.png", dpi=1000)
    plt.figure()


def _sort(index_group_auc: list, item: int) -> list:
    for i in range(len(index_group_auc) - 1):
        for j in range(len(index_group_auc) - i - 1):
            if index_group_auc[j][item] < index_group_auc[j+1][item]:
                tmp = index_group_auc[j]
                index_group_auc[j] = index_group_auc[j+1]
                index_group_auc[j+1] = tmp
    return index_group_auc
    

if __name__ == "__main__":
    _main()
                
