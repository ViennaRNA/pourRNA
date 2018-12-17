#!/usr/bin/python3
import sys
import os
import argparse
from min_dijkstra import *
from operator import itemgetter
#from scipy.sparse.csgraph._shortest_path import shortest_path
import math
import re
import matplotlib.pyplot as plt
import RNA
from math import ceil

"""
format: id_from, structure from, energy from, id_to, str. to, energy to, saddle str, saddle energy. 
"""
def read_saddle_csv_file(file):
    paths_dictionary = {}
    struct_id_energy_dict = {}
    struct_id_structure_dict = {}
    min_saddle = 1000000000
    with open(file, 'r') as f:
        new_index = 0
        for line in f:
            #print(line)
            match = re.match("\s*(\d+),\s*([\.\(\)]+),\s*(\-?\d*\.?\d*)," + \
                             "\s*(\d+),\s*([\.\(\)]+),\s*(\-?\d*\.?\d*)," + \
                             "\s*([\.\(\)\~]+),\s*(\-?\d*\.?\d*)\s*", line)
            if match:
                #print("match")
                #structure_count +=1
                #if filter.MaxMinima != None and structure_count > filter.MaxMinima:
                #    break
                
                index_from = match.group(1)
                structure_from = match.group(2)
                energy_from = match.group(3)
                #energy_from = float(energy_from)
                
                index_to = match.group(4)
                structure_to = match.group(5)
                energy_to = match.group(6)
                #energy_to = float(energy_to)
                
                structure_saddle = match.group(7)
                saddle_energy_kcal = match.group(8)
                saddle_energy_kcal = float(saddle_energy_kcal)

                try:
                    index_from = int(index_from)
                    index_to = int(index_to)
                except Exception as e:
                    print("Error: could not parse index from input file!")
                    exit()
                
                if not index_from in struct_id_energy_dict:
                    struct_id_energy_dict[index_from] = energy_from
                if not index_to in struct_id_energy_dict:
                    struct_id_energy_dict[index_to] = energy_to

                if not index_from in paths_dictionary:
                    paths_dictionary[index_from] = {}
                if not index_to in paths_dictionary:
                    paths_dictionary[index_to] = {}
                
                paths_dictionary[index_from][index_to] = saddle_energy_kcal
                paths_dictionary[index_to][index_from] = saddle_energy_kcal
                if saddle_energy_kcal < min_saddle:
                    min_saddle = saddle_energy_kcal

                struct_id_structure_dict[(index_to,index_to)] = structure_to
                struct_id_structure_dict[(index_from,index_from)] = structure_from
                struct_id_structure_dict[(index_from,index_to)] = structure_saddle
                struct_id_structure_dict[(index_to,index_from)] = structure_saddle
                
    return paths_dictionary, min_saddle, struct_id_energy_dict, struct_id_structure_dict

if __name__ == "__main__":
    
    # arguments:
    parser = argparse.ArgumentParser(description='myBot', conflict_handler='resolve')
    parser.add_argument("-f", "--file", action="store", type=str, required=True, help="csv file with saddles")
    parser.add_argument("-i", "--index_from", action="store", type=int, required=False, help="Index from.")
    parser.add_argument("-j", "--index_to", action="store", type=int, required=False, help="Index to.")
    parser.add_argument("-d", "--structure_index_file", type=str, required=False, help="File with two structures, for which the min max saddle is computed.")
    parser.add_argument("-p", "--plot", action='store_true', required=False, help="Create Path Plot.")
    parser.add_argument("-z", "--bpDist_to_final_str", action='store_true', required=False, help="Distance to ground state, else distance to next structure.")
    args = parser.parse_args()
    paths_graph, min_saddle, struct_id_energy_dict, struct_id_structure_dict = read_saddle_csv_file(args.file)
    
    print(paths_graph)
    
    ifrom = args.index_from
    ito = args.index_to
    if args.structure_index_file and ifrom == None and ito == None:
        with open(args.structure_index_file, 'r') as f:
            lines = f.readlines()
            s1 = lines[2].strip()
            s2 = lines[3].strip()
            with open(args.file, 'r') as f:
                for line in f:
                    if ifrom == None and s1 in line:
                        parts = line.split(',')
                        for i,p in enumerate(parts):
                            if s1 in p:
                                ifrom = int(parts[i-1])
                                break
                    if ito == None and s2 in line:
                        parts = line.split(',')
                        for i,p in enumerate(parts):
                            if s2 in p:
                                ito = int(parts[i-1])
                                break
                    if ito != None and ifrom != None:
                        break
                    
    for k in paths_graph.keys():
        for v in paths_graph[k].keys():
            paths_graph[k][v] = paths_graph[k][v] - min_saddle + 1
    
    
    print(ifrom, ito)
    sp, distances = shortestPath(paths_graph, ifrom, ito)
    #sp.reverse()
    #distances.reverse()
    print(sp)
    print(min_saddle)
    """
    #print([ distances[x] + min_saddle - 1 for x in sp])
    #sp.reverse()
    path_from_saddle_to = []
    for i in range(len(sp)):
        id_from = sp[i]
        e_from = struct_id_energy_dict[id_from]
        path_from_saddle_to.append(float(e_from))
        if (i+1 < len(sp)):
            id_to = sp[i+1]
            e_saddle = paths_graph[id_from][id_to] + min_saddle - 1
            path_from_saddle_to.append(float(e_saddle))
    print(path_from_saddle_to)
    """ 
    bn = os.path.basename(args.file)
    bn = os.path.splitext(bn)[0]
    ### print path with energies ####
    text_file_name = bn+"_refolding_path.txt"
    f_path = open(text_file_name,'w')
    #RNA.bp_distance("....","(())")
    path_from_saddle_to = []
    for i in range(len(sp)):
        id_from = sp[i]
        e_from = struct_id_energy_dict[id_from]
        structure_from = struct_id_structure_dict[(id_from,id_from)]
        path_from_saddle_to.append([structure_from, float(e_from)])
        f_path.write(structure_from +" {:.2f}".format(float(e_from)) +"\n")
        if (i+1 < len(sp)):
            id_to = sp[i+1]
            e_saddle = paths_graph[id_from][id_to] + min_saddle - 1
            structure_saddle = struct_id_structure_dict[(id_from,id_to)]
            path_from_saddle_to.append([structure_saddle, float(e_saddle)])
            f_path.write(structure_saddle +" {:.2f}".format(float(e_saddle)) +"\n")
    #print(path_from_saddle_to)
    f_path.close()

    if (args.plot):
        path_indices = [0]
        path_energies = [ path_from_saddle_to[x][1] for x in range(0, len(path_from_saddle_to))]
        if(args.bpDist_to_final_str):
           final_str = path_from_saddle_to[-1][0]
           path_indices = [ RNA.bp_distance(path_from_saddle_to[x][0], final_str) for x in range(0, len(path_from_saddle_to))]
        else: # bp-dist to next structure
           path_indices = [0]
           total_bp_dist = 0
           for x in range(0, len(path_from_saddle_to)):
               if(x+1 < len(path_from_saddle_to)):
                   x_next = x+1
                   bp_dist = RNA.bp_distance(path_from_saddle_to[x][0], path_from_saddle_to[x+1][0])
                   total_bp_dist += bp_dist
                   path_indices.append(total_bp_dist)

        all_lines = []
        lines, = plt.plot(path_indices, path_energies, label="Structure Energies", linewidth=1, color = "black", marker=".")
        all_lines.append(lines)
        path_labels = []
        for x in path_indices[1:-1]:
            if x % 200 == 0:
                path_labels.append(x)
            else:
                path_labels.append('')
        #path_labels = [path_indices[0]] + path_labels + [path_indices[-1]]
        path_labels = [path_indices[0]] + [path_indices[-1]]
        path_indices = [path_indices[0]] + [path_indices[-1]]
        plt.xticks(path_indices, path_labels)
        e_min = ceil(min(path_energies))
        while(e_min % 5 != 0):
            e_min -= 1
        e_min = int(e_min)
        e_max = int(-55)
        e_ticks = [ x for x in range(e_min, e_max+1, 5)]
        plt.yticks(e_ticks)
        #leg = plt.legend(handles=all_lines, fontsize = 'medium')
        #leg_texts = leg.get_texts()
    
        #plt.setp(lines, color='r', linewidth=2.0)
        #ax = plt.axes()
        #ax.set_xscale("log", nonposx='clip')
        #ax.set_yscale("log")
        plt.xlabel('Structures along Optimal Path')
        plt.ylabel('Energy [kcal/mol]') 
        #plt.grid(True)
        res_file=bn+ "_refolding_path.pdf"
        plt.savefig(res_file)
        plt.clf()
        
    
    



