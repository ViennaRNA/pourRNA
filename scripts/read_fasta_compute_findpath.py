#!/usr/bin/python3.5
import RNA
import argparse
import subprocess

def find_saddle(fold_compound, s1, s2, width):
    path = fold_compound.path_findpath(s1, s2, width)
    #path = RNA.get_path(sequence, s1, s2, width)
    my_path=[]
    try:
        i=0
        while(path[i]):
            #print(path[i].s[:len(sequence)-1], path[i].en)
            structure = str(path[i].s[:len(sequence)])
            energy = float(path[i].en)
            my_path.append((structure, energy))
            i += 1
    except Exception as e:
        #print(e)
        pass
    saddle_energy_dcal = fold_compound.path_findpath_saddle(s1, s2, width)
    saddle_energy_kcal = saddle_energy_dcal / 100.0
    return saddle_energy_kcal, my_path

sequence = ""    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read structures and compute saddle.')
    parser.add_argument("-f", "--structure_file", type=str, required=True, help="File with Sequence and RNA secondary structures.")
    parser.add_argument("-m", "--width", type=int, required=False, help="Width for paths to follow.")
    parser.add_argument("-P", "--energy_parameter", type=str, required=False, help="Vienna RNA Parameter File.")
    args = parser.parse_args()

    if args.energy_parameter:
        RNA.read_parameter_file(args.energy_parameter)

    width = 1
    if(args.width):
        width = args.width
    
    s1 = ""
    s2 = ""
    with open(args.structure_file, "r") as f:
        lines = f.readlines()
        sequence = str(lines[1]).strip()
        s1 = str(lines[2]).strip()
        s2 = str(lines[3]).strip()
        #sequence = sequence.strip()
    print(sequence)
    print(s1)
    print(s2)
    fold_compound = RNA.fold_compound(sequence)
    saddle_energy_kcal, my_path = find_saddle(fold_compound, s1, s2, width)
    
    #s = my_path[0][0]
    #pt = RNA.ptable_from_string(s)
    #pt = [ int(x) for x in list(pt)][1:]
    #print(pt)
    #pathg = fold_compound.path_gradient(tuple(pt)) #, RNA.PATH_STEEPEST_DESCENT)
    #print(pathg)
    #exit()
    sequence = sequence.strip()
    print(sequence)
    
    min_path = []
    try:
        for m in my_path:
            s = m[0]
            toCall = "echo -e \"" + sequence + "\n" +s+"\n\" | RNAlocmin"
            #print(toCall)
            output = subprocess.check_output(toCall, shell=True)
            lines = output.split()
            #print(lines)
            seq, ind, struc, en, ind2 = lines
            ind = int(ind)
            struc = str(struc.strip())
            en = float(en)
            min_path.append((struc, en))
    except Exception as e:
        print("ERROR: could not find minima with RNAlocmin!")
        pass
        
    max_diff_s = 0
    max_diff_m = 0
    for i, s_e in enumerate(my_path):
        print(s_e[0], s_e[1], min_path[i][0], min_path[i][1])
        if i+1 < len(my_path):
            s_e2 = my_path[i+1]	
            diff_s = abs(s_e[1] - s_e2[1])
            diff_m = abs(min_path[i][1] - min_path[i+1][1])
            if diff_s > max_diff_s:
                max_diff_s = diff_s
            if diff_m > max_diff_m:
                max_diff_m = diff_m
    print(saddle_energy_kcal)
    print(max_diff_s, max_diff_m)
            
            
            
            
        
        
   
