#This script is fully functional. Needs annotation (6/30/22)

import csv
import pandas as pd

working_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\DDA_results"
MS2_path = r"D:\HyPep1.0\HyPep_Simple_ASMS_Results\Raw_Files\RawConverter_output\2021_0817_Brain_1.ms2"
tissue_type = 'Brain1_3'

with open(MS2_path) as input:
    lst = [line.strip() for line in input]

new_list= []
final_lst = []
final_lst.append(['m/z', 'resolution', 'charge', 'intensity', 'MS2','scan_number','precursor_charge'])
ms2_list = []

new = lst

for i in new:
    new_list.append(i.split())
    if '@' in i:
        x = i.split()
        for y in x:
            if '@' in y:
                ms2 = y[0:y.index('@')]
                ms2_list.append(str(ms2))

header_list = new_list[0:26]
new_list = new_list[26:] # starts from line 26 to remove the first few header lines so that program could proceed
seperation_list = []
scan_number_list = []    
precursor_charge_list = []

for i in header_list:
    if 'S' in i:
        scan_number_list.append(i[1])
    if 'Z' in i:
        precursor_charge_list.append(i[1])

for i in range(len(new_list)):
    #print(i, new_list[i])
    if 'RetTime' in new_list[i]:
        seperation_list.append(i-1)
    if 'PrecursorInt' in new_list[i]:
        seperation_list.append(i+2)
    if 'S' in new_list[i]:
        scan_number_list.append(new_list[i][1])
    if 'Z' in new_list[i]:
        precursor_charge_list.append(new_list[i][1])

seperation_pairs = []
start = 0
for i in range(int(len(seperation_list)/2)):
    seperation_pairs.append((seperation_list[i+start],seperation_list[i+start+1]))
    start +=1 
 
update_index = 0
for start,end in seperation_pairs:
    start += update_index
    end += update_index
    new_list[start:end] = '-'
    update_index -= (end-start-1)

ms2_list_index = 0
scan_number_index = 0
precursor_charge_index = 0

for element in new_list:
    if element == '-':
        ms2_list_index+=1
        scan_number_index+=1
        precursor_charge_index+=1
        continue   
    element.append(ms2_list[ms2_list_index])
    element.append(scan_number_list[scan_number_index])
    element.append(precursor_charge_list[precursor_charge_index])
    final_lst.append(element)

out_name = working_directory + '\\'+tissue_type+'_formatted.txt'
with open(out_name,'w') as output:
    for i in final_lst:
        for j in i:
            output.write(str(j + ','))
        output.write('\n')
# =============================================================================
# csv_name =  rawconverter_formatted_out_file_name_path_csv
# read_txt = pd.read_csv(out_name)
# read_txt.to_csv(csv_name, index=None)
# print('#9: ',datetime.now())
# read_csv_path = csv_name
# SIM_result_imp = pd.read_csv(read_csv_path, sep = ' ')
# SIM_result = pd.DataFrame()
# SIM_result['m/z'] = SIM_result_imp['m/z']
# SIM_result['resolution'] = SIM_result_imp['resolution']
# SIM_result['charge'] = SIM_result_imp['charge']
# SIM_result['intensity'] = SIM_result_imp['intensity']
# SIM_result['MS2'] = SIM_result_imp['MS2']
# SIM_result['Scan #'] = SIM_result_imp['scan_number']
# SIM_result['Precursor Charge'] = SIM_result_imp['precursor_charge']
# csv_path = rawconverter_formatted_out_file_name_path_csv
# with open(csv_path,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         SIM_result.to_csv(filec,index=False)
# =============================================================================
