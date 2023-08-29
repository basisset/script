#!/usr/bin/env python3
import numpy as np

def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        values = []

        for line in lines:
            line = line.strip()
            try:
                value = float(line)
                values.append(value)
            except ValueError:
                pass        
    return values

def find_zero_indexes(lst):
    return [index for index, value in enumerate(lst) if value == 0]

# Read Mz data from input files
file_path_list1 = 'UiCA-3500eV.txt'
file_path_list2 = 'UiCA-3630eV.txt'
file_path_list3 = 'UiCA-3900eV.txt'

list1 = read_list_from_file(file_path_list1)
list2 = read_list_from_file(file_path_list2)
list3 = read_list_from_file(file_path_list3)

# Read Subtracted and integrated data from input files
sub_file1 = 'peakareas_UiCA_3500eV.txt'
sub_file2 = 'peakareas_UiCA_3630eV.txt'
sub_file3 = 'peakareas_UiCA_3900eV.txt'

sub_data1 = read_list_from_file(sub_file1)
sub_data2 = read_list_from_file(sub_file2)
sub_data3 = read_list_from_file(sub_file3)

# Combine lists and get unique values
combined_list = list1 + list2 + list3
unique_values = set(combined_list)

# Convert set to sorted list and print
sorted_unique_values = sorted(unique_values)

print(sorted_unique_values)

# Read values from list1.txt again and replace missing values with 0
list1_updated = list1.copy()
list2_updated = list2.copy()
list3_updated = list3.copy()

for i, value in enumerate(sorted_unique_values):
    if value not in list1_updated:
        list1_updated.insert(i, 0)
    if value not in list2_updated:
        list2_updated.insert(i, 0)
    if value not in list3_updated:
        list3_updated.insert(i, 0)

print('\n')
print(f'Mz for data1: {list1_updated}')
print(f'Mz for data2: {list2_updated}')
print(f'Mz for data3: {list3_updated}')
print('\n')

zero_indexes1 = find_zero_indexes(list1_updated)
zero_indexes2 = find_zero_indexes(list2_updated)
zero_indexes3 = find_zero_indexes(list3_updated)

for index in zero_indexes1:
    sub_data1.insert(index,0)
for index in zero_indexes2:
    sub_data2.insert(index,0)
for index in zero_indexes3:
    sub_data3.insert(index,0)

sub_data1=np.array(sub_data1)
sub_data2=np.array(sub_data2)
sub_data3=np.array(sub_data3)

np.savetxt('generated_pa1.txt',sub_data1)

print(f'Integrated A for data1: {sub_data1}')
print(f'Integrated A for data2: {sub_data2}')
print(f'Integrated A for data3: {sub_data3}')

print('Done!')




