"""
data processing of lab 15 from the course 'experiments in polymer science'
(Fudan University). data files are assumed to be named as [student id]-
[temperature in celsius degree]-1_DTA.txt. call analyze_all(student id)
to process the data. the outputs are in .eps and .txt formats.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.stats import linregress

def analyze(student_id, temperature, key_crysts):
    
    # read data from txt file
    data_file_path = str(student_id) + '-' + str(temperature) + '-1_DTA.txt'
    column_names = ['time', 'temp_rate', 'I']
    data = pd.read_csv(data_file_path, names=column_names, header=None)
    data['time'] = data['time'] * 10.0
    
    # calculate crystallinity
    I_0   = data['I'].min()
    I_inf = data['I'].max()
    data['cryst'] = 1 - ( (I_inf-data['I']) / (I_inf-I_0))
 
    # drop rows where crystallinity is lesser than 0.01 or larger than 0.99
    # (drop data when crystallinity has not begun or has finished)
    head = data[data['cryst'] > 0.01].index[0]
    tail = data[data['cryst'] > 0.99].index[0]
    data = data[head:tail]
    data['time'] = data['time'] - data['time'].min()
    
    plt.cla()
    plt.scatter(data['time'], data['cryst'])
    plt.title('crystallinity over time at ' + str(temperature) + '$^\circ$C')
    plt.xlabel('time / s')
    plt.ylabel('crystallinity')
    plt.savefig('part1_' + str(temperature) + '_degree_celsius.eps', 
                format='eps', dpi=1000)
    
    # return a crystallinity-time dictionary
    cryst_time_dict = {} 
    for cryst in key_crysts:
        row_id = data[data['cryst'] > cryst].index[0]
        cryst_time_dict[cryst] = round(data.at[row_id, 'time'], 2)
    return cryst_time_dict
    

def analyze_all(student_id, temps=[120, 110, 100, 90, 80],
                key_crysts = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):
    
    # page86 data processing: part 1
    part1_summary = open('part1_summary.txt', 'w')
    part1_summary.write('------------------part1------------------')
    part1_summary.write('\ncrystallinity')
    for cryst in key_crysts:
        part1_summary.write('\t   ' + str(cryst))
    part1_summary.write('\ntime / s\ntemperature')
    
    temp_thalf_dict = {}
    for temperature in temps:
        part1_summary.write('\n' + str(temperature) + '         ')
        cryst_time_dict = analyze(student_id, temperature, key_crysts)
        temp_thalf_dict[temperature] = cryst_time_dict[0.5]
        for cryst in cryst_time_dict.keys():
            part1_summary.write('\t   ' + str(cryst_time_dict[cryst]))
    part1_summary.close()
    
    # page87 data processing: part2
    part2_summary = open('part2_summary.txt', 'w')
    part2_summary.write('------------------part2------------------')
    part2_summary.write('\ntemperature / degree celsius')
    for temp in temp_thalf_dict.keys():
        part2_summary.write('\t   ' + str(temp))
    part2_summary.write('\nt_1/2 / s                   ')
    for temp in temp_thalf_dict.keys():
        value = temp_thalf_dict[temp]
        part2_summary.write('\t   ' + str(value))
    part2_summary.write('\n1/t_1/2 / s^-1              ')
    for temp in temp_thalf_dict.keys():
        value = round(1.0 / temp_thalf_dict[temp], 2)
        part2_summary.write('\t   ' + str(value))
    part2_summary.close()
          
    x, y = zip(* temp_thalf_dict.items())
    y = [1.0 / thalf for thalf in y]
    plt.cla()
    plt.scatter(x, y)
    plt.title('reciprocal of crystallization half-time against temperature')
    plt.xlabel('T')
    plt.ylabel('1/$t_{1/2}$ / $s^{-1}$')
    plt.savefig('part2.eps', format='eps', dpi=1000)
    
    # page87 data processing: part3
    part3_summary = open('part3_summary.txt', 'w')
    crysts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    cryst_time_dict = analyze(student_id, 120, key_crysts)
    x = [math.log10(cryst_time_dict[cryst]) for cryst in crysts]
    y = [math.log10(-math.log(1.0-cryst))   for cryst in crysts]
    slope, intercept, r_value, p_value, stderr = linregress(x, y)
    x_fitting = np.linspace(x[0], x[-1], 1000)
    y_fitting = [(slope * v + intercept) for v in x_fitting]
    
    part3_summary.write('------------------part3------------------')
    part3_summary.write('\nlog(-ln(crystallinity))')
    for v in y:
        part3_summary.write('\t   ' + str(round(v, 3)))
    part3_summary.write('\nlogt                   ')
    for v in x:
        part3_summary.write('\t   ' + str(round(v, 3)))
    part3_summary.write('\n\nlinear fitting: \nslope = ' + str(round(slope, 2)) 
                     + '\nintercept = ' + str(round(intercept, 2))
                     + '\nR^2 = ' + str(round(r_value**2, 5)))
    part3_summary.close()
    
    plt.cla()
    plt.scatter(x, y)
    plt.title('Avrami double-log plot')
    plt.xlabel('logt')
    plt.ylabel('log(-ln( ($I_\infty-I_t$) / ($I_\infty-I_0$) ))')
    plt.plot(x_fitting, y_fitting, '--r')
    plt.savefig('part3.eps', format='eps', dpi=1000)
    
    
    
        
        
        
        