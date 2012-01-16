"""
calculate the mean area proportional value in each ring zone after running a 
a certain amount of the loops.
"""
import numpy as np
import scipy as sp
import pylab

def plot_ratio_function(zone_position, each_time_ratio, type_fiber, 
                        probability_function):
    iteration_times = len(each_time_ratio) / 2
    each_kind_ratio = [0] * type_fiber
    each_coefficients = [0.0173,  0.01975]
    for i_type in sp.arange(type_fiber):
        each_kind_ratio[i_type] = []
        for i_iteration in sp.arange(iteration_times):
            index_same_kind = type_fiber * i_iteration + i_type
            each_kind_ratio[i_type].append(each_time_ratio[index_same_kind])
        each_kind_ratio[i_type] = sp.array(each_kind_ratio[i_type])
        print 'each_kind_ratio', each_kind_ratio[i_type]
        mean_ratio = sp.zeros(len(zone_position), float)
        value_in_iteration = sp.zeros(iteration_times, float)
        for i_position in sp.arange(len(zone_position) - 1):
            for i_iteration in sp.arange(iteration_times):
                print each_kind_ratio[i_type][i_iteration][i_position]
                value_in_iteration[i_iteration] = each_kind_ratio[i_type][i_iteration][i_position]
            mean_ratio[i_position] = np.mean(value_in_iteration)
        print 'mean_ratio', mean_ratio
        position_for_function = sp.linspace(0, zone_position[-1], 50)
        each_probability_function = probability_function[i_type]
        coefficient = each_coefficients[i_type]
        pylab.figure()
        pylab.plot(position_for_function, each_probability_function(position_for_function) / coefficient, '--')
        pylab.plot(zone_position, mean_ratio, 'o')
        pylab.plot()
        pylab.xlabel('Relative position in the yarn domain')
        pylab.ylabel('Probability value')
        pylab.xlim(0., 1.05)
        pylab.ylim(0., 0.90)
        pylab.title('Area probability for %d th kind of fiber'% (i_type+1))
        
        
        #print 'each_kind_ratio value', each_kind_ratio