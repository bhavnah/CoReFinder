#!/usr/bin/env python

################################################
# Script: plotMeanCoverage.py
# Author: Bhavna Hurgobin
# Email:  b.hurgobin@uq.edu.au
# Usage: plotMeanCoverage.py test.cov 1000 1 1000000
################################################

import sys
import numpy as np
import math
import itertools
import matplotlib
import matplotlib.pyplot as plt

file1=open(sys.argv[1],'r')
bin_size=sys.argv[2]
start=sys.argv[3]
end=sys.argv[4]

r0=[]
r1=[]
r2=[]
pos=[]

for line in file1:
    if not line.startswith('chrom'):
       l=line.rstrip().split('\t')
       pos.append(int(l[1]))
       r0.append(float(l[2]))
       r1.append(float(l[3]))
       r2.append(float(l[4]))

plot_title="Mean coverage for " + str(bin_size) + " bp bins for " +str(start) + " to " + str(end) + " bp"

new_pos=pos[int(start)-1:int(end)-1]
new_r0=r0[int(start)-1:int(end)-1]
new_r1=r1[int(start)-1:int(end)-1]
new_r2=r2[int(start)-1:int(end)-1]

bin=int(bin_size)
counter=new_pos[0]

counter_list=[]
while counter <= new_pos[-1]:
      counter_list.append(int(counter))
      counter += bin

#if counter_list[-1] < new_pos[-1]:
#   counter_list.append(new_pos[-1])

def partition(nplist):
    partitions = [nplist[i:i+bin] for i in xrange(int(start), int(end), bin)]
    return partitions

r0_bin=partition(new_r0)
r1_bin=partition(new_r1)
r2_bin=partition(new_r2)


def mean_list(nplist,partitions):
    mean_list=[sum(x)/float(bin) for x in partitions]
    return mean_list

mean_r0=mean_list(r0,r0_bin)
mean_r1=mean_list(r1,r1_bin)
mean_r2=mean_list(r2,r2_bin)

#mean_r0_norm=[]
#mean_r1_norm=[]
#mean_r2_norm=[]

#max_mean_r0=max(mean_r0)
#max_mean_r1=max(mean_r1)
#max_mean_r2=max(mean_r2)

'''
for i in mean_r0:
    mean_r0_norm.append(i/max_mean_r0)

for i in mean_r1:
    mean_r1_norm.append(i/max_mean_r1)

for i in mean_r2:
    mean_r2_norm.append(i/max_mean_r2)

#define plot size in inches (width, height) & resolution(DPI)
plt.figure(1,figsize=(20, 8))
'''

x=counter_list
y_r0=mean_r0
y_r1=mean_r1
y_r2=mean_r2

plt.plot(x,y_r2,'r', label="-r2")
plt.plot(x,y_r1,'g', label="-r1")
plt.plot(x,y_r0,'b', label="-r0")
#plt.xticks(np.arange(int(counter_list[0]),int(counter_list[-1]),1000000))
plt.xlabel("Distance along chromsome (bp)")
plt.ylabel("Mean Coverage")
plt.legend(loc='upper right')

plt.savefig(plot_title)
