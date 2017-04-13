#!/usr/bin/python3.5

print("-- Starting trace analysis!")

print("-- Opening trace file")
traceFile = open("tmp_diffsel_result.trace")

print("-- Listing categories: ",end='')
def strip(str):
    if str[0]=='#':
        return str[1:]
    elif str[-1:]=='\n':
        return str[:-1]
    else:
        return str
categories = list(map(strip,traceFile.readline().split('\t')))
for cat in categories:
    print(cat, end='')
    if cat != categories[-1]:
        print(", ",end='')
    else:
        print("")

print("-- Building data map")
mymap = {}
for cat in categories:
    mymap[cat] = []
for line in traceFile:
    lineData = list(map(float,map(strip,line.split('\t'))))
    for cat in range(len(categories)):
        mymap[categories[cat]].append(lineData[cat])

print("-- Computing mean and variance for every category")
from statistics import mean, variance
for cat in categories:
    print("  * "+cat+": \t"+str(round(mean(mymap[cat]), 2))+"\t"+str(round(variance(mymap[cat]), 3)))
