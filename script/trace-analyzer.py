#!/usr/bin/python3.5

print("-- Parsing command line arguments")
from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='A small script that takes a coevol trace file and computes the mean and variance for every column.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the trace file')
parser.add_argument('-b', '--burnin', type=float, default=0, help="the burn-in percentage (a float bewteen 0 and 100)")
args = parser.parse_args()
traceFile = args.inputFile[0]
print("-- Trace file is "+traceFile.name)
burnin = args.burnin
print("-- Burn-in is "+str(burnin)+"%")

print("-- Starting trace analysis!")

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
entries = len(mymap[categories[0]])
print("-- Trace file contains "+str(entries)+" entries")
burninEnd = int(burnin*entries/100)
print("-- The end of burn-in is "+str(burninEnd))

print("-- Computing mean and variance for every category")
from statistics import mean, variance
for cat in categories:
    print("  * "+cat+": \t"+str(round(mean(mymap[cat][burninEnd:]), 2))+"\t"+str(round(variance(mymap[cat][burninEnd:]), 3)))
