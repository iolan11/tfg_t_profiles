# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:17:19 2021

@author: leoel
"""
import csv


def readMyFile(filename):
    y = []
    x = []

    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        for row in csvReader:
            y.append(row[0])
            x.append(row[1])

    return y, x

def readMyList(filename):
    x = []

    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        for row in csvReader:
            x.append(row[0])

    return x