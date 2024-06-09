# Exact algorithms for the moving firefighter problem on trees

This repository contains the experiments developed for the article "Exact algorithms for the moving firefighter problem on trees".

## Abstract

The Moving Firefighter Problem (MFP) is a generalization of the Firefighter Problem (FP) that models more realistic situations where firefighters have to spend time defending locations but also traveling between the locations they intend to defend. Unfortunately, the only known exact solution does not scale. In this paper, we establish that the MFP is NP-complete on trees of maximum degree three and present three alternative methods to find exact solutions for the case of arbitrary trees, a single initial fire, and a firefighter. The first method is a Dynamic Programming algorithm, while the other two methods are based on an Integer Quadratically Constrained Program (IQCP) and an Integer Linear Program (ILP), respectively. The latter two methods leverage tree properties to reduce the number of constraints to improve scalability. We present the results of a series of experiments that assess the performance of these proposed solutions.

## Algorithms

The algorithms developed through the research can be found under the folders named dp, greedy, ilp, iqcp, and miqcp.

The file containing the execution of the experiments is run_test.py.

## Results

The description of the instances used and the executions' results described in the article are included as a set of JSON files under the results folder.

