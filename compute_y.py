#!/usr/bin/python 
# coding: utf-8
import math


def do_one(k):
    v = pow(k - 1, -0.5)
    r = [v for i in xrange(k - 1)] # l in k - 1 dimensional
    return r

def bigger_than_one(_i, k):
    v = (1 + pow(k, 0.5)) / (pow(k - 1, (3 / 2))) * -1
    r = [v for i in xrange(k - 1)] # l in k - 1 dimensional
    t = [0 for i in xrange(k - 1)] # e in k - 1 dimensional
    t[_i - 2] = pow((k / (k - 1)), 0.5)
    for i in xrange(k - 1):
        r[i] = r[i] + t[i]
    return r

def compute_y(k):
    v = []
    v.append(do_one(k))
    for i in xrange(2, k + 1): # 2, 3 ... k
        v.append(bigger_than_one(i, k))
    return v 

if __name__ == '__main__':
    for i in xrange(2, 5):
        print compute_y(i)



