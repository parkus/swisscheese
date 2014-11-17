# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 15:38:42 2014

@author: Parke

RUN IN INTERACTIVE SCIPY SESSON
"""
import swisscheese as sc
from numpy.random import normal,random_sample
from numpy import array,arange,sqrt,vstack,linspace,cos,sin,pi
from matplotlib.pyplot import plot,text

circ_scale = 1.0
field_scale = 10.0
Ngrid=1e5
precision = 1./sqrt(Ngrid)

def makecircs(N,M):
    circs = []
    for i in range(M):
        r = abs(normal(1,1.0,N))*circ_scale
        x,y = random_sample([2,N])*field_scale
        circs.append(array([x,y,r]).T)
    return circs
    
def plotem(circs):
    allcircs = vstack(circs)
    divider = len(circs[0])
    ta = linspace(0,2*pi,200)
    for i,circ in enumerate(allcircs):
        xc,yc,r = circ
        x,y = r*cos(ta) + xc, r*sin(ta) + yc
        color = 'b' if i < divider else 'r'
        plot(x,y,color)
        text(xc,yc,i,ha='center',va='center')

#%% TEST UNION
def union(Ncircs, Ntrials):
    for j in range(Ntrials):
        circs, = makecircs(Ncircs,1)
        exact = sc.circles_area_union(circs)
        brute = sc.circles_area_union(circs,brute=Ngrid)
        if abs(exact - brute) > brute/precision:
            for i in arange(1,10):
                newbrute = sc.circles_area_union(circs,brute=Ngrid)
                brute = (brute*i + newbrute)/(i+1)
                if abs(exact - brute) < brute/precision/sqrt(i+1):
                    break
            if i == 9:
                print 'broke it'
                break
    if j == Ntrials - 1:
        print 'so far so good!'
        
#%% TEST SUBTRACTION
def subtraction(Ncircs,Ntrials):
    for j in range(Ntrials):
        circs0,circs1 = makecircs(Ncircs,2)
        exact = sc.circleset_area_subtract(circs0,circs1)
        brute = sc.circleset_area_subtract(circs0,circs1,brute=Ngrid)
        if abs(exact - brute) > brute/precision:
            for i in arange(1,10):
                newbrute = sc.circleset_area_subtract(circs0,circs1,brute=Ngrid)
                brute = (brute*i + newbrute)/(i+1)
                if abs(exact - brute) < brute/precision/sqrt(i+1):
                    break
            if i == 9:
                print 'broke it'
                break
    if j == Ntrials - 1:
        print 'so far so good!'
        
#%% TEST DIFFERENCE
def difference(Ncircs,Ntrials):
    for j in range(Ntrials):
        circs0,circs1 = makecircs(Ncircs,2)
        exact = sc.circleset_area_difference(circs0,circs1)
        brute = sc.circleset_area_difference(circs0,circs1,brute=Ngrid)
        if abs(exact - brute) > brute/precision:
            for i in arange(1,10):
                newbrute = sc.circleset_area_difference(circs0,circs1,brute=Ngrid)
                brute = (brute*i + newbrute)/(i+1)
                if abs(exact - brute) < brute/precision/sqrt(i+1):
                    break
            if i == 9:
                print 'broke it'
                break
    if j == Ntrials - 1:
        print 'so far so good!'