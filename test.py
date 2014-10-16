# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 15:38:42 2014

@author: Parke

RUN IN INTERACTIVE SCIPY SESSON
"""

def makecircs(N,M):
    circs = []
    for i in range(M):
        r = abs(normal(1,1.0,N))
        x,y = random([2,N])*10
        circs.append(array([x,y,r]).T)
    return circs
    
def plotem(circs):
    allcircs = np.vstack(circs)
    divider = len(circs[0])
    ta = np.linspace(0,2*pi,200)
    for i,circ in enumerate(allcircs):
        xc,yc,r = circ
        x,y = r*cos(ta) + xc, r*sin(ta) + yc
        color = 'b' if i < divider else 'r'
        plot(x,y,color)
        text(xc,yc,i,ha='center',va='center')

#%% TEST UNION
for j in range(100):
    circs, = makecircs(10,1)
    exact = sc.circles_area_union(circs)
    brute = sc.circles_area_union(circs,brutegrid=1e5)
    if abs(exact - brute) > 0.4:
        for i in arange(1,10):
            newbrute = sc.circles_area_union(circs,brutegrid=1e5)
            brute = (brute*i + newbrute)/(i+1)
            if abs(exact - brute) < 0.4/sqrt(i+1):
                break
        if i == 9:
            print 'broke it'
            break
            
#%% TEST SUBTRACTION

for j in range(100):
    circs0,circs1 = makecircs(10,2)
    exact = sc.circleset_area_subtract(circs0,circs1)
    brute = sc.circleset_area_subtract(circs0,circs1,brutegrid=1e5)
    if abs(exact - brute) > 0.4:
        for i in arange(1,10):
            newbrute = sc.circleset_area_subtract(circs0,circs1,brutegrid=1e5)
            brute = (brute*i + newbrute)/(i+1)
            if abs(exact - brute) < 0.4/sqrt(i+1):
                break
        if i == 9:
            print 'broke it'
            break
        
#%% TEST DIFFERENCE
for j in range(100):
    circs0,circs1 = makecircs(10,2)
    exact = sc.circleset_area_difference(circs0,circs1)
    brute = sc.circleset_area_difference(circs0,circs1,brutegrid=1e5)
    if abs(exact - brute) > 0.4:
        for i in arange(1,10):
            newbrute = sc.circleset_area_difference(circs0,circs1,brutegrid=1e5)
            brute = (brute*i + newbrute)/(i+1)
            if abs(exact - brute) < 0.4/sqrt(i+1):
                break
        if i == 9:
            print 'broke it'
            break