# -*- coding: utf-8 -*-
"""
Created on Sun Sep 07 13:39:08 2014

@author: Parke
"""
import numpy as np
import math

def dist(x0,y0,x1,y1):
    """Compute the distance between two points."""
    return np.sqrt((x1-x0)**2 + (y1-y0)**2)
    
def polygon_area(x,y):
    """
    Compute the area of a polygon given the coordinates of the vertices.
    
    Points given in counterclockwise order result in a positive area. 
    Points in CW order result in negative area.
    """
    return 0.5*(np.sum(x[:-1]*y[1:]) - np.sum(x[1:]*y[:-1]))

def brute_area(xrange,yrange,in_or_out,Nsamples=1e4):
    """
    Compute the area of arbitrary shapes using a Monte-Carlo point-sampling
    method.
    
    Parameters
    ----------
    xrange : 2-element array-like
        The limits in the x coordinate beyond which the shapes do not extend.
    yrange : 2-element array-like
        As above, for the y coordinate.
    in_or_out: function
        A function that returns True or False depending on whether a provided
        point is inside or outside of the shapes. The function must accept
        either the coordinates of a single point (function(x,y)) or vectors of 
        the coordinates (function(xvec, yvec)). As always, a vectorized
        function is faster.
    Nsamples: integer, optional
        The number of random points to sample within the xrange,yrange region.
        
    Returns
    -------
    area : float
        The estimated area of the shapes.
    """
    dx, dy = xrange[1] - xrange[0], yrange[1] - yrange[0]
    area = dx*dy
    xvec = np.random.uniform(xrange[0],xrange[1],Nsamples)
    yvec = np.random.uniform(yrange[0],yrange[1],Nsamples)
    
    try:
        cnt = np.sum(in_or_out(xvec,yvec))
    except:
        cnt = 0
        for x,y in zip(xvec,yvec):
            cnt += in_or_out(x,y)
        
    return area*float(cnt)/Nsamples

def circle_area_segment(*args):
    """Compute the area of the segment of a circle.
    
    Parameters
    ----------
    args : one of either
        r,x0,y0,x1,y1 : floats
            Radius of the circle and coordinates of the chord endpoints.
        r,h : floats
            Radius of the circle and the center-to-chord distance.
    
    Returns
    -------
    area : float
        Area of the circle segment.
    """
    N = len(args)
    if N not in [2,5]:
        raise ValueError('Not a valid number of arguments. See docstring.')
    if N == 5:
        r = args[0]
        chordlen = dist(*args[1:]) #length of the chord
        sectionangle = 2.0*math.asin(chordlen/2.0/r)
        h = math.sqrt(r**2 - (chordlen/2.0)**2) #center-chord distance
    elif N == 2:
        r,h = args
        sectionangle = 2*math.acos(h/r)
        chordlen = 2*math.sqrt(r**2 - h**2)
        
    sectionarea = sectionangle/2.0*r**2 #area of the "pie slice"
    triarea = 0.5*h*chordlen #area of the triangular part of the pie slice
    segmentarea = sectionarea - triarea
        
    return segmentarea

def circle_circle_pts(circles):
    """
    Compute the points where two circles interesect.
    
    Parameters
    ----------
    circles : 2x3 array-like
        [[x0,y0,r0],[x1,y1,r1]] - the radius and center coordinates of the two
        circles
        
    Returns
    -------
    xpts : list or array
        If the circles do not interset, an empty list is returned. If one 
        circle is enclosed in another, the index ([0] or [1]) of that circle is
        returned. Otherwise, the intersection points are returned as the array 
        [[xi0,yi0], [xi1,yi1]] with xi0 >= xi1. In the rare case that the
        just touch, the one intersection point will be returned twice.
    """
    [x0,y0,r0],[x1,y1,r1] = circles       
    d = dist(x0,y0,x1,y1)
    if d >= (r0+r1):
        return []
    elif d < abs(r1-r0):
        return [1] if r1 < r0 else [0]
    else:
        #compute intersection. some variables are arbitrarly assigned just to
        #make the math more concise. the math may be worked out by
        #simultaneously solving the equations for two circles
        q = (r0**2 - r1**2 + x1**2 - x0**2 + y1**2 - y0**2)/2.0
        dx, dy = (x1 - x0), (y1 - y0)
        a = 1 + dx**2/dy**2
        b = -2*x0 - 2*q*dx/dy**2 + 2*y0*dx/dy
        c = x0**2 + y0**2 + q**2/dy**2 - 2*q*y0/dy - r0**2
        xi0 = (-b + math.sqrt(b**2 - 4*a*c))/2/a
        xi1 = (-b - math.sqrt(b**2 - 4*a*c))/2/a
        yi0, yi1 = (q - xi0*dx)/dy, (q - xi1*dx)/dy
        return np.array([[xi0,yi0], [xi1,yi1]])

def circle_circle_area(circles,method='union',xpts=None):
    """
    Find the area of the union, intersection, or (symmetric) difference of two 
    circles.
    
    Parameters
    ----------
    circles : 2x3 array-like
        [[x0,y0,r0],[x1,y1,r1]] - the radius and center coordinates of the two
        circles
    method : {'union'|'u'|'intersection'|'x'|'difference'|'d'}, optional
        'union' or 'u' : computes all of the area covered by the two circles
        'intersection' or 'x' : computes only the area covered by both circles
        'difference' or 'd': computes the area covered by only one (not both)
            circles
    xpts : 2x2 array-like, optional
        The intersection points of the two cirlces as [[x0,y0],[x1,y1]]. This 
        is provided so that the user can speed up the computation if the
        intersection points are already known.
    
    Returns
    -------
    area : float
        The computed union, intersection, or difference area. Returns 0.0 if
        the circles do not intersect.
    """
    circs = np.asarray(circles)
    method = method.lower()
    
    #find intersection pts if not provided
    x = circle_circle_pts(circs) if xpts is None else np.asanyarray(xpts)
    
    #return 0 if the circles do not intersect
    if len(x) == 0: return 0.0
    
    #otherwise, we'll need the areas of the entire circles
    cAreas = __cArea(circs) 
    
    #if one circle is enclosed in another
    if len(x) == 1:
        if method in ['union', 'u']: return max(cAreas)
        if method in ['intersection','x']: return min(cAreas)    
        if method in ['difference','d']: return max(cAreas) - min(cAreas)
    
    #if the circles intersect...
    #compute areas of the segments defined by the chord of intersection
    seginput = [[c[-1]]+list(x.ravel()) for c in circs]
    segAreas = [circle_area_segment(*s) for s in seginput]
    #and the "opposite" areas
    oppAreas = cAreas - segAreas
    
    #check whether the circle centers lie to the same side of the chord. This
    #will be true if the dot product of the chord-center to circle-center vectors
    #is positive
    chordcen = np.sum(x,0)/2.0
    cenvecs = circs[:,:2] - chordcen[np.newaxis,:]
    sameside = (np.dot(*cenvecs) > 0)
    
    #compute and return the appropriate area
    if method in ['union', 'u']:
        area = max(oppAreas) + max(segAreas) if sameside else sum(oppAreas)
    if method in ['intersection', 'x']:
        area = min(oppAreas) + min(segAreas) if sameside else sum(segAreas)
    if method in ['difference', 'd']:
        if sameside:
            area = max(oppAreas) + max(segAreas) - min(oppAreas) - min(segAreas)
        else:
            area = sum(oppAreas) - sum(segAreas)
    return area
        
def circles_area_union(circles, intersections=None, brute=False):
    """
    Compute the area of the union of an arbitrary number of circles.
    
    The function works by identifying points where circles intersect that are
    not inside any other circles. These define polygons with a circle segment
    extending beyond the end of each leg of the polygon. The polygon area and
    the areas of all the circular segments are summed. Isolated pairs of
    intersecting circles, isolated circles are dealt with separately.
    
    This algorithm is similar to finding the area of intersection of N circles,
    as illustrated elegantly by Ben Frekerickson: 
    http://www.benfrederickson.com/calculating-the-intersection-of-3-or-more-circles/
    
    This algortihm is not efficient. It is bottlenecked by finding all of the
    interesction points, which is and order N**2 process.
    
    Parameters
    ----------
    circles : Nx3 array-like
        The radius and cetner coordinates of N circles, given as [[x0,y0,r0],
        [x1,y1,r1],...,[xN,yN,rN]]
    intersections : list, optional
        The output from cricle_intersection_pts. This parameter is used to
        avoid multiple calls to the computationally intense
        circle_intersection_pts if different operations are carried out with
        the same set of cirlces. If not provided, this function will call
        circle_intersection_pts on its own.
    brute : {False|int}, optional
        Use a brute-force Monte-Carlo method to compute the union area. If
        not False, a number will specify the number of random points to use
        in estimating the area.
        
    Returns
    -------
    area : float
        The area of the union of the circles.
    """
    circs = np.asarray(circles)
    
    #he brute force method
    if brute:
        xc,yc,r = circs.T
        
        #function to decide whether points are in or out of the circles
        def in_circs(x,y):
            d = dist(x,y,xc,yc)
            return any(d < r)
        
        #edges of a rectangle containing all the cricles
        xrange = [np.min(xc-r), np.max(xc+r)]
        yrange = [np.min(yc-r), np.max(yc+r)]
        
        return brute_area(xrange, yrange, in_circs, brute)
    
    #function to decide whether an intersection point is on the border of a
    #group of intersecting circles
    def onborder(pt):
        xi,yi,c0,c1 = pt
        others = circs[[i for i in range(len(circs)) if i not in [c0,c1]]]
        d = dist(xi,yi,others[:,0],others[:,1])
        return not any(d < others[:,2])
    
    #get intersection points and keep only those on borders
    if intersections is  None: intersections = circle_intersection_pts(circs)
    xpts,xflags,incircs = intersections
    xpts = np.array(filter(onborder, xpts))
    
    #compute group and pair areas
    grouparea = __circle_group_areas(circs,xpts,'union')
    
    #compute area of isolated circles
    inflags = np.array([len(inc) > 0 for inc in incircs], bool)
    loners = np.logical_not(np.logical_or(inflags,xflags))
    lonerareas = list(__cArea(circs[loners]))
    
    return grouparea + sum(lonerareas)
        
def circleset_area_intersect(circles0, circles1, intersections=None,
                             brute=False):
    """
    Compute the area of interestion of two sets of circles, each set made up
    of an arbitrary number of circles.
    
    Note that this is different than the intersection of all of the circles
    the two sets combined. Instead, this is the intersection of two shapes,
    each of which are defined as the union of some number of circles.
    
    The algorithm operates by identifying intersection points between circles
    that lie on the boundary of the area overlapping area(s) of the two sets.
    The area(s) of the polygon(s) defined by these points are computed, along
    with the circular segments protuding beyond each edge of the polygon(s).
    Pairs of circles and circles within circles are dealt with separately.
    
    This function was created to compute the area of a complex aperture in
    astronomy that consists of one or more circular signal areas with one or
    more circular areas (such as covering a contaminating source) masked out.
    
    Parameters
    ----------
    circles0, circles1 : Nx3 and Mx3, array like
        The circles that make up each set, specified by their radii and
        the coordinates of their centers, i.e. [[x0,y0,r0],...,[xN,yN,rN]].
    intersections : list, optional
        The output from cricle_intersection_pts when provided with all of the
        circles from both sets. This parameter is used to
        avoid multiple calls to the computationally intense
        circle_intersection_pts if different operations are carried out with
        the same set(s) of cirlces. If not provided, this function will call
        circle_intersection_pts on its own.
    brute : {False|int}, optional
        Use a brute-force Monte-Carlo method to compute the union area. If
        not False, a number will specify the number of random points to use
        in estimating the area.
        
    Returns
    -------
    area : float
        The area of the intersection of the two shapes defined by the two
        sets of circles.
    
    """
    circles = np.vstack([circles0,circles1])
    N = len(circles)
    divider = len(circles0) #slice index separating the two sets of circles
    
    if brute:
        xc,yc,r = circles.T
        def in_intersection(x,y):
            d = dist(x,y,xc,yc)
            incircs = d < r
            return any(incircs[:divider]) and any(incircs[divider:])
        
        xrange = [np.min(xc-r), np.max(xc+r)]
        yrange = [np.min(yc-r), np.max(yc+r)]
        
        return brute_area(xrange, yrange, in_intersection, brute)
    
    #get intersection points
    if intersections is None: intersections = circle_intersection_pts(circles)
    xpts,xflags,incircs = intersections
    
    #find the points on the edge of the intersecting regions
    #NOTE: I tried preallocating an array. It is not any faster than just
    #appending to a list (at least not without using Cython)
    xselfflags = np.zeros(N, bool) #track whether a circle ever intersects its own set   
    edgepts = []
    for pt in xpts:
        xi,yi,c0,c1 = pt
        
        #figure out which circles the point is in
        others = circles[[i for i in range(len(circles)) if i not in [c0,c1]]]
        d = dist(xi,yi,others[:,0],others[:,1])
        newdiv = divider - int(c0 < divider) - int(c1 < divider)
        inside = d < others[:,2]
        inset0 = any(inside[:newdiv])
        inset1 = any(inside[newdiv:])
        
        #for intersections between circles of different sets, keep only points
        #that are in one set or no sets
        if (c0 < divider and c1 >= divider) or (c0 >= divider and c1 < divider):
            if inset0 + inset1 < 2:
                edgepts.append(pt)
        #for intersections between circles of the same set, keep only points
        #that are ONLU in the opposing set
        else:
            xselfflags[[c0,c1]] = True
            if c0 < divider and not inset0 and inset1:
                edgepts.append(pt)
            if c0 >= divider and not inset1 and inset0:
                edgepts.append(pt)
    edgepts = np.array(edgepts)
    
    #identify lone circles fully within just the opposing set
    inotherset = np.zeros(N, bool)
    for i,inc in enumerate(incircs):
        if len(inc) == 0: continue
        inc = np.asarray(inc)
        if i < divider:
            if any(inc >= divider) and (not any(inc < divider)): 
                inotherset[i] = True
        else:
            if any(inc < divider) and (not any(inc >= divider)):
                inotherset[i] = True
    #and that do not intersect other circles in their set
    loners = np.logical_and(np.logical_not(xselfflags), inotherset)
    lonerAreas = __cArea(circles[loners]) #compute their area
    
    #discard any xpts on loner circles. these circles can be fully within
    #a circle of another set and intersect other circles in that set, so if we
    #don't remove these points they'll get counted twice
    if len(edgepts):
        keep = np.ones(len(edgepts), bool)
        for c in np.nonzero(loners)[0]:
            keep = np.logical_and(edgepts[:,2] != c, edgepts[:,3] != c)
            edgepts = edgepts[keep]
    
    #compute group and pair areas
    groupareas = __circle_group_areas(circles,edgepts,'intersection',divider)

    return np.sum(lonerAreas) + np.sum(groupareas)
    
def circleset_area_difference(circles0,circles1,intersections=None,
                              brute=False):
    """
    Compute the area of the (symmetric) difference of two shapes defined by
    two sets of circles, each of an arbitrary number of circles.
    
    This function simply uses the circles_area_union to find the area of both
    sets of circles combined, then subtracts the intersection area found with
    circleset_area_intersect.
    
    Parameters
    ----------
    circles0, circles1 : Nx3 and Mx3, array like
        The circles that make up each set, specified by their radii and
        the coordinates of their centers, i.e. [[x0,y0,r0],...,[xN,yN,rN]].
    intersections : list, optional
        The output from cricle_intersection_pts when provided with all of the
        circles from both sets. This parameter is used to
        avoid multiple calls to the computationally intense
        circle_intersection_pts if different operations are carried out with
        the same set(s) of cirlces. If not provided, this function will call
        circle_intersection_pts on its own.
    brute : {False|int}, optional
        Use a brute-force Monte-Carlo method to compute the union area. If
        not False, a number will specify the number of random points to use
        in estimating the area.
        
    Returns
    -------
    area : float
        The area of the difference of the two shapes defined by the two
        sets of circles.
    """
    allcircles = np.vstack([circles0,circles1])
    divider = len(circles0) #slice index separating the two sets of circles
    
    if brute:
        xc,yc,r = allcircles.T
        
        #function to decide if a point is in one set but not the other
        def in_diff(x,y):
            d = dist(x,y,xc,yc)
            incircs = d < r
            in0 = any(incircs[:divider])
            in1 = any(incircs[divider:])
            return (in0 and not in1) or (in1 and not in0)
        
        #x and y limits of the smallest rectangle contianing all circles
        xrange = [np.min(xc-r), np.max(xc+r)]
        yrange = [np.min(yc-r), np.max(yc+r)]
        
        return brute_area(xrange, yrange, in_diff, brute)
    
    if intersections is None: xpts = circle_intersection_pts(allcircles)
    
    Aunion = circles_area_union(allcircles, xpts)
    Ax = circleset_area_intersect(circles0,circles1,xpts)
    return Aunion - Ax
    
def circleset_area_subtract(circles, subcircles, intersections=None, 
                            brute=False):
    """
    Compute the area of the union of a set of circles and subtract any area that
    overlaps with another set of circles.
    
    Parameters
    ----------
    circles : Nx3 array-like
        The circles that make up the base set, specified by their radii and
        the coordinates of their centers, i.e. [[x0,y0,r0],...,[xN,yN,rN]].
    subcircles : Nx3 array-like
        The circles to be subtracted from the base set, specified in the same
        way.
    intersections : list, optional
        The output from cricle_intersection_pts when provided with all of the
        circles from both sets. This parameter is used to
        avoid multiple calls to the computationally intense
        circle_intersection_pts if different operations are carried out with
        the same set(s) of cirlces. If not provided, this function will call
        circle_intersection_pts on its own.
    brute : {False|int}, optional
        Use a brute-force Monte-Carlo method to compute the union area. If
        not False, specifies the number of random points to use
        in estimating the area.
        
    Returns
    -------
    area : float
        The area of the union of the first set of circles, excluding any area
        covered by the subcircles set.
    """
    allcircles = np.vstack([circles,subcircles])
    divider = len(circles) #slice index for the division between the two sets
    
    if brute:
        xc,yc,r = allcircles.T
        
        #function to decide if points are in circles but not subcircles
        def good(x,y):
            d = dist(x,y,xc,yc)
            incircs = d < r
            in0 = any(incircs[:divider])
            in1 = any(incircs[divider:])
            return (in0 and not in1)
        
        #bounds of the rectangle including all of the circles
        xrange = [np.min(xc-r), np.max(xc+r)]
        yrange = [np.min(yc-r), np.max(yc+r)]
        
        return brute_area(xrange, yrange, good, brute)
    
    #get intersection points, if not provided
    if intersections is None: ints = circle_intersection_pts(allcircles)
    xpts,xflags,incircs = ints
    
    #compute the area of the intersection of the two sets
    Aint = circleset_area_intersect(circles,subcircles,ints)
    
    #compute the area of the union of the first set
    #for this, we only need the intersection information for the first set
    #I'll filter out the relevant stuff, but I'm not actually sure this is
    #faster than just running circle_intersection_pts again for circles
    onlyin0 = np.logical_and(xpts[:,2] < divider, xpts[:,3] < divider)
    xpts0 = xpts[onlyin0]
    xflags0 = np.array([i in xpts0[:,2:] for i in range(divider)], bool)
    incircs0 = [filter(lambda i: i < divider, ic) for ic in incircs[:divider]]
    ints0 = [xpts0,xflags0,incircs0]
    Aunion = circles_area_union(circles,ints0)
    
    return Aunion - Aint

def __circpolyarea(circles,polygon,method='union',divider=None):
    """
    Compute the area of a polygon defined by the intersection points between
    circles and the circle segments protruding beyond the polygon edges.
    """
    n = len(polygon)    
    if n <= 3:
        raise ValueError('Polygon with 3 ot fewer vertices encountered.')
    
    #area of the polygon itself
    x,y = polygon.T[:2]
    polyarea = polygon_area(x,y)
    
    #check if a segment line is in a circle of the apposing set
    def checkline(pt0,pt1,c):
        slc = slice(divider,None) if c < divider else slice(divider)
        return __line_in_a_circle(pt0,pt1,circles[slc])
    
    #function to determine how far a circle protrudes from a polygon edge
    def protrude(chordlen,radius,outside):
        dist_cord_cen = np.sqrt(radius**2 - (chordlen/2.0)**2)
        dist_prot = radius + dist_cord_cen if outside else radius - dist_cord_cen
        return dist_prot
    
    #compute the areas of the protruding segments
    borderareas = np.zeros(len(polygon)-1)
    for j in range(n-1): #loop through legs of the polygon
        x0,y0,c00,c01 = polygon[j]
        x1,y1,c10,c11 = polygon[j+1]
        chord = np.array([x1-x0, y1-y0])
        
        #figure out whether each circle is inside or outside of the polygon
        #by checking where a vector to their center (nadir) points relative to
        #the vector defined by the polygon leg (chord)
        cens = circles[[c00,c01],:2]
        xpt0 = np.array([x0,y0])
        nadirs = cens - xpt0[np.newaxis,:]
        outside = __CCWmeasure(chord,nadirs) > 0.0
        
        #if the points are on the same two circles
        if (c00 == c10 and c01 == c11) or (c00 == c11 and c01 == c10):
            #compute how much each circle protrudes from the polygon
            chordlen = dist(x0,y0,x1,y1)
            radii = circles[[c00,c01],2]
            protrusions = np.array(map(protrude, [chordlen]*2, radii, outside))
            
            #now select which circle's protruding area to use
            if method == 'union':
                #use whichever segment sticks out more
                winner = np.argmax(protrusions)
            elif method == 'intersection':
                #use whichever segment sticks out more but is still entirely
                #within another circle
            
                #construct lines from the center of the polygon side (chord)
                #to the far point of the protruding circles
                if x1 == x0:
                    dx, dy = protrusions, np.array([0.0, 0.0])
                else:
                    chordslope = (y1-y0)/(x1-x0)
                    dy = np.sqrt(protrusions**2/(1 + chordslope**2))
                    dx = abs(chordslope*dy)
                if x1 > x0: dy = -dy
                if y1 < y0: dx = -dx
                chordcen = np.array([x0+x1, y0+y1])/2.0
                farpoints = chordcen[np.newaxis,:] + np.vstack([dx,dy]).T
                
                #now see if these lines are in other circles
                incircs = [checkline(chordcen,fp,c) for fp,c in 
                           zip(farpoints,[c00,c01])]
                if all(incircs):
                    #if segments are both within another circle, pick the largest
                    winner = np.argmax(protrusions)
                else:
                    #pick whichever is in a circle
                    winner = 0 if incircs[0] else 1
        else:
            #only one circle "defines" that leg of the polygon. pick it.
            #(both enpoints of the polygon leg are on the circle)
            winner = 0 if c00 in [c10,c11] else 1
        
        #finally actually compute the area of the protrusion
        c = [c00,c01][winner]
        xc,yc,r = circles[c]
        segarea = circle_area_segment(r,x0,y0,x1,y1)
        protarea = __cArea(circles[[c]]) - segarea if outside[winner] else segarea
        borderareas[j] = protarea
    
    #now just sum the polygon area and the protruding areas
    return polyarea + np.sum(borderareas)
    
def __circle_group_areas(circles,xpts,method='union',divider=None):
    """
    Take the circles and their intersetion points, arrange the intersection
    points into polygons (or lines if there are just two intersecting circles)
    for each group of connected circles, and compute the area of it all.
    
    only method = 'union' and 'intersection' are supported. If 'intersection' is
    specified, the divider keyword must also be supplied to let the function
    know the slice index for circles that divides the circles of the first and
    second sets.
    """
    def getconnector(c0,c1,c):
        if c0 in c:
            if c1 in c: return [c0,c1]
            return [c0]
        if c1 in c: return [c1]
        return []
        
    #find nearest counterclockwise intersection point starting from pt on one
    #of the two circles that intersect at pt
    def nextpt(pt,pts):
        c0,c1 = pt[2:] #the two intersecting circles
        cs = pts[:,2:] #the circles involved in all other intersection pts
        
        #figure out which of the other intersection pts are connected
        connectors = [getconnector(c0,c1,c) for c in cs]
        connected = np.array([len(c) for c in connectors], bool)
        
        if sum(connected) == 0:
            raise ValueError('Well, that shouldn\'t have happened. Couldn\'t find a next point.')
        
        #if only oine is connected, the job is done
        elif sum(connected) == 1:
            return np.nonzero(connected)[0][0]
            
        #otherwise, we need to figure out which is the nerest in a
        #counterclockwise sense
        if sum(connected) > 1:
            chords = pts[connected,:2] - pt[np.newaxis,:2] #the connecting lines
            args = np.nonzero(connected)[0] #indices of the connected pts
            
            #if we are computing intersection areas, we must get rid of any
            #connecting lines that fall outside of one or more of the circle sets
            if method == 'intersection':
                keep = np.ones(len(chords), bool)
                for i,arg in enumerate(args):
                    in0 = __line_in_a_circle(pts[arg,:2], pt[:2], circles[:divider])
                    in1 = __line_in_a_circle(pts[arg,:2], pt[:2], circles[divider:])
                    if not (in0 and in1): keep[i] = False
                chords, args = chords[keep], args[keep]
                
            #make a vector that points otward from the intersection of the two
            #circle. this is the reference vector.
            cen1 = circles[c1,:2]
            rad1 = pt[:2] - cen1
            cen0 = circles[c0,:2]
            rad0 = pt[:2] - cen0
            outvec = rad0/norm(rad0) + rad1/norm(rad1)
            
            #now see which intersection point is closest to the reference vector
            #in a counterclockwise sense. I think this is called a line sweep
            ccw = __CCWmeasure(outvec,chords)
            arg = np.argmin(ccw)
            
            #return just the index of that point
            return args[arg]
        
    #-----traverse the intersection points to identify polygons-----
    polygons = []
    while len(xpts) > 0:
        #starting pt
        polygon = [xpts[0]]
        
        #find the first next point
        i = nextpt(xpts[0], xpts[1:]) + 1
        
        #now keep seraching for the next points until we get back to the start
        while True:
            pt = xpts[i]
            polygon.append(pt)
            xpts = np.delete(xpts, i, 0) #we can delete points as we go (just not the start)
            if all(polygon[-1] == polygon[0]): #back at the start point
                polygons.append(np.array(polygon))
                break
            i = nextpt(pt, xpts)
    
    #-----sum area of groups of circles (polygons and associated segments)-----
    areas = np.zeros(len(polygons))
    for i,polygon in enumerate(polygons):
        #if there are three points in the polygon, it's a line
        #the last point is always a copy of the first
        if len(polygon) == 3:
            c0,c1 = polygon[0,2:]
            circs = circles[[c0,c1]]
            xpts = polygon[:2,:2]
            
            if method == 'intersection':
                #if the circles are of the same set, compute the union
                #(this can happen if they are both fully inside of another
                #group of the other set)
                if sum([c0 < divider, c1 < divider]) != 1:
                    area = circle_circle_area(circs,method='u',xpts=xpts)
                else:
                    area = circle_circle_area(circs,method='x',xpts=xpts)
            else:
                area = circle_circle_area(circs,method='u',xpts=xpts)
        #if there are 4 or more points in the polygon, compute it's area!
        else:
            area = __circpolyarea(circles, polygon, method=method, divider=divider)
        areas[i] = area
    
    return np.sum(areas)
    
def circle_intersection_pts(circles):
    """Idetifies all intersection points between the provided circles.
    
    Parameters
    ----------
    circles : Nx3 array-like
        The circles in question, specified by their radii and
        the coordinates of their centers, i.e. [[x0,y0,r0],...,[xN,yN,rN]].
    
    Returns
    -------
    xpts : Mx4 array
        The list of intersection points as the x and y coordinates and the indices
        of the two intersecting circles, i.e. [[x0,y0,c00,c01], ...,
        [xM, yM, cM0, cM1]]
    xflags : len(N) 1-D array
        A True|False flag for whether each circle has any intersections.
    incircs : len(N) list
        A list where item i states which other circles fully envelop 
        circles[i]. E.g. if incircs[2] == [0,5], then circles 0 and 5 both
        surround circle 2.
    """
    circs = np.asarray(circles)
    N = len(circs)
    
    #NOTE: appending to a list is just as fast (maybe slightly faster) as
    #preallocating a maximum-size array
    xpts = []
    xflags = np.zeros(N, bool)
    incircs = [[] for i in range(N)]
    for i in range(N-1):
        for j in range(i+1,N):
            pts = circle_circle_pts(circs[[i,j],:3])
            if len(pts) == 1: #if one circle is inside another
                if pts == [0]: incircs[i].append(j)
                else:          incircs[j].append(i)
            if len(pts) > 1: #if the circles intersect
                xflags[[i,j]] = True
                #append the circle indices to the x,y of the intersection pts
                pts = [list(pt)+[i,j] for pt in pts]
                xpts.extend(pts)
    xpts = np.array(xpts)
    
    return xpts, xflags, incircs

def __cArea(circles):
    """Compute the areas of the input circles, vectorized"""
    return math.pi*circles[:,2]**2

def line_circle_pts(xc,yc,r,x0,y0,x1,y1):
    """
    Find the intersection points of a line and a circle.
    
    Parameters
    ----------
    xc, yc : floats
        x,y coordinates of the circle center
    r : float
        radius of the circle
    x0,y0,x1,y1 : floats
        coordiantes of two points on the line
    
    Returns
    -------
    xpts : tuple
        The intersection points of the line with the circle. The tuple is empty
        if they do not intersect. If they intersect, it contains the (x,y)
        coordinates of the points, i.e. ((x0,y0), (x1,y1)) for two intersection
        points or just ((x,y)) if the line is tangent to the circle.
    """
    
    #if the line is vertical, there will be a zero division unless I deal with
    #the case specially
    if x1 == x0:
        x = x0
        dx = abs(x - xc)
        if dx > r:
            return None
        elif dx == r:
            return x,yc
        else:
            A, B, C = 1.0, -2.0*yc, yc**2 + (x - xc)**2 - r**2
            radical = np.sqrt(B**2 - 4*A*C)
            yi0, yi1 = (-B + radical)/2.0/A, (-B - radical)/2.0/A
            return (x,yi0), (x,yi1)
            
    #you can work out the math...
    m = (y1 - y0)/(x1 - x0)
    b = y0 - x0*m
    y = lambda x: x*m + b
    A = m**2 + 1
    B = 2.0*(b - yc)*m - 2*xc
    C = xc**2 + (b - yc)**2 - r**2
    radicand = B**2 - 4*A*C
    if radicand < 0:
        return ()
    elif radicand == 0:
        xi = -B/2.0/A
        yi = y(xi)
        return (xi,yi),
    else:
        radical = np.sqrt(radicand)
        xi0, xi1 = (-B + radical)/2.0/A, (-B - radical)/2.0/A
        yi0, yi1 = map(y, [xi0, xi1])
        return (xi0, yi0), (xi1, yi1)

def __line_in_a_circle(pt0,pt1,circles):
    """Figure out if the line segmenet defined by the x,y endpoint coordinates
    in pt0 and pt1 falls totally within any of the circles."""
    #if both endpoints are in the a circle, then the line segment is in the
    #circle
    (x0,y0),(x1,y1) = pt0,pt1
    x,y,r = circles.T
    d0, d1 = dist(x0, y0, x, y), dist(x1, y1, x, y)
    D0, D1 = d0 - r, d1 - r
    incirc = np.logical_and(D0 < r/1e6, D1 < r/1e6)
    return any(incirc)
        
def __CCWmeasure(a,bs):
    """Gives a measure of how closely aligned the vectors in b are to a in a CCW sense. It
    does not compute the actual angle or even cos(angle) to save computation
    time. Instead, it returns a number that is greater the more closely aligned
    the b vectors are with a. The result is in the domain [-2,2], where greater
    numbers mean b is more CCW from a and 0 is at 180 deg.  bs is an Nx2 array.
    """
    ccw90 = np.array([-a[1], a[0]])
    bnorms2 = np.sum(bs**2, axis=1)
    anorm2 = np.sum(a**2)
    dots = np.sum(bs*a[np.newaxis,:], axis=1)
    cos_angle = dots/np.sqrt(bnorms2*anorm2)
    sin_sign = np.sum(bs*ccw90[np.newaxis,:], axis=1)
    past180 = (sin_sign < 0)
    before180 = np.logical_not(past180)
    cos_angle[before180] = cos_angle[before180] + 1.0
    cos_angle[past180] = -cos_angle[past180] - 1.0
    return -cos_angle

def norm(vec):
    """
    Computes the euclidean norm (length) of an array (typically a vector).
    
    Parameters
    ----------
    vec : array-like
    
    Returns
    -------
    norm : float
        The Euclidean norm of the array, defined as the square root of the sum
        of the squares of all of the array elements.
    """
    vec2 = np.asarray(vec)
    return np.sqrt(np.sum(vec2**2))