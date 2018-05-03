import svgwrite
import numpy as np
# from svgwrite import cm, mm
import svgpathtools
from svgpathtools import parse_path, Line, Path, wsvg
from svgpathtools import kinks, smoothed_path
from IPython.display import SVG
from fitpath import *
import math

# constants
cm = 35.43307/1.25

# Colors:
col_pt     = svgwrite.rgb(70, 70, 70, '%')    # grey
col_gds    = svgwrite.rgb(0, 0, 100, '%')     # blue
col_sew    = svgwrite.rgb(0, 0, 0, '%') # black
col_hem    = svgwrite.rgb(100, 0, 0, '%') # red

def gshift(coord_x, coord_y, shift_x=20, shift_y=5):
    
    coord_x = coord_x + shift_x
    coord_y = coord_y + shift_y

    return coord_x, coord_y

def pt2cm(coord):

    return (coord[0]*cm, coord[1]*cm)

def pt2ph(coord):
    
    return '%s,%s' % (coord[0]*cm, coord[1]*cm)

def convert2Line(p1, p2):
    
    return Line(complex(p1[0]*cm,p1[1]*cm),complex(p2[0]*cm,p2[1]*cm))

def offset_curve(path, offset_distance, steps=200):
    """Takes in a Path object, `path`, and a distance,
    `offset_distance`, and outputs an piecewise-linear approximation
    of the 'parallel' offset curve."""
    nls = []
    if 'intersect' in locals():
        del intersect 
    
    '''    
    for n, segX in enumerate(path):

        u1 = segX.unit_tangent(0.0)
        u2 = segX.unit_tangent(1.0)
        mag = 1.0*cm
        tan1 = Line(segX.point(0.0), segX.point(0.0) + mag*u1*-1).reversed()
        tan2 = Line(segX.point(1.0), segX.point(1.0) + mag*u2)
        for seg in tan1, segX, tan2:
            for k in range(steps):
                t = k / float(steps)
                offset_vector = offset_distance * seg.normal(t)
                nl = Line(seg.point(t), seg.point(t) + offset_vector)
                nls.append(nl)
    '''
                
    for n, segX in enumerate(path):
        if n == len(path)-1:
            n = -1
        if 'intersect' in locals():
            segX = segX.cropped(intersect[1],1)
            del intersect
        segXN9 = segX.normal(0.9) * offset_distance
        segXN1 = segX.normal(1.0) * offset_distance 
        segYN0 = path[n+1].normal(0.0) * offset_distance
        segYN1 = path[n+1].normal(0.1) * offset_distance
        nlX = Line(segX.point(0.9) + segXN9, segX.point(1) + segXN1)
        nlY = Line(path[n+1].point(0) + segYN0, path[n+1].point(0.1) + segYN1)

        #print nlX.intersect(nlY)
        if nlX.intersect(nlY):            
            intersect = nlX.intersect(nlY)[0]
            seg = segX.cropped(0,intersect[0])                    
            for k in range(steps):
                t = k / float(steps)
                offset_vector = offset_distance * seg.normal(t)
                nl = Line(seg.point(t), seg.point(t) + offset_vector)
                nls.append(nl)        
        else:
            nls_X = []
            nls_Y = []
            for k in range(50):
                t = k / float(50)
                offset_vector_X = offset_distance * segX.normal(t)
                nl = Line(segX.point(t), segX.point(t) + offset_vector_X)
                nls_X.append(nl)

                offset_vector_Y = offset_distance * path[n+1].normal(t)
                nl = Line(path[n+1].point(t), path[n+1].point(t) + offset_vector_Y)
                nls_Y.append(nl)
                
            connect_the_dots_X = [Line(nls_X[k].end, nls_X[k+1].end) for k in range(len(nls_X)-1)]
            connect_the_dots_Y = [Line(nls_Y[k].end, nls_Y[k+1].end) for k in range(len(nls_Y)-1)]

            for n1, x in enumerate(connect_the_dots_X):
                for m, y in enumerate(connect_the_dots_Y):
                    if x.intersect(y):
                        intersect = [n1/50.,m/50.] #x.intersect(y)
                        
            if 'intersect' in locals():                
                seg = segX.cropped(0,intersect[0])
                for k in range(steps):
                    t = k / float(steps)
                    offset_vector = offset_distance * seg.normal(t)
                    nl = Line(seg.point(t), seg.point(t) + offset_vector)
                    nls.append(nl) 
            
            else:
                u1 = segX.unit_tangent(1.0)
                u2 = path[n+1].unit_tangent(0.0)        
                mag = 3.0*cm # to ensure it will intersect
                tan1 = Line(segX.point(1.0), segX.point(1.0) + mag*u1)
                tan2 = Line(path[n+1].point(0.0), path[n+1].point(0.0) + mag*u2*-1).reversed()
                tan1N0 = tan1.normal(0) * offset_distance 
                tan1N1 = tan1.normal(1.0) * offset_distance
                tan2N1 = tan2.normal(1) * offset_distance
                tan2N0 = tan2.normal(0) * offset_distance
                nlX = Line(tan1.point(0) + tan1N0, tan1.point(1) + tan1N1)
                nlY = Line(tan2.point(0) + tan2N0, tan2.point(1) + tan2N1)

                # calc angle
                a = np.array([np.real(tan1[1]), np.imag(tan1[1])])
                b = np.array([np.real(tan1[0]), np.imag(tan1[0])])
                c = np.array([np.real(tan2[0]), np.imag(tan2[0])]) # tan2[0]

                ba = a - b
                bc = c - b

                from math import acos
                from math import sqrt
                from math import pi

                def length(v):
                    return sqrt(v[0]**2+v[1]**2)
                def dot_product(v,w):
                   return v[0]*w[0]+v[1]*w[1]
                def determinant(v,w):
                   return v[0]*w[1]-v[1]*w[0]
                def inner_angle(v,w):
                   cosx=dot_product(v,w)/(length(v)*length(w))
                   rad=acos(cosx) # in radians
                   return rad*180/pi # returns degrees
                def angle_clockwise(A, B):
                    inner=inner_angle(A,B)
                    det = determinant(A,B)
                    if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
                        return inner
                    else: # if the det > 0 then A is immediately clockwise of B
                        return 360-inner


                #cosine_angle = np.dot(bc, ba) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                #angle = np.arccos(cosine_angle)
                # print angle_clockwise(bc, ba)
                if angle_clockwise(bc, ba) < 90.0: 
                    seg = segX
                    for k in range(steps):
                        t = k / float(steps)
                        offset_vector = offset_distance * seg.normal(t)
                        nl = Line(seg.point(t), seg.point(t) + offset_vector)
                        nls.append(nl)

                elif len(nlX.intersect(nlY)) > 0:

                    intersectT = nlX.intersect(nlY)              
                    #print intersectT
                    tan1 = tan1.cropped(0,intersectT[0][0])
                    tan2 = tan2.cropped(intersectT[0][1],1)
                    #print intersectT

                    for seg in segX, tan1, tan2:            
                        for k in range(steps):
                            t = k / float(steps)
                            offset_vector = offset_distance * seg.normal(t)
                            nl = Line(seg.point(t), seg.point(t) + offset_vector)
                            nls.append(nl)
                else:
                    seg = segX
                    for k in range(steps):
                        t = k / float(steps)
                        offset_vector = offset_distance * seg.normal(t)
                        nl = Line(seg.point(t), seg.point(t) + offset_vector)
                        nls.append(nl)

    connect_the_dots = [Line(nls[k].end, nls[k+1].end) for k in range(len(nls)-1)]
    if path.isclosed():
        connect_the_dots.append(Line(nls[-1].end, nls[0].end))
    offset_path = Path(*connect_the_dots)
    return offset_path


