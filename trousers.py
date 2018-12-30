#######################################################################################
# ToDO:
#######################################################################################
'''
- insert text lengths & print checks on pattern
'''

#######################################################################################
# Measurements
#######################################################################################
seat        = 87
waist       = 80
rise        = 25.5
inside      = 76.5
bottom      = 20
waistband   = 3.5

fly_length  = 17
dart_back   = 1.50 # *2
include_hem = True
#######################################################################################
# Begin:
#######################################################################################
import svgwrite
import numpy as np
# from svgwrite import cm, mm
import svgpathtools
from svgpathtools import parse_path, Line, Path, wsvg, kinks, smoothed_path, svg2paths
from IPython.display import SVG

# import our scripts:
from fitpath import *
from utils import *

#######################################################################################
# move back when ready:
trans2 = 35*cm
#######################################################################################

#######################################################################################
# Scripts:
#######################################################################################
def pt(name, coord_x, coord_y, shift=None, cross_sz=0.5):

    # draw cross and label:
    if shift:
        coord_x, coord_y = gshift(coord_x, coord_y)

    dwg.add(dwg.line(((coord_x-cross_sz)*cm, coord_y*cm), ((coord_x+cross_sz)*cm, coord_y*cm), stroke=col_pt))

    dwg.add(dwg.line((coord_x*cm, (coord_y+cross_sz)*cm), (coord_x*cm, (coord_y-cross_sz)*cm), stroke=col_pt))

    dwg.add(dwg.text(name, insert=((coord_x+cross_sz/5)*cm,(coord_y-cross_sz/5)*cm), fill=col_pt))

    if trans2 > 0:

        dwg.add(dwg.line((((coord_x+trans2/cm)-cross_sz)*cm, coord_y*cm), (((coord_x+trans2/cm)+cross_sz)*cm, coord_y*cm), stroke=col_pt))

        dwg.add(dwg.line(((coord_x+trans2/cm)*cm, (coord_y+cross_sz)*cm), ((coord_x+trans2/cm)*cm, (coord_y-cross_sz)*cm), stroke=col_pt))

        dwg.add(dwg.text(name, insert=(((coord_x+trans2/cm)+cross_sz/5)*cm,(coord_y-cross_sz/5)*cm), fill=col_pt))

    return coord_x, coord_y

def ln(point_1, point_2, stroke=col_sew):

    ln_out = dwg.line(pt2cm(point_1), pt2cm(point_2), stroke=stroke)

    return ln_out

def addGrain(pathX, shift=0):

    # Add grain arrows
    # len_x = abs(pathX.bbox()[0]-pathX.bbox()[1])
    len_y = abs(pt0[1]*cm-pt2[1]*cm) # abs(pathX.bbox()[2]-pathX.bbox()[3])
    # center_x = min(pathX.bbox()[0:2]) + len_x/2
    center_x = pt0[0]*cm + shift
    center_y = pt0[1]*cm + len_y/2
    dwg.add(dwg.line((center_x, center_y-(len_y/4)), (center_x, center_y+(len_y/4)), stroke=col_pt, stroke_width=5))
    peak_y = center_y-(len_y/4)
    scale_arrowhead = 100
    dwg.add(dwg.path('M %d,%d,%d,%d L %d,%d Z' % (center_x-(len_y/scale_arrowhead), peak_y+(len_y/scale_arrowhead),
                                                  center_x, peak_y-(len_y/scale_arrowhead), center_x+(len_y/scale_arrowhead),
                                                  peak_y+(len_y/scale_arrowhead)), fill=col_pt))

#######################################################################################
# Begin SVG file
#######################################################################################
h = 250*cm
v = 150*cm
dwg = svgwrite.Drawing('trousers.svg', height=h, width=v)

# Points
pt0 = pt("0", 0, 0, shift=True)
pt1 = pt("1", pt0[0], pt0[1]+rise-waistband)
pt2 = pt("2", pt1[0], pt1[1]+inside)
pt3 = pt("3", pt2[0], pt2[1]-(abs(pt1[1]-pt2[1])/2)-5)
pt4 = pt("4", pt1[0], pt1[1]-rise/4)
pt5 = pt("5", pt1[0]-(seat/8)+1, pt1[1])
pt6 = pt("6", pt5[0], pt4[1])
pt7 = pt("7", pt5[0], pt0[1])
pt8 = pt("8", pt6[0]+(seat/4), pt6[1])
pt9 = pt("9", pt5[0]-(seat/16), pt5[1] + 0.2) # reduced to match
pt10 = pt("10", pt7[0]+1, pt7[1])
pt11 = pt("11", pt10[0]+(waist/4), pt10[1])
pt12 = pt("12", pt2[0]+(bottom/2)-1, pt2[1])
pt13 = pt("13", pt2[0]-(bottom/2)+1, pt2[1])
pt14 = pt("14", pt3[0]+(bottom/2)+0.5, pt3[1])
pt15 = pt("15", pt3[0]-(bottom/2)-0.5, pt3[1])

# back points:
pt16 = pt("16", pt5[0]+(abs(pt1[0]-pt5[0])/4), pt5[1])
pt17 = pt("17", pt16[0], pt4[1])
pt18 = pt("18", pt16[0], pt0[1])
pt19 = pt("19", pt18[0], pt18[1]+(abs(pt18[1]-pt16[1])/2))
pt20 = pt("20", pt18[0]+2, pt18[1])
pt21 = pt("21", pt20[0], pt20[1]-1)
pt22 = pt("22", pt9[0]-(abs(pt5[0]-pt9[0])/2), pt9[1] - 0.2) # increased to match - see pt9 above
pt23 = pt("23", pt22[0], pt22[1]+0.5)
pt24 = pt("24", pt21[0] + (np.sqrt((np.square((waist/4) + 3)) - 1)), pt0[1])
# pt25 = pt("25", # ADD IN 2 darts
pt26 = pt("26", pt17[0] + (seat/4) + 2, pt17[1])
pt27 = pt("27", pt12[0] + 2, pt12[1])
pt28 = pt("28", pt13[0] - 2, pt13[1])
pt29 = pt("29", pt14[0] + 2, pt14[1])
pt30 = pt("30", pt15[0] - 2, pt15[1])


# Define guide lines

ln_0_2 = dwg.add(ln(pt0, pt2, stroke=col_gds))
ln_1_22 = dwg.add(ln(pt22, (pt1[0]+15,pt1[1]), stroke=col_gds))
ln_6_8 = dwg.add(ln(pt6, pt8, stroke=col_gds))
ln_15_14 = dwg.add(ln(pt15, pt14, stroke=col_gds))

ln_9_15 = dwg.add(ln(pt9, pt15, stroke=col_pt))

ln_1_22_back = dwg.add(ln(pt22, (pt1[0]+15,pt1[1]), stroke=col_gds))
ln_1_22_back.translate(trans2)
ln_0_2_back = ln(pt0, pt2, stroke=col_gds)
ln_0_2_back.translate(trans2)
ln_6_26_back = ln(pt6, pt26, stroke=col_gds)
ln_6_26_back.translate(trans2)
ln_30_29_back = ln(pt30, pt29, stroke=col_gds)
ln_30_29_back.translate(trans2)
dwg.add(ln_0_2_back)
dwg.add(ln_6_26_back)
dwg.add(ln_30_29_back)

# Draw paths

# Path1
seg = Line(complex(pt10[0]*cm,pt10[1]*cm),complex(pt6[0]*cm,pt6[1]*cm))
t = seg.ilength(1.0*cm)
pt10A = pt("10A", np.real(seg.point(t))/cm, np.imag(seg.point(t))/cm)

p = (pt2cm(pt10A),pt2cm((pt0[0],pt0[1] + 0.75)),pt2cm(pt11))
pf_10A_11 = fitpath(p,10e-1)
sp_10A_11 = pathtosvg(pf_10A_11)
seg_10A_11 = parse_path(sp_10A_11)
t = seg_10A_11.ilength(4.0*cm)
pt10B = pt("10B", pt10A[0] - abs(pt10A[0] - np.real(seg_10A_11.point(t))/cm), (np.imag(seg_10A_11.point(t))/cm))
#t = seg_10A_11.ilength(1.0*cm)
#pt10C = pt("10C", np.real(seg.point(t))/cm, np.imag(seg.point(t))/cm)

path1 = svgwrite.path.Path('%s' % sp_10A_11, fill="none", stroke=col_sew) # (pt2ph(pt10B),

# only curve 0.25, not 0.5
hipline_curve = 0.2 # not 0.5
seg_11_8 = Line(complex(pt11[0]*cm,pt11[1]*cm),complex(pt8[0]*cm,pt8[1]*cm))
seg_length = seg_11_8.length()
t_mid = seg_11_8.ilength(seg_length/2)
geo_mid = seg_11_8.point(t_mid)
n = seg_11_8.normal(t_mid)
normal_line = Line(seg_11_8.point(t_mid), seg_11_8.point(t_mid) + hipline_curve*cm*n)
ctrl_pt = (np.real(normal_line.point(1)), np.imag(normal_line.point(1)))
ctrl_pt_str = str(np.real(normal_line.point(1))) + ", " + str(np.imag(normal_line.point(1)))

p_11_8 = (pt2cm(pt11),ctrl_pt,pt2cm(pt8))
pf_11_8 = fitpath(p_11_8,10e-1)
sp_11_8 = pathtosvg(pf_11_8).replace('M','L')[17:]

path1.push('%s' % sp_11_8) # curve
# path1.push('L %s' % (pt2ph(pt8))) # curve
path1.push('L %s' % (pt2ph(pt14)))
path1.push('L %s' % (pt2ph(pt12)))
path1.push('L %s' % (pt2ph(pt13)))
path1.push('L %s' % (pt2ph(pt15)))

seg = Line(complex(pt15[0]*cm,pt15[1]*cm),complex(pt9[0]*cm,pt9[1]*cm))
seg_length = seg.length()
t_mid = seg.ilength(seg_length/2)
geo_mid = seg.point(t_mid)
n = seg.normal(t_mid)
normal_line = Line(seg.point(t_mid), seg.point(t_mid) + 0.5*cm*n*-1) # reduced from 0.6
ctrl_pt = (np.real(normal_line.point(1)), np.imag(normal_line.point(1)))
ctrl_pt_str = str(np.real(normal_line.point(1))) + ", " + str(np.imag(normal_line.point(1)))

p_15_9 = (pt2cm(pt15),ctrl_pt,pt2cm(pt9))
pf_15_9 = fitpath(p_15_9,10e-1)
sp_15_9 = pathtosvg(pf_15_9).replace('M','L')
path1.push('%s' % sp_15_9) # curve

##########################################################################
# fly:
# Draw line through pt10A and 6, extrapolating to length of fly.
# Then connect to ctrl_pt5 and 9 to finish fly

ctrl_pt5 = pt("ctrl5", pt5[0] - np.sqrt(1.3), pt5[1] - np.sqrt(1.3)) # orig was 1.5
p = (pt2cm(pt10A),pt2cm(pt6),pt2cm(ctrl_pt5),pt2cm(pt9))
pf = fitpath(p,10e-6)
sp = pathtosvg(pf)
seg = Path(parse_path(sp))
t = seg.ilength((fly_length+1.0)*cm) # hem is necessary for maintaining straight line
ptF1 = pt("F1", np.real(seg.point(t))/cm, np.imag(seg.point(t))/cm)

seg = Line(complex(pt10A[0]*cm,pt10A[1]*cm),complex(ptF1[0]*cm,ptF1[1]*cm))
seg_length = seg.length()
t = seg.ilength(seg_length - 3.5*cm)
n = seg.normal(t)
F2 = seg.point(t) + 4.0*cm*n*-1 # no hem version
ptF2 = pt("F2", np.real(F2)/cm, np.imag(F2)/cm)
t = seg.ilength(seg_length)
n = seg.normal(t)
#F3 = seg.point(t) + -1.0*cm*n
#ptF3 = pt("F3", np.real(F3)/cm, np.imag(F3)/cm)

#seg = Line(complex(pt10[0]*cm,pt10[1]*cm),complex(pt6[0]*cm,pt6[1]*cm))
#t = seg.ilength(1.0*cm)
#pt10A = pt("10A", np.real(seg.point(t))/cm, np.imag(seg.point(t))/cm)

### Next curve to F1
p_9_5_F1 = (pt2cm(pt9),pt2cm(((pt9[0]+1.5), pt9[1])), pt2cm(ctrl_pt5), pt2cm(ptF1))
pf_9_5_F1 = fitpath(p_9_5_F1,10e-6)
sp_9_5_F1 = ' '.join(pathtosvg(pf_9_5_F1).split(' ')[2:])
path1.push('%s' % sp_9_5_F1) # curve

#dwg.add(pathFly)

# check fly line
# dwg.add(ln(pt10A, ptF1, stroke=col_sew))

### finish path 1
#p = (pt2cm(F1),pt2cm(pt6),pt2cm(pt10A))
#pf = fitpath(p,10e-6)
#sp = pathtosvg(pf)[17:]
#path1.push('%s' % sp) # curve
#path1.push('L %s' % pt2ph(pt10A)) # skip pt6 for strainght line F1 - 10A

## Make additional path for fly tracing hem

path1.push('Z') # close path
#path1.push(' Q %s %s L %s Z' % (pt2ph((ptF2[0], ptF1[1])), pt2ph(ptF2), pt2ph(pt10B)))
dwg.add(path1)

#pathFly = path1.copy()
#pathFly.attribs['stroke'] = col_gds #"none"
#pathFly.tostring()
#pathFly.attribs['d'] = pathFly.attribs['d'][:-1] + ' Q %s %s L %s Z' % (pt2ph((ptF2[0], ptF1[1])), pt2ph(ptF2), pt2ph(pt10B))

######################################################################
# Path2
######################################################################
# begin by outlining, then tranlate to x + 5cm
#path2 = svgwrite.path.Path('M %s,%s' % (pt2ph(pt21), pt2ph(pt24)), fill="none", stroke=col_sew)
#path2.push('L %s' % (pt2ph(pt26)))

seg = Line(complex(pt26[0]*cm,pt26[1]*cm),complex(pt29[0]*cm,pt29[1]*cm))
seg_length = seg.length()
t_mid = seg.ilength(seg_length/2)
geo_mid = seg.point(t_mid)
n = seg.normal(t_mid)
normal_line = Line(seg.point(t_mid), seg.point(t_mid) + 0.3*cm*n*-1)
ctrl_pt = (np.real(normal_line.point(1)), np.imag(normal_line.point(1)))
ctrl_pt_str = str(np.real(normal_line.point(1))) + ", " + str(np.imag(normal_line.point(1)))

p_26_29 = (pt2cm(pt26),ctrl_pt,pt2cm(pt29))
pf_26_29 = fitpath(p_26_29,10e-1)
sp_26_29 = pathtosvg(pf_26_29)#.replace('M','L')

path2 = svgwrite.path.Path('%s' % sp_26_29, fill="none", stroke=col_sew)
#path2.push('%s' % sp) # curve
#path2.push('L %s' % (pt2ph(pt29))) # curve

path2.push('L %s' % (pt2ph(pt27)))

seg = Line(complex(pt27[0]*cm,pt27[1]*cm),complex(pt28[0]*cm,pt28[1]*cm))
seg_length = seg.length()
t_mid = seg.ilength(seg_length/2)
geo_mid = seg.point(t_mid)
n = seg.normal(t_mid)
normal_line = Line(seg.point(t_mid), seg.point(t_mid) + 1.0*cm*n)
ctrl_pt = (np.real(normal_line.point(1)), np.imag(normal_line.point(1)))
ctrl_pt_str = str(np.real(normal_line.point(1))) + ", " + str(np.imag(normal_line.point(1)))

p = (pt2cm(pt27),ctrl_pt,pt2cm(pt28))
pf = fitpath(p,10e-6)
sp = pathtosvg(pf).replace('M','L')

path2.push('%s' % sp[17:]) # curve
#path2.push('L %s' % (pt2ph(pt28))) # curve down 1 cm

path2.push('L %s' % (pt2ph(pt30)))

seg = Line(complex(pt30[0]*cm,pt30[1]*cm),complex(pt23[0]*cm,pt23[1]*cm))
seg_length = seg.length()
t_mid = seg.ilength(seg_length/2)
geo_mid = seg.point(t_mid)
n = seg.normal(t_mid)
normal_line = Line(seg.point(t_mid), seg.point(t_mid) + 1.2*cm*n*-1)
ctrl_pt = (np.real(normal_line.point(1)), np.imag(normal_line.point(1)))
ctrl_pt_str = str(np.real(normal_line.point(1))) + ", " + str(np.imag(normal_line.point(1)))

p_30_23 = (pt2cm(pt30),ctrl_pt,pt2cm(pt23))
pf_30_23 = fitpath(p_30_23,10e-4)
sp_30_23 = pathtosvg(pf_30_23)
path2.push('%s' % ('C' + sp_30_23.split('C')[1])) # curve
#path2.push('L %s' % (pt2ph(pt23))) # curve


ctrl_pt16 = pt("ctrl16", pt16[0] - np.sqrt(2.25), pt16[1] - np.sqrt(2.25))
p = (pt2cm(pt23),pt2cm(ctrl_pt16),pt2cm(pt19)) #,pt2cm(pt21))
# p = (pt2cm(pt23),pt2cm(ctrl_pt16),pt2cm(pt21)) # skip point pt19
pf = fitpath(p,10e-1)
sp = pathtosvg(pf)[17:]
# path2.push('%s' % sp) # curve

#path2.push('L %s' % (pt2ph(pt19))) # curve
# path2.push('L %s' % (pt2ph(pt21)))


# rotate back crotch ease
d = svgwrite.Drawing('test.svg', height='200cm', width='100cm')

p = (pt2cm(pt23),pt2cm(ctrl_pt16),pt2cm(pt21))
pf = fitpath(p,10e-6)
sp = pathtosvg(pf)
path5 = svgwrite.path.Path('%s' % sp)
d.add(path5)

path4 = svgwrite.path.Path('M %s L %s' % (pt2ph(pt6), pt2ph(pt26)))
d.add(path4)

d.save()

paths, attributes = svg2paths('test.svg')

for (T1, seg1, t1), (T2, seg2, t2) in paths[0].intersect(paths[1]):
    pX = paths[0].point(T1)
ptB = pt("B", np.real(pX)/cm, np.imag(pX)/cm)
c_shift = 0.5 # from 1
ptC = pt("C", ptB[0] - c_shift, ptB[1])

# rotate
path6 = svgwrite.path.Path('M %s L %s L %s L %s' % (pt2ph(pt26), pt2ph(ptB), pt2ph(pt21), pt2ph(pt24)), fill="none", stroke=col_pt)
path6.tostring()
gap_sz = 1.5
angle = np.rad2deg(np.arcsin(gap_sz / abs(pt26[0] - ptB[0])))
Path6 = parse_path(path6.attribs['d'].encode('utf-8'))
Path6_rot = Path6.rotated(angle, origin=(complex(pt26[0]*cm,pt26[1]*cm)))
Path6_rot_coord = Path6_rot.d().split(" L ")

ptD = pt("D",     float(Path6_rot_coord[1].split(',')[0])/cm,
                  float(Path6_rot_coord[1].split(',')[1])/cm)
pt21A = pt("21A", float(Path6_rot_coord[2].split(',')[0])/cm,
                  float(Path6_rot_coord[2].split(',')[1])/cm)
pt24A = pt("24A", (float(Path6_rot_coord[3].split(',')[0])/cm),
                  (float(Path6_rot_coord[3].split(',')[1])/cm) + 0.5) # reduced to match front outer

# Force to be quadratic bezier curve
p_23_16_21A = (pt2cm(pt23),pt2cm((pt23[0]+5.0,pt23[1])), pt2cm((ctrl_pt16[0]+1.0,ctrl_pt16[1])), pt2cm((ptD[0]+0.5,ptD[1])), pt2cm(pt21A)) # pt2cm(ptC)
pf_23_16_21A = fitpath(p_23_16_21A,10e-6)
sp_23_16_21A = pathtosvg(pf_23_16_21A)
path2.push('%s' % ' '.join(sp_23_16_21A.split(' ')[2:])) #, fill="none", stroke=col_sew)
path2.push('L %s' % pt2ph(pt24A))
#path2.push('L %s' % pt2ph(pt26))
#path7.translate(trans2)
#dwg.add(path7)

## Wrap up path2
path2.push('Z') # close path
path2.translate(trans2) # translate
dwg.add(path2)

######################################################################
# Add darts:
# also extend to hem line... if hem:...
######################################################################
path_darts = svgwrite.path.Path('M %s L %s' % (pt2ph(pt21A), pt2ph(pt24A)), fill="none", stroke=col_sew)
# path_darts.translate(trans2) # translate
path_darts.tostring()
Path_darts = parse_path(path_darts.attribs['d'].encode('utf-8'))

dist_bw_darts = 6.5
d_darts = (dist_bw_darts + dart_back)*cm/2
t_1 = Path_darts[0].ilength((Path_darts[0].length()/2)-d_darts)
t_2 = Path_darts[0].ilength((Path_darts[0].length()/2)+d_darts)
n1 = Path_darts[0].normal(t_1)
n2 = Path_darts[0].normal(t_2)
n_1 = Line(Path_darts[0].point(t_1), Path_darts[0].point(t_1) + 5*cm*n1*-1)
n_2 = Line(Path_darts[0].point(t_2), Path_darts[0].point(t_2) + 5*cm*n2*-1)

d_1 = svgwrite.path.Path(Path(n_1).d(), fill="none", stroke=col_sew)
d_2 = svgwrite.path.Path(Path(n_2).d(), fill="none", stroke=col_sew)
d_1.translate(trans2)
d_2.translate(trans2)
dwg.add(d_1)
dwg.add(d_2)

# draw triangles
# triangle 1
t_1 = Path_darts[0].ilength((Path_darts[0].length()/2)-d_darts-(dart_back/2.*cm))
t_2 = Path_darts[0].ilength((Path_darts[0].length()/2)-d_darts+(dart_back/2.*cm))
n_1 = Line(Path_darts[0].point(t_1), Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)-d_darts)) + 5*cm*n1*-1)
n_2 = Line(Path_darts[0].point(t_2), Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)-d_darts)) + 5*cm*n1*-1)

d_1_1 = svgwrite.path.Path(Path(n_1).d(), fill="none", stroke=col_sew)
d_1_2 = svgwrite.path.Path(Path(n_2).d(), fill="none", stroke=col_sew)
d_1_1.translate(trans2)
d_1_2.translate(trans2)
dwg.add(d_1_1)
dwg.add(d_1_2)

# triangle 2
t_1 = Path_darts[0].ilength((Path_darts[0].length()/2)+d_darts-(dart_back/2.*cm))
t_2 = Path_darts[0].ilength((Path_darts[0].length()/2)+d_darts+(dart_back/2.*cm))
n_1 = Line(Path_darts[0].point(t_1), Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)+d_darts)) + 5*cm*n1*-1)
n_2 = Line(Path_darts[0].point(t_2), Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)+d_darts)) + 5*cm*n1*-1)

d_2_1 = svgwrite.path.Path(Path(n_1).d(), fill="none", stroke=col_sew)
d_2_2 = svgwrite.path.Path(Path(n_2).d(), fill="none", stroke=col_sew)
d_2_1.translate(trans2)
d_2_2.translate(trans2)
dwg.add(d_2_1)
dwg.add(d_2_2)

# also add pocket position
pocket_width = 13
pocket_w = pocket_width*cm/2
n_1 = Line(Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)-pocket_w)) + 5*cm*n1*-1, Path_darts[0].point(Path_darts[0].ilength((Path_darts[0].length()/2)+pocket_w)) + 5*cm*n1*-1)
d_pocket = svgwrite.path.Path(Path(n_1).d(), fill="none", stroke=col_sew)
d_pocket.translate(trans2)
dwg.add(d_pocket)

######################################################################
# Add grain:
######################################################################
addGrain(path1, shift=0)
addGrain(path2, shift=trans2)
dwg.tostring()

#######################################################################################
# Checks
#######################################################################################
print 'Checks:'
print ''

# waist (subtract darts)
ch_front = seg_10A_11.length() / cm
ch_back = (parse_path('M %s L %s' % (pt2cm(pt21A), pt2cm(pt24A))).length() / cm) - (dart_back * 2)
ch_waist = (ch_front + ch_back) * 2
print 'Item\t\tMeasured\tPattern\t\tError'
print 'Waist\t\t%3.2f\t\t%3.2f\t\t%3.2f' % (waist, ch_waist, waist-ch_waist)

# crotch length
ch_back_crotch  = parse_path('%s' % sp_23_16_21A).length() / cm
ch_front_crotch = parse_path('%s L %s' % (pathtosvg(pf_9_5_F1), pt2ph(pt10A))).length() / cm
print ''
print 'Item\t\tBack\t\tFront\t\tTotal + waistband'
print 'Crotch\t\t%3.2f\t\t%3.2f\t\t%3.2f' % (ch_back_crotch, ch_front_crotch, ch_back_crotch+ch_front_crotch+(waistband*2))

# inner leg match
ch_back_leg_inner  = parse_path('M %s C %s' % (pt2ph(pt28), sp_30_23.split(' C ')[1])).length() / cm
ch_front_leg_inner = parse_path('M %s %s' % (pt2ph(pt13), sp_15_9)).length() / cm
print ''
print 'Item\t\tBack\t\tFront\t\tError\t\tMeasured'
print 'Leg inner\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f' % (ch_back_leg_inner, ch_front_leg_inner, ch_back_leg_inner-ch_front_leg_inner, inside)

# outer leg length
ch_back_leg_outer  = parse_path('M %s C %s L %s' % (pt2ph(pt24A), sp_26_29.split(' C ')[1], pt2ph(pt27))).length() / cm
ch_front_leg_outer = parse_path('%s L %s L %s' % (pathtosvg(pf_11_8), pt2ph(pt14), pt2ph(pt12))).length() / cm
print ''
print 'Item\t\tBack\t\tFront\t\tError\t\tMeasured'
print 'Leg outer\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f' % (ch_back_leg_outer + waistband, ch_front_leg_outer + waistband, ch_back_leg_outer-ch_front_leg_outer, inside + rise)

# seat
# Fix to find intersection point
ch_front_seat  = parse_path('M %s L %s' % (pt2ph(pt6), pt2ph(pt8))).length() / cm
ch_back_seat  = parse_path('M %s L %s' % (pt2ph(pt26), pt2ph(ptB))).length() / cm
print ''
print 'Item\t\tBack\t\tFront\t\tTotal\t\tMeasured'
print 'Seat\t\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f' % (ch_back_seat, ch_front_seat, (ch_back_seat+ch_front_seat)*2, seat)

######################################################################
# Add labels
######################################################################


######################################################################
# Outline with hem
######################################################################
hem = 1.0 * cm # to force hem to outside of path
if include_hem:
    
    paths = [parse_path(path1.attribs['d'].encode("utf-8")[:-1] + ' Q %s %s L %s Z' % (pt2ph((ptF2[0], ptF1[1])), pt2ph(ptF2), pt2ph(pt10B))),
             parse_path(path2.attribs['d'].encode("utf-8")).translated(trans2)]

    for path in paths:
        path_hem = offset_curve(path, hem)
        dwg.add(svgwrite.path.Path(path_hem.d(), stroke=col_hem, fill='none'))

######################################################################
# Save file
######################################################################
dwg.save()

