#!/usr/bin/env python

import nebulous
import math
import random
import scipy
import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
import copy



from shapely.geometry import Polygon
from shapely.geometry import Point
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd


class SaveObject:
    '''Create smaller class to save RJMCMC details if do not want to repeat fit'''
    def __init__(self,likearr=None, shapescoordsave=None,probarr=None,thetasave=None,CoMsave=None,ratiosave=None):
        self.likearr=likearr
        self.shapescoordsave = shapescoordsave
        self.probarr = probarr
        self.thetasave = thetasave
        self.CoMsave = CoMsave
        self.ratiosave = ratiosave



parser = argparse.ArgumentParser(description="This program tries to find the optimal subset of phi and hden values that will be integrated over to get observed line ratios, if fit has been run will plot optimal shape and print model ratios.",
usage = "python RJMCMC_code_test.py --fit True")

parser.add_argument("-f","--fit",dest="fit",default="True",help="if you want to fit for optimal set as True, if only want to analyse end result/pickle set False")


args = parser.parse_args()
global sim
sim = 'gridbasic_3solar_grain'#'grid_3solar_200turbdisp'#'grid_3solar_10000turb'#'grid_3solar_100turbdisp'#'grid_3solar_turb'#'grid_3solar_Nh21'#'gridbasic_3solar_SED14'#
global gridPhiSize
gridPhiSize=33#17#
global gridHdenSize
gridHdenSize=29#15#
linetype = 'emerg'#'intrin'#
global ext
if linetype == 'intrin':
        ext = '.lin'
else:
        ext = '.line'
dphi = 0.25#0.5#
dnh = 0.25#0.5#
phi = np.arange(16,24+dphi,dphi)
phired= phi[range(0,len(phi),4)]
philabel = [str(x) for x in phired]
hdenmin = 7
hdenmax=14
hden = np.arange(hdenmin,hdenmax+dnh,dnh)
hdenred= hden[range(0,len(hden),4)]
hdenlabel = [str(x) for x in hdenred]



#open file in pandas
df =pd.read_csv("FakeData.txt",delimiter='\t')
# extract list of line names
linelabels = df.columns
# if the continuum flux value at 1215A, extract that value and remove it from linelabels
if any('Cont 1215' in s for s in linelabels):
    index = np.where(linelabels == 'Cont 1215')
    cont_1215 = df['Cont 1215'][0]
    cont_1215_err = df['Cont 1215'][1]
    linelabels=np.delete(linelabels,index[0])
else:
    cont_1215 = None
    cont_1215_err = None

flux={}
flux_err={}
# save flux and flux errors to associated line
for line in linelabels:
    flux[line] =df[line][0]
    flux_err[line] =df[line][1]


ratio_data = []
ratio_err_data = []
ratio_label = []

# calculate all possible line ratios, using Monte Carlo method
linelabels_red = linelabels
for line1 in linelabels:
    linelabels_red = linelabels_red[1:]
    for line2 in linelabels_red:
        saveArray = np.zeros(1000)
        for ii in range(1000):
            flux1Trial = flux[line1]+np.random.randn()*flux_err[line1]
            flux2Trial = flux[line2]+np.random.randn()*flux_err[line2]
            saveArray[ii] = flux1Trial/flux2Trial
        sorted = np.sort(saveArray)

        ratio_data.append(sorted[499])
        ratio_err_data.append(0.5*(sorted[839]-sorted[159]))
        ratio_label.append(line1+line2)

# print ratio values. Thes values will be used in RJMCMC fit
print('ratios',ratio_data)
print('err', ratio_err_data)

# If have both Lyman alpha flux and continuum value can fit for the covering factor and
# energy budget using the Lyman alpha EW. EW is calculated using a using Monte Carlo method.
if any('H  1  1216A' in s for s in linelabels) and cont_1215:
    saveArray = np.zeros(1000)
    for ii in range(1000):
        LyaTrial = flux['H  1  1216A']+np.random.randn()*flux_err['H  1  1216A']
        cont1215Trial = cont_1215+np.random.randn()*cont_1215_err
        saveArray[ii] = LyaTrial/cont1215Trial
    sorted = np.sort(saveArray)

    LyaEW = sorted[499]
    LyaEW_err = 0.5*(sorted[839]-sorted[159])
else:
    LyaEW = False
    LyaEW_err = False

print('LyaEW',LyaEW)

# convert argparse entry from string to booleen
Fit = args.fit
if Fit == "False":
    Fit = False
if Fit == "True":
    Fit = True
print('fit',Fit)

# define number of RJMCMC steps and burn in steps
steps = 100000
nburn = 50000

# If have Lyman alpha EW value, will fit for the covering fraction and inclination angle (ndim=2), else only fit for inclination angle (ndim=1)
if LyaEW:
    ndim = 2
else:
    ndim = 1

#setup RJMCMC object
RJMCMCobject = nebulous.RJMCMC(ratio_data=ratio_data, ratio_data_err=ratio_err_data,LyaEW=LyaEW,LyaEW_err=LyaEW_err,linelist=linelabels,ndim=ndim)
RJMCMCobject.loadEQWGrid(sim = sim,phi=phi,hden=hden, ext=ext)

if Fit:
    RJMCMCobject.doRJMCMC( cheat=False,n_iterations = steps)
    saveobject = SaveObject(likearr=RJMCMCobject.likearr, shapescoordsave=RJMCMCobject.shapescoordsave,probarr=RJMCMCobject.probarr,thetasave=RJMCMCobject.thetasave,CoMsave=RJMCMCobject.CoMsave,ratiosave=RJMCMCobject.ratiosave)

    print('acceptedfraction',RJMCMCobject.acceptedsteps/(steps))
    # Save fit to file
    f_myfile = open('test'+sim+'RJMCMC.pickle', 'wb')
    pickle.dump(saveobject, f_myfile)
    f_myfile.close()

else:
    # Read fit from file
    f_myfile = open('test'+sim+'RJMCMC.pickle', 'rb')
    RJMCMCobject = pickle.load(f_myfile)
    f_myfile.close()

# print and plot results

index = np.nonzero(RJMCMCobject.likearr == np.max(RJMCMCobject.likearr))

shapeselect= RJMCMCobject.shapescoordsave[index[0][0]]


coords_true = pickle.load( open( "TestShapeForFakeData.p", "rb" ) )
coords_true.append(coords_true[0])

coords_x,coords_y = nebulous.fndshapeindices(shapevertices=shapeselect)
shapeselect.append(shapeselect[0])
x_best,y_best = zip(*shapeselect)
x_true,y_true = zip(*coords_true)
plt.plot(x_best,y_best)
plt.axis([0,33,0,29])
plt.plot(x_true,y_true,color="black")
plt.scatter(coords_x,coords_y)
plt.xlabel(r"$\Phi_H$")
plt.ylabel(r"$n_H$")
plt.xticks(range(0,len(phi),4),philabel)
plt.yticks(range(0,len(hden),4),hdenlabel)
for ii in range(1,6):
    shapeselect= RJMCMCobject.shapescoordsave[-ii*10]
    shapeselect.append(shapeselect[0])
    x,y = zip(*shapeselect)
    plt.plot(x,y)
plt.savefig('test'+sim+"polygon.eps", format='eps', dpi=1000)
plt.close()


plt.plot(RJMCMCobject.probarr)
plt.savefig('test'+sim+"prob.eps", format='eps', dpi=1000)
plt.close()
plt.plot(RJMCMCobject.likearr)
plt.savefig('test'+sim+"like.eps", format='eps', dpi=1000)
plt.close()
if ndim == 2:
    incl, cf =zip(*RJMCMCobject.thetasave)
    plt.plot(cf)
    plt.savefig('test'+sim+"CF.eps", format='eps', dpi=1000)
    plt.close()
else:
    incl =RJMCMCobject.thetasave
plt.plot(incl)
plt.savefig('test'+sim+"incl.eps", format='eps', dpi=1000)
plt.close()

xvals = []
yvals = []

for point in RJMCMCobject.CoMsave:
    xvals.append(point.x)
    yvals.append(point.y)
n_points_toPlot = 100
#interval = int(np.ceil(len(xvals)/n_points_toPlot))
x = xvals[0::5] ### change back to 50 or [0::interval]
y = yvals[0::5]
c = np.arange(np.size(x))
fig = plt.figure(figsize=(5,5))
ax1 = plt.subplot(111)
cm = plt.get_cmap('viridis')

no_points = len(c)
ax1.set_color_cycle([cm(1.*i/(no_points-1))
                     for i in range(no_points-1)])

for i in range(no_points-1):
    bar = ax1.plot(x[i:i+2],y[i:i+2])

plt.savefig('test'+sim+"CoM.eps", format='eps', dpi=1000)
plt.close()


plt.subplot(2, 1, 1)
plt.plot(x)
plt.ylabel('CoM x-coord')
plt.subplot(2, 1, 2)
plt.plot(y)
plt.xlabel('Iterations')
plt.ylabel('CoM y-coord')
plt.savefig('test'+sim+"CoM_xy.eps", format='eps', dpi=1000)
plt.close()

modelsize = []
for shape in RJMCMCobject.shapescoordsave:
    modelsize.append(np.size(shape)/2)
plt.hist(modelsize)
plt.savefig('test'+sim+"modelsize.eps", format='eps', dpi=1000)
plt.close()


zipped = zip(*RJMCMCobject.ratiosave[nburn:])
zipped=list(zipped)

for jj in range(len(ratio_data)):
    plt.plot(np.sort(zipped[jj]))
    #plt.axhline(y=np.median(zipped[jj]),color='red')
    plt.axhline(y=ratio_data[jj],color='black')
    plt.axhline(y=ratio_data[jj]+ratio_err_data[jj],color='black',linestyle="--")
    plt.axhline(y=ratio_data[jj]-ratio_err_data[jj],color='black',linestyle="--")
    plt.savefig('test'+sim+ratio_label[jj]+"ratio_dist.eps", format='eps', dpi=1000)
    plt.close()



shapeselect.append(shapeselect[0])
cm = plt.get_cmap('viridis')
cNorm  = colors.Normalize(vmin=0, vmax=np.size(RJMCMCobject.shapescoordsave))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
x,y = zip(*shapeselect)
plt.plot(x,y,color='red')
plt.axis([0,33,0,29])
plt.xlabel(r"$\Phi_H$")
plt.ylabel(r"$n_H$")

no_points = int(np.size(RJMCMCobject.shapescoordsave)/1000)
#plt.xticks(range(0,len(phi),4),philabel)
#plt.yticks(range(0,len(hden),4),hdenlabel)
for ii in range(0,np.size(RJMCMCobject.shapescoordsave),100):
    shapeselect= RJMCMCobject.shapescoordsave[ii]
    shapeselect.append(shapeselect[0])
    colorVal = scalarMap.to_rgba(ii)
    x,y = zip(*shapeselect)
    plt.plot(x,y,color=colorVal)

plt.plot(x_true,y_true,color="black")
plt.savefig('test'+sim+"polygon_evolve.eps", format='eps', dpi=1000)
plt.close()


coords_x_save = np.array([])
coords_y_save = np.array([])
for ii in range(nburn,steps):
    shapeselect= RJMCMCobject.shapescoordsave[ii]
    if np.size(shapeselect) < 2:
        print(ii)
        print(shapeselect)
    xx,yy = zip(*shapeselect)
    coords_x,coords_y = nebulous.fndshapeindices(shapevertices=shapeselect)
    if np.size(coords_x):
        if np.max(coords_x) > 33:
            print(ii)
            print('x',xx)
    if np.size(coords_y):
        if np.max(coords_y) > 29:
            print(ii)
            print('y',yy)
    coords_x_save=np.append(coords_x_save,coords_x)
    coords_y_save=np.append(coords_y_save,coords_y)

contours = np.zeros((np.size(phi),np.size(hden)))
for jj in range(np.size(coords_x_save)):
    #print(coords_x_save[jj]-1)
    #print(coords_y_save[jj]-1)
    contours[int(coords_x_save[jj]-1),int(coords_y_save[jj]-1)] += 1

hden_true = 7+0.25*np.array(y_true)
phi_true = 16+0.25*np.array(x_true)

plt.contourf(hden,phi,contours)
plt.plot(hden_true,phi_true,color="black")
plt.ylabel(r"$\Phi_H$")
plt.xlabel(r"$n_H$")
plt.savefig('test'+sim+"polygon_contour.eps", format='eps', dpi=1000)
plt.close()
