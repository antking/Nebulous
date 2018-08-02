#!/usr/bin/env python
import numpy as np
import nebulous
import pickle

sim = 'gridbasic_3solar_grain'
gridPhiSize=33
gridHdenSize=29
linetype = 'emerg'#'intrin'#

if linetype == 'intrin':
        ext = '.lin'
else:
        ext = '.line'

dphi = 0.25
dnh = 0.25
phi = np.arange(16,24.1,dphi)
phired= phi[range(0,len(phi),4)]
philabel = [str(x) for x in phired]
hdenmin = 7
hdenmax=14
hden = np.arange(hdenmin,hdenmax+0.1,dnh)
hdenred= hden[range(0,len(hden),4)]
hdenlabel = [str(x) for x in hdenred]


linelabels = ['H  1  6563A','H  1  4861A', 'TOTL  2798A', 'TOTL  1549A']
inclination = np.pi
coveringFrac = 0.2

xlist = [2,8,12,11,5,6] #phi
ylist = [19,17,19,27,26,22] #hden
coords = []
N = len(xlist)
for w in range(N):
    x = xlist[w]
    y = ylist[w]
    coords.append([x,y])


RJMCMCobject = nebulous.RJMCMC(linelist=linelabels)
RJMCMCobject.loadEQWGrid(sim = sim,phi=phi,hden=hden, ext=ext )
RJMCMCobject.shapeTrial = coords
linevalues = RJMCMCobject.line_fluxes(linelist=linelabels,incl=inclination,cf = coveringFrac)

errors = [0.01*x for x in linevalues.values()]
import csv

with open('FakeData.txt','w') as f:
    w = csv.writer(f,delimiter='\t')
    w.writerow(linevalues.keys())
    w.writerow(linevalues.values())
    w.writerow(errors)

with open('TestShapeForFakeData.p', 'wb') as f:
            pickle.dump(coords, f, protocol=2)
