import math
import random
import scipy
import numpy as np
import pickle
import argparse
import copy
import re
from tqdm import trange
import logging
logger = logging.getLogger()
logger.setLevel(logging.ERROR)
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.prepared import prep
import itertools
import shapely.speedups

def initial_shape(N=3,phi_size=None,hden_size=None):
    '''Set initial shape. Default is triangle randomly placed in parameter space.
      Kwargs:

            N (int): Number of vertices of intial shape. Default = 3.

            phi_size (int): set limits of phi parameter space, for which to place vertices.

            hden_size (int): set limits of hden parameter space, for which to place vertices.

        Returns:

            coords (list): Coordinates of the initial shape in form, [[x1,y1],[x2,y2],[x3,y3]] for N=3.

    '''
    coords = []
    for w in range(N):
        x = np.random.randint(0,phi_size-1)
        y = np.random.randint(0,hden_size-1)
        coords.append([x,y])
    return coords

def fn(r,sigma_d):
    '''Calculate the  normal  density :math:`f_n(r|0,\sigma_a)`  defined in Luo (2010) Geophysical
       Journal International, Vol:180, Page:1067.

       Args:

            r (float): The distance the vertex is shifted in a within model-like move.

            sigma_d (float) : variance of the normal distribtion from which the distance was chosen.

       Returns:

            fn (float): Normal  density :math:`f_n(r|0,\sigma_a)`.

    '''
    if sigma_d == 0:
        print('r,sigma_d',r,sigma_d)
        return 0.0
    fn = np.exp(-r**2/(2*sigma_d**2))/(sigma_d*np.sqrt(2*np.pi))
    return fn

def checkduplicates(shapevertices=None):
    '''Check that no vertex points of trial shape have duplicate values.

    Kwargs:

        shapevertices (list) : list of vertex coordinates for trial shape.

    Returns:

         (float): returns log prior value of -np.inf if duplicate values are found, else 0.
    '''
    seen = []
    for coords in shapevertices:
        if coords in seen:
            #print('shapecoords',shapevertices)
            return -np.inf
        seen.append(coords)
    return 0

def checkNearNeighbour(shapevertices=None):
    '''Check that no two vertex points are within 0.0001 of each other.

    Kwargs:

        shapevertices (list) : list of vertex coordinates for trial shape.

    Returns:

         (float): returns log prior value of -np.inf if near neighbours are found, else 0.
    '''
    seen = shapevertices[-1]
    for coords in shapevertices:
        if np.sqrt((coords[0]-seen[0])**2+(coords[1]-seen[1])**2)<0.0001:
            #print('shapecoords',shapevertices)
            return -np.inf
        seen = copy.deepcopy(coords)
    return 0.0

def interiorangle(shapevertices=None):
    '''Calculates the angle interior to each vertices, using dot products of the vectors
       before and after each vertex, as :math: `a \cdot b = |a||b| \cos (\\theta_{interior})`

       Kwargs:

        shapevertices (list) : list of vertex coordinates for trial shape ("shape").

       Returns:

        int_angle (list): list of interior angles associated with the given vertices.
       '''

    size = len(shapevertices)
    int_angle = np.zeros(size)
    a = np.array([shapevertices[-1][0]-shapevertices[0][0],shapevertices[-1][1]-shapevertices[0][1]])
    b = np.array([shapevertices[1][0]-shapevertices[0][0],shapevertices[1][1]-shapevertices[0][1]])
    int_angle[0] =  np.arccos(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))
    a = np.array([shapevertices[-2][0]-shapevertices[-1][0],shapevertices[-2][1]-shapevertices[-1][1]])
    b = np.array([shapevertices[0][0]-shapevertices[-1][0],shapevertices[0][1]-shapevertices[-1][1]])

    int_angle[-1] = np.arccos(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))
    for ii in range(1,size-1):
        a = [shapevertices[ii-1][0]-shapevertices[ii][0],shapevertices[ii-1][1]-shapevertices[ii][1]]
        b = [shapevertices[ii+1][0]-shapevertices[ii][0],shapevertices[ii+1][1]-shapevertices[ii][1]]
        int_angle[ii] = np.arccos(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))
    return int_angle

def fndshapeindices(shapevertices=None):
    '''Define the gridpoints that lie within the 'shape' or on its boundary.
       These points will be included in the sum of emissivities.

       Kwargs:

        shapevertices (list) : list of vertex coordinates for trial shape ("shape").

       Returns:

        pCoords (list): coordinates of grid points within or on boundary of shape.
       '''
    polygon = Polygon(shapevertices)
    prepared_polygon = prep(polygon)
    bounds = polygon.bounds
    ll = bounds[:2]
    ur = bounds[2:]
    xs=[]
    ys=[]

    if (ur[0] -ll[0] < 1)  or (ur[1] - ll[1] < 1):
        return None

    xx =range(int(np.ceil(ll[0])), int(np.floor(ur[0]))+1, 1)
    yy = range(int(np.ceil(ll[1])), int(np.floor(ur[1]))+1, 1)

    listpoints = list(itertools.product(xx, yy))
    points = MultiPoint(listpoints)

    try:
        xvals, yvals = zip(*listpoints)
    except:
        print('ll',ll,'ur',ur)
        print('vertices',shapevertices)
        print('xx',xx,'yy',yy)

    for ii in range(len(xvals)):
        if prepared_polygon.intersects(points[ii]):
            xs.append(xvals[ii])
            ys.append(yvals[ii])
    ##plot to check
    #x,y = polygon.exterior.xy
    #plt.plot(x,y)
    #plt.show()
    pCoords = [xs],[ys]
    return pCoords




class RJMCMC:

    '''
    Name: RJMCMC

    Purpose: Find parameter space in hydrogen density and ionising flux that recreates the observed line ratios using a reversible jump Monte Carlo Markov chain algorithm.

    Explanation:

    **Nebulous** is a reversible jump Markov Chain Monte Carlo that takes Cloudy outputs and fits line ratios to find the optimal :math:`\Phi` and :math:`n_H` parameter space the broad line region occupies under a similar assumption of the Locally Optimally Emitting Cloud (LOC) model (Baldwin et al., 1995).

    The programs first calculates line EW values (with respect to continuum flux at 1215 Å) given a line intensity grid from Cloudy and :math:`\Phi` and :math:`n_H` parameter space, described by a shape (that is fit using the reversible jump MCMC). Then uses these EW values to create line ratios. The shape can have :math:`n>3` vertices and occupy any region in the allowed parameter space (example shown in Figure 1).The RJMCMC fit boths the number of vertices, the position of the vertices in :math:`\Phi` and :math:`n_H` space, and the orientation of the BLR.

    Calling Sequence::

        RJMCMCobject = nebulous.RJMCMC(ratio_data=ratio_data, ratio_data_err=ratio_err_data,LyaEW=LyaEW,LyaEW_err=LyaEW_err,linelist=linelabels,ndim=ndim)
        RJMCMCobject.loadEQWGrid(sim = sim,phi=phi,hden=hden )
        RJMCMCobject.doRJMCMC(n_iterations = n_iterations)

    Required inputs:

         ratio_data (float list):

         ratio_data_err (float list):

         LyaEW (float):

         Lya_err (float):

         linelist (str list):

         ndim (int):
    '''
    def __init__(self, ratio_data=None, ratio_data_err=None,LyaEW=None,LyaEW_err=None,linelist=None,ndim=1):

        self.ratio_data = np.array(ratio_data)
        self.ratio_data_err= np.array(ratio_data_err)
        self.LyaEW = LyaEW
        self.LyaEW_err = LyaEW_err
        self.acceptedsteps = 0
        self.linelist = linelist
        self.ndim = ndim


    def loadEQWGrid(self,sim = None, ext=None,phi=None,hden=None):
        '''
        Define the EQW grid parameters and load total and inward relative line intensity (~EQW) grids for each line in line list.
        Kwargs:

            phi (array): array of log(hydrogen ionising flux) vales covered in Cloudy simulation.

            hden (array): array of log(hydrogen density) values covered in Cloudy simulation.

            sim (str): Filename the relative line intensity (~EQW) values calculated in Cloudy simulation are stored in.

            ext (str): Line extension for file name, either "lin" for intrinsic line intensities or "line" for emergent line intensities.

        '''
        self.phi =phi
        self.hden = hden
        self.ext = ext
        self.sim = sim

        self.dphi =phi[1]-phi[0]
        self.dnh =hden[1]=hden[0]
        self.gridPhiSize = len(phi)
        self.gridHdenSize = len(hden)


        self.EQWttl = {}
        self.EQWinwd = {}
        for line in self.linelist:
            split = line.split()
            inwd = "INWD  "+split[-1]
            self.EQWttl[line] = self.EQWgrid(line=line)
            self.EQWinwd[line] = self.EQWgrid(line=inwd)

        return


    def EQWgrid(self,line=None):
        '''
        Read the EQW grid from cloudy outputs for given line, cloudy simulation and file extension.
        Program opens file, find column header that matches line name, extracts the stored values
        and reshapes vector into appropriate grid size.

        Kwargs:

            line (str): Name of line,

        Returns:

            grid (array) : Array of line intensities, for given line, calculated from Cloudy (Ferland et al., 2013) for clouds with a range of hydrogen gas densities and incident by a range of hydrogen ionisation fluxes.

        '''
        EWFileName = self.sim+self.ext
        with open(EWFileName) as f:
            first_line = f.readline()

        first_line = re.sub('\n', '', first_line)

        header = first_line.split("\t")
        linelist = header[1:]
        symbol =np.where(np.array(linelist)==line)
        symbol = symbol[0][0]

        EW = np.loadtxt(EWFileName,dtype='float',skiprows=1,usecols=(symbol+1,)) # open lin file, read out EW/EW(Lyalpha) values for line of interest
        grid =np.resize(EW,(self.gridPhiSize,self.gridHdenSize))

        return grid
#
    def sum_shape(self,alpha=(-1.0,-1.0),contained_indices=None,EW=None):
        '''
        Programs calculates line EW values (with respect to continuum flux at 1215 Å) given a EQW grid and the shape indices.
        In the LOC approach, the observed emission line spectrum is the sum of emission line contributions from a weighted distribution of 'clouds' with a range of gas densities at a range of radii. The resulting line strength is given as,

        .. math:: L_{line} \propto \int_{r_{min}}^{r_{max}} \int_{n_{H,min}}^{n_{H,max}} W_{1215}(r,n_H) f(r) g(n_H) dn_H dr,

        where :math:`f(r) = r^{\Gamma}`, :math:`g(n_H) = n_H^{\\beta}` are the weighting functions or 'covering fraction'
        of the various clouds (can be thought of as the an number density of clouds at radius :math:`r` and density :math:`n_H`),
        :math:`r_{ min}` and :math:`r_{ max}` are the minimum and maximum radii of the BLR, and :math:`n_{ H,min}` and :math:`n_{ H,max}`
        are the minimum and maximum cloud densities considered. :math:`W_{1215}` is the equivalent width of the line referenced
        to the incident continuum at 1215Å.
        This expression can be rewritten in term of :math:`\log \Phi_H` and :math:`\log n_H`.

        .. math:: L_{line} \propto  \int_{\log \Phi_{H,min}}^{\log \Phi_{vmax}} \int_{\log n_{H,min}}^{\log n_{H,max}} W_{1215}(\log \Phi,\log n_H) 10^{(\\beta+1)\log n_H-0.5(\Gamma+1)\log \Phi} d \log n_H d \log \Phi.

        and the :math:`EW_{line}` is given as

        .. math:: EW_{line} \propto  \\frac{\int_{\log \Phi_{H,min}}^{\log \Phi_{vmax}} \int_{\log n_{H,min}}^{\log n_{H,max}} W_{1215}(\log \Phi,\log n_H) 10^{(\\beta+1)\log n_H-0.5(\Gamma+1)\log \Phi} d \log n_H d \log \Phi}{\int_{\log \Phi_{H,min}}^{\log \Phi_{vmax}} \int_{\log n_{ H,min}}^{\log n_{ H,max}} 10^{(\\beta+1)\log n_H-0.5(\Gamma+1)\log \Phi} d \log n_H d \log \Phi}.


        Kwargs:

            alpha (float list): (:math:`gamma`, :math:`beta`) power law indices for weighted distribution of 'clouds' of various gas densities and ionising fluxes used in LOC sum. Set at (-1.0,-1.0) but could be changes manually or eventually fit with RJMCMC.

            contained_indices (int list): coordinates of grid points within or on boundary of shape.

            EW (array): Array of line intensities calculated from Cloudy (Ferland et al., 2013) for clouds with a range of hydrogen gas densities and incident by a range of hydrogen ionisation fluxes.

        Returns:

            Lline1 (float) : Predicted EW (relative to 1215Å incident light) of given line using LOC-like methodology within prescribed shape.
        '''
        gamma,  beta = alpha
        if contained_indices[0]:
            l1_w1215 = np.zeros(len(contained_indices[0]))
            cf_1 = np.zeros(len(contained_indices[0]))
            for ii in range(len(contained_indices[0])):
                scale =     self.phi[contained_indices[0][ii]]**(-0.5*(gamma+1))*self.hden[contained_indices[1][ii]]**(beta+1)*self.dphi*self.dnh
                l1_w1215[ii]=EW[contained_indices[0][ii],contained_indices[1][ii]]*scale
                cf_1[ii] =  scale
            Lline1 = np.sum(l1_w1215)/np.sum(cf_1)
            return Lline1
        else:
            return np.inf

    def line_fluxes(self, linelist=None,incl=None,cf=None):
        '''
        Calculate the predicted EW values (with respect to ionizing flux at 1215Å) for the lines in the line list.
        Kwargs:

            linelist (list str) : List of line names

            incl (float): inclination angle of the broad line region

            cf (float): covering fraction of the gas.

        Returns:

            model_data (dict) : EW values (with respect to ionizing flux at 1215Å) for the lines in the line list

        '''
        contained1 = fndshapeindices(shapevertices=self.shapeTrial)

        contained= [[],[]]
        gamma = -1.
        beta = -1.
        contained[0] = contained1[0][0]
        contained[1] = contained1[1][0]
        EQW_model = {}
        model = {}
        model_data = {}
        for line in self.linelist:
            EQW_model[line] = 0.5*self.EQWttl[line]*(1+np.cos(incl))-self.EQWinwd[line]*np.cos(incl)

            model_data[line] = self.sum_shape(alpha = (gamma,beta),contained_indices=contained,EW=EQW_model[line])*cf
        return model_data

    def movepoint(self,index=None,ca=None):
        '''
        Move the position of the index-th vertex to create trial shape shapeTrial. The shift is determined by drawing two random numbers
        for distance :math:`r` and direction :math:`\\Theta`. To help occupancy levels of trial jumps, the step size of the
        random walk is restricted by a variance tied to the length of the intersecting sides of the polygon.
        The random distance :math:`r` is drawn from a normal distribution, :math:`f_N(0,\sigma_a)` with a mean of zero and
        variance :math:`\sigma_d = min(d_i^-,d_i^+)c_a` where :math:`d_i^-` and :math:`d_i^+` are the lengths of the two polygon
        lines intersecting at chosen vertex and :math:`c_a` is a constant.
        Kwargs:

         index (int) : Index of vertice to move

         ca (float) : Value controlling the distance a vertices can be moved. The default value of ca is set to 0.25 to minimise self intersecting polygons.

        Returns:

         self.shapetrial (float list) : Trial shape.

         r (float) : Distance the vertex is shifted. Used calculating the acceptance probability.

         sigma_d (float) : Variance of the normal distribtion from which the distance was chosen. Used calculating the acceptance probability.
        '''

        vertex_0 =self.shapeTrial[index]
        if index != len(self.shapeTrial)-1:
            vertex_plus =self.shapeTrial[index+1]
        else:
            vertex_plus =self.shapeTrial[0]
        if index != 0:
            vertex_minus =self.shapeTrial[index-1]
        else:
            vertex_minus =self.shapeTrial[-1]
        d_plus = np.sqrt((vertex_0[0]-vertex_plus[0])**2+(vertex_0[1]-vertex_plus[1])**2)
        d_minus = np.sqrt((vertex_0[0]-vertex_minus[0])**2+(vertex_0[1]-vertex_minus[1])**2)
        sigma_d = np.min([d_plus,d_minus])*ca
        r = np.random.randn(1)*sigma_d
        theta= np.random.random(1)*2*np.pi
        dx = r*np.sin(theta)
        dy = r*np.cos(theta)
        self.shapeTrial[index][0] = float(self.shapeTrial[index][0]+dx)
        self.shapeTrial[index][1] = float(self.shapeTrial[index][1]+dy)
        return {'r':r,'sigma_d':sigma_d}

    def RJMCMC_deathmove(self,shapevertices=None,ca=0.25):
        '''
        Perform death-model move - - proposed model change from :math:`\mathcal{M}_k` to :math:`\mathcal{M}_{k-1}`,
        allowable for :math:`k>3`. The program deletes a randomly chosen vertex and forms a new side by
        joining the two neighbouring vertices.

        Kwargs:

         ca (float) : In birth and within model moves this value controls the distance a vertices can be moved.
         The default value of c_a is set to 0.25 to minimise self intersecting polygons. However in death move influences value of sigma_d.

        Returns:

         self.shapetrial (float list) : List of vertices of the trial shape.

         r (float) : The outputed radius :math:`r` in this case is the distance from the deleted
         vertex to the middle point ofthe new polygon side. Used calculating the acceptance probability.

         sigma_d (float) : :math:`min(d_i^-,d_i^+)c_a` where :math:`d_i^-` and :math:`d_i^+` are the distances between the
         deleted vertex and its preceding and proceding vertex, respectively. Used calculating the acceptance probability.
        '''
        self.shapeTrial = copy.deepcopy(self.shapeKeep)
        index =np.random.randint(len(self.shapeTrial))
        vertex_0 =self.shapeTrial[index]
        if index != len(self.shapeTrial)-1:
            vertex_plus =self.shapeTrial[index+1]
        else:
            vertex_plus =self.shapeTrial[0]
        if index != 0:
            vertex_minus =self.shapeTrial[index-1]
        else:
            vertex_minus =self.shapeTrial[-1]
        d_plus = np.sqrt((vertex_0[0]-vertex_plus[0])**2+(vertex_0[1]-vertex_plus[1])**2)
        d_minus = np.sqrt((vertex_0[0]-vertex_minus[0])**2+(vertex_0[1]-vertex_minus[1])**2)
        sigma_d = np.min([d_plus,d_minus])*ca
        midpoint_x = (vertex_minus[0]+vertex_plus[0])*0.5
        midpoint_y = (vertex_minus[1]+vertex_plus[1])*0.5
        midpoint = [midpoint_x,midpoint_y]
        r = np.sqrt((vertex_0[0]-midpoint[0])**2+(vertex_0[1]-midpoint[1])**2)
        del(self.shapeTrial[index])
        return {'r':r,'sigma_d':sigma_d}


    def RJMCMC_deathmove_accept(self,movedetails=None,bk=None, dk1=None):
        '''
        Determines whether move will be accepted and updates appropriate values (self.thetaKeep,
        self.probKeep, self.shapeKeep,self.ratiosKeep, self.likeKeep, self.acceptedsteps) accordingly.

        The final acceptance probability a death move is

        .. math:: \\alpha_{ birth~move} = \min\{1,\\frac{\pi(\\theta_k,\mathcal{M}_k)}{\pi(\\theta_{k+1},\mathcal{M}_{k+1})}\\frac{\pi(y|\\theta_k,\mathcal{M}_k)}{\pi(y|\\theta_{k+1},\mathcal{M}_{k+1})}\\frac{b_k(k+1)f_n(r|0,\sigma_a)}{2\pi d_{k+1}kr}\},


        Kwargs:

        movedetails (dict): contains value *r* and *sigma_d* of move. See *RJMCMC_deathmove* docstring for detailed explanation.

        bk (float): :math:`b_k` is the probabilty of attempting the birth move from :math:`\mathcal{M}_k` to :math:`\mathcal{M}_{k+1}`

        dk1 (float): :math:`d_{k+1}`is the probability of attempting a death move from :math:`\mathcal{M}_{k+1}` to :math:`\mathcal{M}_k`.

        Returns:

            output (booleen): True if trial move was accepted, False if trial move was rejected.
        '''
        r = movedetails['r']
        sigma_d =movedetails['sigma_d']
        if self.ndim==2:
            sigma_theta = np.array([0.1,0.05])
        else:
            sigma_theta = np.array([0.1])
        Z = np.random.normal(size=np.size(self.thetaKeep))
        self.thetaTrial = Z*sigma_theta+self.thetaKeep
        lnprob_trial = self.lnprob_RJMCMC()
        proposalratio = bk*(len(self.shapeTrial)+1)*fn(r,sigma_d)/(2*np.pi*dk1*len(self.shapeTrial)*r)
        rand = random.random()
        if rand == 0:
            rand = random.random()

        if np.log(rand)<(lnprob_trial-self.probKeep+np.log(proposalratio)):#if np.log(rand)<(lnprob_trial-init_prob):#
            outcome =True
            self.thetaKeep = copy.deepcopy(self.thetaTrial)
            self.probKeep =copy.deepcopy(lnprob_trial)
            self.shapeKeep =copy.deepcopy(self.shapeTrial)
            self.ratiosKeep = copy.deepcopy(self.ratiosTrial)
            self.likeKeep = copy.deepcopy(self.likeTrial)
            self.acceptedsteps +=1.
        else:
            outcome =False

        return outcome

    def RJMCMC_birthmove(self,ca=0.25):
        '''
        Perform birth-model move - a proposed model change from :math:`\mathcal{M}_k`
        to :math:`\mathcal{M}_{k+1}`, allowable for :math:`k<k_{max}`. The program splits a randomly
        chosen side into two at the middle, and make a within-model move from the middle
        point (using the movepoint function).

        Kwargs:

         ca (float) : In birth and within model moves this value controls the distance a vertices can be moved.
         The default value of c_a is set to 0.25 to minimise self intersecting polygons.
        '''
        index =np.random.randint(len(self.shapeKeep))
        vertex_minus =self.shapeKeep[index]
        if index != len(self.shapeKeep)-1:
            vertex_plus =self.shapeKeep[index+1]
        else:
            vertex_plus =self.shapeKeep[0]
        midpoint_x = (vertex_minus[0]+vertex_plus[0])*0.5
        midpoint_y = (vertex_minus[1]+vertex_plus[1])*0.5
        vertex_0 = [midpoint_x,midpoint_y]
        self.shapeTrial = []
        for ii in range(len(self.shapeKeep)+1):
           if ii <= index:
               self.shapeTrial.append(self.shapeKeep[ii])
           if ii ==index+1:
               self.shapeTrial.append(vertex_0)
           if ii>index+1:
               self.shapeTrial.append(self.shapeKeep[ii-1])
        return self.movepoint(index=index+1,ca=ca)


    def RJMCMC_birthmove_accept(self,movedetails=None,bk=None, dk1=None):
        '''
        Determines whether move will be accepted and updates appropriate values (self.thetaKeep,
        self.probKeep, self.shapeKeep,self.ratiosKeep, self.likeKeep, self.acceptedsteps) accordingly.

        The final acceptance probability a birth move is

        .. math:: \\alpha_{birth~move} = \min\{1,\\frac{\pi(\\theta_{k+1},\mathcal{M}_{k+1})}{\pi(\\theta_k,\mathcal{M}_k)}\\frac{\pi(y|\\theta_{k+1},\mathcal{M}_{k+1})}{\pi(y|\\theta_k,\mathcal{M}_k)}\\frac{2\pi d_{k+1}kr}{b_k(k+1)f_n(r|0,\sigma_a)}\},


        Kwargs:

        movedetails (dict): contains value *r* and *sigma_d* of move. See *movepoint* docstring for detailed explanation.

        bk (float): :math:`b_k` is the probabilty of attempting the birth move from :math:`\mathcal{M}_k` to :math:`\mathcal{M}_{k+1}`

        dk1 (float): :math:`d_{k+1}`is the probability of attempting a death move from :math:`\mathcal{M}_{k+1}` to :math:`\mathcal{M}_k`.

        Returns:

            output (booleen): True if trial move was accepted, False if trial move was rejected.
        '''
        r = abs(float(movedetails['r']))
        sigma_d =float(movedetails['sigma_d'])
        if self.ndim==2:
            sigma_theta = np.array([0.1,0.05])
        else:
            sigma_theta = np.array([0.1])
        Z = np.random.normal(size=np.size(self.thetaKeep))
        self.thetaTrial = Z*sigma_theta+self.thetaKeep
        lnprob_trial = self.lnprob_RJMCMC( )
        proposalratio = (2*np.pi*dk1*(len(self.shapeTrial)-1)*r)/bk*(len(self.shapeTrial))*fn(r,sigma_d)
        rand = random.random()

        if np.log(rand)<(lnprob_trial-self.probKeep+np.log(proposalratio)):
            outcome =True
            self.thetaKeep = copy.deepcopy(self.thetaTrial)
            self.probKeep =copy.deepcopy(lnprob_trial)
            self.shapeKeep = copy.deepcopy(self.shapeTrial)
            self.ratiosKeep = copy.deepcopy(self.ratiosTrial)
            self.likeKeep = copy.deepcopy(self.likeTrial)
            self.acceptedsteps +=1.
        else:
            outcome =False

        return outcome

    def RJMCMC_withinmove(self,ca=0.25):
        '''
        Perform Within-model move - a proposed move of one vertex to a nearby location.
        The program selects which vertex to move then performs move using movepoint function.

        Kwargs:

         ca (float) : In birth and within model moves this value controls the distance a vertices can be moved.
         The default value of c_a is set to 0.25 to minimise self intersecting polygons.
        '''
        self.shapeTrial = copy.deepcopy(self.shapeKeep)
        index =np.random.randint(len(self.shapeTrial))
        return self.movepoint(index=index,ca=ca)



    def RJMCMC_withinmove_accept(self,movedetails=None):
        '''
        Determines whether move will be accepted and updates appropriate values (self.thetaKeep,
        self.probKeep, self.shapeKeep,self.ratiosKeep, self.likeKeep, self.acceptedsteps) accordingly.

        The final acceptance probability a within model move is

        .. math:: \\alpha_{within-model} = \min\{1,\\frac{\pi(\\theta_j,\mathcal{M}_j)}{\pi(\\theta_i,\mathcal{M}_i)}\\frac{\pi(y|\\theta_k',\mathcal{M}_k)}{\pi(y|\\theta_k,\mathcal{M}_k)}\},


        Kwargs:

            movedetails (dict): contains value *r* and *sigma_d* of move. See *movepoint* docstring for detailed explanation.

        Returns:

            output (booleen): True if trial move was accepted, False if trial move was rejected.
        '''
        r = movedetails['r']
        sigma_d =movedetails['sigma_d']
        if self.ndim==2:
            sigma_theta = np.array([0.1,0.05])
        else:
            sigma_theta = np.array([0.1])
        Z = np.random.normal(size=np.size(self.thetaKeep))
        self.thetaTrial = Z*sigma_theta+self.thetaKeep
        lnprob_trial = self.lnprob_RJMCMC()
        rand = random.random()
        if np.log(rand)<(lnprob_trial-self.probKeep):
            outcome =True
            self.thetaKeep = copy.deepcopy(self.thetaTrial)
            self.probKeep =copy.deepcopy(lnprob_trial)
            self.shapeKeep = copy.deepcopy(self.shapeTrial)
            self.ratiosKeep = copy.deepcopy(self.ratiosTrial)
            self.likeKeep = copy.deepcopy(self.likeTrial)
            self.acceptedsteps +=1.
        else:
            outcome =False
        #    newtheta =copy.deepcopy(thetaarr[:])
        #    self.propKeep =copy.deepcopy(init_prob)

        return outcome

    def lnlike_RJMCMC(self):
        """
        Calculate the ln(likelihood) of the given model line ratios.

        .. math:: ln(L) = -0.5 \sum_n {[ \\frac{(ratios_{obs}-ratios_{model})^2}{error_{ratios}^2} + ln(2 \pi \sigma^2) ]}

        Returns:

            ln_l (float): Sum of ln(likelihood) values
        """
        if self.ndim==2:
            incl, cf = self.thetaTrial
        else:
            incl = self.thetaTrial
            cf = 0.2 #default value


        model_data = self.line_fluxes(linelist=self.linelist,incl=incl,cf=cf)

        ratio_model = []
        if model_data[self.linelist[0]] ==0.0:
            print("model_data",model_data)

        linelist_red = self.linelist
        for line1 in self.linelist:
            linelist_red = linelist_red[1:]
            for line2 in linelist_red:
                ratio_model.append(model_data[line1]/model_data[line2])


        ratio_model = np.array(ratio_model)
        ratio_model = np.nan_to_num(ratio_model)


        chi2ratios = sum((self.ratio_data-ratio_model)**2/self.ratio_data_err**2)

        self.ratiosTrial = ratio_model
        if np.isnan(chi2ratios):
            print('ratio_model',ratio_model)
        ## need to fix error estimation
        if  (self.LyaEW):
            chi2EW = (self.LyaEW-model_data['H  1  1216A'])**2/self.LyaEW_err**2

        if (self.LyaEW):
            return -0.5*(chi2ratios+chi2EW)
        else:
            return -0.5*(chi2ratios)

    def lnprior_RJMCMC(self):
        '''
        Calculate the ln(priors) of the given model.

        The priors include:
         * Flat prior on the extent of the vertices that restricts them with the :math:`\phi` and :math:`n_H` parameter space.
           The prior does not impose any constraints on the length scale and location of the polygon.
           However, we require that the polygon is contained in the limits of our Cloudy simulations:
           :math:`16\geq\log\Phi\geq24` and :math:`7\geq\log n_H\geq14`.
         * Prior to stop duplicate or vertices that are too close to each other
         * Flat prior on inclination angle betweeo :math:`0<\\theta<\pi`
         * Flat prior on covering fraction between :math:`0<cf<1.0`. A covering fraction greater than 1 breaks the energy budget.
         * Prior on the polygon shape. We adapt the prior model similar to that proposed by Pievatolo & Green (1998)

           .. math:: \pi(\mathcal{M}_k,\\theta_k) \propto \exp(-\\alpha k^{\gamma}-\\frac{\\beta}{k}\sum_{i=1}^k[\phi_i(\\theta_k)-\omega_k]^2), k\geq 3

           where :math:`\gamma\geq1` and penalises more complex models (i.e. models with a higher number of
           vertices), :math:`\phi_i(\omega_k)` is the angle in radians interior to the i-th vertex of polygon
           define by parameter :math:`\\theta_k` and  :math:`\omega_k = (k-2)\pi/k` and is the interior angle of a
           regular-sided polygon. Therefore this prior penalises more complex models with non-regular shapes.
           :math:`\\alpha` and :math:`\\beta` can be altered depending if you want to penalise more vertices
           or irregular shapes. Defaults are :math:`\\alpha=3` and :math:`\\beta = 0.9`.

        Returns:

            ln_p (float): Sum of ln(prior) values
        '''
        polygon = Polygon(self.shapeTrial)
        if self.ndim==2:
            incl, cf = self.thetaTrial
        else:
            incl = self.thetaTrial
            cf = 0.1 # arbitrary value
        k = len(self.shapeTrial)

        #### These values can be tweaked
        gamma_fit = 1.6
        alpha= 10.
        beta = 1.
        xx, yy = zip(*self.shapeTrial)
        if  0<=incl<=1.5*math.pi and np.min(xx)>=0 and np.max(xx)<=33 and np.min(yy)>=0 and np.max(yy)<=29:
            p_theta = 0.0
        else:
            p_theta =-np.inf

        if  0<cf<=1.0:
            p_cf = 0.0
        else:
            p_cf =-np.inf

        test = fndshapeindices(shapevertices=self.shapeTrial)
        if test:
            p_size = 0.0
        else:
            p_size = -np.inf

        omega = (k-2)*np.pi/k
        p_model = -beta*k**gamma_fit-(alpha/k)*np.sum((interiorangle(shapevertices=self.shapeTrial)-omega)**2)

        if np.size(xx)<2:
            p_bounds =-np.inf
        else:
            p_bounds=0.0
        if polygon.is_valid:#.is_simple:
            p_polygon = 0.0
        else:
            p_polygon =-np.inf
        p_duplicate = checkduplicates(shapevertices=self.shapeTrial)
        p_nearneighbour =checkNearNeighbour(shapevertices=self.shapeTrial)

        return p_theta+p_model+p_bounds+p_duplicate+p_nearneighbour+p_cf+p_size+p_polygon


    def lnprob_RJMCMC(self):
        """
        Return the logarithm of the posterior function, to be passed to the RJMCMC
        sampler.

        Returns:

            posterior (float): log posterior which is a combination of the model likelihood and priors.
        """
        lp = self.lnprior_RJMCMC()
        #print('lp',lp)
        if not np.isfinite(lp):
            return -np.inf

        like = self.lnlike_RJMCMC()
        self.likeTrial = copy.deepcopy(like)
        #print('like',like,self.likeTrial)
        #print('prior',lp)
        if not np.isnan(lp + like):
            return lp + like
        else:
            print('lp',lp)
            print('like',like)
            print('Error: nan probability value')
            return -np.inf
            #exit()



    def RJMCMCstepMH(self):
        '''
        Perform reversible jump MCMC step. Three options: within move, birth move and death move. Death move can only happen when there is
        more that 3 vertices. The probability of a within-model move is :math:`w_k= 1-(b_k+d_k)` and we have chosen :math:`b_k=d_{k+1}=d_k=0.1` and :math:`w_k=0.8`
        for our investigations except in the case where :math:`k=3`, then :math:`b_k=0.1`, :math:`d_{k}=0`, :math:`w_k =0.9` and :math:`d_{k+1}=0.1`.
        '''


        Z = np.random.random(1)
        shapesize=len(self.shapeKeep)

        if shapesize >3:
            if Z<0.8:
                newshapedetails = self.RJMCMC_withinmove(ca=0.25)
                accept = self.RJMCMC_withinmove_accept(movedetails=newshapedetails)
            if 0.8<=Z<0.9:
                newshapedetails = self.RJMCMC_birthmove(ca=0.25)
                accept = self.RJMCMC_birthmove_accept(movedetails=newshapedetails,bk=0.1, dk1=0.1)
            if 0.9<=Z<1:
                newshapedetails = self.RJMCMC_deathmove(ca=0.25)
                accept = self.RJMCMC_deathmove_accept(movedetails=newshapedetails,bk=0.1, dk1=0.1)
        else:
            if Z<0.9:
                newshapedetails = self.RJMCMC_withinmove(ca=0.25)
                accept = self.RJMCMC_withinmove_accept(movedetails=newshapedetails)
            if 0.9<=Z<1:
                newshapedetails = self.RJMCMC_birthmove(ca=0.25)
                accept = self.RJMCMC_birthmove_accept(movedetails=newshapedetails,bk=0.1, dk1=0.1)

        return

    def doRJMCMC(self, cheat=False,n_iterations = 100000):
            """
            Run RJMCMC.

            Args:

               cheat (bool): If true, will trial 1000 random inital triangles across the parameter space
                to find one with highest posterior value. This will minimise burnin time. However, I have found it to get stuck in small parameter spaces and without it gets more consistent fits. Default is False.

                n_iteratins (int): Number of iterations to pass to the RJMCMC. Default is 100000.
            """
            shapely.speedups.enable()
            self.probarr = []
            self.likearr = []
            self.shapescoordsave = []
            self.ratiosave = []
            self.thetasave = []
            self.CoMsave = []


            # set up initial guesses
            initialshape = initial_shape(N=3,phi_size=self.gridPhiSize,hden_size=self.gridHdenSize)
            pos = np.random.random(self.ndim)
            if self.ndim==2:
                pos[0] = math.pi
                pos[1] = pos[1]*0.5
            else:
                pos[0] = math.pi
            self.thetaKeep = pos
            self.thetaTrial = pos
            self.shapeKeep = copy.deepcopy(initialshape)
            self.shapeTrial = copy.deepcopy(initialshape)
            # find best triangle to start fit from to help reduce burnin time
            if cheat:
                prob_init = self.lnprob_RJMCMC()
                for nn in range(1000):
                    initialshape_trial = initial_shape(N=3,phi_size=self.gridPhiSize,hden_size=self.gridHdenSize)
                    self.shapeTrial = copy.deepcopy(initialshape_trial)
                    prob_init_trial = self.lnprob_RJMCMC()
                    if prob_init_trial > prob_init:
                        prob_init = copy.deepcopy(prob_init_trial)
                        initialshape = copy.deepcopy(initialshape_trial)
            else:
                initialshape = [[1,1],[self.gridPhiSize-2,1],[self.gridPhiSize-2,self.gridHdenSize-2],[1,self.gridHdenSize-2]]


            pos = np.random.rand(self.ndim)
            if self.ndim==2:
                pos[0] = math.pi
                pos[1] = pos[1]*0.5
            else:
                pos[0] = math.pi
            self.shapeKeep = copy.deepcopy(initialshape)
            self.shapeTrial = copy.deepcopy(initialshape)
            self.probKeep = self.lnprob_RJMCMC()
            self.ratiosKeep = self.ratiosTrial
            self.likeKeep = self.likeTrial
            self.likearr.append(self.likeKeep)
            self.ratiosave.append(self.ratiosKeep)
            self.probarr.append(self.probKeep)
            self.shapescoordsave.append(copy.deepcopy(self.shapeKeep))
            polygon = Polygon(self.shapeKeep)
            self.CoMsave.append(polygon.centroid)
            self.thetasave.append(self.thetaKeep)


            for kk in trange(1,n_iterations):
                newvals = self.RJMCMCstepMH()
                self.thetasave.append(self.thetaKeep)
                self.shapescoordsave.append(copy.deepcopy(self.shapeKeep))
                self.probarr.append(self.probKeep)
                polygon = Polygon(self.shapeKeep)
                self.CoMsave.append(polygon.centroid)
                self.likearr.append(self.likeKeep)
                self.ratiosave.append(self.ratiosKeep)


            return
