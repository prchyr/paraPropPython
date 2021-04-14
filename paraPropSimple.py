# paraPropPython
# s. prohira

#This package is a simple parabolic equation solver for EM fields beneath the surface of the ice.
#References for PE: Levy (2000), Tappert (1977), Feit, Fleck (1976), Brock et.al (1977), Thompson et. al (1983), Apaydin et.al (2017)
#References for RF ice properties: Besson et. al (2007), Barwick et. al (2018)
#Data citations for ice density data: SPICE-Winski et. al (2019), WAIS-Kreutz et. al (2011) 

# GPL v3, see LICENSE

import util
import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.constants as constant
import scipy.io.wavfile as wav
import scipy.signal as sig
import scipy.interpolate as interp
from scipy.signal import butter, lfilter
from numpy import linalg as la
import pickle
import csv
import os

class paraPropSimple:
    """define a paraPropSimple object. 
    used for propagation beneath the ice surface from a dipole source. 
    """
    def __init__(self, freq, xMax, zMax, sourceDepth=30, polarization=0, site="SP", method="II", nProfile="functional", dataType="", debug=0):
        """freq is frequency
    xMax is the maximum x dimension in meters
    zMax is the maximum z depth in meters (positive value for depth. larger value is deeper)
    sourceDepth is the depth of the source in meters
    polarization is the polarization of the source. 0 is vertical (along z) 1 is horizontal (perpindicular to z). currently, the polarization is only vertical (for now)
    site is either \"TD\" for Taylor Dome or \"SP\" for the south pole
    method is the method of solving. "II" is the in-ice model (default, recommended) "WA" is the wide-angle method, "FF" is Feit and Fleck (another wide-angle method)
    nProfile is either \"functional\" or \"data\" for a smooth n(z) fit or a data-driven profile
    debug is a setting that will print out all sorts of stuff."""

        self.debug=debug
        self.path=os.path.dirname(os.path.abspath(__file__))
        self.setVals(freq, xMax, zMax, sourceDepth, polarization,  site, method,nProfile, dataType)
        

    def setVals(self, freq, xMax, zMax, sourceDepth=30, polarization=0,  site="SP", method="II",nProfile="functional", dataType=""):
        """set the values. can be used for updating once a paraPropSimple object has been created"""

        if(freq==0):
            freq=0.001;
        self.freq=freq
        self.xMax=xMax
        self.zMax=zMax
        self.site=site
        nZ=int(zMax*20)
        self.nStepsPerMeter=20
        nX=500;
        if (self.debug>0):#equal resolution in x and z DO NOT USE, will get big and slow and doesn't increase resolution
            nZ=int(zMax*10*self.debug)#int((zMax*10.)/lam);
            nX=int(xMax*10*self.debug)
            self.nStepsPerMeter=10*self.debug
            
        self.fLen=int(nZ/6);#this is the size of the filter **in samples** at the top and bottom of the domain to eliminate artificial reflections
        self.fLenM=int(zMax/6)#this is the size of the filter in meters
        self.nZ=nZ
        self.nX=nX
        self.dz=zMax/(nZ);
        self.dx = xMax/(nX);

        pol=polarization;#1 for h, 0 for v
        z0=sourceDepth;
        self.nProfile=nProfile
        self.dataType=dataType
        self.sourceDepth=sourceDepth

        self.method=method
        
        self.polarization=pol
        z=np.linspace(-zMax-self.fLenM,(zMax+self.fLenM), ((2*self.fLen) + (2*nZ)));
        self.z=z
        self.fullZ=int(len(z))
        self.halfZ=int(len(z)/2)
        self.fullDepthM=zMax+self.fLenM

        self.zPlot=z[self.halfZ:self.fullZ]#useful for plotting

        self.x=np.linspace(0, xMax, nX);
        self.n2Mat=np.zeros((self.nX, self.fullZ), dtype='complex');#matrix to hold the n(x,z) profile

        #z space wavenumber
        self.kz=np.zeros(self.fullZ)
        self.kz[0:self.halfZ]=np.linspace(0,np.pi/self.dz, self.halfZ)
        self.kz[self.halfZ:self.fullZ]=np.linspace(-np.pi/self.dz,0, self.halfZ)

        #set the base index of refraction profile for this object. This is always the functional form, even if the user specifies data-driven, and is used to set values for calculations used.
        self.setBaseIndexOfRefractionProfile(site, nProfile)        
        #get index of refraction correction
        depthIndex=util.findNearest(z, sourceDepth)
        self.debugPrint("index : "+str(depthIndex))

        #various values of n
        nAtDepth=np.sqrt(self.n2Vec[depthIndex])
        nAtSurface=np.sqrt(self.n2Vec[util.findNearest(z, 2)])
        nAtBottom=1.78
        self.nAtDepth=nAtDepth
        self.nAtBottom=nAtBottom
        self.nMean=(nAtBottom+nAtSurface)/2.
        self.nAtSurface=nAtSurface
        self.z0=sourceDepth
        #reference wavelength. uses the depth of the TX to set n. found empirically to give the best agreement.
        lam =  util.c_light/(nAtDepth*abs(freq));
        lam0=util.c_light/abs(freq)

        self.debugPrint((depthIndex, z[depthIndex], nAtDepth, np.sqrt(self.n2Vec[depthIndex])))
        self.debugPrint(util.findNearest(-z, sourceDepth))
        self.lam=lam
        self.omega=2.*np.pi*freq
        #reference wave number
        k0 = 2.*np.pi/lam+0.j;

        self.k0=k0
        
        self.dipoleBeam(self.z)

        

    def dipoleBeam(self, z, A=1+0j):
        """Define a half-wave dipole source. 
        dimensions are defined by the wavelength. The dipole is vertical (axis along z) for now.
        """
        z0=self.z0
        z0Index=util.findNearest(z, z0)
        z0Index1=util.findNearest(z, -z0)

        halfLength=self.lam/2.
        self.debugPrint("halflength: "+str(halfLength))

        npoints=int(halfLength/self.dz)
        self.debugPrint("npoints: "+str(npoints))
        if npoints>len(z)/4:
            npoints=int(len(z)/4)

        zRange=np.linspace(-halfLength, halfLength, 2*npoints, dtype='complex')
        self.debugPrint(zRange.shape)
        zR1=np.linspace(0, 1, npoints, dtype='complex')
        zR2=np.linspace(1, 0, npoints, dtype='complex')
        zRange=np.append(zR1, zR2)

        n_x=np.pi*zRange
        e=[0., 0., 1.]
        beam=np.zeros(len(n_x), dtype='complex');
        f0=np.zeros(len(n_x), dtype='complex');
        for i in range(len(n_x)):
            n=[n_x[i], 0, 0]

            val=(np.cross(np.cross(n, e), n)[2])

            beam[i]=complex(val, val)

        f0=(A*(beam/(np.max(beam))))

        field=np.zeros(self.fullZ, dtype='complex')
        field[z0Index-npoints+1:z0Index+npoints+1]=f0


        self.initField=field;
        return field
        
    def setBaseIndexOfRefractionProfile(self,site, nProfile):
        """sets a base profile for calculations"""
        if site=="TD":
            self.n2TD(abs(self.z), 0)
        else:
            self.n2SP(abs(self.z), 0)

        
        
        
    def makeRandomDensityVector(self, z):
        """make a vector of random density fluctuations. This is currently used with the Taylor Dome n(z) profile.
        the density fluctuations are on the order of a few percent of the density."""
        dz=abs(z[1]-z[0])
        ranVec=util.lowpassFilter(dz, .05, (.6/(1.+(z*.01)))*(np.random.random_sample(len(z))-.5))
        return ranVec


    def n2TD(self, z, nProfile="functional"):
        """return the functional form of the index of refraction for Taylor Dome, with optional random density fluctuations set by ranFlag"""
        rhosurf=.46;
        rhomax=.928;
        const=-.02;
        n=np.full(self.fullZ, 1.0003, dtype='complex')
        n[self.halfZ:self.fullZ]=1.+ .86*(rhosurf +(rhomax-rhosurf)*(1.-np.exp(const*z[self.halfZ:self.fullZ])))
        rann=np.zeros(self.fullZ)

        if nProfile=="data":
            #add random, lowpass filtered noise on the level of a few percent of the density. make it inversely proportional to the depth, such that the variations disappear at large z         
            rann[self.halfZ:self.fullZ]=self.makeRandomDensityVector(z[self.halfZ:self.fullZ])
        self.nFunc=n
        self.n2Func=(n*n)
        self.n2Vec=(n*n)+rann
        for n in range(len(self.n2Mat)):
            self.n2Mat[n]=self.n2Vec
        return self.n2Vec



    
    def n2SP(self, z, nProfile="functional", dataType=""):
        """set the n(z) profile for the south pole, using the functional form. Optional to use the SPICE density profile for depths shallower than 100m"""
        if nProfile=="functional":
            #This is a parameterization based on the SPICE data
            A=1.78
            B=-.43
            C=-.0132

            func1=np.full(self.halfZ, 1.0003)

            func2=A+B*np.exp(C*z[self.halfZ:self.fullZ])

            func3=np.append(func1, func2);
            
            n2Vec=func3*func3
        
            self.n2Vec=n2Vec

            return n2Vec
        
        else:
            if(dataType=="phased"):
                return self.n2Phased()
            elif(dataType=="core2"):
                zz, n=np.loadtxt(self.path+'/share/spice2019_indOfRef_core2_5cm.dat', unpack=True)
            else:
                zz, n=np.loadtxt(self.path+'/share/spice2019_indOfRef_core1_5cm.dat', unpack=True)

                            
                
            nInterp=n #for 5cm data, which is default
            if self.nStepsPerMeter!=20:
                nInterp=util.interpolate(n, int(self.nStepsPerMeter/20)) #for anything coarser than 5cm data

            
            nFlip=np.flip(nInterp, axis=0);


            A=1.78
            B=-.43
            C=-.0132



            func2=A+B*np.exp(C*z[int(self.halfZ)+len(nInterp):])

            tmp0=np.full(int(len(z)/2),1.0003)
            tmp1=np.append(tmp0, nInterp)
            nOut=np.append(tmp1, func2)
            self.n2Vec=nOut*nOut
            return self.n2Vec

    def southpoleFit(self, z):
        """piecewise function that matches the spice data pretty OK, can be used elsewhere"""
        A = 1.78
        B = -0.43
        C = -0.0132
        if z >= 0 and z<30: return 1.361 +z/220
        elif z >=30 and z<60: return 1.4115 + z/365
        elif z >=60 and z < 107: return 1.45 + z/475
        elif z >= 107: return (A+B*(np.exp(C*z)))
        else: return 1.0003


    def n2Phased(self):
        """smooth interpolation between the two spice cores at either end of the simulation domain"""
        zz, n1=np.loadtxt(self.path+'/share/spice2019_indOfRef_core1_5cm.dat', unpack=True)
        zz, n2=np.loadtxt(self.path+'/share/spice2019_indOfRef_core2_5cm.dat', unpack=True)        
        nMat=np.zeros((self.nX, len(n1)))

        n2Mat=np.zeros((self.nX,self.fullZ))
        n2Vec=self.n2SP(self.z, 0)
        
        for step in range(self.nX):
            nVecStep=(n2/self.nX)*step + (n1/self.nX)*(self.nX-step)
            n2Vec[self.halfZ:(self.halfZ)+len(nVecStep)]=nVecStep**2
            n2Mat[step]=n2Vec
        self.n2Mat=n2Mat
        return n2Mat

    
    
        
    def filt(self, z):
        """simple window filter to keep reflections from happening off the artificial boundary at the very top and very bottom of the domain"""
        flen=self.fLen
        win=np.blackman(flen*2)
        zFilt=np.ones(len(z))
        zFilt[0:flen]=win[0:flen]
        zFilt[len(z)-flen:len(z)]=win[flen:len(win)]
        self.zFilt=zFilt
        return zFilt
        
       
    def doSolver(self):
        """This solves the parabolic equation.""" 
        
        
        #initialize the field
        self.field=self.initField;
        #initialize the index of refraction field
        if self.site=="TD":
            self.n2TD(abs(self.z), self.nProfile)
        else:
            self.n2SP(abs(self.z), self.nProfile, self.dataType)



        #first coefficient is common to methods "II", and "WA"
        self.alpha=np.exp(1.j*self.dx*self.k0*(np.sqrt(1.-(self.kz**2/self.k0**2))-1))
        #second coefficient is used for "WA"
        self.beta=np.exp(1.j*self.dx*self.k0*(np.sqrt((self.n2Vec)-1./2)))
        #a third coeff is needed if using the feit and fleck splitting, not used otherwise
        self.gamma=1.
        
        if(self.method=="FF"):
            #Feit and Fleck
            self.alpha=np.exp(-1.j*self.dx*self.kz**2/(2.*self.k0))
            self.beta=np.exp(1.j*self.dx*self.k0*(self.n2Vec-1.)/4.)
            self.gamma=np.exp(1.j*self.dx*self.k0*(self.n2Vec-1.)/4.)
        
        #for the in-ice case, assuming that n!=n(x,z)
        if(self.method=="II"):
            B=self.n2Vec-1
            Y=np.sqrt(1.+self.n2Vec/(self.nAtDepth**2))
            self.beta=np.exp(1.j*self.dx*self.k0*(np.sqrt(B+Y**2) - Y))


        #intializing outputs

        self.u=np.zeros((self.nX, self.nZ), dtype='complex');
        self.psi=np.zeros((self.nX, self.nZ), dtype='complex');
        self.uFull=np.zeros((self.nX, self.fullZ), dtype='complex');
        self.psiFull=np.zeros((self.nX, self.fullZ), dtype='complex');

        self.u[0][:]=self.field[self.halfZ:self.halfZ+self.nZ];
        #filter to remove non-physical boundary reflections from top and bottom of domain
        zFilt=self.filt(abs(self.z))
        
        #do the solver
        for n in range (1,self.nX):
            #if this is 1, n=n(x, z), not just n=n(z), define a n(x,z) matrix
            if self.nProfile=="data" and self.dataType=="phased":

                self.n2Vec=self.n2Mat[n]
                #first coefficient is common to methods "II", and "WA"
                self.alpha=np.exp(1.j*self.dx*self.k0*(np.sqrt(1.-(self.kz**2/self.k0**2))+1))

                #second coefficient is used for "WA"
                self.beta=np.exp(1.j*self.dx*self.k0*(np.sqrt((self.n2Vec)-1./2)))
                #a third coeff is needed if using the feit and fleck splitting, not used otherwise
                self.gamma=1.
                
                if(self.method=="FF"):
                    #Feit and Fleck
                    self.alpha=np.exp(-1.j*self.dx*self.kz**2/(2.*self.k0))
                    self.beta=np.exp(1.j*self.dx*self.k0*(self.n2Vec-1.)/4.)
                    self.gamma=np.exp(1.j*self.dx*self.k0*(self.n2Vec-1.)/4.)
                    

                if(self.method=="II"):
                    B=self.n2Vec-1
                    Y=np.sqrt(1.+self.n2Vec/(self.nAtDepth**2))
                    self.beta=np.exp(1.j*self.dx*self.k0*(np.sqrt(B+Y**2) - Y))

                
            #solve the fields            
            self.field=self.alpha*(util.doFFT(self.gamma*(self.field)));
            self.field=self.beta*(util.doIFFT((self.field)));

            self.field=self.field*zFilt#filter artificial reflections
            
            self.u[n][:]=self.field[self.halfZ:self.halfZ+self.nZ]

            self.psi[n][:]=self.u[n][:]/(np.sqrt(self.dx*n)*np.exp(-1.j*self.k0*self.dx*n))

            self.uFull[n][:]=self.field
            self.psiFull[n][:]=self.field/(np.sqrt(self.dx*n)*np.exp(-1.j*self.k0*self.dx*n))

        self.uFull=np.transpose(self.uFull)
        self.psiFull=np.transpose(self.psiFull)
        self.psi=np.transpose(self.psi)
        self.u=np.transpose(self.u)
        return self.psi


        
    def debugPrint(self, message):
        if self.debug>0:
            print(message)

