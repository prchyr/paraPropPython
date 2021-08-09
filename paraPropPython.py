# paraPropPython
# c. sbrocco, s. prohira

import util
import math as m
import numpy as np
from inspect import signature

class receiver:
    """
    Parameters
    ----------
    x : float
        x position (m)
    z : float
        z position (m)
    """
    def __init__(self, x, z):
        self.x = x
        self.z = z
    
    def setup(self, freq, dt):
        """
        further setup of receiver using simulation parameters
        
        Parameters
        ----------
        freq : float array
            frequencies (GHz)
        dt : float
            time step (ns)
        """
        self.freq = freq
        self.spectrum = np.zeros(len(freq), dtype='complex')
        self.time = np.arange(0, dt*len(freq), dt)
    
    def add_spectrum_component(self, f, A):
        """
        adds the contribution of a frequency to the received signal spectrum
        
        Parameters
        ----------
        f : float
            corresponding frequencie (GHz)
        A : complex float
            complex amplitude of received siganl (V/m???)
        """
        i = util.findNearest(self.freq, f)
        self.spectrum[i] = A
        
    def get_spectrum(self):
        """
        gets received signal spectrum
        
        Returns
        -------
        1-d comlplex float array
        """
        return self.spectrum[:int(len(self.freq)/2)]
    
    def get_signal(self):
        """
        gets received signal
        
        Returns
        -------
        1-d comlplex float array
        """
        return np.flip(util.doIFFT(self.spectrum))
    
    def get_frequency(self):
        """
        gets frequency array
        
        Returns
        -------
        1-d float array
        """
        return abs(self.freq)[:int(len(self.freq)/2)]
    
    def get_time(self):
        """
        gets time array
        
        Returns
        -------
        1-d float array
        """
        return self.time
     
        

class paraProp:
    """
    Parameters
    ----------
    iceLength : float
        length of the simulation (m)
    iceDepth : float
        depth of the ice simulated (m)
    dx : float
        grid spacing in the x direction (m)
    dz : float
        grid spacing in the z direction (m)
    airHeight : float
        amount of air to be simulated above ice (m). Initialized to 25 m
    filterDepth : float
        size of the filtered reason above and below simulated region (m). Initialized to 100 m
    refDepth : float
        reference depth for simulation (m). Initialized to 1 m below surface
    """
    def __init__(self, iceLength, iceDepth, dx, dz, airHeight=25, filterDepth=100, refDepth=1):       
        ### spatial parameters ### 
        # x #
        self.x = np.arange(0, iceLength+dx, dx)
        self.xNum = len(self.x)
        self.dx = dx
        
        # z #
        self.iceDepth = iceDepth
        self.airHeight = airHeight
        self.z = np.arange(-airHeight, iceDepth + dz, dz)
        self.zNum = len(self.z)
        self.dz = dz
        self.refDepth = refDepth            
        
        ### other simulation variables ###       
        # filter information #
        self.fNum0 = int(filterDepth / dz)
        
        self.fNum1, self.fNum2 = self.optimize_filt_size(self.zNum, self.fNum0)
        self.zFull = np.arange(-(airHeight + self.fNum1*dz), iceDepth + self.fNum2*dz + dz, dz)
        self.zNumFull = len(self.zFull)
        win = np.blackman(self.fNum1 + self.fNum2)
        filt = np.ones(self.zNumFull)
        filt[:self.fNum1] = win[:self.fNum1]
        filt[-self.fNum2:] = win[self.fNum1:]
        self.filt = filt
       
        # z wavenumber #
        self.kz = np.fft.fftfreq(self.zNumFull)*2*np.pi/self.dz
        
        # index of refraction array #
        self.n = np.ones(self.zNumFull)
        
        # source array #
        self.source = np.zeros(self.zNumFull, dtype='complex')
        
        # 2d field array #
        self.field = np.zeros((self.xNum, self.zNum), dtype='complex')
        
    def optimize_filt_size(self, zNum, fNum):
        zNumFull = 2*fNum + zNum
        p2 = 2**m.ceil(m.log(zNumFull, 2))
        p3 = 3**m.ceil(m.log(zNumFull, 3))
        p5 = 5**m.ceil(m.log(zNumFull, 5))
        p7 = 7**m.ceil(m.log(zNumFull, 7))
        p = min([p2, p3, p5, p7])
        fNum = p - zNum
        fNum1 = int(fNum/2)
        fNum2 = fNum - fNum1
        return fNum1, fNum2
        
    def get_x(self):
        """
        gets x grid of simulation
        
        Returns
        -------
        1-d float array
        """
        return self.x
    
    def get_z(self):
        """
        gets z grid of simulation
        
        Returns
        -------
        1-d float array
        """
        return self.z
      
    
    ### ice profile functions ###
    def set_n(self, nVal=None, nVec=None, nFunc=None, nAir=1.0003):
        """
        set the index of refraction profile of the simualtion
        
        future implementation plans:
            - complex index of refraction
        
        Parameters
        ----------
        nVal : float
            Postcondition: n(z>=0, x>=0) = nVal
        nVec : array
            1-d or 2-d array of float values
            Precondition: spacing between rows is dz, spacing between columns is dx
            Postcondition: n(z=0,x=0) = nVec[0,0], n(z=dz,x=dx) = nVec[1,1], ..., n(z>=len(nVec[:,0])*dz,x>=len(nVec[0,:])*dx) = nVec[-1,-1]
        nFunc : function
            Precondition: nFunc is a function of one or two variables, z and x, and returns a float value
            Postcondition: n(z>=0,x>=0) = nFunc(z,x)
        nAir : float
            index of refraction of air
            Postcondition: n(z<0) = nAir
        """    
        self.n = np.ones((self.zNumFull, self.xNum))
        
        if nVal != None:
            for i in range(self.zNumFull):
                if self.zFull[i] >= 0:
                    self.n[i,:] = nVal
                else:
                    self.n[i,:] = nAir
        
        elif nVec != None:             
            if len(nVec.shape) == 1:
                a = 0
                nNum = len(nVec)
                for i in range(self.zNumFull):
                    if self.zFull[i] >= 0:
                        ai = a if a < nzNum else -1
                        self.n[i,:] = nVec[ai]
                        a += 1
                    else:
                        self.n[i,:] = nAir
            elif len(nVec.shape) == 2: 
                a = 0
                b = 0
                nzNum = len(nVec[:,0])
                nxNum = len(nVec[0,:])
                for i in range(self.zNumFull):
                    for j in range(self.xNum):
                        if self.zFull[i] >= 0:
                            ai = a if a < nzNum else -1
                            bi = b if b < nxNum else -1
                            self.n[i,j] = nVec[ai,bi]
                            a += 1
                            b += 1
                        else:
                            self.n[i,j] = nAir
         
        elif nFunc != None:
            sig = signature(nFunc)
            numParams = len(sig.parameters)
            if numParams == 1:
                for i in range(self.zNumFull):
                    if self.zFull[i] >= 0:
                        z = self.zFull[i] if self.zFull[i] <= self.iceDepth else self.iceDepth
                        self.n[i,:] = nFunc(z)
                    else:
                        self.n[i,:] = nAir
            elif numParams == 2:
                for i in range(self.zNumFull):
                    for j in range(self.xNum):
                        if self.zFull[i] >= 0:
                            z = self.zFull[i] if self.zFull[i] <= self.iceDepth else self.iceDepth
                            x = self.x[j]
                            self.n[i,j] = nFunc(z,x)
                        else:
                            self.n[i,j] = nAir
                            
        ### set reference index of refraction ###
        self.n0 = self.at_depth(self.n[:,0], self.refDepth)
        
        self.n = np.transpose(self.n) 
        
    def get_n(self):
        """
        gets index of refraction profile of simulation
        
        Returns
        -------
        2-d float array
        """
        return np.transpose(self.n[:,self.fNum1:-self.fNum2])
   
    
    ### source functions ###
    def set_user_source_profile(self, method, z0=0, sVec=None, sFunc=None):
        """
        set the spatial source profile explicitly (no frequency / signal information)
        Precondition: index of refraction profile is already set
        
        Parameters
        ----------   
        method : string
            'vector' for vector defined profile
            'func' for function defined profile
        z0 : float
            Precondition: z0>=0
            reference starting point for sVec (m). Initialized to 0 m
        sVec : array
            if method=='vector', defines the source profile as an array
            Precondition: spacing between elements is dz
            Postcondition: E(z=z0) = sVec[0], E(z=z0+dz) = sVec[1], ... , E(z>=z0+len(sVec)*dz) = sVec[-1], TODO
        sFunc : function
            if method=='func', defines the source profile as a function
            Precondition: sFunc is a function of one variable, z, and returns a float value
            Postcondition: E(z>=0) = sFunc(z)
        """    
        self.source = np.zeros(self.zNumFull, dtype='complex')
        
        ### vector method ###
        if method == 'vector':
            sNum = len(sVec)
            j = 0
            for i in range(self.zNumFull):
                if self.zFull[i] >= z0:
                    if j < sNum:
                        self.source[i] = sVec[j]
                    else:
                        self.source[i] = 0
                    j += 1
                else:
                    self.source[i] = 0
        
        ### functional method ###
        if method == 'func':
            for i in range(self.zNumFull):
                if self.zFull[i] >= 0:
                    self.source[i] = sFunc(self.zFull[i])
                else:
                    self.source[i] = 0      
        
    def set_dipole_source_profile(self, centerFreq, depth, A=1+0.j):
        """
        set the source profile to be a half-wave dipole sized to center frequency
        Precondition: index of refraction profile is already set
        
        Parameters
        ----------  
        centerFreq : float
            center frequency of to model dipole around (GHz)
        depth : float
            Precondition: depth>=0
            depth of middle point of dipole (m)
        A : complex float
            complex amplitude of dipole. Initialized to 1 + 0j
        """
        ### frequency and wavelength in freespace ###
        self.source = np.zeros(self.zNumFull, dtype='complex')
        centerLmbda = util.c_light/centerFreq
        
        ### wavelength at reference depth ###
        centerLmbda0 = centerLmbda/self.n0
        
        ### create dipole ###
        z0 = depth
        z0Index = util.findNearest(self.zFull, z0)
        
        nPoints = int((centerLmbda0/2) / self.dz)
        ZR1 = np.linspace(0,1, nPoints, dtype='complex')
        ZR2 = np.linspace(1,0, nPoints, dtype='complex')
        zRange = np.append(ZR1, ZR2)
        
        n_x = np.pi*zRange
        e = [0., 0., 1.]
        beam = np.zeros(len(n_x), dtype='complex')
        f0 = np.zeros(len(n_x), dtype='complex')
        
        for i in range(len(n_x)):
            n=[n_x[i], 0, 0]
            val = np.cross(np.cross(n,e),n)[2]
            beam[i] = complex(val, val)
        f0 = A*(beam/(np.max(beam)))
        
        self.source[z0Index-nPoints+1:z0Index+nPoints+1]=f0
        
    def get_source_profile(self):
        """
        gets source profile of simulation
        
        Returns
        -------
        1-d comlplex float array
        """
        return self.source[self.fNum1:-self.fNum2]
   
    
    ### signal functions ###
    def set_cw_source_signal(self, freq):
        """
        set a continuous wave signal at a specified frequency
        
        Parameters
        ----------
        freq : float
            frequency of source (GHz) 
        """
        ### frequency ###
        self.freq = np.array([freq], dtype='complex')
        self.freqNum = len(self.freq)
        
        ### wavenumber at reference depth ###
        self.k0 = 2.*np.pi*self.freq*self.n0/util.c_light 
        
        ### coefficient ###
        self.A = np.array([1], dtype='complex')
        
    def set_td_source_signal(self, sigVec, dt):
        ### save input ###
        self.dt = dt
        self.sigVec = sigVec
        
        ### frequencies ###
        df = 1/(len(sigVec)*dt)
        self.freq = np.arange(0, 1/dt, df, dtype='complex')
        self.freqNum = len(self.freq)
        
        ### wavenumbers at reference depth ###
        self.k0 = 2.*np.pi*self.freq*self.n0/util.c_light 
        
        ### coefficient ###
        self.A = util.doFFT(np.flip(sigVec))
        
        # to ignore the DC component #
        self.A[0] = self.k0[0] = 0

        
    def get_spectrum(self):
        """
        gets transmitted signal spectrum
        
        Returns
        -------
        1-d comlplex float array
        """
        return self.A[:int(self.freqNum/2)]
    
    def get_frequency(self):
        """
        gets frequency array
        
        Returns
        -------
        1-d float array
        """
        return abs(self.freq)[:int(self.freqNum/2)]
    
    def get_signal(self):
        """
        gets transmitted signal
        
        Returns
        -------
        1-d comlplex float array
        """
        return self.sigVec
    
    def get_time(self):
        """
        gets time array
        
        Returns
        -------
        1-d float array
        """
        return np.arange(0, self.dt*len(self.sigVec), self.dt)
               
        
    ### field functions ###    
    def do_solver(self, rxList=np.array([])):
        """
        calculates field at points in the simulation
        Precondition: index of refraction and source profiles are set

        future implementation plans:
            - different method options
            - only store last range step option
            
        Parameters
        ----------
        rxList : array of Receiver objects
            optional for cw signal simulation
            required for non cw signal simulation
        """ 
        
        if (self.freqNum != 1):
            ### check for Receivers ###
            if (len(rxList) == 0):
                print("Warning: Running time-domain simulation with no receivers. Field will not be saved.")
            for rx in rxList:
                rx.setup(self.freq, self.dt)
                
        for j in np.arange(0, int(self.freqNum/2)+self.freqNum%2, 1, dtype='int'):
            if (self.freq[j] == 0): continue
            u = 2 * self.A[j] * self.source * self.filt * self.freq[j]
            self.field[0] = u[self.fNum1:-self.fNum2]
            
            alpha = np.exp(1.j * self.dx * self.k0[j] * (np.sqrt(1. - (self.kz**2 / self.k0[j]**2))- 1.))
            B = self.n**2-1
            Y = np.sqrt(1.+(self.n/self.n0)**2)
            beta = np.exp(1.j * self.dx * self.k0[j] * (np.sqrt(B+Y**2)-Y))
            
            for i in range(1, self.xNum):               
                u = alpha * (util.doFFT(u))
                u = beta[i] * (util.doIFFT(u))
                u = self.filt * u

                self.field[i] = u[self.fNum1:-self.fNum2]/(np.sqrt(self.dx*i) * np.exp(-1.j * self.k0[j] * self.dx * i))
                
            if (len(rxList) != 0):
                for rx in rxList:
                    rx.add_spectrum_component(self.freq[j], self.get_field(x0=rx.x, z0=rx.z))
                self.field.fill(0)

          
    def get_field(self, x0=None, z0=None):
        """
        gets field calculated by simulation
        
        future implementation plans:
            - interpolation option
            - specify complex, absolute, real, or imaginary field
            
        Parameters
        ----------
        x0 : float
            position of interest in x-dimension (m). optional
        z0 : float
            position of interest in z-dimension (m). optional
        
        Returns
        -------
        if both x0 and z0 are supplied
            complex float
        if only one of x0 or z0 is supplied
            1-d complex float array
        if neither x0 or z0 are supplied
            2-d complex float array
        """
        if (x0!=None and z0!=None):
            return self.field[util.findNearest(self.x, x0),util.findNearest(self.z,z0)]
        if (x0!=None and z0==None): 
            return self.field[util.findNearest(self.x, x0),:]     
        if (x0==None and z0!=None):
            return self.field[:,util.findNearest(self.z,z0)]
        return self.field
                                   

    ### misc. functions ###
    def at_depth(self, vec, depth):
        """
        find value of vector at specified depth.
        future implementation plans:
            - interpolation option
            - 2D array seraching. paraProp.at_depth() -> paraProp.at()
        
        Parameters
        ----------
        vec : array
            vector of values
            Precondition: len(vec) = len(z)
        depth : float
            depth of interest (m)
        
        Returns
        -------
        base type of vec
        """  
        ### error if depth is out of bounds of simulation ###
        if (depth > self.iceDepth or depth < -self.airHeight):
                print("Error: Looking at z-position of out bounds")
                return np.NaN
         
        # find closest index #
        dIndex = round((depth + self.fNum1*self.dz + self.airHeight) / self.dz)
        
        return vec[dIndex]
    
    