# paraPropPython
# s. prohira, c. sbrocco

### TODO: Time Dependent Signals and Receiver class ###

import util
import numpy as np

class receiver:
    """
    TODO: fully implement
    """
    def __init__(self, x, z):
        self.x = x
        self.z = z

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
        self.zFull = np.arange(-(airHeight + filterDepth), iceDepth + filterDepth + dz, dz)
        self.zNum = len(self.z)
        self.zNumFull = len(self.zFull)
        self.dz = dz
        self.refDepth = refDepth            
        
        ### other simulation variables ###       
        # filter information #
        self.filterDepth = filterDepth
        self.fNum = int(filterDepth / dz)
        win = np.blackman(2*self.fNum)
        filt = np.ones(self.zNumFull)
        filt[:self.fNum] = win[:self.fNum]
        filt[-self.fNum:] = win[self.fNum:]
        self.filt = filt
       
        # z wavenumber #
        self.kz = np.zeros(self.zNumFull)
        self.kz[:int(self.zNumFull/2)] = np.linspace(0, np.pi/self.dz, int(self.zNumFull/2))
        self.kz[-int(self.zNumFull/2):] = np.linspace(-np.pi/self.dz, 0, int(self.zNumFull/2))
        
        # index of refraction array #
        self.n = np.ones(self.zNumFull)
        
        # source array #
        self.source = np.zeros(self.zNumFull, dtype='complex')
        
        # 2d field array #
        self.field = np.zeros((self.xNum, self.zNum), dtype='complex')
        
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
    def set_n(self, method, nVec=None, nFunc=None, nAir=1.0003):
        """
        set the index of refraction profile of the simualtion
        
        future implementation plans:
            - 2-d profiles
            - complex index of refraction
        
        Parameters
        ----------
        method : string
            'vector' for vector defined profile
            'func' for function defined profile
        nVec : array
            if method=='vector', defines the index of refraction profile of ice as an array
            Precondition: spacing between elements is dz
            Postcondition: n(z=0) = nVec[0], n(z=dz) = nVec[1], ... , n(z>=len(nVec)*dz) = nVec[-1]
        nFunc : function
            if method=='func', defines the index of refraction profile of ice as a function
            Precondition: nFunc is a function of one variable, z, and returns a float value
            Postcondition: n(z>=0) = nFunc(z)
        nAir : float
            index of refraction of air
            Postcondition: n(z<0) = nAir
        """    
        self.n = np.ones(self.zNumFull)
        
        ### vector method ###
        if method == 'vector': 
            nNum = len(nVec)
            j = 0
            for i in range(self.zNumFull):
                if self.zFull[i] >= 0:
                    if j < nNum:
                        self.n[i] = nVec[j]
                    else:
                        self.n[i] = nVec[-1]
                    j += 1
                else:
                    self.n[i] = nAir
         
        ### functional method ###
        if method == 'func':
            for i in range(self.zNumFull):
                if self.zFull[i] >= 0:
                    if self.zFull[i] <= self.iceDepth:
                        self.n[i] = nFunc(self.zFull[i])
                    else:
                        self.n[i] = nFunc(self.iceDepth)
                else:
                    self.n[i] = nAir
                    
        ### set reference index of refraction ###
        self.n0 = self.at_depth(self.n, self.refDepth)
        
    def get_n(self):
        """
        gets index of refraction profile of simulation
        
        Returns
        -------
        1-d float array
        """
        return self.n[self.fNum:-self.fNum]
    
    
    
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
            Postcondition: E(z>=0) = nsFunc(z)
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
        return self.source[self.fNum:-self.fNum]

    
    
    ### signal functions ###
    def set_cw_source_signal(self, freq):
        """
        set a continuous wave signal at a specified frequency
        
        Parameters
        ----------
        freq : float
            frequency of source (GHz) 
        """
        ### frequency and wavelength in free space ###
        self.freq = np.array([freq], dtype='complex')
        self.freqNum = len(self.freq)
        self.lmbda = util.c_light/self.freq
        
        ### wavelength and wavenumber at reference depth ###
        self.lmbda0 = self.lmbda/self.n0
        self.k0 = 2.*np.pi/self.lmbda0 
        
        ### coefficient ###
        self.A = np.array([1], dtype='complex')
        
        
        
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
                
        for j in range(self.freqNum):
            u = self.A[j] * self.source * self.filt
            self.field[0,:] = u[self.fNum:-self.fNum]

            ### method II ###
            alpha = np.exp(1.j * self.dx * self.k0[j] * (np.sqrt(1. - (self.kz**2 / self.k0[j]**2))- 1.))
            B = (self.n)**2-1
            Y = np.sqrt(1.+(self.n/self.n0)**2)
            beta = np.exp(1.j * self.dx * self.k0[j] * (np.sqrt(B+Y**2)-Y))

            for i in range(1, self.xNum):           
                u = alpha * (util.doFFT(u))
                u = beta * (util.doIFFT(u))

                u = self.filt * u

                self.field[i,:] = u[self.fNum:-self.fNum]/(np.sqrt(self.dx*i) * np.exp(-1.j * self.k0[j] * self.dx * i))

          
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
        if (self.freqNum != 1):
            print("TODO: Non-CW signal warning")
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
        dIndex = round((depth + self.filterDepth + self.airHeight) / self.dz)
        
        return vec[dIndex]
                     
     

    
        