# paraPropPython
# c. sbrocco, s. prohira

import util
import math as m
import numpy as np
from inspect import signature
from permittivity import eps2m, m2eps
from receiver import receiver

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
        # Alex: Change to an array of complex numbers -> to account for attenuation
        self.n = np.ones((self.zNumFull, self.xNum), dtype='complex')
        self.epsilon_r = np.ones((self.zNumFull, self.xNum), dtype='complex')
        #TODO: Should I define it (zNum, xNum) or (xNum, zNum)? Because it's defined (zNum, xNum) in set_n but then is transformed by np.transpose to (xNum, zNum)

        # source array #
        self.source = np.zeros(self.zNumFull, dtype='complex')
        
        # 2d field array #
        self.field = np.zeros((self.xNum, self.zNum), dtype='complex')
        # 2d field array for reflected signals
        self.field_plus = np.zeros((self.xNum, self.zNum), dtype='complex')
        self.field_minus = np.zeros((self.xNum, self.zNum), dtype='complex')

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
        self.n = np.ones((self.zNumFull, self.xNum), dtype='complex')
        
        if nVal != None:
            for i in range(self.zNumFull):
                if self.zFull[i] >= 0:
                    self.n[i,:] = nVal
                else:
                    self.n[i,:] = nAir
        
        elif nVal == None and nFunc == None:
            if len(nVec.shape) == 1:
                a = 0
                nzNum = len(nVec) #TODO: was originally nNum -> changed to nzNum
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
        self.epsilon_r = m2eps(self.n)

        self.n = np.transpose(self.n)
        self.epsilon_r = np.transpose(self.epsilon_r)

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
    def set_cw_source_signal(self, freq): #TODO: Consider changing the amplitude
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
        self.kp0 = 2.*np.pi*self.freq*self.n0/util.c_light
        self.k0 = 2. * np.pi * self.freq / util.c_light
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
        self.kp0 = 2.*np.pi*self.freq*self.n0/util.c_light
        self.k0 = 2. * np.pi * self.freq / util.c_light

        ### coefficient ###
        self.A = util.doFFT(np.flip(sigVec))
        
        # to ignore the DC component #
        self.A[0] = self.kp0[0] = 0

        
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
    def do_solver(self, rxList=np.array([]), freqMin=None, freqMax=None, solver_mode = 'one-way', refl_threshold=1e-10):
        """
        calculate field across the entire geometry (fd mode) or at receivers (td mode)
        field can be estimate in the forwards or backwards direction or in both directions
        -> modified from do_solver()
        -> calculates forwards field
        -> if dn/dx > 0 -> save position of reflector
        -> use as a source
        -> calculate an ensemble of u_minus

        Precondition: index of refraction and source profiles are set

        future implementation plans:
            - different method options
            - only store last range step option

        Parameters
        ----------
        Optional:
        -rxList : array of Receiver objects
            optional for cw signal simulation
            required for non cw signal (td) simulation
        -freqMin : float (must be less than nyquist frequnecy)
            defines minimum cutoff frequnecy for TD evalaution
        -freqMax : float (must be less than nyquist frequnecy)
            defines maximum cutoff frequuncy for TD evaluation
        -solver_mode : string
            defines the simulation mode
            must be one of three options:
                one-way : only evaluates in the forwards direction (+)
                two-way : evaluates in forwards (+) and backwards direction (-)
                minus : only evaluates in the backwards (-) direction
        -refl_threshold : float
            sets minimum reflection power to be simulated (anything less will be neglected)

        Output:
        FD mode: self.field has solution of E field across the simualtion geometry for the inputs : n, f and z_tx
        TD mode: rxList contains the signal and spectra for the array of receivers
        """
        if solver_mode != 'one-way' and solver_mode != 'two-way':
            print('warning! solver mode must be given as: one-way or two-way')
            print('will default to one way')
            solver_mode = 'one-way'
        if (self.freqNum != 1):
            ### check for Receivers ###
            if (len(rxList) == 0):
                print("Warning: Running time-domain simulation with no receivers. Field will not be saved.")
            for rx in rxList:
                rx.setup(self.freq, self.dt)

        #Check if solving for TD signal or in FD
        if freqMin == None and freqMax == None:
            freq_ints = np.arange(0, int(self.freqNum / 2) + self.freqNum % 2, 1, dtype='int')
        elif freqMin == None and freqMax != None:
            ii_min = util.findNearest(self.freq, freqMin)
            freq_ints = np.arange(ii_min, int(self.freqNum / 2) + self.freqNum % 2, 1, dtype='int')
        elif freqMin != None and freqMax == None:
            ii_max = util.findNearest(self.freq, freqMax)
            freq_ints = np.arange(0, ii_max, 1, dtype='int')
        else:
            ii_min = util.findNearest(self.freq, freqMin)
            ii_max = util.findNearest(self.freq, freqMax)
            freq_ints = np.arange(ii_min, ii_max, 1, dtype='int')

        for j in freq_ints:
            if (self.freq[j] == 0): continue
            u_plus = 2 * self.A[j] * self.source * self.filt * self.freq[j]
            self.field[0] = u_plus[self.fNum1:-self.fNum2]

            alpha_plus = np.exp(1.j * self.dx * self.kp0[j] * (np.sqrt(1. - (self.kz ** 2 / self.kp0[j] ** 2)) - 1.))
            B_plus = self.n ** 2 - 1
            Y_plus = np.sqrt(1. + (self.n / self.n0) ** 2)
            beta_plus = np.exp(1.j * self.dx * self.kp0[j] * (np.sqrt(B_plus + Y_plus ** 2) - Y_plus))

            if solver_mode == 'two-way':
                refl_source_list = []
                x_refl = []
                ix_refl = []
                nRefl = 0

            for i in range(1, self.xNum):
                u_plus_i = u_plus # Record starting reduced field -> in case we need to calculate
                dn = self.x[i] - self.x[i - 1]
                if dn.any() > 0:
                    u_plus *= util.transmission_coefficient(self.n[i], self.n[i-1])
                u_plus = alpha_plus * (util.doFFT(u_plus))
                u_plus = beta_plus[i] * (util.doIFFT(u_plus))
                u_plus = self.filt * u_plus

                delta_x_plus = self.dx * i
                self.field_plus[i] = u_plus[self.fNum1:-self.fNum2] / (
                        np.sqrt(delta_x_plus) * np.exp(-1.j * self.k0[j] * delta_x_plus))

                if solver_mode == 'two-way': #set to reflection modes
                    if dn.any() > 0: #check if the ref index changes in x direction
                        refl_source = util.reflection_coefficient(self.n[i], self.n[i-1]) * u_plus_i
                        if (refl_source**2).any() > refl_threshold: #check if reflected power is above threshold
                            x_refl.append(self.x[i])
                            refl_source_list.append(refl_source)
                            ix_refl.append(i)
                            nRefl = len(refl_source_list)
            if solver_mode == 'two-way' or solver_mode == 'minus':  # set to reflection modes
                if nRefl > 0:
                    u_minus_3arr = np.zeros((self.zNumFull, nRefl), dtype='complex')
                    field_minus_3arr = np.zeros((self.xNum, self.zNum, nRefl), dtype='complex')
                    for l in range(nRefl):
                        ix = ix_refl[l]
                        u_minus_3arr[:, l] = refl_source_list[l]
                        field_minus_3arr[ix, :, l] = u_minus_3arr[self.fNum1:-self.fNum2, l]
                        alpha_minus = np.exp(
                            1.j * self.dx * self.kp0[j] * (np.sqrt(1. - (self.kz ** 2 / self.kp0[j] ** 2)) - 1.))
                        B_minus = self.n ** 2 - 1
                        Y_minus = np.sqrt(1. + (self.n / self.n0) ** 2)
                        beta_minus = np.exp(1.j * self.dx * self.kp0[j] * (np.sqrt(B_minus + Y_minus ** 2) - Y_minus))
                        ix_last = ix_refl[l]
                        for k in np.flip(np.arange(0, ix_last, 1, dtype='int')):
                            x_minus = self.x[k]
                            dx_minus = abs(x_minus - self.x[ix_last])

                            u_minus_3arr[:, l] = alpha_minus * (util.doFFT(u_minus_3arr[:, l]))  # ????
                            u_minus_3arr[:, l] = beta_minus[k] * (util.doIFFT(u_minus_3arr[:, l]))
                            u_minus_3arr[:, l] = self.filt * u_minus_3arr[:, l]
                            field_minus_3arr[k, :, l] = np.transpose(
                                (u_minus_3arr[self.fNum1:-self.fNum2, l] / np.sqrt(dx_minus)) * np.exp(
                                    1j * dx_minus * self.k0[j]))
                        self.field_minus[:,:] += field_minus_3arr[:,:,l]

            if solver_mode == 'one-way':
                self.field[:,:] = self.field_plus[:,:]
                if (len(rxList) != 0):
                    for rx in rxList:
                        rx.add_spectrum_component(self.freq[j], self.get_field(x0=rx.x, z0=rx.z))
                    self.field.fill(0)
            elif solver_mode == 'two-way':
                if (len(rxList) != 0):
                    for rx in rxList:
                        rx.add_spectrum_component_plus(self.freq[j], self.get_field_plus(x0=rx.x, z0=rx.z))
                        rx.add_spectrum_component_minus(self.freq[j], self.get_field_minus(x0=rx.x, z0=rx.z))
                        rx.spectrum = rx.spectrum_plus + rx.spectrum_minus
                    self.field.fill(0)
                    self.field_plus.fill(0)
                    self.field_minus.fill(0)
                else:
                    self.field[:,:] += self.field_plus[:,:]
                    self.field[:,:] += self.field_minus[:,:]

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

    def get_field_plus(self, x0=None, z0=None):
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
        if (x0 != None and z0 != None):
            return self.field_plus[util.findNearest(self.x, x0), util.findNearest(self.z, z0)]
        if (x0 != None and z0 == None):
            return self.field_plus[util.findNearest(self.x, x0), :]
        if (x0 == None and z0 != None):
            return self.field_plus[:, util.findNearest(self.z, z0)]
        return self.field_plus

    def get_field_minus(self, x0=None, z0=None):
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
        if (x0 != None and z0 != None):
            return self.field_minus[util.findNearest(self.x, x0), util.findNearest(self.z, z0)]
        if (x0 != None and z0 == None):
            return self.field_minus[util.findNearest(self.x, x0), :]
        if (x0 == None and z0 != None):
            return self.field_minus[:, util.findNearest(self.z, z0)]
        return self.field_minus

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
    
    