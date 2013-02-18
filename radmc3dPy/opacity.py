# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class radmc3dDustOpac():
    """
    Class for dust opacities

    ATTRIBUTES:
    -----------
        wav   - wavelength grid
        freq  - frequency grid
        nwav  - number of wavelengths
        kabs  - absorption coefficient per unit mass
        ksca  - scattering coefficient per unit mass
        phase_g - phase function
    
    FUNCTIONS:
    ----------
        readopac()
        read_masteropac() 
        write_masteropac() 
    """
# --------------------------------------------------------------------------------------------------
    def __init__(self):

        self.wav     = -1.
        self.freq    = -1.
        self.nwav    = -1.
        self.nfreq   = -1.
        self.kabs    = -1.
        self.ksca    = -1.
        self.phase_g = -1.
        self.ext     = ''
        self.idust   = -1
        self.therm   = True
         
# --------------------------------------------------------------------------------------------------
    def  readopac(self, ext='', used=False):
        """
        Function to read the dust opacity files

        INPUT:
        ------
            ext  : file name extension (file names should look like 'dustkappa_ext.inp')
            used : if set to True the used dust opacity file ('dustkappa_ext.inp.used') file is read, 
                    which is interpolated to the frequency grid of the RADMC3D run
        """

        if not used:
            try:
                rfile = open('dustkappa_'+ext+'.inp', 'r')
            except:
                print 'ERROR'
                print ' No dustkappa_'+ext+'.inp file was found'
                return -1
        else:
            try:
                rfile = open('dustkappa_'+ext+'.inp.used', 'r')
            except:
                print 'ERROR'
                print ' No kappa_'+ext+'.inp.used file was found'
                return -1

        # Read two comment lines if the .used files are to be read
        if used:
            dum = rfile.readline()
            dum = rfile.readline()

        # Read the file format
        iformat = int(rfile.readline())
        if (iformat<1)|(iformat>3):
            print 'ERROR'
            print 'Unknown file format in the dust opacity file'
            rfile.close()
            return -1


        # Read the number of wavelengths in the file
        dum = rfile.readline()
        self.nwav = int(dum)
        self.nfreq = int(dum)
        
        # If only the absorption coefficients are specified
        if iformat==1:
            self.wav = zeros(self.nwav, dtype=float64)
            self.kabs = zeros(self.nwav, dtype=float64)
            for ilam in range(self.nwav):
                dum = rfile.readline().split()
                self.wav[ilam] = dum[0] 
                self.kabs[ilam] = dum[1] 

        # If the absorption and scattering coefficients are specified
        elif iformat==2:
            self.wav = zeros(self.nwav, dtype=float64)
            self.kabs = zeros(self.nwav, dtype=float64)
            self.ksca = zeros(self.nwav, dtype=float64)
            for ilam in range(self.nwav):
                dum = rfile.readline().split()
                self.wav[ilam] = dum[0] 
                self.kabs[ilam] = dum[1] 
                self.ksca[ilam] = dum[2] 
        
        # If the absorption and scattering coefficients and also the scattering phase function are specified
        elif iformat==3:
            self.wav = zeros(self.nwav, dtype=float64)
            self.kabs = zeros(self.nwav, dtype=float64)
            self.ksca = zeros(self.nwav, dtype=float64)
            self.phase_g = zeros(self.nwav, dtype=float64)
            for ilam in range(self.nwav):
                dum = rfile.readline().split()
                self.wav[ilam] = dum[0] 
                self.kabs[ilam] = dum[1] 
                self.ksca[ilam] = dum[2] 
                self.phase_g[ilam] = dum[3] 
   
        rfile.close()
# --------------------------------------------------------------------------------------------------
    def  read_masteropac():
        """
        Function to read the master opacity file 'dustopac.inp' 
        it reads the dustkappa filename extensions (dustkappa_ext.inp) corresponding to dust species indices

        OUTPUT:
        -------
        Returns a dictionary with the following keys:
            'ext'   - list of dustkappa file name extensions
            'therm' - a list of integers specifying whether the dust grain is thermal or quantum heated 
                      (0 - thermal, 1 - quantum heated)
        """
        
        try: 
            rfile = open('dustopac.inp', 'r')
        except:
            print 'Error'
            print ' No dustopac.inp file was found'
            return -1

       
        # file format
        dum = rfile.readline()
        # nr of dust species
        ndust = int(rfile.readline().split()[0])
        # Comment line
        dum = rfile.readline()

        ext = []
        therm= []
        for idust in range(ndust):
            dum = rfile.readline()
            # Check if the dust grain is thermal or quantum heated
            dum = int(rfile.readline().split()[0])
            if dum==0:
                therm.append(True)
            else:
                therm.append(False)
            # Dustkappa filename extension
            dum = rfile.readline().split()[0]
            ext.append(dum)
            #Comment line
            dum = rfile.readline()
        rfile.close()

        return {'ext':ext, 'therm':therm}


# --------------------------------------------------------------------------------------------------
    def  write_masteropac(ext=None, therm=None):
        """
        Function to write the master opacity file 'dustopac.inp' 

        INPUT:
        ------
            ext : list of dustkappa file name extensions
            therm : list of integers specifying whether the dust grain is thermal or quantum heated
                    (0-thermal, 1-quantum)
        """

        print 'Writing dustopac.inp'
       
        if not ext:
            print 'ERROR'
            print 'No file name extension is specified. Without it dustopac.inp cannot be written'
            return -1
        else:
            if (type(ext).__name__=='str'):  ext = [ext]

        if therm:
            if (type(therm).__name__=='int'): therm = [therm]
            if (len(ext)!=len(therm)):
                print 'ERROR'
                print ' The number of dust species in ext and in therm are different'
                return -1

        wfile = open('dustopac.inp', 'w')

        # File format
        wfile.write('%-15s %s\n'%('2', 'Format number of this file'))
        # Number of dust species
        wfile.write('%-15s %s\n'%(str(len(ext)), 'Nr of dust species'))
        # Separator
        wfile.write('%s\n'%'============================================================================')


        for idust in range(len(ext)):
            # Dust opacity will be read from a file
            wfile.write('%-15s %s\n'%('1', 'Way in which this dust species is read'))
            # Check if the dust grain is thermal or quantum heated
            if therm:
                if therm[idust]:
                    wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))
                if therm[idust]:
                    wfile.write('%-15s %s\n'%('1', '0=Thermal grain, 1=Quantum heated'))
            else:
                wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))

            # Dustkappa filename extension
            wfile.write('%s %s %s\n'%(ext[idust], '    ', 'Extension of name of dustkappa_***.inp file'))
            # Separator
            wfile.write('%s\n'%'----------------------------------------------------------------------------')
            
        wfile.close()

# ------------------------------------------------------------------------------------------------------------------------
def run_makedustopac(freq=None, gmin=None, gmax=None, ngs=None, lnk_fname=None, gdens=None):
    """
    Interface function to the F77 code makedust to calculate mass absorption
    coefficients from the optical constants using Mie-theory

    INPUT:
    ------
        freq       - numpy.ndarray containing the frequency grid on which the opacities should be calculated
        gmin       - minimum grain size
        gmax       - maximum grain size
        ngs        - number of grain sizes
        gdens      - density of the dust grain in g/cm^3
        lnk_faname - name of the file in which the optical constants are stored

    OUTPUT:
    -------
        result         - numpy.ndarray[nfreq,ngs] containing the resulting opacities

    FILE OUTPUT:
    ------------
        dustopac_i.inp - Contains the dust opacities in radmc3d format
        dustopac.inp   - Master dust opacity file

    """

#
# Calculate the grain sizes
#
    if ngs>1:
        gsize = gmin * (gmax/gmin)**(arange(ngs, dtype=float64)/(float(ngs)-1.))
    else:
        gsize = [gmin]

#
# Write the frequency.inp file
#
    wfile = open('frequency.inp', 'w')
    wfile.write("%d\n"%freq.shape[0])
    wfile.write("  \n")
    for i in range(freq.shape[0]):
        wfile.write("%.10e\n"%freq[i])
    wfile.close()

#
# Write the dust.inp file (makedust main control file)
#
    wfile = open('dust.inp', 'w')
    for igs in range(ngs):
        wfile.write("%s\n"%lnk_fname)
        wfile.write("%s\n"%"MIE")
        wfile.write("%d %f %f %f %d %f %f %f\n"%(1,0.0,log10(gsize[igs]), log10(gsize[igs]),1.,-3.5,gdens,-2.0))
    wfile.close()

#
# Run the Mie-code
#
    dum = Popen('makedust', shell=True).wait()


#--------------------------------------------------------------------------------------------------------------------
def makeopac(ppar=None, wav=None):
    """
    Function to create dust opacities for RADMC3D using either simple MIE calculation or 
    mixing dust species within a grain using effective medium theory (EMT) with the APMR mixing rule (Min et al. 2008)
    
    INPUT:
    ------
        ppar  - dictionary containing all parameter of the simulation
    
    OPTIONS:
    --------
        wav  - numpy.ndarray containing the wavelength grid on which the mass absorption coefficients should be calculated
    """

#
# Create the wavelength grid if it is not specified
#
    if wav==None:
        grid = radmc3dGrid()
        grid.make_wav_grid(ppar=ppar)
        wav = grid.wav

#
# Do we need to mix the opacities?
#
    if type(ppar['lnk_fname']).__name__=='str':
        ppar['lnk_fname'] = [ppar['lnk_fname']]

    if len(ppar['lnk_fname'])>1:
        for idust in range(len(ppar['lnk_fname'])):
            run_makedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                    lnk_fname=ppar['lnk_fname'][i], gdens=ppar['gdens'])
            for igs in range(ppar['ngs']):
                dum = Popen('cp dustkappa_'+str(igs+1)+'.inp dustkappa_idust_'+str(idust+1)+'_igsize_'+str(igs+1)+'.inp', shell=True).wait()
    else:
        run_makedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                lnk_fname=ppar['lnk_fname'][0], gdens=ppar['gdens'])
                dum = Popen('cp dustkappa_1.inp dustkappa_idust_1_igsize_1.inp', shell=True).wait()

#    #
#    # Loop over all grain sizes
#    #
#            gsize =  ppar['gsmin'] * (ppar['gsmax']/ppar['gsmin'])**(arange(ppar['ngs'], dtype=float64) / (float(ppar['ngs'])-1.))
#            for igs in range(ppar['ngs']):
#
#    #
#    # Loop over the crystallographic axes if we have crystals
#    #
#                lnk = []
#                ncax = len(ppar['lnk_fname'])
#                for iax in range(ncax):
#    #
#    # Calculate the effective refractie index for each dust species
#    #
#                    dum = get_nk_apmr(wav=wav, lnk_fname=ppar['lnk_fname'][iax], mfrac=ppar['mfrac'], gdens=ppar['gdens'], \
#                            wavenum=ppar['wavenum'], gsize=gsize[igs], cpol_type=ppar['cpol_type'])  
#                    lnk.append(dum)
#            
#                neff = lnk[0][1]
#                
#                keff = lnk[0][2]
#                for iax in range(1, ncax):
#                    neff = neff + lnk[iax][1]
#                    keff = keff + lnk[iax][2]
#                
#                neff = neff / float(ncax)
#                keff = keff / float(ncax)
#
#    # 
#    # Write out the refractive indices to a file and call the mie code
#    #
#                write_effnk(wav = lnk[0][0], n = neff, k = keff, fname='optconst_apmr.lnk')
#                run_makedustopac(freq=cc/wav*1e4, gmin=gsize[igs], ngs=1, lnk_fname='optconst_apmr.lnk', gdens=lnk[0][3])
#                dum = Popen('cp dustkappa_1.inp dustkappa_apmr_'+str(igs)+'.inp', shell=True).wait()
#
#    #
#    # Write the master opacity file
#    #
#
#            wfile = open('dustopac.inp', 'w')
#            wfile.write("%4d %s\n"%(2,'             Format number of this file'))
#            wfile.write("%4d %s\n"%(ppar['ngs'],'             Nr of dust species'))
#            wfile.write("%s \n"%('==================================================='))
#            for igs in range(ppar['ngs']):
#                wfile.write("%4d %s\n"%(1, '             Input style for species (1=dust kappa file)'))
#                wfile.write("%4d %s\n"%(0, '             0=Thermal grain, 1=Quantum heated'))
#                wfile.write("%4d %s\n"%(igs,'             Extension of file dustkappa_*.inp'))
#                wfile.write("%s\n"%'---------------------------------------------------')
#            wfile.close()
#
#

# --------------------------------------------------------------------------------------------------
def mixopac(self, ppar=None, mixnames=[], mixspecs=[], mixabun=[], readfile=False, writefile=True):
    """
    Function to mix opacities

    INPUT:
    ------
        ppar     - A dictionary containing all parameters of the actual model setup
                    If any keyword is set besides ppar, the value of the separate keyword
                    will be taken instead of that in ppar.

    OPTIONS:
    --------
        mixnames  - Names of the files into which the mixed dust opacities will be written (not needed if writefile=False)
        mixspecs  - Names of the files from which the dust opacities are read (not needed if readfile=False)
        mixabun   - Abundances of different dust species
        readfile  - If True the opacities to be mixed will be read from files given in mixspecs.
                    If False the opacity values present in the base class' atrributes will be used.
        writefile - If False the mixed opacities will not be written out to files given in mixnames.  
        
    """

    if writefile:
        if len(mixnames)==0:
            if ppar!=None:
                mixnames = ppar['mixnames']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixnames are set in mixopac '
                return
    if readfile:
        if len(mixspecs)==0:
            if ppar!=None:
                mixspecs = ppar['mixspecs']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixspecs are set in mixopac '
                return
    else:
        
    if len(mixabun)==0:
        if ppar!=None:
            mixabun = ppar['mixabun']
        else:
            print 'ERROR'
            print ' Neither ppar nor mixabun are set in mixopac '
            return


# --------------------------------------------------------------------------------------------------
def  read_masteropac():
    """
    Function to read the master opacity file 'dustopac.inp' 
    it reads the dustkappa filename extensions (dustkappa_ext.inp) corresponding to dust species indices

    (This function is an interface function that calls radmc3dDustOpac.read_masteropac())
    """

    res = radmc3dDustOpac.read_masteropac()
    return res
    
# --------------------------------------------------------------------------------------------------
def write_masteropac(ext=None, therm=None):
    """
    Function to write the opacity file 'dustopac.inp'
    (This function is an interface function that calls radmc3dDustOpac.write_masteropac())
    INPUT:
    ------
        ext : list of dustkappa file name extensions
        therm : list of integers specifying whether the dust grain is thermal or quantum heated
                (0-thermal, 1-quantum)

    """
    
    radmc3dDustOpac.write_masteropac(ext=ext, therm=therm)

