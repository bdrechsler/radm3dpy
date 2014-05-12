from radmc3d.analyze import *
from numpy import *
from subprocess import Popen

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def write_effnk(wav=None, n=None, k=None, fname=None):
    """
    Function to write the effective refractive indices to a file

    INPUT:
    ------
        wav   : numpy.ndarray containing the wavelength on which the optical constants are defined
        n     : numpy.ndarray containing the real part of the refractive index
        k     : numpy.ndarray containing the imaginary part of the refractive index

    OPTIONS:
    --------
        fname : str containing the file name into which the refractive indices should be written
                if not specified fname='optconst_emt.lnk' is used
    """

    if not fname:
        fname = 'optconst_emt.lnk'

    wfile = open(fname, 'w')
    for i in range(wav.shape[0]):
        wfile.write("%.10e %.10e %.10e\n"%(wav[i], n[i], k[i]))
    wfile.close()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_nk_apmr(wav=None, lnk_fname=[], wavenum=[], mfrac=[], gdens=[], vfrac=[], gsize=1.0, cpol_type='DHS'):
    """
    Function to calculate the effective refractive indices of a fractal aggregate using 
    effective medium theory with the APMR mixing rule (Min et al. 2008)
    The aggregate is assumed to contain solid constituents and vacuum. 

    INPUT:
    ------
        gsize     : size of the grains 
        wav       : wavelength grid onto which the optical constants should be interpolated
        lnk_fname : list of files containing the refractive indices of the solid constituents
        wavenum   : list containing as man elements as lnk_fname and it specifies if the lnk file 
                    contains wavenumbers (1) instead of wavelengths (0)
        mfrac     : mass fraction of the solid constituents
        gdens     : bulk volume density of the material in g/cm^3 (only needed if mfrac is given instead of vfrac)
        volfrac   : volume fraction of the solid constituents
        cpol_type : scattering theory for the constituents ('DHS', 'CDE', 'MIE')

        NOTE: The mass and volume fractions should be normalized to the mass and volume of the 
        solids!!! The vacuum volume fraction is treated separately

    OUTPUT:
    -------
        res       : two dimensional numpy ndarray
        res[0,:]  : real part of the effective refractive index
        res[1,:]  : imaginary part of the effective refractive index
    """

# First check if the wavelength grid is given and raise and error if not

    if wav==None:
        print 'ERROR'
        print ' No wavelength grid is give on which the effective refractive indicies should be calculated'
        return
    else:
        wav = array(wav)
        nwav = wav.shape[0]

# Convert mass fraction to volume fraction if mfrac is specified instead of vfrac

    if mfrac!=[]:
        vfrac = array(mfrac)/gdens
        vfrac = vfrac/vfrac.sum()

# Now read the optical constants
    ncomp = len(lnk_fname)

    wav_inp = []
    n_inp   = []
    k_inp   = []
    for i in range(ncomp):
        try:
            rfile = open(lnk_fname[i], 'r')
            dum   = rfile.readlines()
            rfile.close()

            dum_w = []
            dum_n = []
            dum_k = []
            for j in range(len(dum)):
                sdum = dum[j].split()
                dum_w.append(float(sdum[0]))
                dum_n.append(float(sdum[1]))
                dum_k.append(float(sdum[2]))
            
            if wavenum[i]==1:
                if dum_w[1]>dum_w[0]:
                    wav_inp.append(1e4/array(dum_w)[::-1])
                    n_inp.append(array(dum_n)[::-1])
                    k_inp.append(array(dum_k)[::-1])
                else:
                    wav_inp.append(1e4/array(dum_w))
                    n_inp.append(array(dum_n))
                    k_inp.append(array(dum_k))
            else:
                wav_inp.append(array(dum_w))
                n_inp.append(array(dum_n))
                k_inp.append(array(dum_k))

            print lnk_fname[i] + ' is read'
        except:
            print 'ERROR'
            print lnk_fname[i]+' cannot be read'
            return

# Interpolate all optical constants to the master wavelength grid
    m = zeros([ncomp,nwav], dtype=complex128)


    lnk_wgrid_min = 0.
    lnk_wgrid_max = 99999999.
    for i in range(ncomp):    
        dn = 10.**interp(log10(wav), log10(wav_inp[i]), log10(n_inp[i]))
        dk = 10.**interp(log10(wav), log10(wav_inp[i]), log10(k_inp[i]))
     
        if wav_inp[i].min()>lnk_wgrid_min:
            lnk_wgrid_min = wav_inp[i].min()
        if wav_inp[i].max()<lnk_wgrid_max:
            lnk_wgrid_max = wav_inp[i].max()

        m[i,:] = vectorize(complex)(dn, dk)
#
# Do the black magic : Calculate the effective refractive index with the APMR mixing rule
#

    gamma = 2.44 # Fractal pre-factor
    Df    = 2.82 # Fractal dimension
    rg    = gamma * gsize**(3./Df)

    vol   = 4e0 * gsize**3. * pi / 3e0
    ffil  = (gsize/rg)**3.

    mf_abs = vol*ffil*vfrac*gdens
    print '************************************************************'
    print ' Volume fractions : ', vfrac
    print ' Average density  : ', mf_abs.sum()/vol 
    print '************************************************************'
# 
# Calculate the polarizability of a unit volume for the constituents
#
    a = zeros(nwav, dtype=complex128)
    for i in range(ncomp):
        if cpol_type == 'DHS': 
            dum1 = (6.0 * m[i,:]**2 + 3.0) / (2.0 * m[i,:]**2 - 2.0)
            dum2 = (2.0 * m[i,:]**2 + 1.0) * (m[i,:]**2 + 2.0) / (9.0 * m[i,:]**2)
            a = a + ffil * vfrac[i] * dum1 * log(dum2) * vol
        
        if cpol_type == 'CDE': 
            a = a + ffil * vfrac[i] * 2.0 * vol * ( (m[i,:]**2 / (m[i,:]**2 - 1.0)) * log(m[i,:]) - 1.0)
        
        if cpol_type == 'MIE':
            a = a + ffil * vfrac[i] * 3.0 * vol * (m[i,:]**2 -1.0) / (m[i,:]**2 + 2.0)

# 
# Calculate the effective refractive index for the aggregate as a whole using MIE theory
#

    meff = sqrt( (3.0 * vol + 2.0 * a) / (3.0 * vol - a))
    neff = meff.real
    keff = meff.imag

    ii = wav[wav<lnk_wgrid_max]>lnk_wgrid_min
    res  = [wav[ii], neff[ii], keff[ii], mf_abs.sum()/vol]
    return res

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def run_makedust(freq=None, gmin=None, gmax=None, ngs=None, lnk_fname=None, gdens=None):
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        grid.makeWavelengthGrid(ppar=ppar)
        wav = grid.wav

#
# Do we need to mix the opacities?
#
    if type(ppar['lnk_fname']).__name__=='str':
        ppar['lnk_fname'] = [ppar['lnk_fname']]

    if len(ppar['lnk_fname'])>1:

#
# Loop over all grain sizes
#
        gsize =  ppar['gsmin'] * (ppar['gsmax']/ppar['gsmin'])**(arange(ppar['ngs'], dtype=float64) / (float(ppar['ngs'])-1.))
        for igs in range(ppar['ngs']):

#
# Loop over the crystallographic axes if we have crystals
#
            lnk = []
            ncax = len(ppar['lnk_fname'])
            for iax in range(ncax):
#
# Calculate the effective refractie index for each dust species
#
                dum = get_nk_apmr(wav=wav, lnk_fname=ppar['lnk_fname'][iax], mfrac=ppar['mfrac'], gdens=ppar['gdens'], \
                        wavenum=ppar['wavenum'], gsize=gsize[igs], cpol_type=ppar['cpol_type'])  
                lnk.append(dum)
        
            neff = lnk[0][1]
            
            keff = lnk[0][2]
            for iax in range(1, ncax):
                neff = neff + lnk[iax][1]
                keff = keff + lnk[iax][2]
            
            neff = neff / float(ncax)
            keff = keff / float(ncax)

# 
# Write out the refractive indices to a file and call the mie code
#
            write_effnk(wav = lnk[0][0], n = neff, k = keff, fname='optconst_apmr.lnk')
            run_makedust(freq=cc/wav*1e4, gmin=gsize[igs], ngs=1, lnk_fname='optconst_apmr.lnk', gdens=lnk[0][3])
            dum = Popen('cp dustkappa_1.inp dustkappa_apmr_'+str(igs)+'.inp', shell=True).wait()

#
# Write the master opacity file
#

        wfile = open('dustopac.inp', 'w')
        wfile.write("%4d %s\n"%(2,'             Format number of this file'))
        wfile.write("%4d %s\n"%(ppar['ngs'],'             Nr of dust species'))
        wfile.write("%s \n"%('==================================================='))
        for igs in range(ppar['ngs']):
            wfile.write("%4d %s\n"%(1, '             Input style for species (1=dust kappa file)'))
            wfile.write("%4d %s\n"%(0, '             0=Thermal grain, 1=Quantum heated'))
            wfile.write("%4d %s\n"%(igs,'             Extension of file dustkappa_*.inp'))
            wfile.write("%s\n"%'---------------------------------------------------')
        wfile.close()


    else:
        run_makedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                lnk_fname=ppar['lnk_fname'][0], gdens=ppar['gdens'])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TEST
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ppar   = {}
#
#ppar['wbound']    = [0.1, 5.5, 30.0, 1e4]
#ppar['nw']        = [30, 300, 40]
#ppar['lnk_fname'] = [['olmg50.lnk', 'pyrmg50.lnk', 'Fo_Suto_hres_295k_b1u_combined.lnk', 'ens_p_combined.lnk'],\
#                    ['olmg50.lnk', 'pyrmg50.lnk', 'Fo_Suto_hres_295k_b2u_combined.lnk', 'ens_s1_combined.lnk'],\
#                    ['olmg50.lnk', 'pyrmg50.lnk', 'Fo_Suto_hres_295k_b3u_combined.lnk', 'ens_s2_combined.lnk']]
#ppar['wavenum']   = [0, 0, 0, 0]
#ppar['mfrac']     = [0.35, 0.45, 0.10, 0.10]
#ppar['gdens']     = [3.71, 3.2, 3.33, 2.71]
#ppar['gsmin']     = 0.1
#ppar['gsmax']     = 100.0
#ppar['ngs']       = 5
#ppar['cpol_type'] = 'DHS'
#makeopac(ppar=ppar)

