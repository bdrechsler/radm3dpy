import numpy as np
import shutil
import os

class StellarAtm():
    """
    """

    def __init__(self):

        self.wav   = None
        self.fnu   = None
        self.bnu   = None
        self.logg  = None
        self.teff  = None
        self.metal = None
        self.fname = None

        self.kuruczDir  = '/disk2/juhasz/Data/kurucz/'
        self.nextgenDir = '/disk2/juhasz/Data/NextGen/SPECTRA/'

    
    def getAtmModel(self, teff=0., logg=None, mstar=None, lstar=None, rstar=None, iwav=None, model='kurucz', wmax=7.):
        """
        """
        ss = 5.6703e-05
        kk = 1.3807e-16
        cc = 29979245800.0
        rs = 69600000000.0
        ls = 3.8525e+33
        ms = 1.99e33
        gg = 6.672e-08
        
        #
        # Get the stellar parameters
        #
        if not rstar:
            rstar = np.sqrt(lstar / (4. * np.pi * ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * ss * teff**4
        
        if logg:
            mstar = 10.**logg * rstar**2/gg
        else:
            logg  = np.log10(gg*mstar/rstar**2)


        #
        # Get the atmosphere model
        #
        if model.strip().lower()=='kurucz':
            sp   = self.getSpectrumKurucz(teff=teff, logg=logg, lstar=lstar)
            wav  = sp['wav'] 
            ilnu = self.rebinSpectrum(wav=sp['wav'], fnu=sp['lnu'], iwav=iwav)
        elif model.strip().lower()=='nextgen':
            sp   = self.getSpectrumNextGen(teff=teff, logg=logg, lstar=lstar)
            wav  = sp['wav'] 
            ilnu = self.rebinSpectrum(wav=sp['wav'], fnu=sp['lnu'], iwav=iwav)
        else:
            print 'ERROR'
            print 'Unknown atmosphere model : ', model
            print 'The model keyword can only be either "Kurucz" or "NextGen"'
            print 'Current model name is : ', model
            
            return

        ilnu = ilnu.clip(1e0, 1e99)
        ii = (iwav>wmax)
        if ii.__contains__(True):
            f0       = 10.**np.interp(np.log10(wmax), np.log10(iwav), np.log10(ilnu))
            ilnu[ii] = f0 * (iwav[ii]/wmax)**(-2.)

        return {'wav':iwav, 'lnu':ilnu}
            

    def getSpectrumKurucz(self, teff=0., logg=0., lstar=None, rstar=None, wmin=None, wmax=None, nwav=None, wav=None):
        """
        Returns the spectrum in lnu in erg/s/cm/cm/Hz
        """
       
        ss = 5.6703e-05
        kk = 1.3807e-16
        cc = 29979245800.0
        rs = 69600000000.0
        ls = 3.8525e+33
        ms = 1.99e33
        gg = 6.672e-08


        # 
        # Get the stellar radius
        # 

        if not rstar:
            rstar = np.sqrt(lstar / (4. * np.pi * ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * ss * teff**4

        mstar = 10.**logg * rstar**2/gg

        print '-------------------------------------------'
        print 'Interpolating in Kurucz model atmospheres'
        print 'Stellar parameters: '
        print '-------------------------------------------'
        print 'Teff [K]          : ', teff
        print 'Radius [Rsun]     : ', rstar/rs
        print 'Luminosity [Lsun] : ', lstar/ls
        print 'Mass [Msun]       : ', mstar/ms 
        print 'logg              : ', logg
        print '-------------------------------------------'


        dum = self.readKuruczGrid(fname=self.kuruczDir+'/fp00k2.pck')
        # 
        # Bracket in Teff
        #
        teff_grid = np.unique(dum['teff'])

        ii   = abs(teff_grid-teff).argmin()
        idt1 = ii
        if teff_grid[ii]>teff:
            idt1 = ii-1
        else:
            idt1 = ii
        idt2 = idt1+1

        # 
        # Bracket in Logg
        #
        
        ii = (dum['teff']==teff_grid[idt1])
        logg_grid_lower = dum['logg'][ii]
        ii = (dum['teff']==teff_grid[idt2])
        logg_grid_upper = dum['logg'][ii]

      
        ii   = abs(logg_grid_lower-logg).argmin()
        if logg<logg_grid_lower[0]:
            idg1 = -1
            idg2 = 0
        elif logg>logg_grid_lower[-1]:
            idg1 = logg_grid_lower.shape[0]-1
            idg2 = -1
        else:
            idg1 = ii
            if logg_grid_lower[ii]>logg:
                idg1 = ii-1
            else:
                idg1 = ii
            idg2 = idg1+1

        idgl1 = idg1
        idgl2 = idg2


        ii   = abs(logg_grid_upper-logg).argmin()
        if logg<logg_grid_upper[0]:
            idg1 = -1
            idg2 = 0
        elif logg>logg_grid_upper[-1]:
            idg1 = logg_grid_upper.shape[0]-1
            idg2 = -1
        else:
            idg1 = ii
            if logg_grid_upper[ii]>logg:
                idg1 = ii-1
            else:
                idg1 = ii
            idg2 = idg1+1
        
        idgu1 = idg1
        idgu2 = idg2

        #
        # Check if we need to do a 3point bilinear interpolation
        #

        if ((idgl1<0)|(idgl2<0)|(idgu1<0)|(idgu2<0)):
            x   = []
            y   = []
            sp  = []
            spc = []
            if idgl1>=0:
                x.append(teff_grid[idt1])
                y.append(logg_grid_lower[idgl1])
                ii = ((dum['teff']==teff_grid[idt1])&(dum['logg']==logg_grid_lower[idgl1]))
                sp.append(squeeze(dum['fnu'][ii,:]))
                spc.append(squeeze(dum['fnucont'][ii,:]))
            if idgl2>=0:
                x.append(teff_grid[idt1])
                y.append(logg_grid_lower[idgl2])
                ii = ((dum['teff']==teff_grid[idt1])&(dum['logg']==logg_grid_lower[idgl2]))
                sp.append(squeeze(dum['fnu'][ii,:]))
                spc.append(squeeze(dum['fnucont'][ii,:]))
            if idgu1>=0:
                x.append(teff_grid[idt2])
                y.append(logg_grid_upper[idgu1])
                ii = ((dum['teff']==teff_grid[idt2])&(dum['logg']==logg_grid_upper[idgu1]))
                sp.append(squeeze(dum['fnu'][ii,:]))
                spc.append(squeeze(dum['fnucont'][ii,:]))
            if idgu2>=0:
                x.append(teff_grid[idt2])
                y.append(logg_grid_upper[idgu2])
                ii = ((dum['teff']==teff_grid[idt2])&(dum['logg']==logg_grid_upper[idgu2]))
                sp.append(squeeze(dum['fnu'][ii,:]))
                spc.append(squeeze(dum['fnucont'][ii,:]))

            if len(x)!=3:
                print 'Something went wrong..'
                print 'Only 3 valid points should have been found and I found '+("%d"%len(x))
                return -1

            else:
                print 'Bracketed spectrum with Teff : ',teff , ' and logg : ',logg
                print 'Teff grid : ', x
                print 'Logg grid : ', y
            
            c1 = ( (y[1]-y[2])*(teff-x[2]) + (x[2]-x[1])*(logg-y[2]) ) / ( (y[1]-y[2])*(x[0]-x[2]) + (x[2]-x[1])*(y[0]-y[2]) )
            c2 = ( (y[2]-y[0])*(teff-x[2]) + (x[0]-x[2])*(logg-y[2]) ) / ( (y[1]-y[2])*(x[0]-x[2]) + (x[2]-x[1])*(y[0]-y[2]) )
            c3 = 1.-c1-c2
            
            fnu = c1*sp[0] + c2*sp[1] + c3*sp[2]
            fnucont = c1*spc[0] + c2*spc[1] + c3*spc[2]


            #fig = figure()
            #loglog(dum['wav'], sp[0], 'b-')
            #loglog(dum['wav'], sp[1], 'r-')
            #loglog(dum['wav'], sp[2], 'g-')
            #loglog(dum['wav'], fnu, 'k-')
            #dum = raw_input() 

        else:
            print 'Bracketed spectrum with Teff : ',teff ,' and logg : ',logg
            print 'Teff grid : ', teff_grid[idt1], teff_grid[idt2]
            print 'Logg grid : ', logg_grid_lower[idgl1], logg_grid_lower[idgl2]
            #
            # Do the standard four point bilinear interpolation
            #
            ii = ((dum['teff']==teff_grid[idt1])&(dum['logg']==logg_grid_lower[idgl1]))
            sp11 = np.squeeze(dum['fnu'][ii,:])
            ii = ((dum['teff']==teff_grid[idt1])&(dum['logg']==logg_grid_lower[idgl2]))
            sp12 = np.squeeze(dum['fnu'][ii,:])
            ii = ((dum['teff']==teff_grid[idt2])&(dum['logg']==logg_grid_upper[idgu1]))
            sp22 = np.squeeze(dum['fnu'][ii,:])
            ii = ((dum['teff']==teff_grid[idt2])&(dum['logg']==logg_grid_upper[idgu2]))
            sp21 = np.squeeze(dum['fnu'][ii,:])

            c11 = (teff_grid[idt2] - teff) * (logg_grid_upper[idgu2]-logg)
            c12 = (teff_grid[idt2] - teff) * (logg-logg_grid_upper[idgu1])
            c22 = (teff-teff_grid[idt1]) * (logg-logg_grid_lower[idgl1])
            c21 = (teff-teff_grid[idt1]) * (logg_grid_lower[idgl2]-logg)  
            c00 = 1./( (teff_grid[idt2]-teff_grid[idt1]) * (logg_grid_lower[idgl2]-logg_grid_lower[idgl1]))

            fnu     = c00 * (c11*sp11 + c12*sp12 + c22*sp22 + c21*sp21)
            fnucont = c00 * (c11*sp11 + c12*sp12 + c22*sp22 + c21*sp21)
                    
      
        nu  = cc/dum['wav']*1e4
        lum = (0.5 * abs(nu[1:] - nu[:-1]) * (fnu[1:] + fnu[:-1])).sum()
        fnu *= lstar / lum

        return {'wav':dum['wav'], 'lnu':fnu, 'fnucont':fnucont}

    def readKuruczGrid(self, fname=''):
        """
        """

        rfile = open(fname, 'r')
        #
        # Skip the program part
        #
        for i in range(22):
            dum = rfile.readline()
        
        # 
        # Read the wavelength grid
        #
        wav     = []
        n       = 10
        for i in range(153):
            dum = rfile.readline().split()
            for j in range(len(dum)):
                wav.append(float(dum[j]))

        #
        # Convert the wavelength in Angstrom to micron
        #
        wav = np.array(wav) * 1e-3
        #
        # Now read the grid of spectra
        #
        nwav         = wav.shape[0]
        tgrid_list   = []
        logg_list    = []
        fnu_list     = []
        fnucont_list = []


        #
        # Read the first section header
        #
        dum = rfile.readline()
        while dum.strip()!='':
            #print '>>>> ', dum, len(dum.strip())
            sdum = dum.split()
            tgrid_list.append(float(sdum[1]))
            logg_list.append(float(sdum[3]))

            # 
            # Read the stellar spectrum
            #
            arr = []
            for i in range(152):
                dum = rfile.readline()
                for j in range(8):
                    arr.append(float(dum[j*n:(j+1)*n]))
            dum = rfile.readline()
            for j in range(5):
                arr.append(float(dum[j*n:(j+1)*n]))
            fnu_list.append(np.array(arr))
            # 
            # Read the continuum spectrum
            #
            arr = []
            for i in range(152):
                dum = rfile.readline()
                for j in range(8):
                    arr.append(float(dum[j*n:(j+1)*n]))
            dum = rfile.readline()
            for j in range(5):
                arr.append(float(dum[j*n:(j+1)*n]))
            fnucont_list.append(np.array(arr))
       
            #
            # Read the next section header
            #
            dum = rfile.readline()

        rfile.close()

        teff_grid  = np.array(tgrid_list)
        logg_grid  = np.array(logg_list)
        fnu        = np.array(fnu_list)
        fnucont    = np.array(fnucont_list)
      

        return {'wav':wav, 'fnu':fnu, 'fnucont':fnucont, 'teff':teff_grid, 'logg':logg_grid, 'nwav':nwav}

    def readNextGenSpectrum(self,fname=''):
        """
        """

        print 'Reading : ', fname
        rfile = open(fname, 'r')
        dum = rfile.readline()
        sdum = dum.split()
        teff = float(sdum[0])
        logg = float(sdum[1])
        mph  = float(sdum[2])
        dum  = rfile.readline()
        nwav = float(dum.split()[0])

        bigline = []
        dum = rfile.readline()
        while dum.strip()!='':
            sdum = dum.split()
            for i in range(len(sdum)):
                bigline.append(float(sdum[i]))
            dum = rfile.readline()
        rfile.close()
        
        bigline = np.array(bigline)
        # Convert wavelength from angstrom to micron
        wav = bigline[:nwav]/1e4
        fnu = bigline[nwav:2*nwav]
        bnu = bigline[nwav*2:nwav*3]

        ii = wav.argsort()
        wav = wav[ii]
        fnu = fnu[ii]*1e-8 * wav * 1e4 /np.pi / (29979245800.0/wav*1e4)
        bnu = bnu[ii]*1e-8 * wav * 1e4 /np.pi / (29979245800.0/wav*1e4)
        
        #
        # The unit is now erg/s/cm/Hz/ster
        #

        return {'teff':teff, 'logg':logg, 'mph':mph, 'nwav':nwav, 'wav':wav, 'fnu':fnu, 'bnu':bnu} 


    def rebinSpectrum(self, wav=None, fnu=None, iwav=None):
        """
        """

        # 
        # Generate the interface grid
        #
        iiwav = np.zeros(iwav.shape[0]+1, dtype=float)
        iiwav[1:-1] = np.sqrt(iwav[1:]*iwav[:-1])
        iiwav[0]    = iwav[0]**2 / iwav[1]
        iiwav[-1]   = iwav[-1]**2 / iwav[-2]


        ifnu = np.zeros(iiwav.shape[0]-1, dtype=float)
        for i in range(iiwav.shape[0]-1):
            print i, iiwav.shape[0]
            ii =( (wav>iiwav[i])&(wav<=iiwav[i+1]) )
            if ii.__contains__(True):
                x = wav[ii]
                y = fnu[ii]
                ifnu[i] = (0.5 * (x[1:]-x[:-1])*(y[1:]+y[:-1])).sum() / (iiwav[i+1]-iiwav[i])

        return ifnu

    def getSpectrumNextGen(self, teff=0., logg=0., lstar=None, rstar=None, wmin=None, wmax=None, nwav=None, wav=None):
        """
        """
        ss = 5.6703e-05
        kk = 1.3807e-16
        cc = 29979245800.0
        rs = 69600000000.0
        ls = 3.8525e+33
        ms = 1.99e33
        gg = 6.672e-08
        
        # 
        # Get the stellar radius
        # 

        if not rstar:
            rstar = np.sqrt(lstar / (4. * np.pi * ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * ss * teff**4

        mstar = 10.**logg * rstar**2/gg

        print '-------------------------------------------'
        print 'Interpolating in NextGen model atmospheres'
        print 'Stellar parameters: '
        print '-------------------------------------------'
        print 'Teff [K]          : ', teff
        print 'Radius [Rsun]     : ', rstar/rs
        print 'Luminosity [Lsun] : ', lstar/ls
        print 'Mass [Msun]       : ', mstar/ms 
        print 'logg              : ', logg
        print '-------------------------------------------'


        teff_grid = np.append( (np.arange(31.)+9.), (np.arange(30.)*2.+40.))
        logg_grid = np.arange(6.)*0.5 + 3.5


        # Bracket the input teff and logg values
        ii   = abs(teff_grid-teff/100.).argmin()
        idt1 = ii
        if teff_grid[ii]>teff/100.:
            idt1 = ii-1
        else:
            idt1 = ii
        idt2 = idt1+1


        ii   = abs(logg_grid-logg).argmin()
        if logg<logg_grid[0]:
            idg1 = 0
            idg2 = 0
        elif logg>logg_grid[-1]:
            idg1 = logg_grid.shape[0]-1
            idg2 = logg_grid.shape[0]-1
        else:
            idg1 = ii
            if logg_grid[ii]>logg:
                idg1 = ii-1
            else:
                idg1 = ii
            idg2 = idg1+1

        
        print 'Bracketing spectrum : ', teff, logg
        print 'Teff  : ', teff_grid[idt1], teff_grid[idt2]
        print 'log(g): ', logg_grid[idg1], logg_grid[idg2]

        # Generate the spectral file names
        mph = '0.0'
        fname_11 = 'lte'+("%d"%teff_grid[idt1])+'-'+("%.1f"%logg_grid[idg1])+'-'+mph+'.NextGen.spec.gz'
        fname_12 = 'lte'+("%d"%teff_grid[idt1])+'-'+("%.1f"%logg_grid[idg2])+'-'+mph+'.NextGen.spec.gz'
        fname_22 = 'lte'+("%d"%teff_grid[idt2])+'-'+("%.1f"%logg_grid[idg2])+'-'+mph+'.NextGen.spec.gz'
        fname_21 = 'lte'+("%d"%teff_grid[idt2])+'-'+("%.1f"%logg_grid[idg1])+'-'+mph+'.NextGen.spec.gz'

        # Create a directory
        if os.path.exists("./tmp")==True:
            shutil.rmtree("./tmp")

        os.system('mkdir tmp')
        os.system('cp -v '+self.nextgenDir+'/'+fname_11+' ./tmp')
        os.system('cp -v '+self.nextgenDir+'/'+fname_12+' ./tmp')
        os.system('cp -v '+self.nextgenDir+'/'+fname_22+' ./tmp')
        os.system('cp -v '+self.nextgenDir+'/'+fname_21+' ./tmp')

        # Unzip the files
        os.chdir('./tmp')
        os.system('gunzip '+fname_11)
        os.system('gunzip '+fname_12)
        os.system('gunzip '+fname_22)
        os.system('gunzip '+fname_21)
        os.chdir('../')

        # Read the spectra
        sp11 = self.readNextGenSpectrum(fname='./tmp/'+fname_11[:-3])
        sp12 = self.readNextGenSpectrum(fname='./tmp/'+fname_12[:-3])
        sp22 = self.readNextGenSpectrum(fname='./tmp/'+fname_22[:-3])
        sp21 = self.readNextGenSpectrum(fname='./tmp/'+fname_21[:-3])


        # Do the interpolation
        c11 = (teff_grid[idt2] - teff/100.)*(logg_grid[idg2]-logg)
        c12 = (teff_grid[idt2] - teff/100.)*(logg-logg_grid[idg1])
        c22 = (teff/100. - teff_grid[idt1])*(logg-logg_grid[idg1])
        c21 = (teff/100. - teff_grid[idt1])*(logg_grid[idg2]-logg)
        c00 = 1./((teff_grid[idt2]-teff_grid[idt1]) * (logg_grid[idg2]-logg_grid[idg1]))
        
        fnu = c00 * (c11*sp11['fnu'] + c12*sp12['fnu'] + c22*sp22['fnu'] + c21*sp21['fnu']) 
        bnu = c00 * (c11*sp11['bnu'] + c12*sp12['bnu'] + c22*sp22['bnu'] + c21*sp21['bnu']) 
        
        shutil.rmtree('./tmp')

        #
        # Scale the spectrum to give the same luminosity as required
        #
        nu = cc/sp11['wav']*1e4
        lum = (0.5 * abs(nu[1:] - nu[:-1]) * (fnu[1:] + fnu[:-1])).sum()
        fnu *= lstar / lum
        
        lum = (0.5 * abs(nu[1:] - nu[:-1]) * (bnu[1:] + bnu[:-1])).sum()
        bnu *= lstar / lum
        
        return {'wav':sp11['wav'], 'lnu':fnu, 'bnu':bnu}



