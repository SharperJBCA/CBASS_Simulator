# Dipole, code for adding a CMB dipole to the data
import numpy as np
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

from CBASS_Simulator_Modules.Tools import Coordinates
from CBASS_Simulator_Modules.Simulators.QU_Simulator_Functions import *

class DipoleModel(BaseModel.BaseModel):

    def __init__(self,dummy=False,**kwargs):
        super().__init__(**kwargs)

        # CMB dipole parameters, from Table 2 of "Planck 2018 results. I. Overview,
        # and the cosmological legacy of Planck": 2020. A&A, 641, p.A1.
        # corresponding direction in celestial coords: right ascension and declination
        self.dipole={'ra':167.942,
                    'dec':  -6.944 ,
                    'amplitude': 0.0033621} # amplitude in K (3.36208  0.00099 mK)
                    

    def __call__(self, filename, tod):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)
        tod = self.setup_tod(tod,hdu)

        theta, phi = (np.pi/2.-hdu[self.ihdu].data['DEC']), \
                     np.mod(hdu[self.ihdu].data['RA'],2*np.pi)


        # Convert spherical angles of to a unit vector in the direction of the dipole
        nside = 64
        npix = 12*nside**2
        pix = np.arange(npix,dtype=int)
        vec = np.array(hp.pixelfunc.pix2vec(nside,pix))
        vec_r = np.array(hp.ang2vec((90-self.dipole['dec'])*np.pi/180.,
                           self.dipole['ra']*np.pi/180.))
        C = vec.dot(vec.T)/vec.shape[1]*3.
        res = np.linalg.solve(C,vec_r[:,None])

        vec_tod = np.array(hp.pixelfunc.ang2vec(theta,phi)).T
        for k in tod.keys():
            if 'I' in k:
                tod[k] += (vec_tod[0]*res[0]+vec_tod[1]*res[1]+vec_tod[2]*res[2])*self.dipole['amplitude']
        

        hdu.close()
        return tod
