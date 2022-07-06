import numpy as np
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

from CBASS_Simulator_Modules.Tools import Coordinates
from CBASS_Simulator_Modules.Simulators.QU_Simulator_Functions import *

class SkyModel(BaseModel.BaseModel):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)

    def __call__(self, filename, tod, maps=None):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)
        tod = self.setup_tod(tod,hdu)

        pa = Coordinates.pa(np.mod(hdu[1].data['RA'],2*np.pi),\
                            hdu[1].data['DEC'],hdu[1].data['MJD'],
                            antenna_observatory['OVRO'][0]*np.pi/180.,
                            antenna_observatory['OVRO'][1]*np.pi/180.,degrees=False)

        theta, phi = (np.pi/2.-hdu[self.ihdu].data['DEC']), \
                     np.mod(hdu[self.ihdu].data['RA'],2*np.pi)
        nside  = hp.npix2nside(maps['I'].size)
        pixels = hp.ang2pix(nside,theta,phi)
        sky_tod = {'{}{}'.format(pol,i+1):m[pixels] for pol, m in maps.items() for i in range(2)}
        # Clockwise rotation matrix: +PA for telescope to sky, -PA for sky to telescope
        tod['I1'] = sky_tod['I1']
        tod['I2'] = sky_tod['I2']
        tod['Q1'],tod['U1'] = rotateQU(sky_tod['Q1'],sky_tod['U1'],-pa)
        tod['Q2'],tod['U2'] = rotateQU(sky_tod['Q2'],sky_tod['U2'],-pa)

        hdu.close()
        return tod
