import numpy as np
from matplotlib import pyplot
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

def load_ground_model(filename):
    d = np.loadtxt(filename)

    keys = ['AZ','I1','I2','Q1','Q2','U1','U2']

    return {k:d[:,i] for i,k in enumerate(keys)}


def ground_signal(az, tods, pa, ground_model_file):
    """
    """

    
    ground_model = load_ground_model(ground_model_file)
    for k in tods.keys():
        tods[k] += np.interp( az, ground_model['AZ'],ground_model[k],period=360)

    return tods


class GroundModel(BaseModel.BaseModel):

    def __init__(self, 
                 ground_model_file='',
                 random_offsets=False,
                 random_slopes=False,
                 **kwargs):
        super().__init__(**kwargs)

        self.random_offsets = random_offsets
        self.random_slopes = random_slopes
        self.ground_model_file = ground_model_file

    def __call__(self, filename, tod):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)
        tod = self.setup_tod(tod,hdu)

        ground_model = load_ground_model(self.ground_model_file)
        for k in tod.keys():
            tod[k] += np.interp( hdu[self.ihdu].data['AZ'], 
                                 ground_model['AZ'],ground_model[k],period=360)

        hdu.close()
        return tod
