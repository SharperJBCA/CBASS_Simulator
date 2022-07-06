import numpy as np
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

from CBASS_Simulator_Modules.Tools import Coordinates
from CBASS_Simulator_Modules.Simulators.QU_Simulator_Functions import *

class LoadModel(BaseModel.BaseModel):

    def __init__(self,load_rms=0,**kwargs):
        super().__init__(**kwargs)
        self.load_rms = load_rms

    def __call__(self, filename, tod):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)
        tod = self.setup_tod(tod,hdu)
        tod = self.basic_load_signal(tod)

        hdu.close()
        return tod

    def basic_load_signal(self,tod):
        """
        Generate a load signal with just a sine wave
        """
        start = 0
        end   = tod['I1'].size

        for i,(s,e) in enumerate(zip([start],[end])):
            sample_rate = 100. # Hz
            t = np.arange(e-s)/sample_rate
            v = 1.2 # Hz
            phase = 0
            for k in tod.keys():
                tod[k][s:e]  += np.sin(2*np.pi*t*v + phase)*self.load_rms            
        return tod
