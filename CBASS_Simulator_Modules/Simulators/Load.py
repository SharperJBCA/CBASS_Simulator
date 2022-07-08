import numpy as np
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

from CBASS_Simulator_Modules.Tools import Coordinates
from CBASS_Simulator_Modules.Simulators.QU_Simulator_Functions import *
import time, os
import pickle

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

def interpolate_frequencies(time_series):

    # Open data

    try:
        delta_t, freq = pickle.load( open( './mains_data.p', 'rb' ) )

    except:
        raise NameError('Could not find mains data: You need the file mains_data.p to be in the same directory as this script!')
        
    # Seed the random number generator based on the current time and the current process number
    seed = int(time.time() + os.getpid())
    np.random.seed(seed)

    # Generate one random integer along the mains data - this is going to be the starting index
    starting_index = int(np.round(np.random.uniform(0,np.shape(freq)[0]-1)))

    # Rotate a month worth's of mains frequency data randomly
    freq = np.roll(freq, starting_index)

    # Interpolate frequencies
    from scipy.interpolate import interp1d
    f = interp1d(delta_t, freq)
    interpolated_frequencies = f(time_series)

    return interpolated_frequencies

def simulate_sinusoidal(amplitude, time_series, varying_frequency=False):

    ''' 1.2 Hz signal sinusoidal simulator. The inputs are defined below:

    * amplitude = amplitude(s) of the sinusoid in the same units as the data. This can be fixed or varying. If you
      want it to vary you should pass an array with a size matching the "time_series".
    * time_series = a series of the times of the data in seconds (the relative start, end times do not matter).
    * varying_frequency = True/False. By default, I use a 1.2 Hz fixed frequency signal, but if you want, you can set this
      to true so UK mains January 2020 data is used instead to simulate variations in the frequency. You will need the
      file with the data (mains_data.fits) to be in the same directory.

    The starting phase of the sinusoid is uniformly randomised every single time you run the function. '''

    if varying_frequency==False:
        frequency = 1.2 # set fixed frequency

    elif varying_frequency==True:
        frequency = interpolate_frequencies(time_series) # simulate varying frequencies

    else:
        raise ValueError('The varying_frequency argument needs to be either True or False')


    # Normalise time series so that the first point is zero
    time_series = np.subtract(time_series, np.nanmin(time_series))

    # Add random starting phase
    starting_phase = np.random.uniform(low=0, high=2*np.pi)

    # Generate sinusoid
    sinusoidal_signal = np.multiply(amplitude, np.sin( np.add( np.multiply(2*np.pi, np.multiply(frequency,time_series) ) , starting_phase) ))


    return sinusoidal_signal
