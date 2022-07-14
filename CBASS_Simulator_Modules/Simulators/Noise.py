import numpy as np
from CBASS_Simulator_Modules.Simulators import BaseModel
from astropy.io import fits

class WnoiseModel(BaseModel.BaseModel):

    def __init__(self,rms=0,use_noise_stats=False, noise_stats_hdu='',
                 **kwargs):
        super().__init__(**kwargs)

        self.rms = rms
        self.use_noise_stats = use_noise_stats
        self.noise_stats_hdu = noise_stats_hdu

    def __call__(self, filename, tod):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)
        tod = self.setup_tod(tod,hdu)
        
        if self.use_noise_stats:
            tod = self.stat_noise(tod,hdu)
        else:
            tod = self.flat_noise(tod)
        hdu.close() 
        return tod

    def stat_noise(self,tod,hdu):
        for k in tod.keys():
            rms = np.nanmean(hdu[self.noise_stats_hdu].data[f'{k}_sigma'])
            tod[k] += np.random.normal(size=tod[k].size,scale=rms)
        return tod

    def flat_noise(self,tod):
        for k in tod.keys():
            tod[k] += np.random.normal(size=tod[k].size,scale=self.rms)
        return tod


class FnoiseModel(BaseModel.BaseModel):

    def __init__(self,rms=0,alpha=0,knee=0,
                 use_noise_stats=False, noise_stats_hdu='',
                 **kwargs):
        super().__init__(**kwargs)

        self.rms  = rms
        self.knee = knee
        self.alpha= alpha
        self.use_noise_stats = use_noise_stats
        self.noise_stats_hdu = noise_stats_hdu

    def __call__(self, filename, tod):
        """
        Filename -- name of the fits file to replace
        tod  -- dictionary containing the output tod's
        """
        hdu = fits.open(filename)

        tod = self.setup_tod(tod,hdu)
        
        if self.use_noise_stats:
            tod = self.stat_noise(tod,hdu)
        else:
            tod = self.flat_noise(tod)
        hdu.close() 

        return tod

    def stat_noise(self,tod,hdu):
        for k in tod.keys():
            rms   = np.nanmean(hdu[self.noise_stats_hdu].data[f'{k}_sigma'])
            knee  = np.nanmean(hdu[self.noise_stats_hdu].data[f'{k}_fknee'])
            alpha = np.nanmean(hdu[self.noise_stats_hdu].data[f'{k}_alpha'])
            fnoise,spectrum = self.Fnoise(tod[k].size,
                                          rms,
                                          knee,
                                          alpha)
            tod[k] += fnoise
        return tod

    def flat_noise(self,tod):
        for k in tod.keys():
            fnoise,spectrum = self.Fnoise(tod[k].size,
                                          self.rms,
                                          self.knee,
                                          self.alpha)
            tod[k] += fnoise

        return tod

    def Fnoise(self, N, rms, fk, alpha, samplerate=100.):

        noise = np.random.normal(scale=rms,size=N)
        
        f = np.fft.fftfreq(N,d = 1./samplerate)
        f[0] = f[1]
        ps =  rms**2*(fk/np.abs(f))**alpha
        ps[0] = 0
        
        spec = (np.random.normal(scale=np.sqrt(ps)) + \
                1j*np.random.normal(scale=np.sqrt(ps)))
        data = np.real(np.fft.ifft(spec))*np.sqrt(N)
        return data, ps
