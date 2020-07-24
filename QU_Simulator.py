import numpy as np
from matplotlib import pyplot
import healpy as hp
from tqdm import tqdm
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
import glob
import importlib
from QU_Simulator_Functions import *
import Coordinates
import click

def populate_hdu(hdu,tods,todvec,rms=None):
    for pol,tod in zip(list(tods.keys()), todvec):
        for ipol in range(2):
            if isinstance(rms,type(None)):
                noise = 0
            else:
                noise = np.random.normal(size=tod.size,scale=rms)

            hdu[1].data['{}{}'.format(pol,ipol+1)] = tod + noise


def generate_simulated_tod(filename, maps,include_noise=False,output_dir=''):
    """
    Generates a simulated C-BASS TOD stream
    
    filename - name of reloaded fits file
    maps     - dictionary containing model maps (I/Q/U) in celestial frame
    """
    hdu = fits.open(filename)
    fstrp = filename.split('/')[-1]
    
    el = hdu[1].data['EL']

    if include_noise:
        rms = Tsys/np.sqrt(bw/sr) # From QU_Simulator_Functions
    else:
        rms = None

    # Calculate parallactic angle
    pa = Coordinates.pa(np.mod(hdu[1].data['RA'],2*np.pi),hdu[1].data['DEC'],hdu[1].data['MJD'],
                        antenna_observatory['OVRO'][0]*np.pi/180.,
                        antenna_observatory['OVRO'][1]*np.pi/180.,degrees=False)
    
    # Interpolate from map:
    theta, phi = (np.pi/2.-hdu[1].data['DEC']), np.mod(hdu[1].data['RA'],2*np.pi)
    tods = {pol:hp.get_interp_val(m,theta,phi) for pol, m in maps.items()}

    # Clockwise rotation matrix: +PA for telescope to sky, -PA for sky to telescope
    tods['Q'],tods['U'] = rotateQU(tods['Q'],tods['U'],-pa)

    todvec = np.array([tods['I'],tods['Q'],tods['U']])
    populate_hdu(hdu,tods,todvec,rms=rms)

    hdu.writeto('{}/{}'.format(output_dir,fstrp),overwrite=True)
    hdu.close()
    return todvec,pa,theta, phi, el

@click.command()
@click.argument('filelist')
@click.option('--wmap_filename',default='wmap_band_iqumap_r9_9yr_K_v5.fits',help='WMAP Filename with I/Q/U in HDUs 0,1,2 respectively.')
@click.option('--include_noise', default=False, help='Include white noise in simulation (Default: False)',type=bool)
@click.option('--output_directory',default='', help='Output directory for simulated TOD files.')
def main_call(filelist, wmap_filename,include_noise,output_directory):
    """
    Arguments:
    filelist - Filelist of C-BASS FITS files
    """
    main(filelist,wmap_filename,include_noise,output_directory)

def main(filelist, wmap_filename= 'wmap_band_iqumap_r9_9yr_K_v5.fits',include_noise=False,output_directory=''):

    # READ THE DATA
    filelist = np.loadtxt(filelist,dtype=str,ndmin=1)
    
    # Read I/Q/U maps from WMAP maps
    maps = {pol:hp.read_map(wmap_filename,ipol) for ipol,pol in enumerate(['I','Q','U'])}
    
    rot = hp.rotator.Rotator(coord=['G','C'])
    
    
    # Rotate maps into IAU Celestial frame
    maps['I'] = rot.rotate_map_pixel(maps['I'])
    maps['Q'],maps['U'] = rotate_QU_frame(maps['Q'],maps['U'])
    maps['U'] *= -1 # CMB to IAU convention

    # Scaling factor for WMAP -> CBASS
    factor = (4.76/22.8)**-2.7/1e3
    for k in maps.keys():
         maps[k] *= factor
        
    # Read in fits files, replace I/Q/U
    for ifile,filename in enumerate(tqdm(filelist)):
         todvec,pa,theta,phi, el = generate_simulated_tod(filename, maps, include_noise=include_noise,output_dir=output_directory)

if __name__ == "__main__":
    main_call()
