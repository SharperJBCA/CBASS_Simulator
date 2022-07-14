import numpy as np
from matplotlib import pyplot
import healpy as hp
from tqdm import tqdm
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
import glob
import importlib
import Coordinates
import click
import Analysis_Funcs
import os

from CBASS_Simulator_Modules import Simulators
from CBASS_Simulator_Modules.Simulators.QU_Simulator_Functions import *
from CBASS_Simulator_Modules.Tools import Parser

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def populate_hdu(hdu,tods):
    for pol,tod in tods.items():
        hdu[1].data['{}'.format(pol)] = tod 

def get_kwargs(parameters,ancil_data):
    if 'ancil_data' in parameters:
        return {k:ancil_data[k] for k in parameters['ancil_data']}
    else:
        return {}

def load_models(model_names, parameters, ancil_data):
    """
    return a list of model functions
    """
    models = {model_name:{'function':Simulators.model_functions[model_name](**parameters[model_name]),
                          'kwargs':get_kwargs(parameters[model_name],ancil_data)} for model_name in model_names}
    return models
    

def generate_simulated_tod(filename, parameters, ancil_data):
    """
    Generates a simulated C-BASS TOD stream
    
    filename - name of reloaded fits file
    maps     - dictionary containing model maps (I/Q/U) in celestial frame
    """
    models = load_models(parameters['Main']['models'],parameters, ancil_data)
    tod = {}
    for model_name, model in models.items():
        tod = model['function'](filename, tod, **model['kwargs'])

    hdu = fits.open(filename)
    # Fix date of filename
    mjd = hdu[1].data['MJD']
    start_date = Time(mjd[0],format='mjd').datetime
    start_date_name = datetime.strftime(start_date, '%d-%b-%Y:%H:%M:%S_reload.fits')
    fstrp = '{}'.format(start_date_name)

    output_directory = parameters['Main']['output_directory']

    populate_hdu(hdu,tod)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    hdu.writeto('{}/{}'.format(output_directory,fstrp),overwrite=True)
    hdu.close() 
    return 

@click.command()
@click.option('--parameter-file',default='none')
@click.option('--filelist',default='none')
@click.option('--modelsky_filename',default='wmap_band_iqumap_r9_9yr_K_v5.fits',help='WMAP Filename with I/Q/U in HDUs 0,1,2 respectively.')
@click.option('--ground_model_file',default='',help='Ground model file to use')
@click.option('--include-noise', default=False, help='Include white noise in simulation (Default: False)',type=bool)
@click.option('--include-sky', default=False, help='Include sky signal (Default: True)',type=bool)
@click.option('--include-load', default=False, help='Include load signal (Default: True)',type=bool)
@click.option('--include-fnoise', default=False, help='Include load signal (Default: True)',type=bool)
@click.option('--include-ground', default=False, help='Include ground signal (Default: True)',type=bool)
@click.option('--load_rms', default=0, help='Set load signal rms (will scale by root(2)) (Default: 0)',type=float)
@click.option('--output-directory',default='', help='Output directory for simulated TOD files.')
@click.option('--wnoise_rms', default=3e-2, help='',type=float)
@click.option('--fnoise_fknee', default=0.1, help='',type=float)
@click.option('--fnoise_alpha', default=1, help='',type=float)
def main_call(**kwargs):#filelist, wmap_filename,include_noise,include_sky,include_load,output_directory):
    """
    Arguments:
    filelist - Filelist of C-BASS FITS files
    """
    main( **kwargs)#filelist,wmap_filename,include_noise,include_sky,include_load,output_directory)

def main(**kwargs):
    np.random.seed(101)

    if not 'none' in kwargs['parameter_file']:
        parameters = Parser.Parser(kwargs['parameter_file'])
    else:
        parameters = {'Main':kwargs}

    # READ THE DATA
    filelist = np.loadtxt(parameters['Main']['filelist'],dtype=str,ndmin=1)
    
    ancil_data = {}
    # Read I/Q/U maps from WMAP maps
    if 'sky' in parameters['Main']['models']:
        ancil_data['maps'] = {pol:hp.read_map(parameters['Main']['modelsky_filename'],ipol) for ipol,pol in enumerate(['I','Q','U'])}
    else:
        maps = {'I':None,'Q':None,'U':None}


    file_idx = np.sort(np.mod(np.arange(len(filelist),dtype=int),size).astype(int))
    file_idx = np.where((file_idx == rank))[0]
    
        
    # Read in fits files, replace I/Q/U
    for ifile,filename in enumerate(tqdm(filelist[file_idx])):
         generate_simulated_tod(filename, parameters, ancil_data)

if __name__ == "__main__":
    main_call()
