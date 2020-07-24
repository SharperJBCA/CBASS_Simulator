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

# SOME CONSTANTS
nodes = {1:'I1',
         2:'I2',
         3:'Q1',
         4:'U1',
         5:'Q2',
         6:'U2'}
Tsys = 30. # K
bw = 0.5e9 #Hz
sr = 100.#Hz
antenna_observatory = {'OVRO':[-118.2822,37.2339,1222.0],
                       'HartRAO':[27.6854,-25.8898,1415.0],
                       'Klerefontein':[21.98037,-30.97144,1371.0]};

def sex2deg(d,m,s,hours=False):
    """
    Convert dms to degrees
    """
    if hours:
        const = 15
    else:
        const = 1
    if d != 0:
        sign = d/np.abs(d)
    else:
        sign = 1
    return sign*(np.abs(d) + m/60. + s/60.**2)*const

def rotateQU(Q,U,pa):
    """
    Clockwise rotation matrix.
    """
    Qp = Q*np.cos(2*pa) - U*np.sin(2*pa)
    Up = Q*np.sin(2*pa) + U*np.cos(2*pa)
    return Qp,Up    

def rotate_QU_frame(q,u, coord=['G','C']):
    """
    Rotate coordinate from of QU angles on a HEALPIX map
    """
    nside = int(np.sqrt(q.size/12.))
    pix = np.arange(12*nside**2).astype(int)
    vecs = np.array(hp.pix2vec(nside, pix))

    rot = hp.rotator.Rotator(coord=coord)
    
    # First rotate the map pixels
    qp = rot.rotate_map_pixel(q)
    up = rot.rotate_map_pixel(u)
    qp[qp < -1e25] = hp.UNSEEN
    up[up < -1e25] = hp.UNSEEN

    # Then caculate the reference angles in the rotated frame
    vecs_r = rot(vecs,inv=True)
    angles = rot.angle_ref(vecs_r)
    
    L_map = (qp + 1j * up)*np.exp(1j*2*angles)

    return np.real(L_map), np.imag(L_map)

def mueller_matrix(alpha, beta):
    """
    Define Mueller matrix with I->Q/U leakage terms.
    """
    M = np.array([
        [1    , alpha, beta],
        [alpha,   1  ,   0 ],
        [beta ,   0  ,   1 ]
    ])
    return M
