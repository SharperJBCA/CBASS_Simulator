import numpy as np

class BaseModel:

    def __init__(self,ihdu=1, 
                 data_keys=['I1','I2','Q1','Q2','U1','U2'],
                 **kwargs):
        self.ihdu = ihdu
        self.data_keys = data_keys

    def setup_tod(self,tod,hdu):
        """
        Makes sure there is data inside the tod
        """
        
        if not tod: # if it is empty populate it
            return {k:np.zeros(hdu[self.ihdu].data[k].size) for k in self.data_keys}
        else:
            return tod
