import numpy as np
import matplotlib.pyplot as plt
import rho


def patch_center(P, hmin, hmax, npatch):
    """Function find centres of the patches
       
       Args:
       -----
       P      : rotatinal period
       hmin   : minimum emission height (in km).
       hmax   : maximum emission height (in km).
       npatch : number of emission patches
       
       
       Returns:
       --------
       patch_centerx : the patch center projection on the x-axis 
       patch_centery : the patch center projection on the y-axis 
    """
#   centers of the patches in the rotational plane (x,y coords):
#    np.random.seed(0)
    patch_center = 2 * np.pi * np.random.random()
    patch_centerx = np.zeros(npatch)
    patch_centery = np.zeros(npatch)
    opa = rho.rho(P, hmin, hmax, npatch)  # opening angle of the beam (rho is degrees)

#   for short periods 
    if P <= 0.15:
        if hmin >= 950 and hmax <= 1000:
            for i in np.arange(npatch):
                patch_centerx[i] = opa * np.sin(patch_center)
                patch_centery[i] = opa * np.cos(patch_center)
                
#   for longer periods
    if P > 0.15:
        if hmin >= 20 and hmax <= 1000:
            for i in np.arange(npatch):
                patch_centerx[i] = opa[i] * np.sin(patch_center)
                patch_centery[i] = opa[i] * np.cos(patch_center)
        
    return patch_centerx, patch_centery


########################### simple test ###########################
if __name__ == "__main__":
    P = 0.16
    hmin = 30
    hmax = 990
    npatch = 3
    cx, cy = patch_center(P, hmin, hmax, npatch)
    print cx, cy 
