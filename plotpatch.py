def plotpatch(P, hmin, hmax, npatch):
    """Function to plot the patches for a given height range.
    
       Args:
       -----
       P       : rotational period (seconds)
       hmin    : minimum emission height (in km)
       hmax    : maximum emission height (in km)
       npatch  : number of emission patches
       
       Returns:
       --------
       A plot of the patches projected on to observational plane.
    
    """
    
    angle = np.linspace(-180, 180, num=50, endpoint=True)
    patchwidths = patch_width(P, hmin, hmax)
    
#   choose random patch widths (wp) depending on how patches specify:
    np.random.shuffle(patchwidths)
    wp = patchwidths[0:npatch]
    
    #x0 = wp[0] * d2r.cosD(angle)
    #y0 = wp[0] * d2r.sinD(angle)
    #x1 = wp[1] * d2r.cosD(angle)
    #y1 = wp[1] * d2r.sinD(angle)
    #for i in np.arange(len(wp)):
    #    x = wp[i] * d2r.cosD(angle)
    #    y = wp[i] * d2r.sinD(angle)  
    X = []
    Y = []
    
    for i in np.arange(len(wp)):
        X.append(wp[i] * d2r.cosD(angle))
        Y.append(wp[i] * d2r.sinD(angle))
    
    x0 = X[0]
    y0 = Y[0]
    x1 = X[1]
    y1 = Y[1]
    x2 = X[2]
    y2 = Y[2] 
    #fig = plt.figure()
    #axs = fig.gca()
    #axs.set_xticks(np.arange(-90, 90, 5))
    #axs.set_yticks(np.arange(-90, 90, 5))
    plt.plot(x0, y0)
    plt.plot(x1, y1)
    plt.plot(x2, y2)
    plt.grid()
    plt.show()
    
    
############################ simple test ##########################
if __name__ == "__main__":
    P = 0.16
    hmin = 30
    hmax = 990
    npatch = 3
    plotpatch(P, hmin, hmax, npatch)
