from multiprocessing import Pool
import os, sys
import numpy as np

# vars=(xi, yi, zi, Dir, Cube)
# retur = CubeII

def Fill_CdP((Dir, Data_Slice, masker, MESHZ, MESHX, xi, zi, o)):
    '''Interpolate discrete values at constant CdP
    To regularize the data'''
    FileName= Dir+os.sep+'CdP_Grid_'+str(o)+'.npy'
    if not os.path.isfile(FileName):
        CdP_mask = np.zeros((6250, len(masker)), dtype=np.int)
        CdP_mask[:,np.where(masker[:,o]>0.1)] +=1
        Z, X = np.where(CdP_mask>0.1)
        Coords = np.vstack((zi[Z],xi[X])).T
        CdP=SC.griddata(Coords,Data_Slice[Z, X], (MESHZ, MESHX), method='linear', fill_value=0).T
        np.save(FileName, CdP)
    return

def Fill_MO((Dir, Data_Slice, masker, Meshz, Meshy, zi, yi, s)):
    '''Interpolate discrete values at constant Offset
    To regularize the data'''
    FileName= Dir+os.sep+'MO_Grid_'+str(s)+'.npy'
    if not os.path.isfile(FileName):
        MO_mask = np.zeros((6250, len(masker[0])))
        MO_mask[:, np.where(masker[s,:]>0.1)] +=1
        Z, Y = np.where(MO_mask>0.1)
        Coords = np.vstack((zi[Z],yi[Y])).T
        MO=SC.griddata(Coords,Data_Slice[Z, Y], (Meshz, Meshy), method='linear', fill_value=0).T
        np.save(FileName, MO)
    return

def SaveCDP((i, CdP)):
    np.savetxt('CDP_R_W8ted/New_CDP_'+str(i)+'.dat', CdP.T, delimiter=',')
    return
########
def Smooth_Data(xi, yi, zi, Dir, Cube):
    main()
    MESHZ, MESHX = np.meshgrid(zi, xi)
    DirCdP = Dir[:-3]+'CdP_filt'
    if __name__ == '__main__':
        pool = Pool(12)
        pool.map(Fill_CdP, [(DirCdP, Cube[:,:,o], masker, MESHZ, MESHX, xi, zi, o) for o in range(len(yi))])
    pool.close()
    pool.join()

####### Where interpolation failed, put the original data##########
#####interpolation fails when there is no gap between traces #####
##################################################################
    for i in range(len(yi)):
        FileName = DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy'
        if os.path.isfile(FileName):
            CubeI[:,:,i] = np.load(FileName)
        else:
            CubeI[:,:,i] = Cube[:,:,i]

######## Filter Offset_wise #######
####################################
    DirMO = Dir[:-3]+'MO_filt'
    Meshz, Meshy = np.meshgrid(zi, yi)
    if __name__ == '__main__':
        pool2 = Pool(12)
        pool2.map(Fill_MO, [(DirMO, CubeI[:,s,:], masker, Meshz, Meshy, zi, yi,s) for s in range(len(xi))])
    pool2.close()
    pool2.join()

    CubeII = np.zeros(shape=(6250, len(xi), len(yi)))
    for i in range(len(xi)):
        FileName = DirMO+os.sep+'MO_Grid_'+str(s)+'.npy'
        if os.path.isfile(FileName):
            CubeII[:,i,:] = np.load(FileName)
        else:
            CubeII[:,i,:] = CubeI[:,i,:]
################## Weighting ##################
    CubeII += Cube
    #np.save(CubeII)
################# writing CdP ##################
    if __name__ == '__main__':
        pool3 = Pool(12)
        pool3.map(SaveCDP, [(i, CubeII[:,:,i].T) for i in range(len(yi))])
    pool3.close()
    pool3.join()

    return CubeII
