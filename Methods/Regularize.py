from multiprocessing import Pool
import os, sys, m8r
import numpy as np
import scipy.interpolate as SC
import matplotlib.pyplot as plt

def LoadRsf (Dir, File):
    F = m8r.Input(Dir+os.sep+File)
    Data = np.array(F)
    return Data

def Fill_CdP((Dir, Data_Slice, masker, MESHZ, MESHX, xi, zi, o, n1)):
    '''Interpolate discrete values at constant CdP
    To regularize the data'''
    FileName= Dir+os.sep+'CdP_Grid_'+str(o)+'.npy'
    if not os.path.isfile(FileName):
        if Data_Slice.var() >0.005:
            print 'Processing '+FileName
            CdP_mask = np.zeros((n1, len(masker)), dtype=np.int)
            CdP_mask[:,np.where(masker[:,o]>0.1)] +=1
            Z, X = np.where(CdP_mask>0.1)
            Coords = np.vstack((zi[Z],xi[X])).T
            CdP=SC.griddata(Coords,Data_Slice[Z, X], (MESHZ, MESHX), method='linear', fill_value=0, rescale=True).T
            np.save(FileName, CdP)
    return

def Fill_MO((Dir, Data_Slice, masker, Meshz, Meshy, zi, yi, s, n1)):
    '''Interpolate discrete values at constant Offset
    To regularize the data'''
    FileName= Dir+os.sep+'MO_Grid_'+str(s)+'.npy'
    if not os.path.isfile(FileName):
        if Data_Slice.var()>0.005:
            print 'Processing '+FileName
            MO_mask = np.zeros((n1, len(masker[0])))
            MO_mask[:, np.where(masker[s,:]>0.1)] +=1
            Z, Y = np.where(MO_mask>0.1)
            Coords = np.vstack((zi[Z],yi[Y])).T
            MO=SC.griddata(Coords,Data_Slice[Z, Y], (Meshz, Meshy), method='linear', fill_value=0, rescale=True).T
            np.save(FileName, MO)
    return

def SaveCDP((i, CdP)):
    np.savetxt('CDP_R_W8ted/New_CDP_'+str(i)+'.dat', CdP.T, delimiter=',')
    return
#
# def Normalise_Stacks(Cube, Stacks):
#     '''normalise amplitude where Stacks of traces in the same bin'''
#     Stacks.sort(key= lambda x: x[0])
#     n= 0
#     while n < len(Stacks):
#         S = 2
#         m=n
#         while Stacks[m] == Stacks[n]:
#             S +=1
#             m +=1
#             if m>len(Stacks):
#                 break
#         Cube[:,Stacks[n][0], Stacks[n][1]] /= S
#         n += S-1
    ###### normalise amplitude where Stacks ########
    ######## of traces in the same bin   ###########
    #
    # return Cube

def Cubize(Dir, Filt_Line, n1):
    Params= np.genfromtxt(Dir+os.sep+'Utils/Parameters.txt', usecols=1, dtype=np.int)
    MinOff, MaxOff, Low_bound, High_bound = Params
    Offsets, Pointer = np.load(Dir+os.sep+'Utils/Offsets.npy').astype(np.int), np.load(Dir+os.sep+'Utils/Pointer.npy').astype(np.int)
    xi = np.arange(MinOff, MaxOff+1)
    yi = np.arange(Low_bound, High_bound)
    zi = np.arange(n1)
    Filtered_line= LoadRsf (Dir,Filt_Line)
    Stacks =[]
    Cube = np.zeros(shape=(n1, len(xi), len(yi)))
    masker = np.zeros(shape=(len(xi), len(yi)), dtype=np.int)
    for i in range(len(Filtered_line)):
        o, c= Offsets[i]-MinOff, Pointer[i]-Low_bound
        if Cube[1000,o, c] != 0:
            Stacks.append((o, c))
        masker[o, c] += 1
        Cube[:,o, c] += Filtered_line[i]
    plt.imshow(masker)
    plt.show()
    ##### Normalize stacked traces
    for i in range(len(masker)):
        for j in range(len(masker[0])):
            if masker[i,j] >1:
                Cube[:,i,j] /=masker[i,j]
    np.save(Dir+os.sep+'Utils/masker.npy', masker)
    np.save(Dir+os.sep+'Utils/Cube.npy', Cube)
    # Cube=np.load(Dir+os.sep+'Utils/Cube.npy')
    # masker=np.load(Dir+os.sep+'Utils/masker.npy')
    Smooth_Data(xi, yi, zi, Dir, Cube, masker, n1)
    return

def Smooth_Data(xi, yi, zi, Dir, Cube, masker, n1):
    ######## Workaround SCons.compat module renaming ########
    import imp
    del sys.modules['pickle']
    del sys.modules['cPickle']
    sys.modules['pickle'] = imp.load_module('pickle', *imp.find_module('pickle'))
    sys.modules['cPickle'] = imp.load_module('cPickle', *imp.find_module('cPickle'))
    import pickle
    import cPickle
    #print "(pickle == cPickle) = ", (pickle == cPickle)

    # MESHZ, MESHX = np.meshgrid(zi, xi)
    # DirCdP = Dir+'/CdP_filt'
    # if __name__ :#== '__main__':
    #     pool = Pool(12)
    #     pool.map(Fill_CdP, [(DirCdP, Cube[:,:,o], masker, MESHZ, MESHX, xi, zi, o, n1) for o in range(len(yi))])
    # pool.close()
    # pool.join()
    # pool.terminate()
    # print 'Done interpolating in Offset direction'
    # ####### Where interpolation failed, put the original data##########
    # #####interpolation fails when there is no gap between traces #####
    # ##################################################################
    # CubeI = np.zeros(shape=(n1, len(xi), len(yi)))
    # for i in range(len(yi)):
    #     FileName = DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy'
    #     if os.path.isfile(FileName):
    #         Slice= np.load(FileName)
    #         if np.shape(Slice) != (n1, len(xi)):
    #             Slice= np.transpose(Slice)
    #         CubeI[:,:,i] = Slice
    #     else:
    #         CubeI[:,:,i] = Cube[:,:,i]
    #
    # ######## Filter Offset_wise #######
    # ####################################
    # DirMO = Dir+'/MO_filt'
    # Meshz, Meshy = np.meshgrid(zi, yi)
    # if __name__ :# == '__main__':
    #     pool2 = Pool(4)
    #     pool2.map(Fill_MO, [(DirMO, CubeI[:,s,:], masker, Meshz, Meshy, zi, yi,s, n1) for s in range(len(xi))])
    # pool2.close()
    # pool2.join()
    # pool2.terminate()
    # print 'Done interpolating in CdP direction'
    # CubeII = np.zeros(shape=(n1, len(xi), len(yi)))
    # for i in range(len(xi)):
    #     FileName = DirMO+os.sep+'MO_Grid_'+str(i)+'.npy'
    #     if os.path.isfile(FileName):
    #         CubeII[:,i,:] = np.load(FileName)
    #     else:
    #         CubeII[:,i,:] = CubeI[:,i,:]
    #
    # ################## Weighting ##################
    # CubeII += Cube
    #np.save(CubeII)
    ################# writing CdP ##################
    if __name__:# == '__main__':
        pool3 = Pool(18)
        pool3.map(SaveCDP, [(i, Cube[:,:,i]) for i in range(len(yi))])
    pool3.close()
    pool3.join()
    pool3.terminate()
    return
######Traces observed to be problematic during QC##########
###########################################################
# ManDeadTr = np.array([[20,39],[20,43],[21,41],[21,110],[23,35],[28,10],[70,48],[81,8],[159,67],[175,26],[187,115],[187,116],[188,111],[188,112],[201,38],[284,55],[285,51],[435,112],[435,116],[435,119],[436,44],[437,108],[437,115],[438,111],[443,55],[443,91],[444,51]])
# for dt in range(len(ManDeadTr)):
#     t= ManDeadTr[0]*120+ManDeadTr[1]
#     masker[Offsets[i], Pointer[i]-Low_bound] = 0 ### ignores all stack in case large noise
