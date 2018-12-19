import sys, os, m8r
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
import scipy.interpolate as SC
from multiprocessing import Pool

def Interpolate(X, Y, Z, Bin, deg, S):
    ''' Interpolate a 3D spline according to distance between points'''
    #find NaN
    Fin = np.isfinite(X)
    X, Y, Z = X[:][Fin], Y[:][Fin], Z[:][Fin]
    Ds = np.zeros(len(X), dtype =np.float)
    # sort by distance from first point
    Ds[1:] += (X[1:]-X[0])**2
    Ds[1:] += (Y[1:]-Y[0])**2
    Ds[1:] += (Z[1:]-Z[0])**2
    Ds = Ds**0.5
    Ind = np.argsort(Ds)
    # Reinitialise coordinates add padding with linear interpolation
    Xs, Ys, Zs = np.zeros(len(X)+2), np.zeros(len(X)+2), np.zeros(len(X)+2)
    Xs[1:-1], Ys[1:-1], Zs[1:-1] = X[Ind], Y[Ind], Z[Ind]
    n = 4
    Xs[0], Ys[0], Zs[0]    = n*X[0]-(n-1)*X[1],   n*Y[0]-(n-1)*Y[1],   n*Z[0]-(n-1)*Z[1]
    Xs[-1], Ys[-1], Zs[-1] = n*X[-1]-(n-1)*X[-2], n*Y[-1]-(n-1)*Y[-2], n*Z[-1]-(n-1)*Z[-2]
    D = np.zeros(len(Xs))
    D[1:] += (Xs[1:] - Xs[:-1])**2
    D[1:] += (Ys[1:] - Ys[:-1])**2
    D[1:] += (Zs[1:] - Zs[:-1])**2
    D = D**0.5
    D = D.cumsum()
    X_Spline = SC.UnivariateSpline( D, Xs, k=deg, s = S)
    Y_Spline = SC.UnivariateSpline( D, Ys, k=deg, s = S)
    Z_Spline = SC.UnivariateSpline( D, Zs, k=deg, s = S)
    New_D = np.arange(D[0], D[-1], Bin)
    Spl_X = X_Spline(New_D)
    Spl_Y = Y_Spline(New_D)
    Spl_Z = Z_Spline(New_D)
    return X_Spline, Y_Spline, Z_Spline, Spl_X, Spl_Y, Spl_Z

def getData (Dir, File, Header):
    F = m8r.Input(Dir+os.sep+File)
    TrH = m8r.Input(Dir+os.sep+Header)
    Data = np.array(F)
    TrHead= np.array(TrH)
    return Data, TrHead

def LoadRsf (Dir, File):
    F = m8r.Input(Dir+os.sep+File)
    Data = np.array(F)
    return Data
# def WriteRsf(Array, Dir, Name, Tsample, Rspacing):
#     '''Writes numpy array and its mask into optimised rsf file'''
#     n1,n2,n3 = np.shape(Array)
#     Out      = m8r.Output(Dir+os.sep+Name)
#     Mask_Out = m8r.Output(Dir+os.sep+'Mask_'+Name)
#     # array gets transposed in rsf
#     axis = [{'n':n3,'d':Tsample,'o':0,'l':'Time','u':'s'},
#             {'n':n2,'d':Rspacing,'o':0,'l':'Offset','u':'m'},
#             {'n':n1,'d':1,'o':0,'l':'Shot','u':''}]
#     for i in range(len(axis)):
#         Out.putaxis(axis[i], i+1)
#         Mask_Out.putaxis(axis[i], i+1)
#     Out.write(Array.data)
#     Mask_Out.write(Array.mask)
    # return
def WriteRsf(Array, Dir, Name, *axis):
    '''Writes numpy array and its mask into optimised rsf file'''
    Out      = m8r.Output(Dir+os.sep+Name)
    Out.type = 'int'
    #Mask_Out = m8r.Output(Dir+os.sep+'Mask_'+Name)
    dim = np.shape(Array)[::-1] #Rsf transposed compared to numpy
    for i in range(len(axis)):
        axis[i]['n']=dim[i]
        Out.putaxis(axis[i], i+1)
        #Mask_Out.putaxis(axis[i], i+1)
    Out.write(Array)#.data
    #Mask_Out.write(Array.mask)
    return

def MakeData(L, H, CorrTr, Cube, TraceH, S, Nh, Receivers):
    n=0
    for f in range(len(L)):
        if CorrTr[f] > 0:  ### correct for missing trace assume tailing
            Buff  = np.vstack((L[f], np.zeros((CorrTr[f],S ))))
            HBuff = np.vstack((H[f], np.zeros((CorrTr[f],Nh))))
        else:
            Buff  = L[f]
            HBuff = H[f]

        for i in range(0,len(Buff),Receivers):
            idx = n + i/Receivers
            Cube  [idx] =  Buff[i:i+Receivers].T
            TraceH[idx] = HBuff[i:i+Receivers].T
        n = idx +1
    return Cube, TraceH

def QCshots(Dir, Cube):
    '''Plot all the shots'''
    print ('saving the shots for QC')
    nrows, ncols = 3,4
    tot = nrows*ncols
    FIGS= len(Cube)//tot
    print FIGS
    for j in range(FIGS+1):
        print j
        fig, axes = plt.subplots(nrows, ncols)
        for i, ax in enumerate(fig.axes):
            idx = i+tot*j
            if idx<len(Cube):
                ax.set_title('ShotPoint '+str(idx))
                ax.imshow(Cube[idx,:2000,:], vmin = -1, vmax = 1, aspect ='auto', origin='upper')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
        fig. tight_layout()
        fig.savefig(Dir+os.sep+"QC/RawShots_p"+str(j)+".pdf", dpi=200 )
        plt.close()
    return

def STD(Cube, Receivers):
    '''test traces with std'''
    STD= np.zeros((len(Cube), Receivers))
    for s in range(len(Cube)):
        SSTD = Cube[s].std()
        for tr in range(Receivers):
            STD [s,tr] = np.std(Cube[s,:, tr:tr+1])
        STD [s] -= SSTD
    STD = np.ma.masked_invalid(STD)
    return

def NullTrace(Cube, Receivers):
    '''Generate a list of traces with insignificant data in window Sample 50-2000
        With DeadTr [shots, TraceNb]'''
    DeadTr =[]
    So, Se = 100, 2000 #window to look for significant data
    for s in range(len(Cube)):
         for tr in range(Receivers):
            a = np.sum(Cube[s,So:Se, tr]**2)
            if a < 0.001:
                     DeadTr.append([s,tr])
    return np.array(DeadTr)

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpassb(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def Frequencyfilt():
    Sf = 250 #Hz 4ms
    LowCut = 30 #Hz
    HighCut = 125 #Hz

    Test = np.zeros(np.shape(Cube[90]))
    #b, a = butter_bandpass(LowCut, HighCut, Sf, order=3)
    for i in range(len(Test[0])):
        Test[:,i]= butter_bandpass_filter(Cube[90,:,i], LowCut, HighCut, Sf, order=3)

    #Wave Nb filtering
    Df = 4 #Geoph every 25m
    LowCutK =1
    HighCutK = 2
    for j in range(len(Test)):
        Test[j,:]= butter_bandpass_filter(Cube[90,j,:], LowCutK, HighCutK, Df, order=3)

    #plot test freq

    fig, axes = plt.subplots(1, 2)
    m= axes[0].imshow(Test[:2000,:], aspect='auto')
    n= axes[1].imshow(Cube[90,:2000,:], aspect='auto')
    plt.colorbar(m)

    plt.show()

def PlotCoord(X,Y,Z):
    fig, ax = plt.subplots(2,2)
    #IDX =np.argsort(YS)
    ax[0,0].plot(X,Y)
    ax[0,1].plot(X,Z)
    ax[1,0].plot(Y,Z)
    return

def closest_node(pt, line):
    dist_2 = np.sum((line-pt)**2, axis=1)
    return np.argmin(dist_2)

def Calc_CMP(TraceH, CoX, CoY):
    '''make CMPs arrays tr X, Y, Z
    Construct a pointer between shortest binned curve and each CMP, extract fold'''
    ScZ, ScX = abs(TraceH[:,19,:]), abs(TraceH[:,20,:])
    XS= TraceH[:,21,:]/ScX+CoX
    YS= TraceH[:,22,:]/ScX+CoY
    ZS= TraceH[:,13,:]/ScZ
    Xr, Yr, Zr    = TraceH[:,23,:]/ScX+CoX, TraceH[:,24,:]/ScX+CoY, TraceH[:,12,:]/ScZ
    CMPo, CMPm, CMPe = np.zeros((len(XS),4)), np.zeros((len(XS),4)),np.zeros((len(XS),4))

    for o in range(len(TraceH)):
            CMPo[o], CMPm[o], CMPe[o] = [o, (Xr[o,0]+XS[o,0])/2, (Yr[o,0]+YS[o,0])/2, (Zr[o,0]+ZS[o,0])/2], \
                                        [o, (Xr[o,60]+XS[o,60])/2, (Yr[o,60]+YS[o,60])/2, (Zr[o,60]+ZS[o,60])/2], \
                                        [o, (Xr[o,119]+XS[o,119])/2, (Yr[o,119]+YS[o,119])/2, (Zr[o,119]+ZS[o,119])/2]
            for i in range(Receivers):
                CMPs[o*Receivers+i,0], CMPs[o*Receivers+i,1], CMPs[o*Receivers+i,2], CMPs[o*Receivers+i,3] = o*Receivers+i, (Xr[o,i]+XS[o,i])/2, (Yr[o,i]+YS[o,i])/2, (Zr[o,i]+ZS[o,i])/2

    ### Prepare spline:
    xyz = np.unique(np.ma.compress_rowcols(np.ma.masked_invalid(CMPm[:,1:]), axis = 0), axis=0)
    xyz = xyz[xyz[:,1].argsort()] #sorting to compensate for weird sorting from np.unique
    ### calculate spline to make binning notice subsampling and intense smoothing
    X_Spline, Y_Spline, Z_Spline, Spl_X, Spl_Y, Spl_Z = Interpolate(xyz[:,0][::6], xyz[:,1][::6], xyz[:,2][::6], Bin=25, deg=3, S=3000)
        #PlotCoord(Spl_X, Spl_Y, Spl_Z)
    Pointer = np.ones(len(CMPs), dtype=np.int)*-1 # contains index of nearest Spline node
    Spline= np.zeros((len(Spl_X),3))
    Spline[:,0], Spline[:,1], Spline[:,2] = Spl_X, Spl_Y, Spl_Z
    for a in range(len(CMPs)):
        if np.isfinite(CMPs[a, 1]):
            Pointer[a] = closest_node(CMPs[a,1:], Spline)
    Pointer = np.ma.masked_where(Pointer<0, Pointer)
    BinP, fold = np.unique(Pointer, return_counts=True)
    ### plot fold along line
    ##plt.plot(b[:-1], label='fold') #remove nans at the end
    ##plt.legend()
    ##plt.show()
    return Xr, Yr, Zr, CMPs, Spline, Pointer, BinP, fold[:-1]

def Calc_CMP2D(TraceH, CoX, CoY, CdP_Bin):
    '''make CMPs arrays tr X, Y, Z
    Construct a pointer between shortest binned curve and each CMP, extract fold'''
    Tr=TraceH[:,0]
    ScZ, ScX = abs(TraceH[:,19]), abs(TraceH[:,20])
    XS= TraceH[:,21]/ScX+CoX
    YS= TraceH[:,22]/ScX+CoY
    ZS= TraceH[:,13]/ScZ
    Xr, Yr, Zr    = TraceH[:,23]/ScX+CoX, TraceH[:,24]/ScX+CoY, TraceH[:,12]/ScZ
    CMPs = np.zeros((len(XS),4))
    for o in range(len(TraceH)):
        CMPs[o,0], CMPs[o,1], CMPs[o,2], CMPs[o,3] = Tr[o], (Xr[o]+XS[o])/2, (Yr[o]+YS[o])/2, (Zr[o]+ZS[o])/2
    CMPm = CMPs[::Receivers/2]
    ### Prepare spline:
    xyz = np.unique(CMPm[:,1:], axis=0)
    xyz = xyz[xyz[:,1].argsort()] #sorting to compensate for weird sorting from np.unique
    ### calculate spline to make binning notice subsampling and intense smoothing
    X_Spline, Y_Spline, Z_Spline, Spl_X, Spl_Y, Spl_Z = Interpolate(xyz[:,0][::6], xyz[:,1][::6], xyz[:,2][::6], Bin=CdP_Bin, deg=3, S=3000)
    Pointer = np.ones(len(CMPs), dtype=np.int)*-1 # contains index of nearest Spline node
    Spline= np.zeros((len(Spl_X),3))
    Spline[:,0], Spline[:,1], Spline[:,2] = Spl_X, Spl_Y, Spl_Z
    for a in range(len(CMPs)):
        if np.isfinite(CMPs[a, 1]):
            Pointer[a] = closest_node(CMPs[a,1:], Spline)
    #Pointer = np.ma.masked_where(Pointer<0, Pointer)
    BinP, fold = np.unique(Pointer, return_counts=True)
    return Xr, Yr, Zr, CMPs, Spline, Pointer, BinP, fold[:-1]

def RemoveGroundRoll(Sh_Rec, V0=30, padding=1000):
    '''Remove Ground Roll in shot receiver space by rotating
    then dip-filter the low velocity component
    and trashing the low frequency, back rotating after
    Adaptated from Fomel 2002 http://ahay.org/RSF/book/sep/pwd/paper_html/paper.html'''
    nj     = (3,2)
    Shot   = m8r.File(Sh_Rec.T)
    Shot   = m8r.put(label1='Time', unit1='Seconds', d1=0.004, o1=0,label2='Receivers', o2=1, n2=120, d2 =1)[Shot]
    Tilted = m8r.pad(beg1=padding)[Shot]
    Tilted = m8r.stretch(rule='l', v0 =V0, half='n', delay=0, verb ='y')[Tilted]
    noiz   = m8r.mutter(v0=V0*2, half='n')[Tilted]
    noiz   = m8r.bandpass(fhi=10) [noiz]
    sign   = m8r.bandpass(flo=10) [Tilted]
    mask   = m8r.math(output='abs(input)')[noiz]
    mask   = m8r.mask(min=0.02)[mask]
    mask   = m8r.dd (type='float')[mask]
    dat2   = m8r.add(mode='p')[Tilted, mask]
    ndip   = m8r.twodip2(order=3, nj1=nj[0], nj2=nj[0], eps=8, verb='n', gauss='n', p0=2, q0=-2)[noiz]
    sdip   = m8r.twodip2(order=3, nj1=nj[1], nj2=nj[1], eps=8, verb='n', gauss='n', p0=-2, q0=2)[sign]
    pq     = m8r.twodip2(order=3, nj1=nj[0], nj2=nj[1], eps=8, verb='n', gauss='n', p0=-2, q0=2, dip1=ndip, dip2=sdip, mask=mask)[dat2]
    sdip2  = m8r.window(n3=1)[pq]
    ndip2  = m8r.window(f3=1)[pq]
    signoi = m8r.planesignoi(sdip=sdip2, ndip=ndip2, order =3, nj1=nj[0], nj2=nj[1], eps=2, verb='n', niter=100)[dat2]
    noiz2  = m8r.window(n3=1)[signoi]
    noiz2  = m8r.stretch(rule='l',inv='y', v0=V0, half='n', delay=0)[noiz2]
    noiz2  = m8r.window(f1=padding)[noiz2]
    signW  = m8r.window(f3=1)[signoi]
    buff   = m8r.add(scale='1,1,-1')[mask,Tilted,dat2]
    buff   = m8r.window(f3=1)[buff]
    sign2  = m8r.add(mode='p')[signW,buff]
    sign2  = m8r.stretch(rule='l',inv='y', v0=V0, half='n', delay=0)[sign2]
    sign2  = m8r.window(f1=padding)[sign2]
    lowf   = m8r.bandpass (fhi=8)[sign2]
    sign3  = m8r.add(scale='1,-1')[sign2,lowf]
    noiz3  = m8r.add(scale='1, 1')[noiz2,lowf] #assuming sfadd without instruction is just an addition...
    return

##def get_trace(ShotPoint, Trace):
##    Data, TrHead = getData (Dir, File, Header)
##    extract trace and headers
##    return trace and header
#### File management
########## Data gathering ##########
Dir ='/media/julien/NuDrive/Himalayas/dummy/SGY/INDEPTH/RSF'
Files   =  ['TIB01.rsf']#,'TIB01_B.rsf','TIB01_C.rsf','TIB01_D.rsf','TIB01_E.rsf']
Headers = Files[:]
L, H, TT    = [], [], []
for a in range(len(Files)):
    Headers[a] = Files[a][:-4]+'_T'+Files[a][-4:]
    DD, h = getData (Dir, Files[a], Headers[a])
    L.append(DD)
    H.append(h)
    TT.append(h[-1,1])

########## Data consolidation ########
Receivers= H[0][0,-13]
Nh = len(H[0][0])
S = H[0][0,38]
CorrTr=[0]
#CorrTr = [0,2,0,0,13] ##L[1] missing 2 traces L[4] missing 13 traces QC might be done directly from trace headers?
Tr = sum(TT+CorrTr)
Geom =(Tr/Receivers,S, Receivers) ## assuming first shot is correct
TraceH, Cube  = np.zeros((Geom[0], Nh, Receivers)), np.zeros(Geom)

########### Make Data Files ###########
#Cube, TraceH = MakeData(L, H, CorrTr, Cube, TraceH, S, Nh, Receivers)

########### save first QC plot ############
#QCshots(Dir, Cube)

######Creating a mask for the Cube #########
## dead shots checked 2 times for possible info (time shift and not scaled values)
# DeadShots = [0, 3, 4, 24, 25, 37, 38, 47, 48, 52, 53, 74, 75, 87, 106, 108, 111, 112,120, 121, 137, 153, 154, 171, 186, 202, 203, 220, 235, 238, 297, 298, 309, 310, 316, 317, 335, 336, 349, 350, 355, 374, 375, 396, 397, 419, 420, 431, 432]

# WeirdShots=[36,88,138,170,201,234,325, 337,338,427, 436]
#
# ### Mask for dead shots
# MASK = np.zeros(Geom, bool)
# MASK[DeadShots] = True
# MASK[WeirdShots]=True
#
# ### Mask for deadTrace
# DeadTr = NullTrace(Cube, Receivers)
# MASK[DeadTr[:,0],:,DeadTr[:,1]]=True
# MASK[ManDeadTr[:,0],:,ManDeadTr[:,1]]=True
#
# ## Masking
# Cube   = np.ma.MaskedArray(Cube, MASK)
# TraceH = np.ma.MaskedArray(TraceH, MASK[:,:len(TraceH[0]),:])

###############Export to rsf###############
###########################################
# axis = [{'d':'','o':'','l':'','u':''},
#        {'d':'','o':'','l':'Traces','u':''}]
#WriteRsf(Cube.transpose(0,2,1), Dir, 'Cube_test', Tsample=0.004, Rspacing=50) #Cube in Shot Time Offset needs ordering before export
#WriteRsf(TraceH.transpose(0,2,1), Dir, 'Trace_test.rsf', Tsample=1, Rspacing=1)

#### Source and Receiver positions in real space, get CMPs etc
CoX, CoY = 700000, 3000000  #Coordinate modifier
TraceH=H[0].astype(np.float)
Xr, Yr, Zr, CMPs, Spline, Pointer, BinP, fold = Calc_CMP2D(TraceH, CoX, CoY, CdP_Bin=100)

########### Create CMP Gather Geometry ###########
Offsets = np.zeros(len(Pointer), dtype=np.int)
for i in range(len(Pointer)):
        #o,r = i/Receivers//1, i%Receivers #pointing inside Offset cube headers
        Offsets[i] =np.around((([Xr[i], Yr[i], Zr[i]]-Spline[Pointer[i]])**2).sum()**.5/50, decimals=0)

########EXPORT HEADER TO ASCII########
# New_Head= np.zeros((3,len(Offsets)), dtype=np.int)
# New_Head[0]=H[0][:,0]
# New_Head[1]=Offsets
# New_Head[2]=Pointer
# np.savetxt('New_Head.dat',New_Head, delimiter=',')



########Interpolate Missing Data ##########
Low_bound, High_bound = BinP[0], BinP[-1]+1
MinOff, MaxOff= int(Offsets.min()),int(Offsets.max()//1)+1
xi = np.arange(MinOff, MaxOff)
yi = np.arange(Low_bound, High_bound)
zi =np.arange(6250)
Filtered_line= LoadRsf (Dir,'Masked_line.rsf')
Stacks =[]
Cube = np.zeros(shape=(6250, len(xi), len(yi)))
masker = np.zeros(shape=(len(xi), len(yi)), dtype=np.int)
for i in range(len(Filtered_line)):
    if Cube[100,Offsets[i], Pointer[i]-Low_bound] != 0:
        Stacks.append((Offsets[i], Pointer[i]-Low_bound))
    masker[Offsets[i], Pointer[i]-Low_bound] += 1
    Cube[:,Offsets[i], Pointer[i]-Low_bound] += Filtered_line[i]

######Traces observed to be problematic during QC##########
###########################################################
ManDeadTr = np.array([[20,39],[20,43],[21,41],[21,110],[23,35],[28,10],[70,48],[81,8],[159,67],[175,26],[187,115],[187,116],[188,111],[188,112],[201,38],[284,55],[285,51],[435,112],[435,116],[435,119],[436,44],[437,108],[437,115],[438,111],[443,55],[443,91],[444,51]])
for dt in range(len(ManDeadTr)):
    t= ManDeadTr[0]*120+ManDeadTr[1]
    masker[Offsets[i], Pointer[i]-Low_bound] = 0 ### ignores all stack in case large noise

###### normalise amplitude where Stacks ########
####### of 2 traces in the same bin   ##########
Stacks.sort(key= lambda x: x[0])
n= 0
while n < len(Stacks):
    S = 2
    m=n
    while Stacks[m] == Stacks[n]:
        S +=1
        m +=1
        if m>len(Stacks):
            break
    Cube[:,Stacks[n][0], Stacks[n][1]] /= S
    n += S-1

####### interpolate traces #####
###### CubeI.npy save ##########

def Fill_Slice((Data_Slice, Meshx, Meshy, x, y, t)):
    '''Interpolate discrete values at constant t
    To regularize the data'''
    Slice=SC.griddata((x,y),Data_Slice, (Meshx, Meshy), method='cubic', fill_value=0)
    np.savetxt(Dir+os.sep+'CdP_Slice_'+str(t)+'.dat', Slice, delimiter=',')
    return

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


MESHZ, MESHX = np.meshgrid(zi, xi)
DirCdP = Dir[:-3]+'CdP_filt'
if __name__ == '__main__':
    pool = Pool()
    pool.map(Fill_CdP, [(DirCdP, Cube[:,:,o], masker, MESHZ, MESHX, xi, zi, o) for o in range(len(yi))])
pool.close()
pool.join()
####### When interpolation failed, put the original data##########
#####interpolation fails when there is no gap between traces #####
##################################################################
 for i in range(551,560):
     np.save(DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy',Cube[:,:,i])
for i in range(594,600):
     np.save(DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy',Cube[:,:,i])
for i in range(611,620):
     np.save(DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy',Cube[:,:,i])
np.save(DirCdP+os.sep+'CdP_Grid_'+str(869)+'.npy',Cube[:,:,i])
CubeI = np.zeros(shape=(6250, len(xi), len(yi)))
for i in range(len(yi)):
    CubeI[:,:,i] = np.load(DirCdP+os.sep+'CdP_Grid_'+str(i)+'.npy')

DirMO = Dir[:-3]+'MO_filt'
Meshz, Meshy = np.meshgrid(zi, yi)
if __name__ == '__main__':
    pool2 = Pool()
    pool2.map(Fill_MO, [(DirMO, CubeI[:,s,:], masker, Meshz, Meshy, zi, yi,s) for s in range(10,50)])
pool2.close()
pool2.join()
CubeII = np.zeros(shape=(6250, len(xi), len(yi)))
for i in range(len(xi)):
    CubeII[:,:,i] = np.load(Dir+os.sep+'MO_Grid_'+str(i))
# ax[0].imshow(masker, vmax=1, aspect='auto', cmap='gray')
# ax[1].imshow(Cube[1000,:,:], vmin=-0.5, vmax=0.5, aspect='auto')
# ax[2].imshow(Slice.T, vmin=-0.5, vmax=0.5, aspect='auto')
# plt.show()
# for i in range(len(yi)):
#     np.savetxt(Dir[:-4]+os.sep+'CDP'+os.sep+'CdP_py_'+str(i+Low_bound)+'.dat', Cube[:,:,i], delimiter=',')
# axis = [{'d':'0.004','o':'0','l':'Traces','u':'s'}, {'d':'50','o':'0','l':'Offsets','u':'m'},     {'d':'100','o':str(Low_bound),'l':'CdP','u':'m'}]
# Out = m8r.Output(Dir+os.sep+'Cube_bined.rsf')
# dim = np.shape(Cube)[::-1]
# axis=[]#Rsf transposed compared to numpy
# for i in range(len(axis)):
#     axis[i]['n']=dim[i]
#     Out.putaxis(axis[i], i+1)
# Out.write(Cube.transpose())



















###### plot Offset vs CMP bin
##plt.scatter(Pointer, Offsets, marker ='+')
##plt.xlabel('CMP Number')
##plt.ylabel('Offset in m')
##plt.show()










######## Figure dump ########

##    ax[0].plot(Xs, Ys, 'x')
##    ax[1].plot(Xr,Yr,'-')
##    ax[2].plot(CMPx, CMPy, 'o')
##    plt.setp(ax, aspect='equal')
##    ax[0].set_title('Source')
##    ax[1].set_title('Receivers')
##    ax[2].set_title('CMP')
##plt.show()


##fig, ax = plt.subplots(2,2)
##ax[0,0].plot(CMPo[:,1],CMPo[:,2], label = '0')
##ax[0,0].plot(CMPm[:,1],CMPm[:,2], label = '60')
##ax[0,0].plot(CMPe[:,1],CMPe[:,2], label = '119')
##ax[0,0].plot(Spl_X, Spl_Y, label ='spline')
##ax[0,1].plot(CMPo[:,1],CMPo[:,3])
##ax[0,1].plot(CMPm[:,1],CMPm[:,3])
##ax[0,1].plot(CMPe[:,1],CMPe[:,3])
##ax[0,1].plot(Spl_X, Spl_Z, label ='spline')
##ax[1,1].plot(CMPo[:,3],CMPo[:,0], label = '1st CMP')
##ax[1,1].plot(CMPm[:,3],CMPm[:,0], label = 'Mid CMP')
##ax[1,1].plot(CMPe[:,3],CMPe[:,0], label = 'Last CMP')
##
##ax[1,0].scatter(CMPs[:,1],CMPs[:,2], marker ='.',label = 'CMPs')
##ax[1,0].plot(Spl_X, Spl_Y, 'r', label ='spline')
##
##ax[0,0].set_xlabel('Easting')
##ax[0,0].set_ylabel('Northing')
##ax[1,0].set_xlabel('Easting')
##ax[1,0].set_ylabel('Shot Nb')
##ax[0,1].set_xlabel('Easting')
##ax[0,1].set_ylabel('Altitude')
##ax[1,1].set_xlabel('Altitude')
##ax[1,1].set_ylabel('Shot Nb')
##plt.legend()
##plt.show()


##fig, ax = plt.subplots(1,2)
##ax[0].imshow(np.array(Shot), aspect='auto')
##ax[1].imshow(np.array(test), aspect='auto')
##plt.show()
##ax[0].imshow(Cube[115,:2000,:].data, extent =[0, 6000, 0, 12000])
##ax[1].imshow(Cube[116,:2000,:].data, extent =[0, 6000, 0, 12000])
##plt.show()
### mayavi plot
#obj=mlab.volume_slice(Cube[:,:2000,:],vmin=-1, vmax=1,plane_orientation='y_axes')
### use scalar cut plane in mayavi pipe
