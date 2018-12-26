import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as SC

dx, dy= 10, 10
data = np.genfromtxt('/media/julien/NuDrive/Himalayas/Tomo.txt', names=True)
### capture topo
TopX, TopY,TopV =[data['Distance'][0]],[data['Z'][0]],[data['VP'][0]]
i = 1
while i < len(data):
    if data['Distance'][i] ==  data['Distance'][i-1]:
        i+=1
    else:
        if  data['VP'][i] <0:
            i+=1
        else:
            TopX.append(data['Distance'][i])
            TopY.append(data['Z'][i])
            TopV.append(data['VP'][i])
            i+=1
# ####Compute distance stretch
# TopX, TopY,TopV =np.array(TopX),np.array(TopY),np.array(TopV)
# d= ((TopX[1:]-TopX[:-1])**2+(TopY[1:]-TopY[:-1])**2)**.5
Vo= 1794 #mean Static Velocity

#Map Velocity on CMP
def closest_node(pt, line):
    dist_2 = np.sum((line-pt)**2, axis=1)
    return np.argmin(dist_2)




# minX, maxX = data['Distance'].min(), data['Distance'].max()
# minY, maxY = data['Z'].min(), data['Z'].max()
# minV, maxV = data['VP'].min(), data['VP'].max()
# xi= np.arange(minX-dx, maxX+dx, dx)
# yi= np.arange(round(minY-dy,dy), round(maxY+dy), dy)
# MeshX, MeshY =np.meshgrid(xi, yi)
# Grid= SC.griddata((data['Distance'],data['Z']), data['VP'], (MeshX, MeshY), method='cubic', fill_value=0)
# fig, ax = plt.subplots(1,1)
# ax.imshow(Grid, aspect='auto', vmin=minV, vmax=4500, extent=[minX, maxX,maxY,minY])
# ax.plot(TopX,TopY, 'r')
# ax.invert_yaxis()
# plt.show()