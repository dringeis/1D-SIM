import numpy as np  
import pylab as plt

from matplotlib.colors import LogNorm,SymLogNorm
from matplotlib import ticker

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

#available data in ../output/
datas={
    'vrs': ['A','div','Er','h','u'],
    'ri' : ['_00600s_010km_solv1_IMEX0_adv3_BDF20_ts','_00600s_010km_solv2_IMEX0_adv3_BDF20_ts'],
    'tps': ['000001','000100','001440'],
    'px' : ['.01'],
}

def fname_ex(pc,datas):
    fname=str(datas['vrs'][pc[0]]+datas['ri'][pc[1]]+datas['tps'][pc[2]]+datas['px'][pc[3]])
    return '../output/'+fname

def loaddata(pc,datas):
    fname=fname_ex(pc,datas)
    data=np.loadtxt(fname)
    return data

def plotdatapc(data_pc,ax,x,lab):
    ax.plot(x,data_pc)
    return None 

def mainplot(datas,pcs):
    
    fig, ax = plt.subplots()
    
    for i in range(len(pcs)):
        data_pc=loaddata(pcs[i],datas)
        x=np.linspace(0,2000,len(data_pc))
        plotdatapc(data_pc,ax,x,lab)

    ax.set_xlabel('x [km]')
    ax.set_ylabel(str(datas['vrs'][pcs[0][0]]))
    
    return None

def plotloop(datas,pcsl):
    for j in range(len(pcsl)):
        mainplot(datas,pcsl[j])
    return None

# plotchoice lists
pcsl=[ [[3,0,0,0],[3,0,1,0],[3,0,2,0]],
       [[0,0,0,0],[0,0,1,0],[0,0,2,0]],
       [[1,0,0,0],[1,0,1,0],[1,0,2,0]],
       [[4,0,0,0],[4,0,1,0],[4,0,2,0]],
     ]

plotloop(datas,pcsl)

plt.show()

