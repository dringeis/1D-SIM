import numpy as np  
import pylab as plt
import xarray as xr
from datetime import datetime

from matplotlib.colors import LogNorm,SymLogNorm
from matplotlib import ticker

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

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

# # plotchoice lists
# pcsl=[ [[3,0,0,0],[3,0,1,0],[3,0,2,0]],
#        [[0,0,0,0],[0,0,1,0],[0,0,2,0]],
#        [[1,0,0,0],[1,0,1,0],[1,0,2,0]],
#        [[4,0,0,0],[4,0,1,0],[4,0,2,0]],
#      ]

# plotloop(datas,pcsl)

# plt.show()

def loadxr_coords(datas,dt=600):
    
    data_pc=loaddata([0,0,0,0],datas)
    x=np.linspace(0,2000,len(data_pc))

    tsps= np.array([int(i) for i in datas['tps']])

    coords = dict(
        timesteps =  (["timestep"], tsps),
        time      = (["time"], tsps*dt) ,
        x         =  (["x"], x)
    ) 
    return coords

def loadxr1D(datas,dt=600):

    coords = loadxr_coords(datas,dt=600)

    d_vars=dict()

    # create coords
    for var in datas['vrs']: 
        for run in datas['ri']:
            for tpst in datas['tps']:
                for px in datas['px']:
                    fname=str('../output/'+var+run+tpst+px)
                    data=np.loadtxt(fname)
                    if var in d_vars:
                        d_vars[var] = np.append(d_vars[var],[data],axis=0)
                    else:
                        d_vars[var] = np.array([data])

    data_vars=dict()
    for var in d_vars:
        data_vars[var]=(['time','x'],d_vars[var])

    
    attrs = {'creation_date':datetime.now(), 
             'author':'Damien Ringeisen', 
             'email':'damien.ringeisen@mcgill.ca'}

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs=attrs,
    )

    return ds

#available data in ../output/
datas={
    # 'vrs': ['Er','u'],
    'vrs': ['A','h','div'],
    'ri' : ['_00600s_010km_solv2_IMEX0_adv3_BDF20_ts'],
    'tps': ['000001','000100','001440'],
    'px' : ['.01'],
}

ds=loadxr1D(datas)
print(ds)

# ds["Er"][0].plot(label=str("t="+str(ds['time'][0].values)+" s"))
# ds["Er"][1].plot(label=str("t="+str(ds['time'][1].values)+" s"))
# ds["Er"][2].plot(label=str("t="+str(ds['time'][2].values)+" s"))

ds["div"][0].plot(label=str("t="+str(ds['time'][0].values)+" s"))
ds["div"][1].plot(label=str("t="+str(ds['time'][1].values)+" s"))
ds["div"][2].plot(label=str("t="+str(ds['time'][2].values)+" s"))

plt.legend()
plt.show()



