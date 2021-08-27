import numpy as np
import pylab as plt
import xarray as xr
from datetime import datetime

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

def fname_ex(pc,datas):
    fname = str(datas['vrs'][pc[0]] + '_' + datas['ri'] + '_ts' + datas['tps'][pc[1]] + datas['px'])
    return datas['fdr']+fname

def loaddata(pc, datas):
    fname = fname_ex(pc, datas)
    data = np.loadtxt(fname)
    return data

def plotdatapc(data_pc, ax, x):
    ax.plot(x, data_pc)
    return None

def mainplot(datas, pcs):
    fig, ax = plt.subplots()
    for i in range(len(pcs)):
        data_pc = loaddata(pcs[i], datas)
        x = np.linspace(0, 2000, len(data_pc))
        plotdatapc(data_pc, ax, x)

    ax.set_xlabel('x [km]')
    ax.set_ylabel(str(datas['vrs'][pcs[0][0]]))
    return None

def plotloop(datas, pcsl):
    for j in range(len(pcsl)):
        mainplot(datas, pcsl[j])
    return None

# # plotchoice lists
# pcsl=[ [[3,0,0,0],[3,0,1,0],[3,0,2,0]],
#        [[0,0,0,0],[0,0,1,0],[0,0,2,0]],
#        [[1,0,0,0],[1,0,1,0],[1,0,2,0]],
#        [[4,0,0,0],[4,0,1,0],[4,0,2,0]],
#      ]

# plotloop(datas,pcsl)

# plt.show()

def loadxr_coords(datas):

    data_pc = loaddata([0, 0], datas)
    xc = np.linspace(0, 2e6, len(data_pc))
    xg = (xc[1:] + xc[:-1]) / 2

    tsps = np.array([int(i) for i in datas['tps']])

    coords = dict(
        timesteps=(["timestep"], tsps),
        time=(["time"], tsps*datas['dt']),
        xc= (["xc"], xc),
        xg= (["xg"], xg)
    )
    return coords

def saveNC(ds,name='', path=''):
    if name == '':
        name = ds.attrs['entry']['ri']+'.nc'
    if path == '':
        path = ds.attrs['entry']['fdr']
        if path[-1] != '/':
            path = path + '/'

    ds.attrs['entry'] = str(ds.attrs['entry'])

    ds.to_netcdf(path=str(path+name))

    return None

def loadxr1D(datas, attrs={}, save=False):

    coords = loadxr_coords(datas)

    d_vars = dict()

    # create coords
    run = datas['ri']
    px = datas['px']
    for var in datas['vrs']:
        for tpst in datas['tps']:
            fname = str('../output/'+var+'_'+run+'_ts'+tpst+px)
            data = np.loadtxt(fname)
            if var in d_vars:
                d_vars[var] = np.append(d_vars[var], [data], axis=0)
            else:
                d_vars[var] = np.array([data])

    data_vars = dict()
    for var in d_vars:
        data_vars[var] = (['time', datas['vx'][var]], d_vars[var])

    attrs['entry'] = datas

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs=attrs,
    )

    if save:
        saveNC(ds)

    return ds

#available data in ../output/
datas = {
    'fdr': '../output/',
    'vrs': ['A', 'h', 'div', 'Er', 'u'],
    'vx': {'A': 'xc', 'h': 'xc', 'div': 'xc', 'Er': 'xg', 'u': 'xg'},
    'ri': '00600s_010km_solv2_IMEX0_adv3_BDF20',
    'tps': ['000001', '000100', '001440'],
    'px': '.01',
    'dt': 600

}

attrs = {'creation_date': datetime.now().strftime("%m/%d/%Y, %H:%M:%S"),
         'author': 'Damien Ringeisen',
         'email': 'damien.ringeisen@mcgill.ca'}

ds = loadxr1D(datas, attrs=attrs, save=True)



# ds["Er"][0].plot(label=str("t="+str(ds['time'][0].values)+" s"))
# ds["Er"][1].plot(label=str("t="+str(ds['time'][1].values)+" s"))
# ds["Er"][2].plot(label=str("t="+str(ds['time'][2].values)+" s"))

# ds["A"][0].plot(label=str("t="+str(ds['time'][0].values)+" s"))
# ds["A"][1].plot(label=str("t="+str(ds['time'][1].values)+" s"))
# ds["A"][2].plot(label=str("t="+str(ds['time'][2].values)+" s"))

# plt.legend()
# plt.show()
