import numpy as np
import pylab as plt
import xarray as xr
from datetime import datetime
import matplotlib as mpl

# Using LaTeX in figures
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans')

def loadxr_sNC(ds,name='', path=''):

    if name == '':
        name = ds.attrs['entry']['ri'] + ds.attrs['entry']['exp'] + '.nc'
    if path == '':
        path = ds.attrs['entry']['fdr']
        if path[-1] != '/':
            path = path + '/'

    ds.attrs['entry'] = str(ds.attrs['entry'])

    ds.to_netcdf(path=str(path+name))

    return None

def loadxr_coords(d_vars, datas):

    # Looking at datas for dimensions
    nx = 0
    for j in d_vars:
        if d_vars[j].shape[1] > nx:
            nx = d_vars[j].shape[1]

    # Cell centers coordinates
    xc = np.linspace(0, 260, nx)

    # Cell Boundaries coordinates
    xg = (xc[1:] + xc[:-1]) / 2

    # Time coordinates
    tsps = np.array([int(i) for i in datas['tps']])

    # Create the coordinates dictionnary
    coords = dict(
        timesteps=(["timestep"], tsps),
        time=(["time"], tsps*datas['dt']),
        xc= (["xc"], xc),
        xg= (["xg"], xg)
    )

    return coords

def loadxr_data(datas):

    # Create a dictionnary with all the variable at all timesteps
    d_vars = dict()
    run = datas['ri']
    exp = datas['exp']
    for var in datas['vrs']:
        print(var)
        for tpst in datas['tps']:
            fname = str(datas['fdr']+var+'_'+run+'_ts'+tpst+exp)
            data = np.loadtxt(fname)
            if var in d_vars:
                d_vars[var] = np.append(d_vars[var], [data], axis=0)
            else:
                d_vars[var] = np.array([data])

    return d_vars

def loadxr_SIM1D(datas, attrs={}, save=False):
    # Extract the data from the files to a dictionnary
    d_vars = loadxr_data(datas)

    # Extract the coordinates for the data dictionnary
    coords = loadxr_coords(d_vars,datas)

    #Create the data array dictionnary with the coordinates
    data_vars = dict()
    for var in d_vars:
        data_vars[var] = (['time', datas['vx'][var]], d_vars[var])

    # Add import dictionnary to the attributes
    attrs['entry'] = datas

    # Create the xarray dataset
    ds = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs=attrs,
    )

    # Save if asked
    if save:
        loadxr_sNC(ds)

    #return the dataset
    return ds

def plot_SIM1D(ds, varn, tslim=[0, 1], step=1):
    for var in varn:
        if var in ds:

            nbrmax=int(tslim[1]*ds[var].shape[0])
            nbrmin=int(tslim[0]*ds[var].shape[0])
            lin_nbrs = range(nbrmin, nbrmax, step)
            lin_nbrs_cb = range(nbrmin, nbrmax, step)
            n_lines = len(lin_nbrs)

            norm = mpl.colors.Normalize(vmin=nbrmin, vmax=nbrmax)
            cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)

            plt.figure()
            for i in lin_nbrs:
                ds[var][i].plot(linestyle='dashed', c=cmap.to_rgba(i + 1))
            plt.grid()
            plt.colorbar(cmap, label='Timesteps')
            plt.title(str('variable ' + var))
        else:
            print(var,' is not in the dataset')
    return None

def plot_SIM1D_PM(ds, varn):
    for var in varn:
        if var in ds:
            plt.figure()
            ds[var].plot()
            plt.title(str('variable ' + var))
        else:
            print(var,' is not in the dataset')
    return None

def timesteps(tb,te,step):
    tps = range(tb, te, step)
    tps = list(map(str, tps))
    tps = [str(item).zfill(6) for item in tps]
    return tps

##########################

if __name__ == "__main__":

    # Exemple of use

    # Creation of the timestep array
    tps = timesteps(10, 1800, 10)

    #data to load in an xarray
    datas = {
        'fdr': '../output.04/',
        'vrs': ['A', 'h', 'div', 'Er', 'u'],
        'vx': {'A': 'xc', 'h': 'xc', 'div': 'xc', 'Er': 'xg', 'u': 'xg'},
        'ri': '00010s_001km_solv1_IMEX0_adv3_BDF20',
        'tps': tps,
        'exp': '.04',
        'dt': 10
    }

    # attribute if the NC array is to be shared (Optional)
    attrs = {'creation_date': datetime.now().strftime("%m/%d/%Y, %H:%M:%S"),
             'author': 'Damien Ringeisen',
             'email': 'damien.ringeisen@mcgill.ca'}

    # Loading all datas specified in "datas" as a xarray dataset
    ds = loadxr_SIM1D(datas, attrs=attrs, save=True)
    print(ds)

    # Plotting for one dataset
    varn = ["u", "h", "A", "div", "Er"]

    # lines plots
    plot_SIM1D(ds, varn, step=10)

    # pcolormesh plots
    plot_SIM1D_PM(ds, varn)

    # ds['A'].plot(vmin=0., cmap='Blues_r')
    # ds['div'].plot(cmap='RdBu_r')

    plt.show()
