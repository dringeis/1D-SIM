import numpy as np
import pylab as plt
import xarray as xr
from datetime import datetime

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

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

def loadxr_coords(d_vars):

    # Looking at datas for dimensions
    nx = 0
    for j in d_vars:
        if d_vars[j].shape[1] > nx:
            nx = d_vars[j].shape[1]

    # Cell centers coordinates
    xc = np.linspace(0, 2e6, nx)

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
        for tpst in datas['tps']:
            fname = str('../output/'+var+'_'+run+'_ts'+tpst+exp)
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
    coords = loadxr_coords(d_vars)

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
        plt.figure()
        for i in range(int(tslim[0]*ds[var].shape[0]), int(tslim[1]*ds[var].shape[0]), step):
            ds[var][i].plot(linestyle='dashed', label=str("t="+str(ds['time'][i].values) + " s"))
        plt.legend()
        plt.grid()
        plt.title(str('variable ' + var))

    plt.show()

    return None

##########################

if __name__ == "__main__":

    # Exemple of use

    #data to load in an xarray
    datas = {
        'fdr': '../output/',
        'vrs': ['A', 'h', 'div', 'Er', 'u'],
        'vx': {'A': 'xc', 'h': 'xc', 'div': 'xc', 'Er': 'xg', 'u': 'xg'},
        'ri': '00600s_010km_solv2_IMEX0_adv3_BDF20',
        'tps': ['000001', '000100', '001440'],
        'exp': '.01',
        'dt': 600
    }

    # attribute if the NC array is to be shared (Optional)
    attrs = {'creation_date': datetime.now().strftime("%m/%d/%Y, %H:%M:%S"),
             'author': 'Damien Ringeisen',
             'email': 'damien.ringeisen@mcgill.ca'}

    # Loading all datas specified in "datas" as a xarray dataset
    ds = loadxr_SIM1D(datas, attrs=attrs, save=True)

    # Plotting for one dataset
    varn = ["u", "h", "A", "div", "Er"]
    plot_SIM1D(ds, varn)