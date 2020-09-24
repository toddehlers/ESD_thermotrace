'''
This file contains all the needed custom functions to run the ESD_thermotrace simulations
Last edited by A. Madella on 17th September 2020
'''

import numpy as np                                 # library for arrays
import pandas as pd                                # library for tables
import geopandas as gpd                            # library for georeferenced tables
from collections import OrderedDict                # ordered dictionary objects
import matplotlib.pyplot as plt                    # plotting library
import matplotlib.gridspec as gs                   # library to make gridded subplots
from matplotlib.pyplot import cm                   # colour management liibrary
import seaborn as sns                              # pretty statistical plotting library
import scipy.interpolate as intr                   # interpolation functions
from sklearn.linear_model import LinearRegression  # linear regression function
import utm                                         # conversion from and to UTM coordinates
from os import mkdir, path                         # operating system utilities

class DEM:
    '''
    class with all needed DEM attributes
    '''
    def __init__(self, name):
        self.name = name

    def from_ascii(self, filename):
        '''
        Function to read an ASCII raster, most commonly exported from ArcGIS or QGIS.
        The first 6 lines must inform:

        ncols
        nrows
        xllcorner
        yllcorner
        cellsize
        NODATA_value

        followed by the raster values separated by spaces

        if verbose=False, metadata is not printed, but only returned
        '''

        # read DEM text file
        fid = open(filename, 'r')

        # make a table of the dem info and convert values to suitable data types (integer, float)
        dem_info = [fid.readline().split() for i in range(6)]
        dem_info = np.array(dem_info).transpose()
        dem_info = dict([(k,v) for k,v in zip(dem_info[0],dem_info[1])])
        dtypes = (int,int,float,float,float,float)
        for i,f in zip(dem_info,dtypes):
            dem_info[i]=f(dem_info[i])

        # get the dem data as a list of strings
        dem_ls_of_str = [fid.readline().split() for i in range(dem_info['nrows'])]

        # then convert all strings to floats
        dem = np.array([[float(i) for i in dem_ls_of_str[j]] for j in range(dem_info['nrows'])])
        if dem.shape != (dem_info['nrows'], dem_info['ncols']):
            print('something went wrong while parsing the DEM, nrows and/or ncols do not match the original input')

        # change NODATA_value to np.nan, unless it equals 0
        if dem_info['NODATA_value'] != 0:
            dem[dem==dem_info['NODATA_value']]=np.nan
            dem_info['NODATA_value']=np.nan

        self.z = dem
        self.info_dict = dem_info
        # specify the figure's geographical extent in lat,lon
        self.xllcorner = dem_info['xllcorner']
        self.yllcorner = dem_info['yllcorner']
        self.ncols = dem_info['ncols']
        self.nrows = dem_info['nrows']
        self.cellsize = dem_info['cellsize']
        self.nodata = dem_info['NODATA_value']
        self.extent84 = (self.xllcorner, self.xllcorner+self.ncols*self.cellsize,
                         self.yllcorner, self.yllcorner+self.nrows*self.cellsize)

        # build coordinate grids and arrays
        # convert llcorner and urcorner coordinates to utm and define extentUTM
        self.xyll = utm.from_latlon(self.extent84[2], self.extent84[0]) #force_zone_number=19
        self.xyur = utm.from_latlon(self.extent84[3], self.extent84[1]) #force_zone_number=19
        self.extentUTM = (self.xyll[0], self.xyur[0], self.xyll[1], self.xyur[1])

        # make easting and northing vectors
        Xi = np.linspace(self.xyll[0], self.xyur[0], self.ncols)
        Yi = np.linspace(self.xyll[1], self.xyur[1], self.nrows)
        self.xi, yi = np.meshgrid(Xi,Yi)
        self.yi = yi[::-1] ################ flipped row order for latitude to decrease from top
        # 1d vectors, needed for linear interpolation
        self.xi_1d = self.xi.reshape(dem.size)
        self.yi_1d = self.yi.reshape(dem.size)
        self.zi_1d = dem.reshape(dem.size)

    def info(self):
        '''
        Prints details of imported DEM, except nodata value
        '''
        print('')
        print('METADATA OF '+self.name)
        print('')
        print('xllcorner = {}'.format(self.xllcorner))
        print('yllcorner = {}'.format(self.yllcorner))
        print('ncols = {}'.format(self.ncols))
        print('nrows = {}'.format(self.nrows))
        print('cellsize [km] = {}'.format(self.cellsize))
        print('cellsize [m] = ~{}'.format(int(np.around(self.cellsize*110000,1)/10)*10))
        print('min value = {}'.format(np.nanmin(self.z)))
        print('max value = {}'.format(np.nanmax(self.z)))
        print('NODATA_value = {}'.format(self.nodata))

    def resample(self, res, xyll=None, xyur=None, extent84=None, method='nearest'):
        '''
        Method to resample the input rasters to desired resolution .
        By default it uses a nearest neighbour algorithm
        res = resolution in meters
        xyll = (x,y) coordinates of lower left corner of target extent of the resampled dem
        xyur = (x,y) coordinates of upper right corner of target extent of the resampled dem
        extent84 = (left,right,bottom,top) in lon-lat, of the target extent of the resampled dem
        method --> https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
        '''
        if xyll == None:
            xyll=self.xyll
        if xyur == None:
            xyur=self.xyur
        if extent84 == None:
            extent84=self.extent84
        # first make resampled easting and northing vectors
        Xi_res = np.arange(xyll[0], xyur[0]+res, res)
        Yi_res = np.arange(xyll[1], xyur[1]+res, res)
        self.xi_res, yi_res = np.meshgrid(Xi_res, Yi_res)
        self.yi_res = yi_res[::-1] ################# flipped row order for latitude to decrease from top
        self.xi_res_1d = self.xi_res.reshape(self.xi_res.size)
        self.yi_res_1d = self.yi_res.reshape(self.xi_res.size)

        # also make lon-lat vectors at same resolution
        Xi_res84 = np.linspace(extent84[0], extent84[1], self.xi_res.shape[1])
        Yi_res84 = np.linspace(extent84[2], extent84[3], self.xi_res.shape[0])
        self.xi_res84, yi_res84 = np.meshgrid(Xi_res84,Yi_res84)
        self.yi_res84 = yi_res84[::-1] ############## flipped row order for latitude to decrease from top

        # resample by interpolating at new grid nodes with resolution "res"
        # input_coords are organized in a 2D array, with columns representing x,y
        input_coords = np.concatenate(([self.xi_1d],[self.yi_1d])).transpose()
        # resampled coords are organized in a 2D array, with columns representing x,y
        self.res_coords = np.concatenate(([self.xi_res_1d],[self.yi_res_1d])).transpose()
        # Now resample: the 'values' variable refers to the known elevations of the input dem
        self.zi_res_1d = intr.griddata(points=input_coords, values=self.zi_1d, xi=self.res_coords, method=method)
        self.zi_res = self.zi_res_1d.reshape(self.xi_res.shape) # reshape from 1D to 2D


def extrapolation(gdop, gdopx, gdopy, gdopz, data, datax, datay, dataz, ext_rad):
    '''
    Extrapolates data within wanted radius
    using an inverse distance weighted linear regression of the available data points
    input arguments:
    gdop: griddata output,
          a 1D-array that contains the interpolated data
          as well as the nans that you want to replace
    gdopx
    gdopy
    gdopz: 1D coordinate arrays of griddata output

    data: 1D-array of known ages

    datax
    datay
    dataz: 1D coordinate arrays of known data points

    ext_rad: maximum extrapolation radius in meters
    '''

    # select nans from the griddata output
    nans = gdop[gdop!=gdop]
    nansx = gdopx[gdop!=gdop]
    nansy = gdopy[gdop!=gdop]
    nansz = gdopz[gdop!=gdop]

    # define distance function
    def dist3D(xyz1, xyz2):
        '''
        Calculates the distance between two points in 3D.
        xyz1 - list or tuple of x,y,z coords for first point
        xyz2 - list or tuple of x,y,z coords for second point
        '''
        return np.sqrt((xyz1[0]-xyz2[0])**2+(xyz1[1]-xyz2[1])**2+(xyz1[2]-xyz2[2])**2)

    # This is the workflow of the extrapolation function:
    # for each of the nans:
    # calculate inverse distance from NaN to all samples, drop samples too far away
    # multiply inverse distances by related age and store in a [1 x M] vector of weighted values
    # summate and divide by M

    for i in np.arange(nans.size):
        # make array of ages divided by distance and number of data points
        dists = np.array([dist3D((nansx[i],nansy[i],nansz[i]),
                                 (datax[j],datay[j],dataz[j])) for j in np.arange(data.size)])
        dists1 = dists[dists < ext_rad] # do not consider points farther than extra_rad
        if dists1.size > 0:
            data1 = data[dists < ext_rad] # select related ages
            dataz1 = dataz[dists < ext_rad]
            dists1 = 1/dists1 # invert distances
            # make linear regression based on age-elevation from data points within extra_rad
            f_z = LinearRegression().fit(dataz1.reshape((-1,1)), data1, sample_weight=dists1)
            nans[i] = f_z.intercept_+f_z.coef_*nansz[i]
        else:
            nans[i] = np.nan

    # now substitute nans with extrapolated values
    gdop[gdop!=gdop] = nans
    return gdop

def clip_to_ws(raster, shp_filename, extent, input_folder, output_folder):
    '''
    This function clips the input raster to the watershed shapefile,
    so that all raster cells can be used to predict detrital distributions
    - raster: 2D np.array
    - shp_filename: string, filename of the watershed shapefile
    - extent: tuple or list, extent of the raster in wgs1984 (west, east, south, north)
    '''
    import fiona
    import rasterio
    from rasterio.plot import show
    from rasterio.mask import mask
    from shapely.geometry import Polygon

    # calculate x and y cellsize in degrees
    xsize = np.abs(extent[0]-extent[1])/raster.shape[1]
    ysize = np.abs(extent[2]-extent[3])/raster.shape[0]
    # define rasterio transform function
    transform = rasterio.transform.from_origin(extent[0], extent[3], xsize, ysize)
    # define coordinate reference system to wgs1984
    crs = rasterio.crs.CRS.from_epsg(4326) # wgs1984: 4326
    # make new raster file based on input, necessary step to use rasterio's functions,
    # and define the metadata
    src = rasterio.open(output_folder+'/raster.tif', 'w', driver='GTiff',
                        height = raster.shape[0], width = raster.shape[1],
                        count = 1, dtype = str(raster.dtype),
                        crs = crs, transform=transform)
    # write and close the new tif file
    src.write(raster, 1)
    src.close()
    # get watershed polygon vertices
    with fiona.open(input_folder+'/'+shp_filename,'r') as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    # read the raster and make masking information
    with rasterio.open(output_folder+'/raster.tif','r') as src:
        out_image, out_transform = mask(src, shapes, nodata=np.nan)
        out_meta = src.meta
    # update the metadata accordingly
    out_meta.update({"driver": "GTiff", "height": out_image.shape[1],
                     "width": out_image.shape[2], "transform": out_transform})
    # write the clipped raster
    with rasterio.open(output_folder+'/raster_clipped.tif', 'w', **out_meta) as dest:
        dest.write(out_image)

    return rasterio.open(output_folder+'/raster_clipped.tif').read(1)

def make_dists(pop):
    '''
    makes a normalized cdf of input population
    returns a dataframe with values, counts, frequency
    - pop: 1D np.array of ages, or list of ages
    '''
    pop = np.array(pop) # make array
    df = pd.DataFrame() # initialise dataframe
    df['vals'], df['valcount'] = np.unique(pop, return_counts=True)
    df.index = df.vals.apply(int) # use ages as index
    df['cdf_y'] = df.valcount.cumsum()/df.valcount.cumsum().iloc[-1] # make cumulative frequency
    return df

def make_spdf(pop,sd=0.05):
    '''
    spdf (synoptic probability density function) equation
    from Ruhl & Hodges 2003, Vermeesch 2007, Riebe et al. 2015

    pop = 1d-array of grain ages
    sd = list/array of sigma1, or scalar informing preassigned 1sigma (e.g. 0.05 = 5% uncertainty)

    returns 1d-array where each element informs the frequency of the respective sorted value from
    '''
    pop = np.array(pop)
    if np.isscalar(sd):
        sigma1 = np.array([t*sd for t in pop])
    else:
        sigma1 = np.array(sd)
    vals = np.unique(pop, return_counts=False)
    spdf_y = [1.0/pop.size*np.sum(np.array([np.sqrt(2.0*np.pi)*u*np.exp((t-ti)**2.0/(2.0*u**2.0))
                                            for ti,u in zip(pop,sigma1)])**-1) for t in vals]
    return np.array(spdf_y)

def get_KS(res_p, d):
    '''
    res_p = resampled population p, with n=k
    d = real distribution of p
    returns KS statistic
    '''
    res_p = np.array(res_p)
    res_p.sort()
    res_p_vals, res_p_valcount = np.unique(res_p, return_counts=True)

    res_d = res_p_valcount.cumsum()/res_p_valcount.sum()
    # return KS statistic between real and interpolated from downsampled
    return max(np.abs(np.interp(x=res_p_vals,xp=d.vals,fp=d.cdf_y)-res_d))

def get_Kui(res_p, d):
    '''
    function to calculate the statistic of Kuiper's test,
    i.e., given a CDF and a number n of observations x_i-n, Kuiper's statistic expresses the sum of the two maxima
    of the vectors CDF-CDF(x_i-n) and CDF(x_i-n)-CDF
    that inform vertical distances both above and below the reference CDF

    res_p = resampled population p, with n=k
    d = real distribution of p
    
    returns Kuiper test's statistic (D)
    '''
    res_p = np.array(res_p)
    res_p.sort()
    Dplus = np.max((np.arange(res_p.size)+1)/res_p.size-np.interp(x=res_p,xp=d.vals,fp=d.cdf_y))
    Dmin = np.max(np.interp(x=res_p,xp=d.vals,fp=d.cdf_y)-np.arange(res_p.size)/res_p.size)
    return Dplus + Dmin

def get_KL(res_p, d):
    '''
    res_p = resampled population p, with n=k
    d = real distribution of p
    bounds = limits of the distribution, defined at the beginning
    returns characteristic KL divergence at the chosen confidence level
    '''
    from scipy.stats import entropy
    res_p = np.array(res_p)
    res_p.sort()
    res_p_vals, res_p_valcount = np.unique(res_p, return_counts=True)
    res_d = res_p_valcount.cumsum()/res_p_valcount.sum()
    # return KL divergence between real and interpolated from downsampled
    return entropy(np.interp(x=res_p_vals,xp=d.vals,fp=d.cdf_y), res_d)

def smooth(y, window_size=3, order=1, deriv=0, rate=1):
    '''
    dirty smoothing function, so that curves can be calculated on the base of ~1000 iterations,
    but still look good, if wanted
    '''
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        print("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y1 = np.concatenate((firstvals, y, lastvals))
    y_new = np.convolve( m[::-1], y1, mode='valid')
    y_new[0] = y[0]
    return y_new

def plot_distributions(pops_dict, dists_dict, ref_scen, detr_labels, saveas, noise_size=50, show_DKW=False, confidence=0.95):
    '''
    plotting function to display all distributions
    pops_dict: dictionary of all simulated populations
    dists_dict: dictionary of distribution dataframes, having "vals" and "cdf_y" columns
    ref_scen: label of the reference scenario
    detr_labels: list of labels to recognize measured detrital destributions
    saveas: path to save to
    noise_size: sample size for which you want to generate noise around the reference scenario
    show_DKW: boolean, if True the DKW confidence interval is calculated and plotted as alternative to the iterations
    confidence: scalar, confidence interval for the DKW
    '''
    fig = plt.figure(figsize=(14,14))
    gspec = gs.GridSpec(2,1,figure=fig)
    ax1 = fig.add_subplot(gspec[0])

    # plot all scenarios spdf
    color=iter(cm.rainbow(np.linspace(0,1,len(dists_dict))))
    xlim = (pops_dict[ref_scen].min(), pops_dict[ref_scen].max())

    for scen,df in dists_dict.items():
        if sum([i==scen for i in detr_labels])>0:
            ls='--'
        else:
            ls='-'
        sns.kdeplot(pops_dict[scen], color=next(color), linewidth=4, linestyle=ls, label='_nolegend_', ax=ax1)
        ax1.set(xlim=xlim, ylabel='frequency', xlabel='', xticks=[])
        ax1.set_title('Modeled vs Detrital: Kernel Density Estimation',
                     fontdict={'weight':'bold'})

    # second row
    ax2 = fig.add_subplot(gspec[1])
    color=iter(cm.rainbow(np.linspace(0,1,len(dists_dict))))
    c_ref = next(color)
    # plot all scenarios cdf
    if not show_DKW:
        for i in np.arange(100): # include 100 random subsamples of reference scenario?
            pop1 = np.random.choice(pops_dict[ref_scen],noise_size)
            if i==0:
                dist1 = make_dists(pop1)
                sns.lineplot(x=dist1.vals, y=dist1.cdf_y,
                             color=c_ref, alpha=0.1, lw=1, ax=ax2, label=ref_scen+', n='+str(noise_size))
            else:
                dist1 = make_dists(pop1)
                sns.lineplot(x=dist1.vals, y=dist1.cdf_y,
                             color=c_ref, alpha=0.1, lw=1, ax=ax2, label='_nolegend_')

    color=iter(cm.rainbow(np.linspace(0,1,len(dists_dict))))
    for scen,df in dists_dict.items():
        if sum([i==scen for i in detr_labels])>0:
            ls = '--'
        else:
            ls = '-'
        if scen == ref_scen and show_DKW:
            DKW = np.sqrt(np.log(2/(1-confidence))/(2*noise_size))
            ax2.fill_between(x=df.vals,
            y1=df.cdf_y-DKW,
            y2=df.cdf_y+DKW,
            color=c_ref,
            alpha=0.3,
            label=str(int(confidence*100))+'% conf of '+ref_scen+', n='+str(noise_size)
            )
        ax2.plot(df.vals, df.cdf_y, color=next(color), linestyle=ls, label=scen+', n='+str(df.valcount.sum()), lw=4)
    ax2.set(xlim=xlim, ylabel='cumulative frequency', xlabel='age [My]')
    ax2.set_title('Modeled vs Detrital: Cumulative Age Distribution',
                     fontdict={'weight':'bold'})
    ax2.legend()
    # save fig
    fig.savefig(saveas, dpi=200)

def plot_confidence(prob_dict, all_k, ref_scen, saveas, num_of_colors):
    '''
    plots the confidence level as function of sample size
    prob_dict: dictionary of arrays, each containing confidence values related to all_k
    all_k: 1D-array of k values (integers)
    ref_scen: the label of the reference scenario, based on which prob_dict has been made
    saveas: path to be saved at
    num_of_colorss: number of colors to iterate through,
                    must equal the number of colors plotted before, to be consistent
    '''
    sns.set_style('white')
    fig,ax = plt.subplots(figsize=(12,8))
    # plot grey fields
    ax.fill_between([all_k[0],all_k[-1]],[95,95],[68,68],color='k',alpha=0.2)
    ax.text(all_k[0]+0.5,96,'95%',fontdict=dict(size=20))
    ax.fill_between([all_k[0],all_k[-1]],[68,68],color='k',alpha=0.4)
    ax.text(all_k[0]+0.5,69,'68%',fontdict=dict(size=20))
    color=iter(cm.rainbow(np.linspace(0,1,num_of_colors)))
    next(color) # skip one color for Euni
    leg=[]
    for key,i in prob_dict.items():
        c = next(color)
        ax.plot(all_k,
                smooth(i),
                c=c,
                alpha=1)
        leg.append(key)
    ax.set(xlim=(all_k[0],all_k[-1]), ylim=(0,101), yticks=[0,20,40,60,80,100])
    ax.set_xlabel('number of grains',fontdict=dict(size=20, weight='bold'))
    ax.set_ylabel('confidence level [%]',fontdict=dict(size=20, weight='bold'))
    ax.set_title('Confidence of discerning from "'+ref_scen+'" as function of sample size',
                 fontdict=dict(size=20,weight='bold'), pad=10)
    ax.legend(leg, loc='lower right')
    fig.savefig(saveas, dpi=200)

def plot_violins(data, label, column, saveas, k_iter, sam_size):
    '''
    plotting function for violins figure
    data: pd.dataframe with column 'scenario' (for the x axis) and another column to be plotted as y
    label: string, the key of the wanted population for the dictionary of detrital data
    column: string, column label of the dataframe to plot
    saveas: path to save to
    k_iter: number of iterations
    sam_size: number of grains used in each iteration
    '''
    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot()
    # also useful: sns.color_palette('colorblind')
    sns.violinplot(data=data, x='scenario', y=column, split=True, hue='metric',
    ax=ax, cut=0, scale='area', palette='rainbow')
    ax.set_xlabel('erosion scenario', fontdict={'weight':'bold'})
    ax.set_ylabel('dissimilarity from '+str(k_iter)+' iterations',fontdict={'weight':'bold'})
    ax.set_title('Difference between subsampled erosion scenarios and '+label+', n='+str(sam_size),
                 pad=10, fontdict={'weight':'bold'})
    # save figure
    fig.savefig(saveas, dpi=200)
