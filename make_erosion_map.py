def make_erosion_map(xi, yi, points, kind='inter'):
    '''
    function to generate a synthetic erosion map, based on Gaussian peaks
    xi,yi - coordinate grids, like outputs of numpy.meshgrid()
    points - list of 3-elements-lists (x,y,w) or array with shape (n,3), 1 row per point.
            If kind='gauss' w = ratio of max/min values within the map
            If kind='inter' w = erosion weigth of point
    kind - 'gauss': makes a gaussian surface for each point, with specified max/min ratio of erosiivity
         - 'inter': linearly interpolates a surface through the points
    '''
    import numpy as np

    if kind == 'gauss':
        from scipy.optimize import fsolve, minimize
        def gau2d(xi, yi, m, s):
            '''
            makes gaussian surface with peak at desired location
            xi,yi - coordinate grids
            m - (x,y) coordinates of peak
            s - standard deviation
            '''
            return 1/(s*np.sqrt(2*np.pi))*np.exp(-1/2*(np.sqrt((xi-m[0])**2+(yi-m[1])**2)/s)**2)


        def make_gau_Eratio(xi, yi, m, Eratio):
            '''
            create a normal distribution of erosion for a scenario with given ratio Emax/Emin 'Eratio',
            given location of erosional focus 'm',
            given positions array 'pos'
            '''

            def eq_zero(sd,Eratio,xi,yi,m):
                '''
                this function provides the zero condition to solve for, in order to obtain the wanted erosion ratio
                sd - standard deviation
                Eratio = ratio between max and min erosion within catchment
                xi,yi = coordinate grids of surface
                m = (x,y) coordinates of peak
                '''
                gau1 = gau2d(xi, yi, m, sd)
                gau1 = gau1/gau1.min() # minimum = 1
                Eratio1 = gau1.max()/gau1.min()
                return (Eratio-Eratio1)**2

            #sd1 = fsolve(eq_zero,100,(Eratio,xi,yi,m))
            sd1 = minimize(eq_zero,np.array(50),(Eratio,xi,yi,m), bounds=[(0,xi.max()-xi.min())], tol=1e-2).x
            gau = gau2d(xi, yi, m, sd1)
            return gau/gau.min()

        gau_list = [make_gau_Eratio(xi, yi, p[:2], p[2]) for p in points] # make all gaussian peaks
        gau_final = np.ones(xi.shape)
        for g in gau_list:
            index = [g>gau_final]
            gau_final[tuple(index)] = g[tuple(index)]
        return gau_final/gau_final.min()

    elif kind=='inter':
        import scipy.interpolate as inter

        # known data points, for which the interpolation function is found
        # they are organized in a 2D array, with columns representing x,y,z
        pts = np.array(points)[:,:2]

        # positions where interpolation needs to be made (all the catchment's gridcells)
        # they are organized in a 2D array, with columns representing x,y,z
        pos = np.concatenate(([xi.reshape(xi.size)],[yi.reshape(yi.size)])).transpose()

        # known weigths, 3rd element of each point
        vals = np.array(points)[:,2]
        return inter.griddata(points=pts, values=vals, xi=pos, method='linear').reshape(xi.shape)
    else:
        print('I do not know what you mean by kind="{}"'.format(kind))
