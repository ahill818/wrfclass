import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
from matplotlib.widgets import Slider, Button, RadioButtons
from pylab import cm
import scipy.ndimage as ndimage

class WRFobject:
    """ Pass in WRF filepath to plot 2d and 3d data fields """

    def __init__(self,filepath):
        self.filepath = filepath
        self.wrf = nc.Dataset(filepath)
        self.lons = self.wrf.variables['XLONG'][0]
        self.lats = self.wrf.variables['XLAT'][0]
        self.lat0 = self.wrf.CEN_LAT
        self.lon0 = self.wrf.CEN_LON
        self.lat1 = self.wrf.TRUELAT1
        self.lat2 = self.wrf.TRUELAT2
        self.gp_width = len(self.wrf.dimensions['west_east'])
        self.gp_height = len(self.wrf.dimensions['south_north'])
        self.dx = self.wrf.DX
        self.dy = self.wrf.DY
        self.high = self.dy*(self.gp_height-1)
        self.wide = self.dx*(self.gp_width-1)
        self.times = self.wrf.variables['T'].shape[0]
        self.m = self.make_basemap()
        self.x,self.y = self.make_xy(self.m)
#        self.save_base = '/Users/severe/Research/plots/'
#        self.datestring = filepath[-19:-15]+filepath[-14:-12]+filepath[-11:-9]+filepath[-8:-6]
        self.save_base = '/lustre/work/aarhill/enkfsrc/plots/'
        self.datestring = '20110522'

    def make_basemap(self):
        """ make basemap figure """
        return Basemap(projection='lcc',width=self.wide,height=self.high,resolution='i',
            lat_0=self.lat0,lon_0=self.lon0,lat_1=self.lat1,lat_2=self.lat2)

    def make_xy(self,m):
        return m(self.lons,self.lats)

    def grab_variable(self,var):
        """ grab variable of interest
        type: ncdump -h wrfile.nc into command line for variable options """
        return self.wrf.variables[var]

    def plot_var(self,var,time,threeD=False,lev=0):
        """ plot variables at specific time (input) and level (input) if 3D. Levels are given 
        as model levels, not in another unit (e.g. pressure, height, etc.) """
        plt.figure(figsize=[10,6])
        self.m.drawcoastlines()
        self.m.drawcountries()
        self.m.drawstates()
        if threeD:
            D=plt.contourf(self.x,self.y,self.grab_variable(var)[time,lev,:,:],cmap=matplotlib.cm.spectral)
        else:
            D=plt.contourf(self.x,self.y,self.grab_variable(var)[time,:,:],cmap=matplotlib.cm.spectral)
        plt.colorbar()
        plt.show() 

    def make_plot(self,data,title,save_string,units,t,c_range=[None],other_plots=[False,False],wind_data=[0,0],c_table=cm.jet):
    #    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
        if c_range[0] is not None:
            cmin,cmax,cstep = c_range[0],c_range[1],c_range[2]
        else:
            cmin,cmax,cstep = np.min(data),np.max(data),(np.max(data)-np.min(data))/20.
        print cmin,cmax,cstep
        print c_range
        plt.figure(figsize=[10,6])
        self.m.drawcoastlines()
        self.m.drawcountries()
        self.m.drawstates()
        cs = plt.contourf(self.x,self.y,data[0],np.arange(cmin,cmax,cstep),cmap=c_table)
        cbar = plt.colorbar(cs,orientation='vertical',shrink=0.875,pad=0.02)
        if other_plots[0]: # plot wind
            self.x_wind,self.y_wind = self.m(self.lons[::20,::20],self.lats[::20,::20])            
            plt.barbs(self.x_wind,self.y_wind,wind_data[0][::20,::20],wind_data[1][::20,::20],
                length=5,pivot='middle')
        if other_plots[1]:  # plot SLP contour
            LEVS_slp = range(990,1020,2)
            cs = plt.contour(self.x,self.y,ndimage.gaussian_filter(data[1],sigma=1.0,order=0),
                LEVS_slp,colors='k',linestyles='solid')
            plt.clabel(cs,inline=1,fontsize=10,fmt='%4.0f')
#        cbar.set_label(units,rotation=0,y=1.05,labelpad=-30.0) 
        plt.title(title)
#        plt.savefig("{0}{1}_fhr{2}_{3}.png".format(self.save_base,save_string,"%02d"%t,self.datestring),
#            bbox_inches='tight',format='PNG',dpi=300)
        
#        plt.show()

    def slp(self,t=0):
        temps = self.grab_variable('T2')
        psfc = self.grab_variable('PSFC')
        stemps = temps[t]+6.5*self.grab_variable('HGT')[t]/1000.
        mslp = psfc[t]*np.exp(9.81/(287.0*stemps)*self.grab_variable('HGT')[t])*0.01
        return mslp

    def dewpoint(self,t=0):
        temps = self.grab_variable('T2')
        psfc = self.grab_variable('PSFC')
        qhum = self.grab_variable('Q2')
        # Convert Surface Pressure to Mean Sea Level Pressure
        stemps = temps[t]+6.5*self.grab_variable('HGT')[t]/1000.
        mslp = psfc[t]*np.exp(9.81/(287.0*stemps)*self.grab_variable('HGT')[t])*0.01
        # Find saturation vapor pressure
        es = 6.112 * np.exp(17.67 * temps[t]/(temps[t] + 243.5))
        w = qhum[t]/(1-qhum[t])
        e = (w * psfc[t] / (.622 + w)) / 100
        Td_C = (243.5 * np.log(e/6.112))/(17.67-np.log(e/6.112))
        return Td_C 

    def comp_reflectivity(self,t=0):
        temps = self.grab_variable('T2')
        psfc = self.grab_variable('PSFC')
        QR = self.grab_variable('QRAIN')
        QS = self.grab_variable('QSNOW')
       # Define 'constant' densities (kg m-3)
        rhor = 1000
        rhos = 100
        rhog = 400
        rhoi = 917
        # Define "fixed intercepts" (m-4)
        Norain = 8.0E6
        #Nosnow = 2.0E7
        Nosnow = 2.0E6*np.exp(-0.12 * (temps[t]-273))
        Nograu = 4.0E6
        # First, find the density at the first sigma level
        # above the surface
        density = np.divide(psfc[t],(287.0 * temps[t]))
        #print "Rho: ", np.mean(density)
        Qra_all = QR[t]
        Qsn_all = QS[t]

        for j in range(len(Qra_all[1,:,1])):
            curcol_r = []
            curcol_s = []
            for i in range(len(Qra_all[1,1,:])):
                maxrval = np.max(Qra_all[:,j,i])
                maxsval = np.max(Qsn_all[:,j,i])
                curcol_r.append(maxrval)        
                curcol_s.append(maxsval)
            np_curcol_r = np.array(curcol_r)
            np_curcol_s = np.array(curcol_s)
            if j == 0:
                Qra = np_curcol_r
                Qsn = np_curcol_s
            else:
                Qra = np.row_stack((Qra, np_curcol_r))
                Qsn = np.row_stack((Qsn, np_curcol_s))

        # Calculate slope factor lambda
        lambr = np.divide((3.14159 * Norain * rhor), np.multiply(density, Qra))
        lambr = lambr ** 0.25
        #lambs = np.divide((3.14159 * Nosnow * rhoi), np.multiply(density, Qsn))
        #lambs = lambs ** 0.25
        lambs = np.exp(-0.0536 * (temps[t] - 273))
        # Calculate equivalent reflectivity factor
        Zer = (720.0 * Norain * (lambr ** -7.0)) * 1E18
        Zes = (0.224 * 720.0 * Nosnow * (lambr ** -7.0) * (rhos/rhoi) ** 2) * 1E18
        Zes_int = np.divide((lambs * Qsn * density), Nosnow)
        Zes = ((0.224 * 720 * 1E18) / (3.14159 * rhor) ** 2) * Zes_int ** 2 
        Ze = np.add(Zer, Zes)
        #Ze = Zer
        # Convert to dBZ
        dBZ = 10 * np.log10(Ze) 
        dBZ = np.nan_to_num(dBZ)
        return dBZ

    def sim_reflectivity(self,t=0):
        temps = self.grab_variable('T2')
        psfc = self.grab_variable('PSFC')
        QR = self.grab_variable('QRAIN')
        QS = self.grab_variable('QSNOW')

        # Define 'constant' densities (kg m-3)
        rhor = 1000
        rhos = 100
        rhog = 400
        rhoi = 917
        # Define "fixed intercepts" (m-4)
        Norain = 8.0E6
        #Nosnow = 2.0E7
        Nosnow = 2.0E6*np.exp(-0.12 * (temps[t]-273))
        Nograu = 4.0E6
        # First, find the density at the first sigma level
        # above the surface
        density = np.divide(psfc[t],(287.0 * temps[t]))
        #print "Rho: ", np.mean(density)
        Qra = QR[t,1]
        Qsn = QS[t,1]
        Qra = np.nan_to_num(Qra)
        Qsn = np.nan_to_num(Qsn)
        # Calculate slope factor lambda
        lambr = np.divide((3.14159 * Norain * rhor), np.multiply(density, Qra))
        lambr = lambr ** 0.25
        #lambs = np.divide((3.14159 * Nosnow * rhoi), np.multiply(density, Qsn))
        #lambs = lambs ** 0.25
        lambs = np.exp(-0.0536 * (temps[t] - 273))
        # Calculate equivalent reflectivity factor
        Zer = (720.0 * Norain * (lambr ** -7.0)) * 1E18
        Zes = (0.224 * 720.0 * Nosnow * (lambr ** -7.0) * (rhos/rhoi) ** 2) * 1E18
        Zes_int = np.divide((lambs * Qsn * density), Nosnow)
        Zes = ((0.224 * 720 * 1E18) / (3.14159 * rhor) ** 2) * Zes_int ** 2 
        Ze = np.add(Zer, Zes)
        #Ze = Zer
        # Convert to dBZ
        dBZ = 10 * np.log10(Ze) 
        dBZ = np.nan_to_num(dBZ)
        return dBZ 

    def plot_cref(self,time):
        cref = self.comp_reflectivity(t=time)
        title='test'
        save_string = 'cref'
        units = 'dBZ'
        self.make_plot([cref],title,save_string,units,time,c_range=[-30,70,5])

    def plot_sref(self,time):
        sref = self.sim_reflectivity(t=time)
        title='test'
        save_string = 'sref'
        units = 'dBZ'
        self.make_plot([sref],title,save_string,units,time,c_range=[-30,70,5])    

    def plot_dewpoint(self,time):
        dew = self.dewpoint(t=time)
        title = 'test'
        save_string = "td2"
        units = 'C'
        self.make_plot([dew],title,save_string,units,time)
 

    def plot_surfacevars(self,t):
        """ Plot the mean surface variables, 2-meter T, 2-meter Td, 
            10-meter U and V wind, surface pressure """
        print ('Plotting surface variables at fhr:%s' % (t))
        t2 = self.grab_variable('T2')[t,:,:]-273.15  # in C
        slp = self.slp(t=t)
        u10 = self.grab_variable('U10')
        v10 = self.grab_variable('V10')
        u10knots = u10[t,:,:] * 1.94384
        v10knots = v10[t,:,:] * 1.94384
        title='Mean 2-meter T (C), MSLP (mb), 10m Wind (knots)\n'
        units='C'
        savetitle = 'meansurfvars'
        self.make_plot([t2,slp],title,savetitle,units,t,
            other_plots=[True,True],wind_data=[u10knots,v10knots])

    def plot_temperature(t):
        """ plot the mean 2-meter temperature """
        print ('Plotting 2-meter temp at fhr:%s' % (t))
        t2 = self.grab_variable('T2')[t,:,:]-273.15
        title='2-meter Temp (F)'
        units='C'
        savetitle = 't2'
        self.make_plot([t2],title,savetitle,units,t)

    def plot_mdbz_td2(self,time):
        print "Plotting reflectivity and dewpoint overlayed at fhr:{0}".format(time)
        mdbz = self.comp_reflectivity(t=time)
        td2 = self.dewpoint(t=time)
        dbz_levs = np.arange(-30,70,3)
        td2_levs = np.arange(-10,30,2)
#        D = plt.contourf(x,y,mdbz,dbz_levs,cmap=coltbls.reflect())
#        C = plt.contour(x,y,td2,td2_levs,colors='k',linewidths=0.5)
        units='dBZ'
        title='MDBZ and TD2 test'
        save_string = 'mdbztd2'
        plt.contour(self.x,self.y,td2,td2_levs,colors='k',linewidths=0.5) 
        self.make_plot([mdbz],title,save_string,units,time,c_range=[-30,70,5])



class EnsMemDiff(WRFobject):
    """ Pass in Ensemble Forecasts of members to plot differences. Forecast files have also been 
    interpolated to levels using p_interp """

    def __init__(self,filebase,conv_list=[1],nonconv_list=[None]):
        """ input lists of convecting and non-convecting members along with a file base string """
        self.c_list = conv_list
        self.nc_list = nonconv_list
        self.file_base = filebase
        self.wrf = WrfClass(self.file_base+str(conv_list[0]))
        self.m = self.wrf.make_basemap()
        self.x,self.y = self.wrf.make_xy(self.m)

    def calc_diff(self,time,function):
        c_data = np.zeros_like(self.wrf.lons)
        nc_data = c_data.copy()
        print function
        for mem in self.c_list:
            print mem
            wrf = WrfClass(self.file_base+str(mem))
            c_data += wrf.slp(t=time)
            print c_data[0:10,0:10]
            wrf.wrf.close()

        for mem2 in self.nc_list:
            print mem2
            wrf = WrfClass(self.file_base+str(mem2))
            nc_data += wrf.slp(t=time)
            print nc_data[0:10,0:10]
            wrf.wrf.close()

        c_data = c_data / len(self.c_list)
        nc_data = nc_data / len(self.nc_list)
        return c_data,nc_data


    def plot_slpdiff(self,time):
        print "Slp Differences Between members at fhr:{0}".format(time)
        slp_c, slp_nc = self.calc_diff(time,self.slp)
#        LEVS = np.arange(-2.,2.,0.1)  
        title='SLP Difference (Convecting - Non-convecting) at fhr%s (mb)'%(time)
        units='mb'
        savetitle = 'slpdiff'
        print np.max(slp_c), np.max(slp_nc)
        print np.min(slp_c), np.min(slp_nc)
        self.make_plot([slp_c-slp_nc],title,savetitle,units,time,c_range=[-1.,1.,0.02],c_table=cm.RdBu_r) 


# def plot_td2diff(t):
#     print("   Td2 Diff")
#     nc1 = netcdf.netcdf_file(diff1)
#     nc2 = netcdf.netcdf_file(diff2)  
#     td_c1 = dewpoint(nc1,t)
#     td_c2 = dewpoint(nc2,t)

#     LEVS = np.arange(-5.,5.,0.2)
#     # Contour and fill the dewpoint temperature     
#     Td=plt.contourf(x,y,td_c1-td_c2,LEVS,cmap=cm.RdBu_r,extend='both')
# #   plt.clabel(Td_lev,inline=1,fontsize=7,fmt='%1.0f',inline_spacing=1)
#     title='Diff'
#     units='C'
#     savetitle = 'difftest'
#     make_plot(Td,title,savetitle,units,t,False) 

# def plot_t2diff(t):
#     print(" Difference in T2")
#     f1 = netcdf.netcdf_file(diff1)
#     f2 = netcdf.netcdf_file(diff2)
#     t2diff = f1.variables['T2'][t,:,:]-f2.variables['T2'][t,:,:]
#     LEVS = np.arange(-5.0,5.0,0.2)
#     D = plt.contourf(x,y,t2diff,LEVS,cmap=cm.RdBu_r,extend='both')
# #    C = plt.contour(x,y,t2diff,LEVS,colors='k',linewidths=.5)
#     units='C'
#     title='T2 Difference (Convecting - Non-Convecting) at fhr%s (%s)'%(t,units)
#     savetitle = 't2diffmems'
#     make_plot(D,title,savetitle,units,t,False) 

# def plot_rainncdiff(t):
#     print(" Difference in Hourly Rainfall")
#     f1 = netcdf.netcdf_file(diff1)
#     f2 = netcdf.netcdf_file(diff2)
#     f1b = netcdf.netcdf_file(diff1)
#     f2b = netcdf.netcdf_file(diff2)
#     fhourly = f1.variables['RAINNC'][t,:,:]-f1b.variables['RAINNC'][t-1,:,:]
#     fhourly2 = f2.variables['RAINNC'][t,:,:]-f2b.variables['RAINNC'][t-1,:,:]
#     LEVS = np.arange(-3.0,3.0,0.1)
#     D = plt.contourf(x,y,fhourly-fhourly2,LEVS,cmap=cm.RdBu_r,extend='both')
#     units='mm'
#     title='Rainfall Difference (Convecting - Non-Convecting) at fhr%s (%s)'%(t,units)
#     savetitle = 'rainncdiffmems'
#     make_plot(D,title,savetitle,units,t,False) 
    

# def plot_mdbzdiff(t):
#     print("Difference in MDBZ")
#     conv = netcdf.netcdf_file(diff1)
#     noconv = netcdf.netcdf_file(diff2)
#     conv_mdbz = comp_reflectivity(conv,t)
#     noconv_mdbz = comp_reflectivity(noconv,t)
#     LEVS = np.arange(-40,40,4.0)
#     diff = conv_mdbz-noconv_mdbz
#     D = plt.contourf(x,y,diff,LEVS,cmap=cm.RdBu_r,extend='both')
#     units = 'dBZ'
#     title = 'test MDBZ difference between 36 and 10'
#     savetitle = 'mdbzdiff'
#     make_plot(D,title,savetitle,units,t,False) 

class Interp(WRFobject):
    """ Takes in a file that has interpolated temp, height, etc. from a WRF file. Inherits the WRF class """
    def __init__(self,interp_path):
        self.interp = nc.Dataset(interp_path)

    def grab_interp(self,t,var,level):
        # var == 'TT','GHT'
        levs = {700:4,
                500:5
        }
        return self.grab_variable(var)[t,levs[level],:,:]

    def plot_ua_var(self,t,var,level):
        data = self.grab_interp(t,var,level)

