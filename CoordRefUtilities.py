import numpy as np

#---------------------------------------------------------------------
# x,y = ll2psn(lat, lon)
# 
# Takes in:
# lat - scalar or 1D array of latitudes
# lon - scalar or 1D array of longitudes
# 
# Returns:
# x - scalar or 1D array of equivalent x coordinates in EPSG3413
# y - scalar or 1D array of equivalent y coordinates in EPSG3413
#----------------------------------------------------------------------
def ll2psn(lat,lon):
    phi_c = 70.0
    a = 6378137.0
    e = 0.08181919
    lambda_0 = -45.0
    #Convert to radians:    
    phi=lat*np.pi/180
    phi_c=phi_c*np.pi/180
    lambda_1=lon*np.pi/180
    lambda_0=lambda_0*np.pi/180
    #this is not commented very well. See Snyder for details.
    t=np.tan(np.pi/4-phi/2)/np.power((1-e*np.sin(phi))/(1+e*np.sin(phi)),e/2)
    t_c=np.tan(np.pi/4 - phi_c/2)/np.power((1-e*np.sin(phi_c))/(1+e*np.sin(phi_c)),e/2)
    m_c=np.cos(phi_c)/np.sqrt(1-np.square(e)*np.square(np.sin(phi_c)))
    rho=a*m_c*t/t_c #true scale at lat phi_c
    x=rho*np.sin(lambda_1-lambda_0)
    y=-rho*np.cos(lambda_1 - lambda_0)
    return x,y



#---------------------------------------------------------------------
# lat, lon = psn2ll(x,y)
# 
# Takes in:
# x - scalar or 1D array of equivalent x coordinates in EPSG3413
# y - scalar or 1D array of equivalent y coordinates in EPSG3413
# 
# Returns:
# lat - scalar or 1D array of latitudes
# lon - scalar or 1D array of longitudes
#----------------------------------------------------------------------
def psn2ll(x,y):
    phi_c = 70.0
    a = 6378137.0
    e = 0.08181919
    lambda_0 = -45.0
    #Convert to radians:    
    phi_c=phi_c*np.pi/180
    lambda_0=lambda_0*np.pi/180
    
    #this is not commented very well. See Snyder for details.
    t_c=np.tan(np.pi/4 - phi_c/2)/np.power((1-e*np.sin(phi_c))/(1+e*np.sin(phi_c)),e/2)
    m_c=np.cos(phi_c)/np.sqrt(1-np.square(e)*np.square(np.sin(phi_c)))
    rho=np.sqrt(x**2+y**2)
    t=rho*t_c/(a*m_c)

    #find phi with a series instead of iterating.
    chi=np.pi/2 - 2*np.arctan(t)
    phi=(chi+(e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360)*np.sin(2*chi)
    + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520)*np.sin(4*chi)
    + (7*e**6/120+81*e**8/1120)*np.sin(6*chi)
    + (4279*e**8/161280)*np.sin(8*chi))

    lambda_2=lambda_0 + np.arctan2(x,-y)

    #correct the signs and phasing
    lambda_2=np.mod(lambda_2+np.pi,2*np.pi)-np.pi    #want longitude in the range -pi to pi

    #convert back to degrees
    lat=phi*180/np.pi
    lon=lambda_2*180/np.pi
    return lat,lon
    

#---------------------------------------------------------------------
# x,y = ll2ps(lat, lon)
# 
# Takes in:
# lat - scalar or 1D array of latitudes
# lon - scalar or 1D array of longitudes
# 
# Returns:
# x - scalar or 1D array of equivalent x coordinates in EPSG3031
# y - scalar or 1D array of equivalent y coordinates in EPSG3031
#----------------------------------------------------------------------
def ll2ps(lat,lon):
    phi_c = -71.     # standard parallel - this is different from Andy Bliss' function, which uses -70! 
    a = 6378137.0    # radius of ellipsoid, WGS84
    e = 0.08181919   # eccentricity, WGS84
    lambda_0 = 0.    # meridian along positive Y axis
    lat = lat*np.pi/180
    # Convert to radians
    phi_c = phi_c*np.pi/180
    lon = lon*np.pi/180;
    lambda_0 = lambda_0*np.pi/180
    #this is not commented very well. See Snyder for details.
    t=np.tan(np.pi/4+lat/2)/np.power((1-e*np.sin(-lat))/(1+e*np.sin(-lat)),e/2)
    t_c=np.tan(np.pi/4 + phi_c/2)/np.power((1-e*np.sin(-phi_c))/(1+e*np.sin(-phi_c)),e/2)
    m_c=np.cos(-phi_c)/np.sqrt(1-np.square(e)*np.square(np.sin(-phi_c)))
    rho=a*m_c*t/t_c #true scale at lat phi_c
    x=-rho*np.sin(-lon+lambda_0)
    y=rho*np.cos(-lon + lambda_0)
    return x,y
 

#---------------------------------------------------------------------
# x,y = ps2ll(lat,lon)
# 
# Takes in:
# x - scalar or 1D array of equivalent x coordinates in EPSG3031
# y - scalar or 1D array of equivalent y coordinates in EPSG3031
# 
# Returns:
# lat - scalar or 1D array of latitudes
# lon - scalar or 1D array of longitudes
#----------------------------------------------------------------------
def ps2ll(x,y):
    phi_c = -71.      #standard parallel - this is different from Andy Bliss' function, which uses -70! 
    a = 6378137.0     #radius of ellipsoid, WGS84
    e = 0.08181919    #eccentricity, WGS84
    lambda_0 = 0.     # meridian along positive Y axis
    
    phi_c = -phi_c*np.pi/180
    lambda_0 = -lambda_0*np.pi/180
    x=-x
    y=-y

    #this is not commented very well. See Snyder for details.
    t_c=np.tan(np.pi/4 - phi_c/2)/np.power((1-e*np.sin(phi_c))/(1+e*np.sin(phi_c)),e/2)
    m_c=np.cos(phi_c)/np.sqrt(1-np.square(e)*np.square(np.sin(phi_c)))
    rho=np.sqrt(x**2+y**2)
    t=rho*t_c/(a*m_c)

    #find phi with a series instead of iterating.
    chi=np.pi/2 - 2 * np.arctan(t)
    lat=(chi+(e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360)*np.sin(2*chi)
    + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520)*np.sin(4*chi)
    + (7*e**6/120+81*e**8/1120)*np.sin(6*chi)
    + (4279*e**8/161280)*np.sin(8*chi))

    lon=lambda_0 + np.arctan2(x,-y)

    #correct the signs and phasing
    lat=-lat
    lon=-lon
    lon=np.mod(lon+np.pi,2*np.pi)-np.pi #want longitude in the range -pi to pi

    #convert back to degrees
    lat=lat*180/np.pi
    lon=lon*180/np.pi
    return lat, lon

#---------------------------------------------------------------------
# at = geodetic2alongtrack(lat,lon,EPSG)
#
# Takes in:
# lat - scalar or 1D array of latitudes
# lon - scalar or 1D array of longitudes
# EPSG - scalar integer of the coordinate reference system code for 
# converting latitude and longitude to x, y coordinates in meters. Note 
# that only 3413 (Greenland) and (3031) Antarctica are supported.
# 
# Returns:
# at - 1D array with the cumulative distance in meters along the flight 
# track
#----------------------------------------------------------------------
def geodetic2alongtrack(lat,lon,EPSG):
    if EPSG == 3413:
        x,y = ll2psn(lat,lon)
    elif EPSG == 3031:
        x,y = ll2ps(lat,lon)
    else:
        print("Coordinate reference system not supported.")
    
    at = np.zeros(lat.shape)
    at[1:] = np.cumsum(np.sqrt(np.abs(np.diff(x))**2+np.abs(np.diff(y))**2))
    return at
