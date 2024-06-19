import numpy as np

#------------------------------------------------------------------------------------------------
# corrected = GeometricPowerCorrection(surface, time, radar, n)
#
# Takes in:
# surface - 1D array vertical pixel index of the surface echo for each trace
# time - 1D array of the corresponding two-way travel time for each vertical pixel location
# radar - 2D array of the radar returns
# n - scalar refractive index of the material (1.78 for solid ice)
#
# Returns:
# corrected - 2D array of geometrically corrected radar echo power (same size/shape as radar)
#------------------------------------------------------------------------------------------------
def GeometricPowerCorrection(surface, time, radar, n):
    c = 299792458.
    corrected = np.empty(radar.shape)
    for i in range(radar.shape[1]):
        for j in range(radar.shape[0]):
            if j <= surface[i]:
                corrected[j,i] = radar[j,i]*np.square(time.item(j)*c)
            else:
                to_surf = time.item(surface[i])  
                under_surf =  (time.item(j) - time.item(surface[i]))/n   
                corrected[j,i] = radar[j,i]*np.square(c*(to_surf + under_surf))
    return corrected

#------------------------------------------------------------------------------------------------
# corrected, new_surf_ind = ElevationCorrection_Time(surface, elevation, time, radar)
#
# Takes in:
# surface - 1D array of the vertical pixel index of the surface echo for each trace
# elevation - 1D array of elevation in meters above the geoid of the aircraft for each trace
# time - 1D array of the corresponding two-way travel time for each vertical pixel location
# radar - 2D array of the radar returns
#
# Returns:
# corrected - 2D interpolated array of radar echo returns (same size/shape as radar) that removes
# trace to trace variance due to changes in aircraft altitude
# new_surf_ind - 1D array of vertical pixel index of surface echo in the new corrected 2D array
#------------------------------------------------------------------------------------------------
def ElevationCorrection_Time(surface, elevation, time, radar):
    c = 299792458.
    delta_time = 2*(elevation - elevation[0])/c
    corrected = np.empty(radar.shape)
    new_surf_ind = np.empty(surface.shape)
    for i in range(radar.shape[1]):
        new_time = time - delta_time[i]
        new_surf_ind[i] = surface[i] - np.rint(delta_time[i]/np.mean(np.diff(time)))
        corrected[:,i] = np.interp(time, new_time, radar[:,i])
    return corrected, new_surf_ind.astype(int)

#------------------------------------------------------------------------------------------------
# corrected, depth, new_surf_elev = ElevationCorrection_Time(surface, elevation, time, radar)
#
# Takes in:
# surface - 1D array of the vertical pixel index of the surface echo for each trace
# elevation - 1D array of elevation in meters above the geoid of the aircraft for each trace
# time - 1D array of the corresponding two-way travel time for each vertical pixel location
# radar - 2D array of the radar returns
# n - scalar refractive index of the material (1.78 for solid ice)
# res - scalar vertical resolution in meters of the corrected radargram
#
# Returns:
# corrected - 2D interpolated array of radar echo returns (same size/shape as radar) with
# trace to trace variance due to changes in aircraft altitude removed and converted from two-way
# travel time to elevation (e.g. real physical distance between points)
# elev - 1D array of the elevation in meters corresponding to each vertical pixel in a trace
# new_surf_elev - 1D array of the elevation in meters of the surface at each trace
#------------------------------------------------------------------------------------------------
def ElevationCorrection_Depth(surface, elevation, time, radar, n, res):
    c = 299792458.
    elev_full = np.empty(radar.shape)
    for i in range(elev_full.shape[1]):
        for j in range(elev_full.shape[0]):
            if j <= surface[i]:
                elev_full[j,i] = elevation[i] - 0.5*time[j]*c
            else:
                elev_full[j,i] = elevation[i] - 0.5*time[surface[i]]*c - 0.5*(time[j] - time[surface[i]])*(c/n)
                
    elev = np.arange(np.min(elev_full.flatten()),np.max(elev_full.flatten()),res)
    new_surf_elev = np.empty(surface.shape)
    corrected = np.empty((elev.shape[0], radar.shape[1]))
    for i in range(corrected.shape[1]):
        new_surf_elev[i] = elev_full[surface[i], i]
        corrected[:,i] = np.interp(elev,np.flip(elev_full[:,i]),np.flip(radar[:,i]))
    return np.flip(corrected,axis=0), np.flip(elev), new_surf_elev
 