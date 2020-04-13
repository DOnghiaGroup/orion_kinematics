import numpy as np
import pandas as pd
from astropy.coordinates import Angle
from astropy.coordinates.builtin_frames import LSR
from astropy.coordinates import (CartesianRepresentation,CartesianDifferential,CylindricalRepresentation,CylindricalDifferential)
from astropy.coordinates import ICRS, Galactic, GalacticLSR, Galactocentric
import astropy.units as u

'''
Functions for examining the kinematics of Orion clusters
'''

def load_data(filenames):
    """ Return list of dataframes to be worked with
    Parameters
    ----------
    filenames : tuple
        tuple of files to open
    Returns
    -------
    list
        return list of pandas dataframes of the data
    """    
    dataframes = []
    for f in filenames:
        dataframes.append(pd.read_csv(f))
    return dataframes

def with_rv_only(df):
    """ Return data that only contains rows with radial velocities
    Parameters
    ----------
    df : Pandas DataFrame
    """    
    return df.loc[df.radial_velocity.notnull() | df.vhelio.notnull()]

def combine_chen_cam(filename_chen,filename_cam, radial_velocity_only=False):
    """ Merges data frames to contain a complete sample of radial
    velocities from gaia dr2 and apogee
    Parameters
    ----------
    filename_chen : string
        Filename of SNN output
    filename_cam : string
        Filename of my apogee data with vhelio
    Returns
    -------
    pd.DataFrame
        Complete sample of data with radial velocities
    """    
   
    chen_cam_data = load_data([filename_chen,filename_cam])
    df_chen = chen_cam_data[0]
    df_cam = chen_cam_data[1] 
    df_cam = df_cam.drop_duplicates(['gaia_source_id']) # some apogee data have duplicate measurements
    df_complete = pd.merge(df_chen, df_cam, how='left', left_on='source_id', right_on='gaia_source_id') 
    if radial_velocity_only == True:
        df_complete = with_rv_only(df_complete)
    return df_complete

def grouping_ICRS(df, group_select,ra_dec_region=None):
    """ Puts data sample into ICRS equitorial frame based on chosen stellar group
    Parameters
    ----------
    df : pd.DataFrame
        Complete sample of data
    group_select : list
        List of group label numbers to include for analysis
    Returns
    -------
    ICRS
        ICRS coordinate frame
    """
    
    df_group = df[df['label'].isin(group_select)]
    df_group = with_rv_only(df_group)
    if ra_dec_region != None:
        ra_select = ra_dec_region[0]
        dec_select = ra_dec_region[1]
        df_group = df_group.loc[(df_group.ra > ra_select[0]) & (df_group.ra < ra_select[1]) & (df_group.dec > dec_select[0]) & (df_group.dec < dec_select[1])]
    rv_apogee = (df_group.vhelio.dropna()).values 
    x_apogee = df_group.loc[df_group.vhelio.notnull()].ra.values
    y_apogee = df_group.loc[df_group.vhelio.notnull()].dec.values
    dx_apogee = df_group.loc[df_group.vhelio.notnull()].pmra.values
    dy_apogee = df_group.loc[df_group.vhelio.notnull()].pmdec.values
    z_apogee = 1000/df_group.loc[df_group.vhelio.notnull()].parallax.values
    
    rv_gaia = (df_group.radial_velocity.dropna()).values
    x_gaia = df_group.loc[df_group.radial_velocity.notnull()].ra.values
    y_gaia = df_group.loc[df_group.radial_velocity.notnull()].dec.values
    dx_gaia = df_group.loc[df_group.radial_velocity.notnull()].pmra.values
    dy_gaia = df_group.loc[df_group.radial_velocity.notnull()].pmdec.values
    z_gaia = 1000/df_group.loc[df_group.radial_velocity.notnull()].parallax.values

    rv = np.concatenate((rv_apogee,rv_gaia))
    x = np.concatenate((x_apogee,x_gaia))
    y = np.concatenate((y_apogee,y_gaia))
    z = np.concatenate((z_apogee,z_gaia))
    dx = np.concatenate((dx_apogee,dx_gaia))
    dy = np.concatenate((dy_apogee,dy_gaia))
        
    return ICRS(ra = x*u.deg, dec = y*u.deg, distance =z*u.parsec, pm_ra_cosdec = dx*(u.mas/u.yr), pm_dec = dy*(u.mas/u.yr), radial_velocity=rv*(u.km/u.s))

'''
Next two functions return a coordinate frame (cartesian and equitorial) in a stellar group's reference frame 
by subtracting off LSR and average motion of group. Method used in Kounkel et al. (2018)
'''
def rf_cartesian(input_coord_sys, average_velocity = None):
    """ Convert to reference frame of stellar group after subtraction of LSR (U,V,W) = (11.1,12.24,7.25) kms^-1
    defined in Schönrich et al. (2010) and leave it in cartesian coordinates.
    Parameters
    ----------
    input_coord_sys : ICRS
        ICRS input of group for analysis
    average_velocity : tuple
        average cartesian velocity to be subtracted if known,
        else computes it here
    Returns
    -------
    ICRS
        ICRS frame
    """
    galactic = input_coord_sys.transform_to(GalacticLSR) # This conversion implicitly subtracts off the LSR
    galactic.representation_type = 'cartesian'
    galactic.differential_type = 'cartesian'

    x_new = galactic.cartesian.x
    y_new = galactic.cartesian.y
    z_new = galactic.cartesian.z
    
    if average_velocity == None:
        dx_median = np.median(galactic.velocity.d_x.value)
        dy_median = np.median(galactic.velocity.d_y.value)
        dz_median = np.median(galactic.velocity.d_z.value)
        print('median(dx) = {} km/s, median(dy) = {} km/s, median(dz) = {} km/s'.format(dx_median, dy_median, dz_median))
    else:
        dx_median = average_velocity[0]
        dy_median = average_velocity[1]
        dz_median = average_velocity[2]
        
    dx_new = (galactic.velocity.d_x.value - dx_median)*(u.km/u.s)
    dy_new = (galactic.velocity.d_y.value - dy_median)*(u.km/u.s)
    dz_new = (galactic.velocity.d_z.value - dz_median)*(u.km/u.s)
    
    cartesian_representation = ICRS(x_new, y_new, z_new, dx_new, dy_new, dz_new, representation_type = CartesianRepresentation, differential_type = CartesianDifferential)
    return cartesian_representation

def rf_equitorial(input_coord_sys, average_velocities=None):
    """ Convert to reference frame of stellar group after subtraction of LSR (U,V,W) = (11.1,12.24,7.25) kms^-1
    defined in Schönrich et al. (2010) and convert it back to equitorial coordinates. 
    Parameters
    ----------
    input_coord_sys : ICRS
        ICRS input of group for analysis
    average_velocity : tuple
        average cartesian velocity to be subtracted if known,
        else computes it here
    Returns
    -------
    ICRS
        ICRS frame
    """    

    cartesian_representation = rf_cartesian(input_coord_sys,average_velocities=None)
    revert = cartesian_representation.transform_to(ICRS)
    crep = revert
    revert.representation_type = 'spherical'
    revert.differential_type = 'spherical'

    galactic_revert = Galactic(l=revert.ra, b=revert.dec, distance=revert.distance, pm_l_cosb=revert.pm_ra, pm_b = revert.pm_dec, radial_velocity = revert.radial_velocity)

    rf_lsr_equitorial = galactic_revert.transform_to(ICRS)

    return rf_lsr_equitorial 




