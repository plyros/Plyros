# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 16:43:57 2017

@author: echo
"""
# Python standard library imports
import tempfile
import shutil
import sys

# Third-party imports
import numpy as np
import healpy as hp

import astropy.units as u
from astropy.table import Table
import astropy.coordinates

import subprocess
import requests

from lxml import etree
import pymysql
import os
import ply_gen

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

def getxmlval(root,tag):
    ''' a simple function to retrieve a value from an XML etree object'''
    vals=[]
    for element in root.iter(tag=tag):
        vals.append(element.text)
    return vals
        
def gets(ra,dec,field):
    '''Just a convenient function that returns the corners of
    a field of view from its center'''
    demichamp = field/2
    racorners = np.zeros(4)
    deccorners = np.zeros(4)
    racorners[0], deccorners[0] = ra + demichamp, dec - demichamp
    racorners[1], deccorners[1] = ra + demichamp, dec + demichamp
    racorners[2], deccorners[2] = ra - demichamp, dec + demichamp
    racorners[3], deccorners[3] = ra - demichamp, dec - demichamp
    return racorners, deccorners


def frange(x, y, jump):
    while x < y:
        yield x
        x += jump
        
'''conn = pymysql.connect(host='tarot9.oca.eu', user='tarot', password=pwd, db='ros')'''
                
def get_obs_info(sitename,pwd):
    '''This function retrieves all the necessary info from ROS telescopes database'''
    conn = pymysql.connect(host='tarot9.oca.eu', user='ros', password='ros123', db='ros')
    error, idtelescope = get_db_info(conn, "telescopes", "idtelescope", "name", sitename)
    error, latitude = get_db_info(conn, "telescopes", "latitude", "name", sitename)
    error, longitude = get_db_info(conn, "telescopes", "longitude", "name", sitename)
    error, altitude = get_db_info(conn, "telescopes", "altitude", "name", sitename)
    error, sens = get_db_info(conn, "telescopes", "sens", "name", sitename)
    error, horizontype = get_db_info(conn, "telescopes", "horizontype", "name", sitename)
    error, horizondef = get_db_info(conn, "telescopes", "horizondef", "name", sitename)
    error, limdecmin = get_db_info(conn, "telescopes", "limdecmin", "name", sitename)
    error, limdecmax = get_db_info(conn, "telescopes", "limdecmax", "name", sitename)
    error, limharise = get_db_info(conn, "telescopes", "limharise", "name", sitename)
    error, limhaset = get_db_info(conn, "telescopes", "limhaset", "name", sitename)
	    
    if error != 0 : 
        print("Error recovering data")
    #Turn location information into an EarthLocation object
    location = transform_location(latitude, longitude, sens, altitude)
    conn.close()
    print(horizontype, horizondef,limharise, limhaset)
    print("converting horizon")
    print (horizondef, horizontype)
    #Convert the horizon information into data Table of angle values
    myhorizon = convert_horizon(horizondef,horizontype)
    print(myhorizon)
    hadeclims = (limdecmin,limdecmax,limharise,limhaset)
    return location, myhorizon, horizontype, hadeclims, idtelescope

def read_xml(ID):
    #Main tree
    try:
        tree = etree.parse("ply_conf.xml")
        root = tree.getroot()
    except:
        print("Error: Fail to read ply_conf.xml")
    siteid = list()
    name = list()
    latitude = list()
    longitude = list()
    altitude = list()
    sens = list()
    horizontype = list()
    horizondef = list()
    limdecmin = list()
    limdecmax = list()
    limharise = list()
    limhaset = list()
    nbfield = list()
    freetime = list()
    for param in tree.xpath("/tel/tele"):
        siteid.append(param.get("id"))
    for param in tree.xpath("/tel/tele/name"):
        name.append(param.text)
    for param in tree.xpath("/tel/tele/latitude"):
        latitude.append(param.text)
    for param in tree.xpath("/tel/tele/longitude"):
        longitude.append(param.text)
    for param in tree.xpath("/tel/tele/altitude"):
        altitude.append(param.text)
    for param in tree.xpath("/tel/tele/sens"):
        sens.append(param.text)
    for param in tree.xpath("/tel/tele/horizontype"):
        horizontype.append(param.text)
    for param in tree.xpath("/tel/tele/horizondef"):
        horizondef.append(param.text)
    for param in tree.xpath("/tel/tele/limdecmin"):
        limdecmin.append(param.text)
    for param in tree.xpath("/tel/tele/limdecmax"):
        limdecmax.append(param.text)
    for param in tree.xpath("/tel/tele/limharise"):
        limharise.append(param.text)
    for param in tree.xpath("/tel/tele/limhaset"):
        limhaset.append(param.text)
    for param in tree.xpath("/tel//opt/nbfield"):
        nbfield.append(param.text)
    for param in tree.xpath("/tel//opt/freetime"):
        freetime.append(param.text)

    nom = ''
    lat = ''
    lon = ''
    alt = ''
    se = ''
    horitype = ''
    horidef = ''
    limin = ''
    limax = ''
    limrise = ''
    limset = ''
    nbfields = ''
    freetim = ''

    if ID == 1:
        nom = name[0]
        lat = float(latitude[0])
        lon = float(longitude[0])
        alt = float(altitude[0])
        se = sens[0]
        horitype = horizontype[0]
        horidef = horizondef[0]
        limin = float(limdecmin[0])
        limax = float(limdecmax[0])
        limrise = float(limharise[0])
        limset = float(limhaset[0])
    elif ID == 2:
        nom = name[1]
        lat = float(latitude[1])
        lon = float(longitude[1])
        alt = float(altitude[1])
        se = sens[1]
        horitype = horizontype[1]
        horidef = horizondef[1]
        limin = float(limdecmin[1])
        limax = float(limdecmax[1])
        limrise = float(limharise[1])
        limset = float(limhaset[1])
    elif ID == 8:
        nom = name[2]
        lat = float(latitude[2])
        lon = float(longitude[2])
        alt = float(altitude[2])
        se = sens[2]
        horitype = horizontype[2]
        horidef = horizondef[2]
        limin = float(limdecmin[2])
        limax = float(limdecmax[2])
        limrise = float(limharise[2])
        limset = float(limhaset[2])
    elif ID == 3:
        nom = name[3]
        lat = float(latitude[3])
        lon = float(longitude[3])
        alt = float(altitude[3])
        se = sens[3]
        horitype = horizontype[3]
        horidef = horizondef[3]
        limin = float(limdecmin[3])
        limax = float(limdecmax[3])
        limrise = float(limharise[3])
        limset = float(limhaset[3])
    else:
        print("Wrong site id!")
    
    nbfields = nbfield[0]
    freetim = freetime[0]

    location = transform_location(lat, lon, se, alt)
    myhorizon = convert_horizon(horidef, horitype)
    hadeclims = (limin, limax, limrise, limset)

    return location, myhorizon, horitype, hadeclims, ID, nom

def transform_location(latitude, longitude, sens, altitude):
    '''This function returns an astropy Earthlocation object based on ROS location values'''
    if sens == 'E':
        location = astropy.coordinates.EarthLocation(longitude*u.deg, latitude*u.deg, altitude*u.m)
    elif sens == 'W':
        location = astropy.coordinates.EarthLocation(-longitude*u.deg, latitude*u.deg, altitude*u.m)
    return location
    
def convert_horizon(horizondef,horizontype):
    '''This converts the string of ROS's horizondef into a suitable table of values'''
    lims = horizondef.split("} {")
    lims[0] = lims[0].replace("{","")
    lims[len(lims)-1] = lims[len(lims)-1].replace("}","")
    lims[:]=[x.split(" ") for x in lims]
    if horizontype == "hadec" :
        limites = Table(names =["dec","limrise","limset"], dtype =("f4","f4","f4"))
        return limites
    elif horizontype == "altaz":
        limites = Table(names =["azim","elev"], dtype =("f4","f4"))
    for item in lims :
        limites.add_row(item)
    limites = limites*u.deg
    return limites

def get_db_info(connection, table, entry, keycolumn, keyvalue):
    '''This utilitarian function retrieves a bit of information from a distant table.
    The connection must first be established with pymysql.connect(), then passed in a 1st argument'''
    error = 0
    query = ""
    query = "SELECT " + entry + " FROM " + table + " WHERE " + keycolumn + "= " + keyvalue + ";"
    mycursor = connection.cursor()
    countrow = mycursor.execute(query)
    if countrow != 1 :
        error =1
        print("error fetching data in database")
        value = ""
    value = mycursor.fetchone()[0]
    return error, value
    
def site_field(site):
    '''This is a crude definition of each site's FoV'''
    for case in switch(site):
        if case("'Tarot_Calern'"):
            return 1.86
            break
        if case("'Tarot_Chili'"):
            return 1.86
            break
        if case("'Tarot_Reunion'"):
            return 4.2
            break
        if case("'Zadko_Australia'"):
            return 0.5
            break

def site_number(site):
    '''This is a crude definition of each site's targeted number of fields'''
    for case in switch(site):
        if case("'Tarot_Calern'"):
            return 5
            break
        if case("'Tarot_Chili'"):
            return 5
            break
        if case("'Tarot_Reunion'"):
            return 10
            break
        if case("'Zadko_Australia'"):
            return 5
            break

def site_timings(site):
    '''This is a crude definition of each site's timing parameters
    as well as exposure scheme'''
    for case in switch(site):
        if case("'Tarot_Calern'"):
            settime = 30     
            readoutTime = 10
            exptime = 120
            exps = [exptime,exptime,0,0,0,0]
            fil = "1"
            filters = [fil,fil,fil,fil,fil,fil]
            return settime,readoutTime,exps,filters
            break
        if case("'Tarot_Chili'"):
            settime = 30     
            readoutTime = 10
            exptime = 120
            exps = [exptime,exptime,0,0,0,0]
            fil = "1"
            filters = [fil,fil,fil,fil,fil,fil]
            return settime,readoutTime,exps,filters
            break
        if case("'Tarot_Reunion'"):
            settime = 30     
            readoutTime = 10
            exptime = 120
            exps = [exptime,exptime,exptime,0,0,0]
            fil = "0"
            filters = [fil,fil,fil,fil,fil,fil]
            return settime,readoutTime,exps,filters
            break
        if case("'Zadko_Australia'"):
            settime = 90     
            readoutTime = 10
            exptime = 120
            exps = [exptime,exptime,0,0,0,0]
            fil = "1"
            filters = [fil,fil,fil,fil,fil,fil]
            return settime,readoutTime,exps,filters
            break

def replica_is_running():
    '''This simple function checks if replica is running on this machine
    It is useless unless replica runs locally'''
    ps_replica = subprocess.check_output(["ps","-edf"])
    if ps_replica.count("replica") > 0:
        return True
    else:
        return False


def download_skymap(myurl,alertname,name):
    '''This function actually downloads the file from an url then moves it to
    the directory of the corresponding name.'''
    #subprocess.check_call(['curl', '-O', '--netrc', myurl])
    subprocess.check_call(['wget','--auth-no-challenge', myurl])
    shutil.move(name,os.path.join(alertname,name))

def load_skymap(myfile):
    '''Loads a healpix skymap'''
    hpx, header = hp.read_map(myfile, h=True, verbose=True)
    
    return hpx, header

def skymap_properties(hpx):
    '''Reads and displays some of the main properties of a skymap'''
    print ("Number of pixels:"); print (len(hpx))

    nside = hp.npix2nside(len(hpx))
    print ("The lateral resolution (nside) is:"); print (nside)
    
    sky_area = 4 * 180**2 / np.pi
    print("pixel per degree:")
    print(len(hpx) / sky_area)
    

def ratophi(ra):
    '''Converts a RA value into a phi angle for healpy to use. Also uses checkphi
    to check for overflows in the operation
    ra*u.deg -> phi*u.rad'''
    myphi = np.deg2rad(ra)
    myphi = checkphi(myphi)
    return myphi
def phitora(phi):
    '''Converts a phi angle from healpy coordinates into RA angle
    phi*u.rad -> ra*u.deg'''
    ra = np.rad2deg(phi)
    return ra
        
    
def dectotheta(dec):
    '''Converts a DEC value into a theta angle for healpy to use. Also uses checktheta
    to check for overflows in the operation
    dec*u.deg -> theta*u.rad'''
    theta = 0.5*np.pi - np.deg2rad(dec)
    theta = checktheta(theta)
    return theta
def thetatodec(theta):
    '''Converts a theta angle from healpy into a DEC angle
    theta*u.rad -> dec*u.deg'''
    dec = np.rad2deg(0.5*np.pi - theta)
    if dec > 180 : print("error"); return -1000;
    if dec < -180 : print("error"); return -1000;
    else : return dec
    
    
def checkphi(phi):
    '''Checks for a possible overflow in the value of a phi angle
    so that it remains consistant with the system used in healpy'''
    phi = phi % (np.pi*2)
    return phi    
    
    
def checktheta(theta):
    '''Checks for a possible overflow in the value of a phi angle
    so that it remains consistant with the system used in healpy'''
    theta = theta % np.pi
    return theta

def checkra(ra):
    '''Checks for a possible overflow in the value of a RA angle
    so that it remains consistant with the system used in healpy'''
    if type(ra) == float : ra = ra % 360
    else : print("error: need float type for checkra"); return;
    return ra

def slicesky(d0, s0,s1,field):
    '''This function builds a list of coordinates that cover the whole
    sky with FoV of a givent instrument. Parameters d0, s0 and s1 are 
    used to give tiling variations'''
#d0, s0, s1 are variables that allow shifted and skew slicing
#d0<field
#s0<field
#s1<s0

#create a declination slicing
    decgrid = list(np.arange(-90 + d0,90,field))
    mylistra, mylistdec = [], []
    n = 0
#for each declination, create a ra slicing
    for dec in decgrid :
        rastep = field/(np.cos(np.deg2rad(dec)))
        ragrid = list(np.arange(s0 + n*s1, s0 + n*s1 + 360, rastep))
        decfields = []
        for myfield in ragrid:
            decfields.append(dec)
        mylistra.extend(ragrid)
        mylistdec.extend(decfields)
        n += 1


    if len(mylistra) != len(mylistdec) : print ("error, couldn't match ra to dec"); return 1;
    return mylistra, mylistdec
    
    
def get_field_value(hpx,ra,dec,field):
    '''This function returns the integral probability inside a disc field of view.
    It is more reliable than the get_fast_field_value, but also slower.'''
    nside = hp.npix2nside(len(hpx))
#    pixel = hp.ang2pix(mynside, dectotheta(dec), ratophi(ra))
    xyz = hp.ang2vec(dectotheta(dec),ratophi(ra))
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(field /2))#* (sq2+1)/4)here radius seems to be a diameter instead * (sq2+1)/4)
    totdisc = hpx[ipix_disc].sum()
    return totdisc

def get_fast_field_value(hpx,ra,dec,field):
    '''This returns the skymap pixel value of the center of a field'''
    nside = hp.npix2nside(len(hpx))
    mypix = hp.ang2pix(nside,dectotheta(dec),ratophi(ra))
    return hpx[mypix]
    
def get_fields_value(hpx,myfieldsra,myfieldsdec, field):
    '''This returns the disc integral values for a whole list of fields'''
    nside = hp.npix2nside(len(hpx))
    prob=[]
    for indec in range(0, len(myfieldsra)):
#        cornerra, cornerdec = getcorners(keptfields["ra"][indec],keptfields["dec"][indec],field=4.2)
#        xyz = hp.ang2vec(checktheta(dectotheta(cornerdec)),checkphi(ratophi(cornerdec)))
#        hp.query_polygon(nside,xyz)
        xyz = hp.ang2vec(dectotheta(myfieldsdec[indec]),ratophi(myfieldsra[indec]))
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(field  /2))#* (sq2+1)/4)here radius seems to be a diameter instead * (sq2+1)/4)
        totdisc = hpx[ipix_disc].sum()
        prob.append(totdisc)
    return prob
    

def calculate_efficiency(hpx, d0, s0, s1, nfields, fieldop):
    '''This function determines the score of a given set of tiling parameters for the optimization.
    The figure of merit is the sum value of the for the <nfields> higher fields of the tiling'''
    nside = hp.npix2nside(len(hpx))

    keptfields = ply_gen.build_fields(hpx, d0,s0,s1,nfields,fieldop)
#    totaldots = np.sum(keptfields["prob"])
    total= 0
    prob_integral = []
    for indec in range(0, len(keptfields)):
#        cornerra, cornerdec = getcorners(keptfields["ra"][indec],keptfields["dec"][indec],field=4.2)
#        xyz = hp.ang2vec(checktheta(dectotheta(cornerdec)),checkphi(ratophi(cornerdec)))
#        hp.query_polygon(nside,xyz)
        xyz = hp.ang2vec(dectotheta(keptfields["dec"][indec]),ratophi(keptfields["ra"][indec]))
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(fieldop))#here radius seems to be a diameter instead * (sq2+1)/4)
        totdisc = hpx[ipix_disc].sum()
        prob_integral.append(totdisc)
        total += totdisc
    
#    efficiency = total / nfields
    #print (total)
    return total, prob_integral
        
        
def optimize_quin(hpx, nfields, fieldaz):
    '''This function returns uptimized tiling parameters based on their calculate_efficiency() score for <nfields>'''
    mymax, d0max, s0max, s1max = 0,0,0,0
    numberofiterations = 0
    for d0 in np.arange(0,fieldaz*3/4,fieldaz/4):
        for s0 in np.arange(0,fieldaz*3/4,fieldaz/4):
            for s1 in np.arange(-fieldaz/3,fieldaz/3,fieldaz/9):
                value, prob_integral = calculate_efficiency(hpx, d0, s0, s1, nfields, fieldaz)
                numberofiterations += 1
                if value > mymax :
                    mymax, d0max, s0max, s1max = value, d0, s0, s1
    print(numberofiterations, mymax)
    return mymax, d0max, s0max, s1max
 
   
'''
total, a,b,c = optimize_quin(hpx, 5, 4.2)
build_fields(hpx, a, b, c, 5, 4.2)
'''


def get_skymap(root):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.
    """
    # Read out URL of sky map.
    # This will be something like
    # https://gracedb.ligo.org/apibasic/events/M131141/files/bayestar.fits.gz
    skymap_url = root.find(
        "./What/Param[@name='SKYMAP_URL_FITS_BASIC']").attrib['value']

    # Send HTTP request for sky map
    response = requests.get(skymap_url, stream=True)

    # Uncomment to save VOEvent payload to file
    # open('example.xml', 'w').write(payload)

    # Raise an exception unless the download succeeded (HTTP 200 OK)
    response.raise_for_status()

    # Create a temporary file to store the downloaded FITS file
    with tempfile.NamedTemporaryFile() as tmpfile:
        # Save the FITS file to the temporary file
        shutil.copyfileobj(response.raw, tmpfile)
        tmpfile.flush()

        # Uncomment to save FITS payload to file
        # shutil.copyfileobj(reponse.raw, open('example.fits.gz', 'wb'))

        # Read HEALPix data from the temporary file
        skymap, header = hp.read_map(tmpfile.name, h=True, verbose=False)
        header = dict(header)

    # Done!
    return skymap, header


def interpolate(dataset, target,horizontype):
    '''This function is used to interpolate an horizon definition.
    It works with altaz horizons but was adapted to ROS hadec horizons (pourly standardized)
    This assumes that the horizondef absissa is monotonous and increasing!!!
    This assumes that the intervals cover ALL the possibilities
    It should check that the target is inside the interval '''
    #print (target)
    if horizontype == "altaz":

        value = 90*u.deg #this default value will prevent observation in case of a problem
        if ((target < dataset["azim"][0])or(target > dataset["azim"][len(dataset["azim"])-1])):
            print ("error, target out of range")
            return value
        for index in np.arange(0,len(dataset["azim"])-1,1):# "-1" allows us to go from 0 to the penultimate
            #print("Between",dataset["azim"][index],"and ",dataset["azim"][index+1])
            if ((dataset["azim"][index] <= target) and (dataset["azim"][index+1]>= target)):
                #Interpolation of target in
                value = dataset["elev"][index] + (target - dataset["azim"][index]) * (dataset["elev"][index+1] - dataset["elev"][index]) / (dataset["azim"][index+1] - dataset["azim"][index])
                break
        if value == 90*u.deg :
            print ("error: Couldn't interpolate horizon. Azimuth, result =============================== ")
            print(target,value)
            return value
        else :
            #print ("horizon elevation at site is")
            #print (value)
            return value
    elif horizontype == "hadec":
        print ("this function doesn't support hadec horizon yet")
    else : print ("error: unknown horizon type")

    return value

def checkhorizon(horizondef,horizontype):
    return

def checksun(coordinates, time, mylocation):
    '''This function checks for the presence of night at the site(-5° of sun elevation)
    and an agular separation to the sun greater than 30°'''
    observable = 1
    astropy.coordinates.solar_system_ephemeris.set('builtin')
    SunObject = astropy.coordinates.get_sun(time)
    sunaltaz = SunObject.transform_to(astropy.coordinates.AltAz(obstime=time, location=mylocation))
    if sunaltaz.alt >= -5*u.deg:
        observable = 0
        return observable
    if SunObject.separation(coordinates) <= 30*u.deg:
        observable = 0
        print("Too close to the sun")
    return observable

def checkmoon(coordinates, time, mylocation):
    '''This checks for a moon distance larger than 30°'''
    observable = 1
    astropy.coordinates.solar_system_ephemeris.set('builtin')
    MoonObject = astropy.coordinates.get_moon(time)
    if MoonObject.separation(coordinates) <= 30*u.deg:
        observable = 0
        print("Object too close to the moon")
    return observable

def checkelev(coordinates, time, mylocation, horizondef, horizontype):
    '''This makes a basic check for an elevation higher tha0 10° AND in case 
    of an altaz horizon definition, checks for an elevation higher than this horizon.'''
    observable = 1
    #objaltaz = coordinates.transform_to(astropy.coordinates.AltAz(time), location=mylocation)
    objaltaz = coordinates.transform_to(astropy.coordinates.AltAz(location=mylocation, obstime=time))
    if objaltaz.alt <= 10*u.deg:
        observable = 0
        #print("Object too low (10°)")
        return observable
    #This condition is assumed sufficient for unsupported hadec horizons
    if horizontype == "hadec":
        return observable
    local_horizon = interpolate(horizondef, objaltaz.az, horizontype)
    if objaltaz.alt - local_horizon <= 0*u.deg:
        observable = 0
        print("Object below site horizon")
    return observable

def checkhadec(coordinates, time, mylocation, hadeclims):
    '''This function checks for compliance to simple hadec mount limits, in 
    a format compliant to the ROS telescopes database table'''
    #hadeclims = (limdecmin,limdecmax,limharise,limhaset)
    observable = 1
    #objaltaz = coordinates.transform_to(astropy.coordinates.AltAz(time), location=mylocation)
    loctime = Time(time,location = mylocation)
    LST = loctime.sidereal_time("apparent")
    
    LHA = coordinates.ra + LST
    #testing for a valid hour angle
    if LHA >= 24 * u.hourangle or LHA < 0 * u.hourangle :
        #print ("hour angle overflow",LHA)
        LHA = LHA - 24*u.hourangle
        #print ("corrected:",LHA)
    #Testing if object is in hadec blind spot
    if LHA >= hadeclims[2]*u.deg and LHA <= hadeclims[3]*u.deg:
        observable = 0
        print("Object HA below limits", LHA)
        return observable
        
    if coordinates.dec <= hadeclims[0]*u.deg:
        observable = 0
        print("Object DEC below limits=================================", coordinates.dec)
        return observable
    if coordinates.dec >= hadeclims[1]*u.deg:
        observable = 0
        print("Object DEC over limits==================================", coordinates.dec)
        return observable
    
    return observable

def check_declination(latitude,declination):
    '''This method is a preliminary check designed to eliminate targets never observable from
    the latitude of the observatory, and reduce calculation steps. This assumes an elevation margin of 10°'''
    observable = 1
    if latitude > 0:
        if declination <  (latitude - 80*u.deg):
            observable = 0
            print("field rejected for declination too low",declination)
    elif latitude < 0:
        if declination > latitude + 80*u.deg:
            observable = 0
            print("field rejected for declination too low",declination)
    return observable
    
def is_observable(coords,time,location,horizondef,horizontype,hadeclims):
    '''Thid function regroups all observability check into one'''
    sunok = checksun(coords, time, location)
    #Here the sunok is checked and return first to save computation time
    #since night/day condition is a major cause of non-observability
    if sunok == 0: return 0
    moonok = checkmoon(coords, time, location)
    elevok = checkelev(coords, time, location, horizondef,horizontype)
    hadecok = checkhadec(coords, time, location, hadeclims)
    #print (sunok, moonok, elevok,hadecok)
    observable = sunok * moonok * elevok * hadecok
    return observable
    

def next_sunset(mylocation, mytime):
    '''This function returns a rough estimation of the next sunset time'''
    astropy.coordinates.solar_system_ephemeris.set('builtin')
    SunObject = astropy.coordinates.get_sun(mytime)
    sunaltaz = SunObject.transform_to(astropy.coordinates.AltAz(obstime=mytime, location=mylocation))
    if sunaltaz.alt > -10*u.deg :
        night = 0
    else :
        night = 1
    set_time=-100*u.hour
    for i in np.arange(0, 24, 0.25):
        if night == 0:
            
            time = mytime + i*u.hour
            sunaltaz = SunObject.transform_to(astropy.coordinates.AltAz(obstime=time, location=mylocation))
            if night == 0 and sunaltaz.alt < -10*u.deg:
                set_time = time
                break
        elif night == 1:
            night = 1
    return set_time

