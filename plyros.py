# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 10:24:13 2017

@author: RLaugier
"""



# Python standard library imports
import tempfile
import shutil
import sys
# Third-party imports
import gcn
import gcn.handlers
import gcn.notice_types
import requests
import healpy as hp
import numpy as np

from datetime import datetime
import astropy.coordinates
from astropy.time import Time
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii

from time import sleep

import os
import ply_utility
import ply_gen
import ply_sav_fields


            
#parsing for global variables
#print ("Parsing user password")
#try : 
#    userpassword = str(sys.argv[1])
#except:
#    print("error: expected the password for user alert inbetween \"s like \"****\" ")
#print ("Parsing database password")
#try: 
    #databasepassword = str(sys.argv[2])
#    databasepassword = "c1974"
#except:
#    print("error: expected the database password for user ros inbetween \"s like \"****\" ")
print ("initiating test counter")
testcounter = 0
testcycle = 100
try:
    runmode = str(sys.argv[3])
    print("mode",runmode)
except:
    runmode = "run"
            

# Function to call every time a GCN is received.
# Run only for notices of type LVC_INITIAL or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_TEST)
#def process_gcn(payload, root):
def process_gcn():
    global testcounter
    
    # Print the alert
    print('Got VOEvent:')
    #print(payload)

    # Respond only to 'test' events.
    # VERY IMPORTANT! Replce with the following line of code
    # to respond to only real 'observation' events.
    # if root.attrib['role'] != 'observation': return
    
    #This is the test sequence procedure, we run global in test mode once every x time (defined by testcycle)
    if root.attrib['role'] == 'test' :
        print("This is a test: checking testcounter"); sys.stdout.flush()
        print(testcounter); sys.stdout.flush()
        if testcounter == 0:
            print("We run the test this time")
            print(payload)
            #skymap_url = root.find("./What/Param[@name='SKYMAP_URL_FITS_BASIC']").attrib['value']
            skymap_url = "./"
            #process_global(skymap_url,userpassword,databasepassword,test=True)
            process_global(skymap_url,test=True)
        else :
            print ("The test counter is at %i: not doing the drill this time.")
        testcounter += 1
        if testcounter >= testcycle:
            testcounter=0
    
    #This is what happens if the alert is NOT A TEST: we run global in non-test mode
    if root.attrib['role'] == 'observation' :
        print(payload)
        #skymap_url = root.find("./What/Param[@name='SKYMAP_URL_FITS_BASIC']").attrib['value']
        skymap_url = "./" 
        process_global(skymap_url,test=False)


    # Respond only to 'CBC' events. Change 'CBC' to "Burst' to respond to only
    # unmodeled burst events.
    if root.find("./What/Param[@name='Group']").attrib['value'] != 'CBC': return

    # Read out integer notice type (note: not doing anythin with this right now)
    notice_type = int(root.find("./What/Param[@name='Packet_Type']").attrib['value'])

    # Read sky map
    #skymap, header = ply_utility.get_skymap(root)

def process_global(url,test=False):
    '''This function is designed to do all the steps required for the scheduling of
    a whole GW event followup for all the telescopes in the network, based solely on the url of a skymap.
    warning: the URL is used assume the id of the event, so let's hope they do not change their url scheme for now
    In case the alert is a test, the scenes will be teleted after 10s'''
    if test == True:
        print ("Ok, let's go through the drill one more time")
    else :
        print ("Ok this is for real now! This is not an exercise!")
    start_time = Time(datetime.utcnow(), scale='utc')
    sitenames = ["'Tarot_Calern'","'Tarot_Chili'","'Tarot_Reunion'"]
    siteids = [1,2,8]
    scenes = Table()
    #name = os.path.basename(url)
    #alertname = url.split("/")[5]
    #if not os.path.isdir(alertname):
    #    os.makedirs(alertname)
    #print("retrieving skymap from %s"% url); sys.stdout.flush()
    #download_skymap(url, alertname, name)

    
    print("loading skymap named %s"% name); sys.stdout.flush()
    #hpx, header = load_skymap(os.path.join(alertname,name))
    hpx, header = ply_utility.load_skymap("LALInference_skymap.fits")
    #for site in sitenames:
    sitescenes = Table()
    siteID = 1 #it will come later from the xml request file
    location, horizondef, horizontype, hadeclims, idtelescope, site = ply_utility.read_xml(siteID)
    fields,mylocation,sitescenes = main(hpx, header, siteID)
    ply_sav_fields.sav(fields,site,location) #sav the fields at your format
    #this fixes problems arising when a table was empty
    if len(sitescenes)>0 and type(sitescenes) == astropy.table.table.Table:
        if len(scenes)>0:
            print ("joining")
            scenes.pprint(max_width=-1)
            sitescenes.pprint(max_width=-1)
            scenes = astropy.table.vstack([scenes,sitescenes])
        else :
            scenes = sitescenes

    print (scenes)
    
    settime,readoutTime,exps,filters = ply_utility.site_timings(site)

    #idreq = cadorrest.post_request(alertname + "autogen","0","90",pwd)
    idscene = []
    ra = []
    dec = []
    timeisot = []

    for index in np.arange(0,len(scenes),1):
        #cadorrest.post_scene(prefix,idreq,entry,exps,filters,pwd):
        #myidscene, myra, mydec, mytimeisot = cadorrest.post_scene(alertname+"_",idreq,scenes[index],exps,filters,pwd)
        idscene.append(myidscene)
        ra.append(myra)
        dec.append(mydec)
        timeisot.append(mytimeisot)
    scenes ["idscene"] = idscene
    scenes ["ra"] = ra
    scenes ["dec"] = dec
    scenes ["timeisot"] = timeisot
        
    end_time = Time(datetime.utcnow(), scale='utc')
    length = end_time - start_time
    print ("Executed in ",length.sec,"seconds")
    print ("Waiting for planification (10s)")
    sys.stdout.flush()
    sleep(10)
    
    #Downloading planification logs
    print("Keeping a little souvenir"); sys.stdout.flush()
    #for idteles in siteids:
    #    planiurl = "http://cador.obs-hp.fr/ros/sequenced" + str(idteles) + ".txt"
    #    rejecurl = "http://cador.obs-hp.fr/ros/rejected" + str(idteles) + ".txt"
    #    try :
    #        download_skymap(planiurl,alertname,"sequenced"+str(idteles)+".txt")
    #        download_skymap(rejecurl,alertname,"rejected"+str(idteles)+".txt")
    #        print ("Succesfully downloaded plani files")
    #    except:
    #        print ("error copying plani")
    sys.stdout.flush()    
    ascii.write(scenes, os.path.join(alertname,'Planification_table.csv'), format='csv', fast_writer=False)
    if test == True:
        print ("This was just an exercise, let's delete the request now")
        #cadorrest.remove_request(idreq,pwd)
    else :
        print ("This was for real: the scenes will be observed")
#        for i in np.arange(0,900,10):
#            if pyrosutilities.replica_is_running():
#                print("Waiting until replica has finished")
#                sleep (10)
#            else:
#                print("Lauching replica to propagate planification")
#                #subprocess.check_call("nohup", "php", "/srv/develop/ros_private_cador/src/replica2/replica_slow.php", "1>", "/srv/www/htdocs/ros/logs/replica/replica_slow.txt")
#                try:
#                    subprocess.check_call("php", "/srv/develop/ros_private_cador/src/replica2/replica_slow.php")
#                    print("Replica finished successfuly, we are done, here")
#                except:
#                    print("Replica was probably crashed: waited for %i and got no response"% i)
#                break
        print("Unfortunately, we cannot run replica from vega yet. We have to be patient and wait for it to run on schedule")
        
    return scenes
        

'''"https://gracedb.ligo.org/apibasic/events/M131141/files/bayestar.fits.gz"'''
'''"https://gracedb.ligo.org/apibasic/events/G277583/files/skyprobcc_cWB.fits"'''
def main(hpx,header,siteID):
    '''This function handles the scenes and planning creation for a single telescope.
    It returns a list of scene entries in a table, ready to be submited to CADOR'''
    current_time = Time(datetime.utcnow(), scale='utc')
    print("Date is");print(current_time.value)

##################################################################################
###Site    
    #site = "'Tarot_Reunion'"
    location, horizondef, horizontype, hadeclims, idtelescope, site = ply_utility.read_xml(siteID)
    nfields = ply_utility.site_number(site)
    print ("working on site: %s"% (site)); sys.stdout.flush()
    #location, horizondef, horizontype, hadeclims, idtelescope = ply_utility.get_obs_info(site,dbpwd)
    total, a,b,c = ply_utility.optimize_quin(hpx, nfields, ply_utility.site_field(site))
    myfields = ply_gen.build_fields(hpx, a, b, c, nfields*2, ply_utility.site_field(site))
    print (myfields); sys.stdout.flush()
    thefields = ply_scheduler.clean_table(myfields)
    print (thefields); sys.stdout.flush()
    print(thefields["coords"][0])
    print ("Observavility from %s, at location %s" % (site, location))
    toremove = []
    for index in np.arange(0,len(thefields),1):
        if ply_utility.check_declination(location.latitude,thefields[index]["coords"].dec) == 0:
            toremove.append(index)
    thefields.remove_rows(toremove)
    mycyclegrid,scenelength = ply_scheduler.build_cyclegrid(15,site,24)
    #Determining observability to remove excess fields
    fieldobservability = []
    for i in np.arange(0,len(thefields),1):
        obsevability = []
        for j in np.arange(0,len(mycyclegrid),1):
            result = ply_utility.is_observable(thefields["coords"][i],current_time + mycyclegrid["date"][j]*u.second,location,horizondef,horizontype,hadeclims)
            obsevability.append(result)
        print (obsevability)
        fieldobservability.append(obsevability)
    #thefields["obsevability"]=fieldobservability
    thefields["observability"] = fieldobservability
    print (thefields,fieldobservability)
    toremove=[]
    for i in np.arange(0,len(thefields),1):
        if np.count_nonzero(thefields["observability"][i])<3:
            toremove.append(i)
    thefields.remove_rows(toremove)
    print(thefields)
    if len(thefields)>nfields:
        thefields=thefields[0:nfields]
    elif len(thefields)==0:
        return thefields,location,fieldobservability
    thefields["index"]=np.arange(0,len(thefields),1)
    print(thefields)
    print("proba covered:",np.sum(thefields["proba"]))
    #unable to sort the fields at this point because of SkyCoord object.
    #thefields["temp"]=thefields["coords"].ra
    mycyclegrid,scenelength = ply_scheduler.build_cyclegrid(len(thefields),site,48)
    
    scenes = []
    for i in np.arange(0,len(thefields),1):
        obsevability = []
        timeindex = 0
        for j in np.arange(0,len(mycyclegrid),1):
            exacttime = ply_scheduler.mycyclegrid["date"][j]*u.second + i*scenelength*u.second + current_time
            result = ply_utility.is_observable(thefields["coords"][i],exacttime,location,horizondef,horizontype,hadeclims)
            if result == 1:
                scenes.append([site,thefields["index"][i],timeindex,exacttime,thefields["coords"][i]])
                timeindex += 1
        
        fieldobservability.append(obsevability)
    revscenes = np.transpose(scenes)
    thescenes = Table()
    thescenes["site"]=revscenes[0]
    thescenes["index"]=revscenes[1]
    thescenes["tindex"]=revscenes[2]
    thescenes["time"]=revscenes[3]
    thescenes["coords"]=revscenes[4]
    thescenes["idtelescope"] = idtelescope

    return thefields,location,thescenes

if runmode == "run":
    print ("the program is now ready: listening for events")
    sys.stdout.flush()
    #gcn.listen(port=8096, handler=process_gcn)
    process_gcn
    
elif runmode =="tools":
    print ("The tools are loaded, ready for interactive mode (use ipython)")
