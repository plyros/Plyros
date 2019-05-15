# -*- coding: utf-8 -*-
"""
Choose your fields format to be saved 
"""

# Python xml libraries
from lxml import etree

def sav(fields,sitename,location): #xml format
    #Main tree
    coord = etree.Element("coord")
    #Sub tree Roots
    ra = etree.SubElement(coord, "ra")
    ra.text = fields["ra"]
    dec = etree.SubElement(coord, "dec")
    dec.text = fields["dec"]
    frame = etree.SubElement(coord, "frame")
    frame.text = fields["frame"]
    loc = etree.SubElement(coord, "location")
    loc.text = location

    file = open("ply_sav_fields_"+sitename+".xml", "w", encoding='UTF-8')
    file.write(str(etree.tostring(coord, pretty_print = True)))

    file.close()

#---- test part ----
#fields = {}
#fields["ra"] = "0.1"
#fields["dec"] = "0.2"
#fields["frame"] = "bla bla"
#site = "pc"
#location = "(4581869.27527051, 556376.42045728, 4389143.90956099) m"

#sav(fields,site,location)
