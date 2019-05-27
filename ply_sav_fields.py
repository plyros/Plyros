# -*- coding: utf-8 -*-
"""
Choose your fields format to be saved 
"""

# Python xml libraries
from lxml import etree

def sav(fields,sitename): #xml format
    #Main tree
    field = etree.Element("field")
    #Sub tree Roots
    for i in range(len(fields)):
        index = etree.SubElement(field, "index")
        index.text = str(i)
        proba = etree.SubElement(index, "proba")
        proba.text = str(fields["proba"][i])
        coords = etree.SubElement(index, "coords")
        coords.text = str(fields["coords"][i])
        obser = etree.SubElement(index, "observability")
        obser.text = str(fields["observability"][i])

    file = open("ply_sav_fields_"+sitename+".xml", "w", encoding='UTF-8')
    file.write(str(etree.tostring(field, pretty_print = True)))

    file.close()

#---- test part ----
#fields = {}
#fields["ra"] = "0.1"
#fields["dec"] = "0.2"
#fields["frame"] = "bla bla"
#site = "pc"
#location = "(4581869.27527051, 556376.42045728, 4389143.90956099) m"

#sav(fields,site,location)
