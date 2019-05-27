#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Creates a XML file containing all Tarot telescope informations
"""

# Python xml libraries
from lxml import etree

# Main tree
tel = etree.Element("tel")

# Sub tree
param = [
    ("1", "Tarot_Calern", "43.75203", "6.92353", "1320.0", "E", "hadec", "{-27 0.0 0.0} {-23.5 330.0 30.0} \
{-22 311.3 48.8} {-21 295.0 65.8} {-19 290.0 70.8} {-15 290.0 73.3} {-10 282.5 78.0} {-5 278.8 83.8} \
{0 272.5 85.5} {5 266.3 90.5} {10 260.8 96.0} {15 257.0 101.3} {20 250.5 107.0} {25 245.0 113.8} \
{30 236.3 120.0} {40 225.0 135.0} {45 233.8 135.0}", "-27.0", "90.0", "180.0", "180.0"),
    ("2", "Tarot_Chili", "-29.259917", "70.7326", "2398.0", "W", "hadec", "{37 19h45m 3h32m} {30 19h15m 4h32m} \
{20 19h10m 5h20m} {10 19h05m 5h30m} {0 18h45m 5h55m} {-20 18h10m 6h55m} {-27 18h30m 7h15m} {-30 18h30m 7h25m} \
{-40 17h55m 8h00m} {-42 17h55m 8h15m} {-45 18h40m 8h30m} {-50 18h40m 8h20m} {-55 17h00m 8h20m} {-60 16h25m 8h20m} \
{-70 16h40m 7h00m}", "-90.0", "37.0", "180.0", "180.0"),
    ("8", "Tarot_Reunion", "-21.198822", "55.410222", "991.0", "E", "altaz", "{0 2} {46 2} {130 14} {146 20} \
{158 26} {197 22} {309 19} {329 9} {349 3} {360 2}", "-89.9", "42.0", "180.5", "179.5"),
    ("3", "Zadko_Australia", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"),
]

# Sub tree Roots
for i in param:
    tele = etree.SubElement(tel, "tele")
    tele.set("id",i[0])
    nom = etree.SubElement(tele, "name")
    nom.text = i[1]
    lat = etree.SubElement(tele, "latitude")
    lat.text = i[2]
    lon = etree.SubElement(tele, "longitude")
    lon.text = i[3]
    alt = etree.SubElement(tele, "altitude")
    alt.text = i[4]
    sen = etree.SubElement(tele, "sens")
    sen.text = i[5]
    hot = etree.SubElement(tele, "horizontype")
    hot.text = i[6]
    hod = etree.SubElement(tele, "horizondef")
    hod.text = i[7]
    lin = etree.SubElement(tele, "limdecmin")
    lin.text = i[8]
    lax = etree.SubElement(tele, "limdecmax")
    lax.text = i[9]
    lis = etree.SubElement(tele, "limharise")
    lis.text = i[10]
    let = etree.SubElement(tele, "limhaset")
    let.text = i[11]

# Sub tree for Plyros optimization options
opt = etree.SubElement(tel, "opt")
etree.SubElement(opt, "nbfield").text = "0"
etree.SubElement(opt, "freetime").text = "0"

file = open("ply_conf.xml", "w", encoding='UTF-8')
file.write(str(etree.tostring(tel, pretty_print = True)))



# End tree
#tree = etree.ElementTree(tel)
# Xml file name
#tree.write("ply_conf.xml")

file.close()
