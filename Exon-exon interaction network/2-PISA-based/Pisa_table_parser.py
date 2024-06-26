#!/usr/bin/python3
"""
Code for parsing pisa xml files and saving data in a csv file.

  How to use
  ----------

First you need to have the interfacetable.xml and each hydrogenbond.xml and saltbridge.xml files
given by PDBePISA web server.

Then you can run the script with the following command :

    python Pisa_xml_parser.py interfacetable.xml

Note that the hydrogen bond and slat bridge files must be on the same directory as interfacetable

  Author
  ------
    Hocine Meraoun --> modified by Khalique Newaz

"""

import argparse
import pandas as pd
import re
import os.path


def xmlbond_parser(xml_file):
    """
    The function to parse the interaction xml files.
    It haven't been tested on disulfide and covalent bonds.

    Parameters
    ----------
    xml_file : string
        the name of the xml file

    Returns
    -------
    list
    """
    value = []
    lst = []

    if os.path.exists(xml_file):
        with open(xml_file, "r") as f_xml:
            for line in f_xml :
                if line.startswith("<STRUCTURE1>"):
                    value.append(line[12:13])
                    value.append(line[14:27])
                elif line.startswith("<DISTANCE>"):
                    value.append(float(re.split('<|>',line[10:17])[0]))
                elif line.startswith("<STRUCTURE2>"):
                    value.append(line[13:14])
                    value.append(line[15:28])
                    lst.append(value)
                    value = []
    # else:
    #     print("No ",xml_file, "found")

    return(lst)


def interfacetable_parse(xml_file):
    """
    Function to parse interfacetable.xml and calls xmlbond_parser().

    Parameters
    ----------
    xml_file : string
        interface table xml file name

    Returns
    -------
    list
    """
    lst = []
    intern_lst = []

    path = '/'.join(xml_file.split('/')[:-1])
    # print(xml_file)

    with open(xml_file, "r") as f_xml:
        for line in f_xml :
            if line.startswith("<INTERFACENO>"):
                i = int(line[13:15].strip("<"))
                intern_lst.append(i)
                intern_lst.append(xmlbond_parser(path+"/hydrogenbond"+str(i-1)+".xml"))
                intern_lst.append(xmlbond_parser(path+"/saltbridge"+str(i-1)+".xml"))
            elif line.startswith("<INTERFACEAREA>"):
                intern_lst.append(float(line.split('>')[1].split('<')[0]))
            elif line.startswith("<INTERFACEDELTAGPVALUE>"):
                intern_lst.append(float(line.split('>')[1].split('<')[0]))
                lst.append(intern_lst)
                intern_lst = []

    return(lst)

def create_df(lst):
    """
    Function that creates the dataframe to save in a csv file.

    Parameters
    ----------
    lst : list
        list of lists given by interfacetable_parse()

    Returns
    -------
    pandas DataFrame
    """
    data = {'chain1': [], 'res1': [], 
    'distance': [], 'chain2': [], 
    'res2': [], 'interaction type': [], 'ΔiG kcal/mol': [], 
    'ΔiG P-value': []}

    for c in lst:
        for i in c[1]:
            data['chain1'].append(i[0])
            data['res1'].append(i[1])
            data['distance'].append(i[2])
            data['chain2'].append(i[3])
            data['res2'].append(i[4])
            data['interaction type'].append("Hydrogen bond")
            data['ΔiG kcal/mol'].append(c[3])
            data['ΔiG P-value'].append(c[4])
        if c[2] != []:
            for j in c[2]:
                data['chain1'].append(j[0])
                data['res1'].append(j[1])
                data['distance'].append(j[2])
                data['chain2'].append(j[3])
                data['res2'].append(j[4])
                data['interaction type'].append("Salt bridge")
                data['ΔiG kcal/mol'].append(c[3])
                data['ΔiG P-value'].append(c[4])

    return(pd.DataFrame.from_dict(data))


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("xml_file", help="the pisa interface table xml file", type=str)
    PARSER.add_argument("store_dir", help="the folder to store results", type=str)
    PARSER.add_argument("pdbid", help="the pdb ID", type=str)

    ARGS = PARSER.parse_args()

    XML_FILE = ARGS.xml_file
    STR_DIR = ARGS.store_dir
    PDB_ID = ARGS.pdbid

    create_df(interfacetable_parse(XML_FILE)).to_csv(STR_DIR+"/"+PDB_ID+".csv")
