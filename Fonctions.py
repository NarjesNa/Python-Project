#!/usr/bin/env python
#-*- coding : utf8 -*-


from math import sqrt
import string
import re
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def ParserPDB(PDBfile):

    """ Fonction qui parse un fichier PDB afin d'en extraire le numero de conformation, la liste des residus(ainsi que leurs noms)
        composant la proteine, la liste des atomes constituant les residus, et les coordonnees x, y et z de chaque atome
        Entree : fichier au format PDB
        Sortie : un dictionnaire
    """

    try:
		PDB_file = open(PDBfile,'r')
    except:
		print "L'ouverture a echoue"
		sys.exit(1)

    lines = PDB_file.readlines()

    dico_PDB = {}
    Flag = False

    for line in lines:

		if line[0:5] == "MODEL":
			num_model = int(string.strip(line[10:14]))
			dico_PDB[num_model] = {}
			num_model = int(num_model)
			dico_PDB[num_model]["Liste_Chaines"] = []  # numero du residu
			dico_PDB[num_model]["ResList"] = []     # noms des residus
			# le dictionnaire a la cle "Liste_Chaines" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM":
			chain = line[24:27]
												   # on ne selectionne que les lignes qui contiennent des ATOM
			if chain not in dico_PDB[num_model]["Liste_Chaines"]:
				dico_PDB[num_model]["Liste_Chaines"].append(chain)
				dico_PDB[num_model]["ResList"].append(string.strip(line[17:20]))
				dico_PDB[num_model][chain] = {}
                                                    # pour la cle number ayant pour cle "resname"
			if dico_PDB[num_model][chain].has_key("resname") == False:
				dico_PDB[num_model][chain]["resname"] = string.strip(line[17:20])
				dico_PDB[num_model][chain]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atom = string.strip(line[12:16])

			dico_PDB[num_model][chain]["atomlist"].append(atom) # ajout de l'atome a la liste

			dico_PDB[num_model][chain][atom] = {}    # cree un dictionnaire dans dicPBD[chain][number]

			dico_PDB[num_model][chain][atom]["x"] = float(line[30:38])
			dico_PDB[num_model][chain][atom]["y"] = float(line[38:46])
			dico_PDB[num_model][chain][atom]["z"] = float(line[46:54])
			dico_PDB[num_model][chain][atom]["id"] = line[7:11].strip()

    PDB_file.close()
    return(dico_PDB)


def centredemasse(dico_PDB) :

    """ Fonction qui calcule les coordonnees du centre de masse d'une proteine a partir du fichier PDB parse
        Entree : un dictionnaire (fichier PDB d'origine parse)
        Sortie : un dictionnaire contenant les coordonnees (x, y, z) du centre de masse de la proteine
    """

    dict_coord = {}

    for chain in range(0, len(dico_PDB)):
        dict_coord[chain]={}
        x=y=z=cpt=0.0
        ResList=dico_PDB[chain]["Liste_Chaines"]
        for res in ResList:
            atomlist = dico_PDB[chain][res]["atomlist"]
            for atom in atomlist:
                x += dico_PDB[chain][res][atom]["x"]
                y += dico_PDB[chain][res][atom]["y"]
                z += dico_PDB[chain][res][atom]["z"]
                cpt+=1
        dict_coord[chain]["x"]= float(x)/cpt
        dict_coord[chain]["y"]= float(y)/cpt
        dict_coord[chain]["z"]= float(z)/cpt

    return dict_coord


def RMSDglob (dico_PDB):

    """Fonction qui calcule le RMSD global d'une proteine
       Entree : un dictionnaire
       Sortie : une liste
    """

    result_glob = []

    for chaine in dico_PDB.keys() :
		count = somme = 0
		reslist = dico_PDB[chaine]["Liste_Chaines"]
		for res in reslist :
			atomlist = dico_PDB[chaine][res]["atomlist"] # boucle sur la liste d'atomes
			for atom in atomlist :
				measure = ((dico_PDB[chaine][res][atom]["x"] - dico_PDB[0][res][atom]["x"]))**2
				+ ((dico_PDB[chaine][res][atom]["y"] - dico_PDB[0][res][atom]["y"]))**2
				+((dico_PDB[chaine][res][atom]["z"] - dico_PDB[0][res][atom]["z"]))**2
				somme+=measure
				count+=1
		RMSD=sqrt(somme/count)
		result_glob.append(RMSD)

    return(result_glob)

def RMSDlocal (dico_PDB):

        """Fonction qui calcule le RMSD local d'une proteine
           Entree : un dictionnaire
           Sortie : une liste
       """
	result_loc = []

	for chaine in dico_PDB.keys() :
		reslist = dico_PDB[chaine]["Liste_Chaines"]
		for res in reslist :
			count = somme = 0
			atomlist = dico_PDB[chaine][res]["atomlist"] # boucle sur la liste d'atomes
			for atom in atomlist :
				measure = ((dico_PDB[chaine][res][atom]["x"] - dico_PDB[0][res][atom]["x"]))**2
				+ ((dico_PDB[chaine][res][atom]["y"] - dico_PDB[0][res][atom]["y"]))**2
				+((dico_PDB[chaine][res][atom]["z"] - dico_PDB[0][res][atom]["z"]))**2
				somme+=measure
				count+=1
			RMSD=sqrt(somme/count)
			result_loc.append(RMSD)

 	return (result_loc)


def rayon1(xA, x0, yA, y0, zA, z0):
    """Fonction qui calcule une distance entre 2 points
    """
    return sqrt((xA-x0)**2+(yA-y0)**2+(zA-z0)**2)

def gir_global(dico_PDB, CDM):
    """Fonction qui calcule le rayon giratoire entre l'atome le plus eloigne et
        par rapport au centre de masse
        Entree : un dictionnaire, centre de masse
        Sortie : une liste
    """

    lst_memoire=[]

    for chain in dico_PDB.keys():
        ResList=dico_PDB[chain]["Liste_Chaines"]
        rayon_max=0
        for res in ResList:
            atomlist = dico_PDB[chain][res]["atomlist"]
            for atom in atomlist:
                xA=dico_PDB[chain][res][atom]["x"]
                yA=dico_PDB[chain][res][atom]["y"]
                zA=dico_PDB[chain][res][atom]["z"]
                rayon=rayon1(xA, CDM[chain]['x'], yA, CDM[chain]['y'], zA, CDM[chain]['z'])
                if rayon>=rayon_max :
                    rayon_max=rayon
        lst_memoire.append(rayon_max)

    return lst_memoire


def distance(dico_PDB, CDM):

    """ Fonction qui calcule une distance entre les residus et le centre de masse
        Entree : un dictionnaire, centre de masse
        Sortie : une liste de distances
    """

    dist =[]

    for chaine in dico_PDB.keys() :
        reslist = dico_PDB[chaine]["Liste_Chaines"]
        for res in reslist :
            atomlist = dico_PDB[chaine][res]["atomlist"]
            for atom in atomlist:
                measure=sqrt(((dico_PDB[chaine][res][atom]["x"] - CDM[chaine]["x"]))**2 + ((dico_PDB[chaine][res][atom]["y"] - CDM[chaine]["y"]))**2 + ((dico_PDB[chaine][res][atom]["z"] - CDM[chaine]["z"]))**2)
            dist.append(measure)
    return dist
