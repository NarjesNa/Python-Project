#!/usr/bin/env python
#-*- coding : utf8 -*-

# Usage : Main.py "fichier_d'entree.pdb"

from math import sqrt
import string
import re
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from Fonctions import *

# On parse notre fichier PDB entre en argument"
dico_PDB = ParserPDB(sys.argv[1])

# variables pour le modele global
nb_modeles = len(dico_PDB) # on recupere le nombre de conformations presentes dans le fichier PDB


CDM = centredemasse(dico_PDB)


# Partie globale
Liste_RMSD=(RMSDglob(dico_PDB))

Liste_RayonGir=(gir_global(dico_PDB,CDM))
#print Liste_RayonGir
#print range (0,len(Liste_RayonGir))

# Partie locale
RMSD_loc = (RMSDlocal(dico_PDB))

# Calcul RMSD moyen sur tous les modeles
resultats=[]
resultats= np.zeros(89)

for i in range(0,89): # boucle sur les residus
    for y in range(0,nb_modeles):
        resultats[i]+= RMSD_loc[i+89*y]
    resultats[i]=resultats[i]/nb_modeles


RMSD_moyen=resultats # moyenne de tous les RMSD de tous modeles pour chaque residu

# Les deux residus particuliers sont
ARG=[]
LYS=[]

for cle in range(0,len(dico_PDB)):
    for i in range(0,89):
        if i==15:
            LYS.append(Liste_RMSD[cle])
print LYS

for cle in range(0,len(dico_PDB)):
    for i in range(0,89):
        if i==17:
            ARG.append(Liste_RMSD[cle])
print ARG

# Distance moyenne
vecteur_distance=distance(dico_PDB, CDM)
dist_moy=np.zeros(89)
for i in range(0,89):
    for y in range(0,nb_modeles):
        dist_moy[i]+=vecteur_distance[i+89*y]
    dist_moy[i]=dist_moy[i]/nb_modeles



# Graphes

# Graphe donnant le RMSD global des modeles avec le rayon giratoire
fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(1,len(Liste_RMSD)),Liste_RMSD[1:])
ax[0].set_ylabel("RMSD")
ax[1].plot(range(1,len(Liste_RMSD)),Liste_RayonGir[1:],range(1,len(Liste_RMSD)),[Liste_RayonGir[0] for x in range(1,len(Liste_RMSD))],'r-')
ax[1].set_ylabel("Rayon giratoire"); ax[1].set_xlabel("Modeles")
fig.suptitle('RMSD', fontsize=12)
fig.text(.5,.5,'Rayon giratoire',fontsize=12,ha='center')


# Graphe donnant les RMSD moyens et les distances moyennes
fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(0,len(RMSD_moyen)),RMSD_moyen)
ax[0].set_ylabel("RMSD")
ax[1].plot(range(0,len(RMSD_moyen)),dist_moy,'r-')
ax[1].set_ylabel("Distance"); ax[1].set_xlabel("AA")
fig.suptitle('RMSD moyen', fontsize=12)
fig.text(.5,.5,'Enfouissement',fontsize=12,ha='center')
plt.show()

# Graphes des acides amines d'interet de l'article
fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(0,len(LYS)), LYS)
ax[0].set_ylabel("RMSD 15 (en Angstrom)")
ax[1].plot(range(0,len(ARG)), ARG)
ax[1].set_ylabel("RMSD 17 (en Angstrom)")
fig.suptitle('15 LYS', fontsize=12)
fig.text(.5,.5,'17 ARG',fontsize=12,ha='center')
plt.show()


##################### Fichiers de sortie

# output pour avoir les RMSD et les rayons de giration

output = open("global_10.txt", "w")
output.write("Rayon giratoire a t = 0 : %f\n\n" % (Liste_RayonGir[0]))
output.write("Modele\tRMSD\t\tRayon giratoire\n\n")

for i in range(1,nb_modeles):
	output.write("%i\t%f\t%f\n" % (i, Liste_RMSD[i],Liste_RayonGir[i]))
output.close()

# output pour avoir la moyenne des RMSD et des distances

output = open ("Moyenne_10.txt", "w")
output.write("AA\tRMSD\tRayon Giratoire\n")

for i in range(len(RMSD_loc[1])):
	output.write("%s\t%s\t%s\n" % (dico2[0][i*5+1], RMSD_moyen[i], dist_moyenne[i]))
output.close()

# output pour avoir toutes les valeurs locales

output = open ("Locales.txt", "w")
output.write("Modele\tAA\tRMSD\tDistance\n")

for i in range(1,len(RMSD_loc)):
	for j in range (len(RMSD_loc[1])):
		output.write("%i\t%i\t%s\t%s\n" % (i,j, RMSD_loc[i][j], RayonGir_loc[i][j]))
output.close()
