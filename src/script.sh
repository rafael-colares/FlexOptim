#!/bin/bash
 
#===============================================================================
# exemples d'options
 
#SBATCH --partition=normal    # choix de la partition où soumettre le job
#SBATCH --nodelist=django     # choix du noeud
#SBATCH --mem=8000            # mémoire nécessaire (par noeud) en Mo
 
#===============================================================================
#exécution du programme (remplacer exe par le nom du programme
# ou la ligne de commande à exécuter)
./exe teste3.txt