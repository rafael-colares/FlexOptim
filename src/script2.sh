#!/bin/bash
 
#===============================================================================
# exemples d'options
 
#SBATCH --partition=court        # choix de la partition où soumettre le job
#SBATCH --nodes=1
#SBATCH --nodelist=kephren        # choix du noeud
#SBATCH --ntasks=1                # nombre processus
#SBATCH --mem=8000                # mémoire nécessaire (par noeud) en Mo
 
#===============================================================================
#exécution du programme (remplacer exe par le nom du programme
# ou la ligne de commande à exécuter)
./exe teste3.txt