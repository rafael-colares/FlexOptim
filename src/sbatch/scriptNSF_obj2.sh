#!/bin/bash
 
#===============================================================================
# exemples d'options
 
#SBATCH --partition=normal         # choix de la partition où soumettre le job
#SBATCH --nodes=1
#SBATCH --nodelist=django          # choix du noeud
#SBATCH --ntasks=1                 # nombre processus
#SBATCH --mem=15000                # mémoire nécessaire (par noeud) en Mo
 
#===============================================================================
#exécution du programme (remplacer exe par le nom du programme
# ou la ligne de commande à exécuter)
cd ..
./exec teste_NSF_0_obj2.txt NSF/ 2 3 8 0