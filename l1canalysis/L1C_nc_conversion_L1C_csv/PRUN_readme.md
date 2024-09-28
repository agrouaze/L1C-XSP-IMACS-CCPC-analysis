Lancer un job PRUN pour génération de dataframe à partir de fichiers Sentinel1 L1C : 


1. Créer et sauvegarder un listing des fichiers à traiter au format .txt. Utiliser le notebook 'lister_L1C.ipynb'


2. Configurer le script 'convert_l1c.py' en changeant le 'root_savepath' et les variables à séléctionner selon le besoin. (appelle le fichier 'l1c_b07_conversion')


3. Configurer le fichier 'convert_l1c.pbs' où il faut mettre l'environnement nécessaire au fonctionnement du script python, et spécifier la ressource (temps, RAM) (appelle le fichier 'convert_l1c.py').


4. rendre le script exécutable : 'chmod +x /home1/datahome/toto/my_awsome_script.pbs'


5. Configurer le fichier 'exe_convert_l1c_bash' (appelle convert_l1c.pbs et le listing): changer le nom du job PRUN et le chemin vers le listing

'--split-max-lines 100' découpe ton listing en sous listing de 100 fichiers que prun traite ensuite en parallèle. 
Pour la ressource en mémoire : 'Diminue la mémoire je pense aussi (mem=3G), passe le à 500m. L'idéal c'est de mettre exactement ce qu'il te faut au max de mémoire pour un fichier + une marge de 10%.' 
Pour savoir combien de mémoire est nécessaire au maximum pour un fichier : 'il y a la lib psutil. Après si tu n'as pas beaucoup modifié ma fonction, 500m ça doit suffire'. 


6. Se rendre dans le répertoire du fichier bash et le lancer: 
'bash script.bash'
ou 
'chmod +x script.bash' puis './script.bash'

Arreter de force : ctrl+D
Sur Datarmor : 
qstat -ftx 424011[].datarmor0 | grep Exit_status : montre le exit status de chaque bash, marche avec used.walltime et used.mem 
qstat -ftx 424011[xx].datarmor0 pour voir les stats du job xx en particulier 