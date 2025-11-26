# RaMiLASS

Ce projet est un assembleur de génome en C++ qui utilise des graphes de De Bruijn.

## Auteurs

- Raphaël Ribes
- Mickael Coquerelle
- Loïk Galtier

Le nom de l'exécutable, `ramilass`, est une combinaison des prénoms des auteurs (Raphaël, Mickael, Loïk) et du mot "Assembleur".

## Structure du projet
![workflow](./imgs/workflow.png)
-   `src/`: Contient le code source de la librairie et de l'application.
    -   `graphdbj.cpp`/`.h`: Implémentation du graphe de De Bruijn.
    -   `bitvector.cpp`/`.h`: Implémentation d'un vecteur de bits.
    -   `compare.cpp`/`.h`: Fonctions de comparaison.
    -   `convert.cpp`/`.h`: Fonctions de conversion.
    -   `main.cpp`: Point d'entrée de l'application `ramilass`.
-   `tests/`: Contient les tests unitaires.
-   `build/`: Contient les fichiers de compilation et les exécutables.
-   `pixi.toml`: Fichier de configuration pour l'environnement et les tâches `pixi`.
-   `CMakeLists.txt`: Fichier de configuration pour `CMake`.

## Dépendances

- C++17
- CMake
- GoogleTest
- ragtag
- minia
- quast

## 1. Installation et Compilation

Ce projet utilise `pixi` pour gérer les dépendances et les tâches.

1.  **Installer pixi**

    ```bash
    curl -fsSL https://pixi.sh/install.sh | sh
    ```
2.  **Installer les dépendances du projet**

    ```bash
    pixi install
    ```

3.  **Compiler le projet**

    ```bash
    pixi run build
    ```

    Cela va créer les exécutables dans le dossier `build/`.


## 2. Prétraitement

Notre logiciel prend en entrée des fichiers FASTA et non pas FASTQ. 
Nous utilisons donc une commande `pixi run tofasta <input.fastq> <output.fasta>` pour convertir les fichiers FASTQ en FASTA avant de les utiliser avec `ramilass`.
Par défaut, `pixi run tofasta` sans argument va convertir `Tests_Et_Ref/reads.fastq.fq` en `/Tests_Et_Ref/reads.fasta`.


## 3. Exécuter l'application

```
██████╗  █████╗ ███╗   ███╗██╗██╗      █████╗ ███████╗███████╗
██╔══██╗██╔══██╗████╗ ████║██║██║     ██╔══██╗██╔════╝██╔════╝
██████╔╝███████║██╔████╔██║██║██║     ███████║███████╗███████╗
██╔══██╗██╔══██║██║╚██╔╝██║██║██║     ██╔══██║╚════██║╚════██║
██║  ██║██║  ██║██║ ╚═╝ ██║██║███████╗██║  ██║███████║███████║
╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝

    
Usage: ./build/ramilass <input.fasta> <output_prefix> [OPTIONS]

Arguments obligatoires:
  <input.fasta>        Fichier contenant les lectures (FASTA)
  <output_prefix>      Prefixe pour les fichiers de sortie

Options Generales:
  -k <int>             Taille des k-mers (defaut: 31)
  --fuse               Activer l'etape de fusion des contigs (defaut: inactif)
  --gfa                Exporter le graphe au format GFA (defaut: non)
  --debug              Afficher les temps d'execution et infos detailles
  --help, -h           Afficher ce message

Options de l'Assembleur (GraphDBJ):
  --simplification-passes <int>       Nb max de passes de simplification (defaut: 50)
  --popping-passes <int>        Nb max de passes de suppression de tips et bulle (defaut: 10)

  --overlap-err <dbl>  % d'erreur autorise pour chevauchement (defaut: 0.05)
  --contained-err <dbl>% d'erreur autorise pour inclusion (defaut: 0.02)

  --cov-ratio <dbl>    Ratio de couverture pour bifurcations (defaut: 1)
  --tip-ratio <dbl>    Ratio couverture ancrage/bout (defaut: 2)

  --search-depth <dbl> Facteur de profondeur de recherche (defaut: 20)
```

Pour lancer l'application principale :

```bash
pixi run ramilass <chemin_vers_fichier_fasta> <nom_fichier_sortie>
```

Par exemple :

```bash
pixi run ramilass ./Tests_Et_Ref/reads.fasta contigs_output
```

Pour avoir une execution complète, nous recommandons d'utiliser `pixi run start` qui revient a taper la commande complète :

```bash
./build/ramilass ./Tests_Et_Ref/reads.fasta ramilass --fuse --gfa
```

## BONUS: Image singularity/apptainer

Nous fournissons une image apptainer/singularity pour exécuter l'assembleur dans un environnement isolé.
Cette dernière est construite avec `pixi run build_container` puis exécutée avec `pixi run run_container`.
Cette denière exécute la commande suivante :

```bash
pixi run apptainer exec --bind $INIT_CWD/data:/data pixigularity.sif /app/build/ramilass /data/reads.fasta output --fuse --gfa --debug --max-contig-len 11000 --popping-passes 0
```

Pour faire tourner l'outil depuis l'image avec vos propres paramètres, vous pouvez utiliser la commande suivante :

```bash
pixi run apptainer exec --bind path/to/your/data:/data pixigularity.sif /app/build/ramilass /data/your_reads.fasta [OPTIONS]
```