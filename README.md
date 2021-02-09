# MOS2.2_Informatique_Graphique

Projet de Ray tracing réalisé pendant le cours MOS 2.2 d'Infromatique Graphique à l'Ecole Centrale de Lyon

Enseignant : [Nicolas Bonnel](https://perso.liris.cnrs.fr/nicolas.bonneel/teaching.html)

## Table of Contents

 * [Description](#Description)
 * [Prerequis](#Prerequis)
 * [Utilisation](#Utilisation)

## Description

## Prerequis

1. Nécessite un compileur C++ type g++ que l'on peut obtenir en installant [MinGW](https://sourceforge.net/projects/mingw-w64/)
2. Ajouter le chemin de mingW au variable d'envirronement windows (C:\MinGW\bin)

## Utilisation

* Executer test.cpp avec openMP (Plus rapide, utilisation de tout les coeurs logiques de votre processeurs)
```sh
g++ -fopenmp test.cpp -o test
```

* Executer test.cpp sans openMP 
```sh
g++ test.cpp -o test
```
