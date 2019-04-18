# REMD
Changements effectués sur le Dockerfile :

1) Toujours préciser la version de Ubuntu (ou autre distribution). Ici, j'ai pris Ubuntu xemial. Sinon, il prendra la "latest" et celle-ci, par définition, va changer au cours du temps;

2) J'ai enlevé "RUN apt-get upgrade -y" qui (a priori) n'est pas utile;

3) Tu peux étaler une commande sur plusieurs ligne avec un antislash (c'est juste esthétique, ça évite d'avoir plein de "RUN apt-get install");

4) J'ai ajouté l'option -q (quiet) pour apt-get (Note: il existe -qq pour very quiet !);

5) **IMPORTANT** : Il faudrait que tu précises ce que font les paquets (bison ??) et qui a besoin de quoi. Ces informations peuvent s'avérer très utiles par la suite. Profite qu'elles sont encore fraîches dans ta tête !;

6) J'ai déplacé l'installation de vim et wget, dans la liste d'outils que j'installe ensuite et qui sont utiles lorsqu'on lance docker en mode interactif (auto-complétion, coloration syntaxique, screen etc.);

7) C'est un détail, mais je pense qu'il est inutile de faire un chmod, vu que dans le conteneur, tu es toujours root (mais je l'ai laissé, on verra quand les scripts bash seront fonctionnels);

8) J'ai réorganisé les actions : installation de GROMACS (copie, puis make etc.), installation de AMBERTOOLS etc.

9) J'ai aussi réorganisé les fichiers, directement sur GitHub.
Tout ce qui est scripts va dans le répertoire du même nom. Tout ce qui est données va dans data/ et logiciel tiers dans src/ (i.e. acpype/). Ceux qui sont trop lourds sont téléchargés depuis le OwnCloud RPBS. Tu peux changer l'organistaion comme tu préfères, car il n'y a en fait aucune règle. Mais il vaut mieux essayer de ranger un peu pour s'y retrouver ensuite...

9) Enfin, je termine par un nettoyage (autoremove, clean etc.).
