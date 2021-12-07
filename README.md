# The Conceptron by Flaviano Morone

Per compilare il codice utilizzare:
gcc -o opp opp_bandwidth.c -lm -O3

Per lanciare il codice utilizzare:
./opp <grafo.txt> <float p>

<grafo.txt> è il file che contiene il tuo grafo (in formato matrice di adiacenza, non lista). Prova con ~ 100 nodi, ma dovrebbe funzionare fino a N ~ 1000. 
 
<float p> è un valore tra 0 e 1 che misura la larghezza della banda del clustering. Metti 0.5 per cominciare. 

In OUTPUT scrivo 4 file:
FILE r:
Per plottare l'energia fai: p 'bandwidth_r0.50_randnet.txt' u 2:3 w dots 
Per plottare l'errore in convergenza: p 'bandwidth_r0.50_randnet.txt' u 2:4 w dots

FILE s:
Per plottare il grafo clusterizzato: 
set view map

sp 'bandwidth_s0.50_randnet.txt' u 1:2:3 w p pt 5 palette


FILE x:
Per plottare la matrice di permutazione:
set view map

sp 'bandwidth_x0.50_randnet.txt' u 1:2:3 w p pt 5 palette



FILE b:

Per plottare la matrice del grafo originale:
set view map

sp 'bandwidth_b0.50_randnet.txt' u 1:2:3 w p pt 5 palette
