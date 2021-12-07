**** BANDWIDTH CLUSTERING ****

To compile the source code download the files "opp.h" and "opp_bandwidth.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_bandwidth.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the bandwidth of the clustering. Start with p = 0.5. 

The program outputs four files: <bandwidth_r0.50_graph.txt>, <bandwidth_s0.50_graph.txt>, <bandwidth_r0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'bandwidth_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'bandwidth_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'bandwidth_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'bandwidth_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'bandwidth_b0.50_graph.txt' u 1:2:3 w p pt 5 palette


**** NESTEDNESS CLUSTERING ****

To compile the source code download the files "opp.h" and "opp_nestedness.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_nestedness.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).

<p> is a float between 0 and 1, measuring the curvature of the nested region of the matrix. The smaller the <p> the higher is the nestedness. Start with p = 0.5. 

The program outputs four files: <nestedness_r0.50_graph.txt>, <nestedness_s0.50_graph.txt>, <nestedness_r0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'nestedness_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'nestedness_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'nestedness_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'nestedness_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'nestedness_b0.50_graph.txt' u 1:2:3 w p pt 5 palette

