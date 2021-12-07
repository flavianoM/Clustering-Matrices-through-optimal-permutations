**** BANDWIDTH CLUSTERING ****

To compile the source code download the files "opp.h" and "opp_bandwidth.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_bandwidth.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the bandwidth of the clustering. Start with p = 0.5. 

The program outputs four files: <bandwidth_r0.50_graph.txt>, <bandwidth_s0.50_graph.txt>, <bandwidth_x0.50_graph.txt>, <bandwidth_b0.50_graph.txt>
  
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

The program outputs four files: <nestedness_r0.50_graph.txt>, <nestedness_s0.50_graph.txt>, <nestedness_x0.50_graph.txt>, <nestedness_b0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'nestedness_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'nestedness_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'nestedness_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'nestedness_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'nestedness_b0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
  
**** BOXED CLUSTERING ****

To compile the source code download the files "opp.h" and "opp_square.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_square.c -lm -O3

To run the executable file type: ./opp <graph.txt> <Q> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).

<Q> is an integer Q = 2,3,4,... counting the number of boxes the matrix will be clustered into. Start with Q = 3. 

The program outputs four files: <square_r3_graph.txt>, <square_s3_graph.txt>, <square_x3_graph.txt>, <square_b3_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'square_r3_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'square_r3_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'square_s3_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'square_x3_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'square_b3_graph.txt' u 1:2:3 w p pt 5 palette

  
**** TRIANGULAR CLUSTERING ****

To compile the source code download the files "opp.h" and "opp_triangular.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_triangular.c -lm -O3

To run the executable file type: ./opp <graph.txt> <Q> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).

<Q> is an integer Q = 2,3,4,... counting the number of triangles the matrix will be clustered into. Start with Q = 3. 

The program outputs four files: <triangular_r3_graph.txt>, <triangular_s3_graph.txt>, <triangular_x3_graph.txt>, <triangular_b3_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'triangular_r3_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'triangular_r3_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'triangular_s3_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'triangular_x3_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'triangular_b3_graph.txt' u 1:2:3 w p pt 5 palette


