//Created on 2021
//Author: Flaviano Morone

**** BANDWIDTH CLUSTERING **** Figure1(a)

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


**** NESTEDNESS CLUSTERING **** Figure1(b)

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
  
  
**** BOXED CLUSTERING **** Figure1(c)

To compile the source code download the files "opp.h" and "opp_square.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_square.c -lm -O3

To run the executable file type: ./opp <graph.txt> <Q> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).

<Q> is an integer Q = 2,3,4,... counting the number of boxes the matrix will be clustered into. Start with Q = 5. 

The program outputs four files: <square_r5_graph.txt>, <square_s5_graph.txt>, <square_x5_graph.txt>, <square_b5_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'square_r5_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'square_r5_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'square_s5_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'square_x5_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'square_b5_graph.txt' u 1:2:3 w p pt 5 palette

  
**** TRIANGULAR CLUSTERING **** Figure1(d)

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
  

**** DIPOLE CLUSTERING **** Figure1(e)

To compile the source code download the files "opp.h" and "opp_dipole.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_dipole.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the transversal width of the dipole. Start with p = 0.5. 

The program outputs four files: <dipole_r0.50_graph.txt>, <dipole_s0.50_graph.txt>, <dipole_x0.50_graph.txt>, <dipole_b0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'dipole_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'dipole_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'dipole_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'dipole_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'dipole_b0.50_graph.txt' u 1:2:3 w p pt 5 palette

  
**** TRIPOLE CLUSTERING **** Figure1(f)

To compile the source code download the files "opp.h" and "opp_tripole.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_tripole.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the transversal width of the tripole. Start with p = 0.5. 

The program outputs four files: <tripole_r0.50_graph.txt>, <tripole_s0.50_graph.txt>, <tripole_x0.50_graph.txt>, <tripole_b0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'tripole_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'tripole_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'tripole_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'tripole_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'tripole_b0.50_graph.txt' u 1:2:3 w p pt 5 palette

  
**** QUADRUPOLE CLUSTERING **** Figure1(g)

To compile the source code download the files "opp.h" and "opp_quadrupole.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_quadrupole.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the transversal width of the quadrupole components. Start with p = 0.5. 

The program outputs four files: <quadrupole_r0.50_graph.txt>, <quadrupole_s0.50_graph.txt>, <quadrupole_x0.50_graph.txt>, <quadrupole_b0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'quadrupole_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'quadrupole_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'quadrupole_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'quadrupole_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'quadrupole_b0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
  
**** PENTAPOLE CLUSTERING **** Figure1(h)

To compile the source code download the files "opp.h" and "opp_pentapole.c" and put them in the same folder. 

Open a terminal and type: gcc -o opp opp_pentapole.c -lm -O3

To run the executable file type: ./opp <graph.txt> <p> 

<graph.txt> is the input file containing the graph formatted as an adjacency matrix. Start with a graph with N < 100 nodes (you should be able to run the algorithm on graphs up to N ~ O(10^3) nodes).
                                                                                                                    
<p> is a float between 0 and 1, measuring the transversal width of the pentapole components. Start with p = 0.5. 

The program outputs four files: <pentapole_r0.50_graph.txt>, <pentapole_s0.50_graph.txt>, <pentapole_x0.50_graph.txt>, <pentapole_b0.50_graph.txt>
  
To plot the energy (cost) function open gnuplot and type: plot 'pentapole_r0.50_graph.txt' u 2:3 w dots

To plot the error function open gnuplot and type: plot 'pentapole_r0.50_graph.txt' u 2:4 w dots
  
To plot the clustered graph open gnuplot and type: set view map, splot 'pentapole_s0.50_graph.txt' u 1:2:3 w p pt 5 palette
  
To plot the optimal permutation matrix open gnuplot and type: set view map, splot 'pentapole_x0.50_graph.txt' u 1:2:3 w p pt 5 palette

To plot the input graph open gnuplot and type: set view map, splot 'pentapole_b0.50_graph.txt' u 1:2:3 w p pt 5 palette


  

  
  
For more information, visit https://theconceptron.com/assets/opp2.pdf


