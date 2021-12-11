//Created in 2020
//Author: Flaviano Morone

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define MAX_DEGREE 1000
#define alpha 1.0
#define M 8

char line[MAX_DEGREE];

int *S, *T, *deg, **blist;
float **x, **y, **b, **s, *R, *L;

/* Max and min functions for integer inputs */
int imax(int i, int j){
    return (i > j ? i :j);
}
int imin(int i, int j){
    return (i < j ? i :j);
}
/* Allot space for intermediate tools I need */
void make_space(int N){
    int i, j;
    S = (int *)calloc(N+1, sizeof(int));
    T = (int *)calloc(N+1, sizeof(int));
    L = (float *)calloc(N+1, sizeof(float));
    R = (float *)calloc(N+1, sizeof(float));
    x = (float **)calloc(N+1, sizeof(float *));
    y = (float **)calloc(N+1, sizeof(float *));
    for(i = 0; i <= N; i++){
        x[i] = (float *)calloc(N+1, sizeof(float));
        y[i] = (float *)calloc(N+1, sizeof(float));
    }
}
/* Allot space for graph matrices*/
int make_network(const char *network){
    int i, j, k, nj, node, N;
    char *start;
    FILE *list;

    list = fopen(network, "r");
    N = 0;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        N++;
    }
    fclose(list);

    deg = (int *)calloc(N+1, sizeof(int));
    b = (float **)calloc(N+1, sizeof(float *));
    s = (float **)calloc(N+1, sizeof(float *));
    for(i = 0; i <= N; i++){
        b[i] = (float *)calloc(N+1, sizeof(float));
        s[i] = (float *)calloc(N+1, sizeof(float));
    }
    list = fopen(network, "r");
    i = 1;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        start = line;
        j = 1;
        while( sscanf(start, "%d%n", &node, &k) == 1) {
            b[i][j] = node;
            deg[i] += (b[i][j] != 0 ? 1 : 0);
            start += k;
			j++;
		}
        deg[0] += deg[i];
		i++;
	}
    fclose(list);

    blist = (int **)calloc(N+1, sizeof(int *));
    for(i = 0; i <= N; i++){
        blist[i] = (int *)calloc(deg[i]+1, sizeof(int));
    }
    for(i = 1; i <= N; i++){
        if(b[i][i]!=0){
            printf("self-loop-detected-on-node: %d\n", i);
        }
    }
    for(i = 1; i <= N; i++){
        nj = 1;
        for(j = 1; j <= N; j++){
            if(b[i][j] != 0){
                blist[i][nj] = j;
                nj++;
            }
        }
    }
    return N;
}
/* Initialize DS-matrix and Sinkhorn Left-vector */
void init_all(int N){
    int i, j;
    float N2;
    srand(time(NULL));
    N2 = 0.5*N; 
    for(i = 1; i <= N; i++){
        L[i] = drand48();
        for(j = 1; j <= N; j++){
            x[i][j] = drand48()/N2;
        }
    }
}
/* Free all space */
void free_space(int N){
    int i, j;
    free(L);
    free(R);
    free(deg);
    for(i = 0; i <= N; i++){
        free(x[i]);
        free(y[i]);
        free(b[i]);
        free(s[i]);
        free(blist[i]);
    }
    free(x);
    free(y);
    free(b);
    free(s);
    free(blist);
}
/* print double-stochastic matrix x */
void print_matrix_on_file(float **A, const char *A_file, int N){
    int i, j;
    FILE *a_file = fopen(A_file, "w");
    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            fprintf(a_file, "%d %d %f\n", i, N-j+1, A[i][j]); 
            fflush(a_file);
        }
    }
    fclose(a_file);
}








