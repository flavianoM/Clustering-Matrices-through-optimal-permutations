//Created on 2020
//Author: Flaviano Morone

#include"opp.h"

/* Normalize matrix A to a doubly-stochastic form */
float row_col_normal(float **A, int N){
    int i, j;
    float sum, diff;

    do{
        for(j = 1; j <= N; j++){
            sum = 0;
            for(i = 1; i <= N; i++){
                sum += (L[i] * A[i][j]);
            }
            R[j] = 1./sum;
        }
        diff = 0;
        for(i = 1; i <= N; i++){
            sum = 0;
            for(j = 1; j <= N; j++){
                sum += (A[i][j] * R[j]);
            }
            sum = 1.0/sum;
            diff = fmaxf(diff, fabs(L[i] - sum));
            L[i] = sum;
        }
    } while(diff > 1e-5);

    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            A[i][j] = L[i] * A[i][j] * R[j];     /* Normalize A into DS-form */
        }
    }
    return diff;
}
void square_exp(int N, int Q, float eps){
    int i, j, k, l, nk, q;
    int lmin, lmax;
    float w, r;
    r = ((float)N)/((float)Q);

    for(i = 1; i <= N; i++){
        for(q = 1; q <= Q; q++){
            for(j = (int)(1 + (q-1)*r); j <= (int)(q*r); j++){
                w = -(eps * x[i][j]);
                lmin = (int)(1 + (q-1)*r);
                lmax = (int)(q*r);
                for(l = lmin; l <= lmax; l++){
                    for(k = 1; k <= deg[i]; k++){
                        nk = blist[i][k];
                        w += x[nk][l];
                    }
                }
                w = (w * 0.01);                   /* multiplied by (beta * eta) = 10 * 0.001 = 0.1 */ 
                w += (0.999 * log(x[i][j]));      /* multiplied by (1-eta) = 1 - 0.001 = 0.999 */
                y[i][j] = exp(w);
            }
        }
    }
}
/* Compute energy (cost-function) */
float square_energy(int N, int Q){
    int i, j, k, l, nk, q;
    int lmin, lmax;
    float w, e = 0, e0 = 0, r;
    r = ((float)N)/((float)Q);

     for(i = 1; i <= N; i++){
        for(q = 1; q <= Q; q++){
            for(j = (int)(1 + (q-1)*r); j <= (int)(q*r); j++){
                w = 0;
                lmin = (int)(1 + (q-1)*r);
                lmax = (int)(q*r);
                for(l = lmin; l <= lmax; l++){
                    for(k = 1; k <= deg[i]; k++){
                        nk = blist[i][k];
                        w += x[nk][l];
                    }
                }
                e += (x[i][j] * w);
            }
        }
     }
    e0 = 0;
    for(j = (int)(1 + (q-1)*r); j <= (int)(q*r); j++){
        lmin = (int)(1 + (q-1)*r);
        lmax = (int)(q*r);
        for(l = lmin; l <= lmax; l++){
            e0 += 1.0;
        }
    }
    return ( (e0 + deg[0] - 2*e)/(N*N) );
}/* One step of the OPP-algorithm */
float onestep(int N, int Q, float eps){
    int i, j;
    float diff;
    
    square_exp(N, Q, eps);
    row_col_normal(y, N);

    diff = 0;
    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            diff = fmaxf(fabs(y[i][j] - x[i][j]), diff);
            x[i][j] = y[i][j];
        }
    }
    return diff;
}
/* Full solution */
void solve(int N, int Q, const char *res_file){
    int i, j, k, l, m, t;
    float fm, eps, DIFF, COST;
    FILE *r_file = fopen(res_file, "w");

    init_all(N);
    t = 0;
    for(m = 0; m < M; m++){
        fm = (1 - ((float)m)/(M-1));
        eps = -(fm * pow(M, alpha));

        do{
            DIFF = onestep(N, Q, eps);
            COST = square_energy(N, Q);
            t++;
            fprintf(r_file, "%f %d %f %f\n", eps, t, COST, DIFF);
            fflush(r_file);
        }while(DIFF > 1e-5);
    }
    for(i = 1; i <= N; i++){
        S[i] = i;  /* Permutation */
        T[i] = i;  /* (Permutation)^-1 */
        for(j = 1; j <= N; j++){
            S[i] = x[i][j] > 0.99 ? j : S[i];
            T[i] = x[j][i] > 0.99 ? j : T[i];
        }
    }
    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            s[i][j] = b[T[i]][T[j]]; /* clustered matrix */
        }
    }
    fclose(r_file);
}

void print_compression_rate(int N, int Q){
    int i, j, jmin, jmax, q;
    float M0 = 0, M1 = 0, Mout, r;

    r = ((float)N)/((float)Q);

    for(q = 1; q <= Q; q++){
        for(i = (int)(1 + (q-1)*r); i <= (int)(q*r); i++){
            jmin = (int)(1 + (q-1)*r);
            jmax = (int)(q*r);
            for(j = jmin; j <= jmax; j++){
                if(s[i][j] == 1){
                    M1 += 1;
                    s[i][j] = 2;
                }
                if(s[i][j] == 0){
                    M0 += 1;
                }
            }
        }
    }
    Mout = deg[0] - M1;

    printf("Q = %d\n", Q);
    printf("Density(filter) = %f\n", (M0+M1)/(N*N));
    printf("Compression rate = %f\n", M1/deg[0]);
}

/* MAIN PROGRAM starts HERE */
int main(int argc, char *argv[]){
    int N, Q;
    const char *network;
    char xfile[128], bfile[128], sfile[128], rfile[128];
    
    network = argv[1];
    Q = atoi(argv[2]);

    sprintf(xfile, "square_x%d_%s", Q, network);
    sprintf(bfile, "square_b%d_%s", Q, network);
    sprintf(sfile, "square_s%d_%s", Q, network);
    sprintf(rfile, "square_r%d_%s", Q, network);

    N = make_network(network);
    printf("N = %d\n", N);
    printf("2M = %d\n\n", deg[0]);

    make_space(N);
    solve(N, Q, rfile);

    print_compression_rate(N, Q);
    print_matrix_on_file(x, xfile, N);
    print_matrix_on_file(b, bfile, N);
    print_matrix_on_file(s, sfile, N);

    printf("\nExperiment on file %s completed.\n", network);
    free_space(N);

    return 0;
}

    
