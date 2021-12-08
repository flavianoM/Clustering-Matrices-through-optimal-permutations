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
int f1j(int I, int N, float p){
    float i = I-1;
    float n = N/3-1;
    float r = i/n;

    return (int)(1 + pow(r,p)*n);
}
int f2j(int I, int N, float p){
    float i = I-N/3;
    float n = N/3-1;
    float r = i/n;

    return (int)(N/3 + pow(r,p)*n);
}
int f3j(int I, int N, float p){
    float i = I-(2*N)/3;
    float n = N/3-1;
    float r = i/n;

    return (int)((2*N)/3 + pow(r,p)*n);
}
void tripole_exp(int N, float p, float eps){
    int i, j, k, l, nk;
    int lmin, lmax;
    float w;

    for(i = 1; i <= N; i++){
        for(j = 1; j < N/3; j++){
            w = -(eps * x[i][j]);
            for(l = f1j(j, N, 1./p); l <= f1j(j, N, p); l++){
                for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            w = (w * 0.01);                   /* multiplied by (beta * eta) = 10 * 0.001 = 0.01 */ 
            w += (0.999 * log(x[i][j]));      /* multiplied by (1-eta) = 1 - 0.001 = 0.999 */
            y[i][j] = exp(w);
        }
        for(j = N/3; j < (2*N)/3; j++){
            w = -(eps * x[i][j]);
            for(l = f2j(j, N, 1./p); l <= f2j(j, N, p); l++){
                 for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            w = (w * 0.01);                   /* multiplied by (beta * eta) = 10 * 0.01 = 0.1 */ 
            w += (0.999 * log(x[i][j]));      /* multiplied by (1-eta) = 1 - 0.01 = 0.99 */
            y[i][j] = exp(w);
        }
        for(j = (2*N)/3; j <= N; j++){
            w = -(eps * x[i][j]);
            for(l = f3j(j, N, 1./p); l <= f3j(j, N, p); l++){
                 for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            w = (w * 0.01);                   /* multiplied by (beta * eta) = 10 * 0.01 = 0.1 */ 
            w += (0.999 * log(x[i][j]));      /* multiplied by (1-eta) = 1 - 0.01 = 0.99 */
            y[i][j] = exp(w);
        }
    }
}
/* Compute energy (cost-function) */
float tripole_energy(int N, float p){
    int i, j, k, l, nk;
    int lmin, lmax;
    float w, e, e0;

    for(i = 1; i <= N; i++){
        for(j = 1; j < N/3; j++){
            w = 0;
            for(l = f1j(j, N, 1./p); l <= f1j(j, N, p); l++){
                for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            e += (x[i][j] * w);
        }
        for(j = N/3; j < (2*N)/3; j++){
             w = 0;
            for(l = f2j(j, N, 1./p); l <= f2j(j, N, p); l++){
                 for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            e += (x[i][j] * w);
        }
        for(j = (2*N)/3; j <= N; j++){
             w = 0;
            for(l = f3j(j, N, 1./p); l <= f3j(j, N, p); l++){
                 for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            e += (x[i][j] * w);
        }
    }

    e0 = 0;
    for(j = 1; j < N/3; j++){
        for(l = f1j(j, N, 1./p); l <= f1j(j, N, p); l++){
            e0 += 1.0;
        }
    }
    for(j = N/3; j < (2*N)/3; j++){
        for(l = f2j(j, N, 1./p); l <= f2j(j, N, p); l++){
            e0 += 1.0;
        }
    }
    for(j = (2*N)/3; j <= N; j++){
        for(l = f3j(j, N, 1./p); l <= f3j(j, N, p); l++){
            e0 += 1.0;
        }
    }
    return ( (e0 + deg[0] - 2*e)/(N*N) );
}
/* One step of the OPP-algorithm */
float onestep(int N, float p, float eps){
    int i, j;
    float diff;
    
    tripole_exp(N, p, eps);
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
void solve(int N, float p, const char *res_file){
    int i, j, k, l, m, t;
    float fm, eps, DIFF, COST;
    FILE *r_file = fopen(res_file, "w");

    init_all(N);
    t = 0;
    for(m = 0; m < M; m++){
        fm = (1 - ((float)m)/(M-1));
        eps = -(fm * pow(M, alpha));

        do{
            DIFF = onestep(N, p, eps);
            COST = tripole_energy(N, p);
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

void print_compression_rate(int N, float p){
    int i, j;
    float M0 = 0, M1 = 0, Mout;

    for(i = 1; i < N/3; i++){
        for(j = f1j(i, N, 1./p); j <= f1j(i, N, p); j++){
            if(s[i][j] == 1){
                M1 += 1;
                s[i][j] = 2;
            }
            if(s[i][j] == 0){
                M0 += 1;
            }
        }
    }
    for(i = N/3; i < (2*N)/3; i++){
        for(j = f2j(i, N, 1./p); j <= f2j(i, N, p); j++){
            if(s[i][j] == 1){
                M1 += 1;
                s[i][j] = 2;
            }
            if(s[i][j] == 0){
                M0 += 1;
            }
        }
    }
    for(i = (2*N)/3; i <= N; i++){
        for(j = f3j(i, N, 1./p); j <= f3j(i, N, p); j++){
            if(s[i][j] == 1){
                M1 += 1;
                s[i][j] = 2;
            }
            if(s[i][j] == 0){
                M0 += 1;
            }
        }
    }
    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            if(s[i][j] != s[j][i]){
                s[i][j] = s[j][i];
            }
        }
    }
    Mout = deg[0] - M1;

    printf("p = %f\n", p);
    printf("Density(filter) = %f\n", (M0+M1)/(N*N));
    printf("Compression rate = %f\n", M1/deg[0]);
}

/* MAIN PROGRAM starts HERE */
int main(int argc, char *argv[]){
    int N;
    float p;
    const char *network;
    char xfile[128], bfile[128], sfile[128], rfile[128];
    
    network = argv[1];
    p = atof(argv[2]);

    sprintf(xfile, "tripole_x%.2f_%s", p, network);
    sprintf(bfile, "tripole_b%.2f_%s", p, network);
    sprintf(sfile, "tripole_s%.2f_%s", p, network);
    sprintf(rfile, "tripole_r%.2f_%s", p, network);

    N = make_network(network);
    printf("N = %d\n", N);
    printf("2M = %d\n\n", deg[0]);

    make_space(N);
    solve(N, p, rfile);

    print_compression_rate(N, p);
    print_matrix_on_file(x, xfile, N);
    print_matrix_on_file(b, bfile, N);
    print_matrix_on_file(s, sfile, N);

    printf("\nExperiment on file %s completed.\n", network);
    free_space(N);

    return 0;
}

    