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

int fj(int I, int N, float p){
    int i = I-1;
    int n = N-1;
    float q = 1./p;

    return (int)(1 + pow(i,q)*pow(n,1-q));
}
void bandwidth_exp(int N, float p, float eps){
    int i, j, k, l, nk;
    int lmin, lmax;
    float w;

    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            w = -(eps * x[i][j]);
            lmin = fj(j, N, p);
            lmax = fj(j, N, 1./p);
            for(l = lmin; l <= lmax; l++){
                for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            w = (w * 0.01);                   /* multiplied by (beta * eta) = 10 * 0.001 = 0.01 */ 
            w += (0.999 * log(x[i][j]));      /* multiplied by (1-eta) = 1 - 0.001 = 0.999 */
            y[i][j] = exp(w);
        }
    }
}
/* Compute energy (cost-function) */
float bandwidth_energy(int N, float p){
    int i, j, k, l, nk;
    int lmin, lmax;
    float w, e, e0;

    e = 0;
    for(i = 1; i <= N; i++){
        for(j = 1; j <= N; j++){
            w = 0;
            lmin = fj(j, N, p);
            lmax = fj(j, N, 1./p);
            for(l = lmin; l <= lmax; l++){
                for(k = 1; k <= deg[i]; k++){
                    nk = blist[i][k];
                    w += x[nk][l];
                }
            }
            e += (x[i][j] * w);
        }
    }
    e0 = 0;
    for(j = 1; j <= N; j++){
        lmin = fj(j, N, p);
        lmax = fj(j, N, 1./p);
        for(l = lmin; l <= lmax; l++){
            e0 += 1.0;
        }
    }
    return ( (e0 + deg[0] - 2*e)/(N*N) );
}
/* One step of the OPP-algorithm */
float onestep(int N, float p, float eps){
    int i, j;
    float diff;
    
    bandwidth_exp(N, p, eps);
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
            COST = bandwidth_energy(N, p);
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
            s[i][j] += b[T[i]][T[j]]; /* clustered matrix */
        }
    }
    fclose(r_file);
}

void print_compression_rate(int N, float p){
    int i, j;
    float M0 = 0, M1 = 0, Mout;

    for(i = 1; i <= N; i++){
        for(j = fj(i, N, p); j <= fj(i, N, 1./p); j++){
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

//    printf("Density(A) = %f\n", (1-p)/(1+p));
}

/* MAIN PROGRAM starts HERE */
int main(int argc, char *argv[]){
    int N;
    float p;
    const char *network;
    char xfile[128], bfile[128], sfile[128], rfile[128];
    
    network = argv[1];
    p = atof(argv[2]);

    sprintf(xfile, "bandwidth_x%.2f_%s", p, network);
    sprintf(bfile, "bandwidth_b%.2f_%s", p, network);
    sprintf(sfile, "bandwidth_s%.2f_%s", p, network);
    sprintf(rfile, "bandwidth_r%.2f_%s", p, network);

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

    