#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <malloc.h>

double cpuSecond()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ((double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9);
}

void matrix_vector_product(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
            c[i] += a[i * n + j] * b[j];
    }
}

/*
    matrix_vector_product_omp: Compute matrix-vector product c[m] = a[m][n] * b[n]
*/
// void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n)
// {
// #pragma omp parallel num_threads(8)
    // {
        // int nthreads = omp_get_num_threads();
        // int threadid = omp_get_thread_num();
        // int items_per_thread = m / nthreads;
        // int lb = threadid * items_per_thread;
        // int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        // for (int i = lb; i <= ub; i++)
        // {
            // c[i] = 0.0;
            // for (int j = 0; j < n; j++)
                // c[i] += a[i * n + j] * b[j];
        // }
    // }
// }

void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n, int numThread)
{
#pragma omp parallel num_threads(numThread)
    {
        //double t = omp_get_wtime();
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++)
        {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
                c[i] += a[i * n + j] * b[j];
        }
        //t = omp_get_wtime() - t;
        //printf("Thread %d items %d [%d - %d], time: %.6f\n", threadid, ub - lb + 1, lb, ub, t);
    }
}


void run_serial(size_t n, size_t m)
{
    double *a, *b, *c;
    a = (double*)malloc(sizeof(double) * m * n);
    b = (double*)malloc(sizeof(double) * n);
    c = (double*)malloc(sizeof(double) * m);

    if (a == NULL || b == NULL || c == NULL)
    {
        free(a);
        free(b);
        free(c);
        printf("Error allocate memory!\n");
        exit(1);
    }
    double t = cpuSecond();
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }

    for (int j = 0; j < n; j++)
        b[j] = j;

    matrix_vector_product(a, b, c, m, n);
    t = cpuSecond() - t;

    printf("Elapsed time (serial): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
}

void run_parallel(size_t n, size_t m, int numThread)
{
    double *a, *b, *c;

    a = (double*)malloc(sizeof(double) * m * n);
    b = (double*)malloc(sizeof(double) * n);
    c = (double*)malloc(sizeof(double) * m);

    if (a == NULL || b == NULL || c == NULL)
    {
        free(a);
        free(b);
        free(c);
        printf("Error allocate memory!\n");
        exit(1);
    }
    double t = cpuSecond();
    #pragma omp parallel num_threads(numThread)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++) {
            for (int j = 0; j < n; j++)
                a[i * n + j] = i + j;
                c[i] = 0.0;
            }
        //printf("!Thread %d items %d [%d - %d], time:\n", threadid, ub - lb + 1, lb, ub);
    }
    for (int j = 0; j < n; j++) b[j] = j;
    matrix_vector_product_omp(a, b, c, m, n, numThread);
    t = cpuSecond() - t;

    printf("Elapsed time (parallel): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
}


int main(int argc, char *argv[])
{
    size_t M = 1000;
    size_t N = 1000;
    int numThread = 4;
    if (argc > 1)
        numThread = atoi(argv[1]);
    if (argc > 2)
        M = atoi(argv[2]);
    if (argc > 3)
        N = atoi(argv[3]);
    //run_serial(M, N);
    run_parallel(M, N, numThread);
    return 0;
}