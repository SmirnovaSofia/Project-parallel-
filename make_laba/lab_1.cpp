#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <malloc.h>

double e = 0.00001;  //эпсила


double *A;  //матрица коэффициентов
double *b;  //вектор ответов на уравнения
int m = 4, n = 4;  //высота на ширину

//std::<std::<double[]>[]> Ab(new std::<double[]>[n]);

char end = 'f';  //флаг, достигнута указанная точность
double t = 0.1;  //коэффициент приближения

double loss = 0;
int parallel_loop = 0;


double cpuSecond()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ((double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9);
}
int t_prov = 0;
void proverka(double *x, double* x_predict, int numThread, double bhg, double bhgv){
    {
        int sm = 0;
        double r = sqrt(bhg)/sqrt(bhgv);
        if (loss < r && loss != 0 && loss > 0.999999999) {
            t = t/10;
            t_prov = 1;
            loss = 0;
            for(int i = 0; i < n; i++){
                x[i] = 0;
            }
        }
        else{
        if(abs(loss - r) < 0.00001 && sm == 0) {
            t = t*-1;
            sm = 1;
        }
        loss = r;
        if(r < e) end = 't';
        }
    }
}

double* matrix_vector_product_omp(double *x, int numThread, double bhgv)
{
    double *x_predict = (double*)malloc(sizeof(double)*n);
    double bhg = 0;
#pragma omp parallel num_threads(numThread)
    {
        double aaaaaa = 0;
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);

        for (int i = lb; i <= ub; i++){
            x_predict[i] = 0;
            for(int j = 0; j < n; j++){
                x_predict[i] = x_predict[i] + A[i*n+j]*x[j];
            }
            x_predict[i] = x_predict[i]-b[i];
            aaaaaa += x_predict[i]*x_predict[i];
            x_predict[i] = x[i] - t*x_predict[i];
        }
        #pragma omp atomic
        bhg += aaaaaa;
    }
    //printf("%f\n", x_predict[0]);
    proverka(x, x_predict, numThread, bhg, bhgv);
    if(t_prov == 1){
        x_predict = x;
        t_prov = 0;
    }
    return x_predict;
}


double* matrix_vector_product_omp_second(double *x, int numThread, double bhgv)
{

    double *x_predict = (double*)malloc(sizeof(double)*n);
    double bhg = 0;
#pragma omp parallel num_threads(numThread)
    {

        #pragma omp for schedule(dynamic, int(m/(numThread*3))) nowait reduction(+:bhg)
        for (int i = 0; i < m; i++){
            //int nthreads = omp_get_num_threads();
            //printf("%d\n", nthreads);
            x_predict[i] = 0;
            for(int j = 0; j < n; j++){
                x_predict[i] = x_predict[i] + A[i*n+j]*x[j];
            }
            x_predict[i] = x_predict[i]-b[i];
            bhg += x_predict[i]*x_predict[i];
            x_predict[i] = x[i] - t*x_predict[i];
        }
    }

    proverka(x, x_predict, numThread, bhg, bhgv);
    if(t_prov == 1){
        x_predict = x;
        t_prov = 0;
    }
    return x_predict;
}

void run_parallel(int numThread)
{
    double *x = (double*)malloc(sizeof(double) * n);
    A = (double*)malloc(sizeof(double) * m * n);
    b =  (double*)malloc(sizeof(double) * m);
    double bhgv = 0;
    if (A == NULL || b == NULL || x == NULL)
    {
        free(A);
        free(b);
        free(x);
        printf("Error allocate memory!\n");
        exit(1);
    }

    double time = cpuSecond();
    #pragma omp parallel num_threads(numThread)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++) {
            for (int j = 0; j < n; j++){
                if(i == j) A[i * n + j] = 2;
                else A[i*n+j] = 1;
            }
            b[i] = n+1;
            x[i] = b[i]/A[i*n+i];
            #pragma omp atomic
            bhgv += b[i]*b[i];
        }
    }

    if(parallel_loop == 0){
        while(end == 'f'){
            x = matrix_vector_product_omp(x, numThread, bhgv);
        }
    }
    else{
        while(end == 'f'){
            x = matrix_vector_product_omp_second(x, numThread, bhgv);
        }
    }
    
    for(int i = 0; i < n; i++){
        printf("%f ", x[i]);    //выводится вектор ответов
    }
    
    time = cpuSecond() - time; // время работы. начиная с инициализации

    printf("Elapsed time (parallel): %.6f sec.\n", time);
        free(A);
        free(b);
        free(x);
}
int main(int argc, char **argv){

    int numThread = 2;
    if (argc > 1)
        numThread = atoi(argv[1]);
    if (argc > 2){
        m = atoi(argv[2]);
        n = atoi(argv[2]);
    }
    if (argc > 3){
        parallel_loop = atoi(argv[3]);
    }
    run_parallel(numThread);
}