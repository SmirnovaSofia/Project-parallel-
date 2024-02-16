#define _CRT_SECURE_NO_WARNINGS
#define M_PI 3.14159265358979323846
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <typeinfo>
#define N 10000000
#ifdef AHAHA
#define T double
#else
#define T float
#endif

int main()
{
    printf("Enter the start of the period on the x-axis:");
    T c;
    scanf("%lf", &c);

    T d = (M_PI * 2 / N);

    printf("%.9lf\n", d);
    T sum = 0;
    T* arr = (T*)malloc(sizeof(T) * (N + 1));

    for (int i = 0; i < N; i++) {
        arr[i] = sin(c + i * d);
        //printf("%.9lf\n", c);
        sum = sum + arr[i];
    }

    for (int i = 0; i < 10; i++) {
        printf("%.9lf\n", arr[i]);
    }

    printf("Sum = %.9lf, Type = %s", sum, typeid(sum).name());
}

