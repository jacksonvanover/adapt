#include <stdio.h>
#include <approx.h> 
#include <math.h>

/*
 * Precimonious output labels located in comments beside variable declarations
 */

#define ABS(x) ( ((x) < 0.0) ? (-(x)) : (x) )

#define PI  3.1415926535897932L
#define ANS 5.795776322412856L // reference answer for endpoint = PI
#define N 1000000

double endpoint;
double h;

double t1 = 0.0;
double t2;
double t3;

double s1 = 0.0;

double d1 = 1.0;
double d2;

double fun (double x)
{
    d2 = d1;    // also d1 in original
    t3 = x;     // also t1 in original

    int k;
    for (k = 1; k <= 5; k+=1)
    {
        d2 = 2.0 * d2;
#pragma approx snapshot(out) out(t3) label("t3")
        t3 = t3 + sin (d2 * x) / d2;
    }
    return t3;
}

void do_fun ()
{
    int i;
    for (i = 1; i <= N; i+=1)
    {
#pragma approx snapshot(out) out(t2) label("t2")
        t2 = fun (i * h);
#pragma approx snapshot(out) out(s1) label("s1")
        s1 = s1 + sqrt (h * h + (t2 - t1) * (t2 - t1));
#pragma approx snapshot(out) out(t1) label("t1")
        t1 = t2;
    }
}

int main (int argc, char **argv)
{
    if (argc > 1){
        endpoint = atof(argv[1]);
    }
    else{
        endpoint = PI;
    }

    h  = endpoint / (double)N;

    HPACRegisterApplicationInput(&endpoint, sizeof(double), "endpoint", HDOUBLE);

    do_fun();

    HPACRegisterApplicationOutput(&s1, sizeof(double), "arclength", HDOUBLE);

    printf("intervals : %d\n", N);
    printf("endpoint  : %f\n", endpoint);
    printf("arclength : %f\n", s1);

    return 0;
}