/**
 * sum2pi_x
 *
 * CRAFT demo app. Calculates pi*x in a computationally-heavy way that
 * demonstrates how to use CRAFT without being too complicated.
 *
 */

#include <stdio.h>
#include <math.h>
#include <approx.h> 

/* macros */
#define ABS(x) ( ((x) < 0.0) ? (-(x)) : (x) )

/* constants */
#define PI     3.1415926535897932384626433832795
#define EPS    5e-7

/* loop  iterations; OUTER is X */
#define INNER    25

int main(int argc, char **argv)
{
    double sum = 0.0;
    double tmp;
    double acc;
    int i, j, OUTER;

    if (argc > 1){
        OUTER = atoi(argv[1]);
    }
    else{
        OUTER = 2000;
    }

  HPACRegisterApplicationInput(&OUTER, sizeof(int), "factor", HINT);

    for (i=0; i<OUTER; i++) {
        acc = 0.0;
        for (j=1; j<INNER; j++) {

            /* accumulatively calculate pi */
#pragma approx snapshot(out) out(tmp) label("tmp")
            tmp = PI / pow(2.0, j);
#pragma approx snapshot(out) out(acc) label("acc")
            acc = acc + tmp;
        }
#pragma approx snapshot(out) out(sum) label("sum")
        sum = sum + acc;
    }

  HPACRegisterApplicationOutput(&sum, sizeof(double), "product", HDOUBLE);

    double answer = (double)OUTER * PI;             /* correct answer */
    double diff = (double)answer-(double)sum;
    double error = ABS(diff);

    if ((double)error < (double)EPS*answer) {
        printf("SUM2PI_X - SUCCESSFUL!\n");
    } else {
        printf("SUM2PI_X - FAILED!!!\n");
    }

    return 0;
}

