#include <time.h>
#include <stdarg.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <approx.h> 

double pi;

double fun(double xarg) {
  double result;

#pragma approx snapshot(out) out(result) label("result")
  result = sin(pi * xarg);

  return result;
}

int main( int argc, char **argv) {


  const int n = 1000000;
  double a; double b;
  double h; double s1; double x;

  if (argc > 1){
      b = atof(argv[1]);
  }
  else{
      b = 1.0;
  }

  HPACRegisterApplicationInput(&b, sizeof(double), "endpoint", HDOUBLE);

  a = 0.0;
  pi = acos(-1.0);
  h = (b - a) / (2.0 * n);
  s1= 0.0;

  x = a;

#pragma approx snapshot(out) out(s1) label("s1_before")
  s1 = (fun(a) + fun(b));

  for(int l = 0; l < n; l++) { // ITERS before

#pragma approx snapshot(out) out(x) label("x_0_during")
    x += h;

#pragma approx snapshot(out) out(s1) label("s1_0_during")
    s1 = s1 + 4.0 * fun(x);

#pragma approx snapshot(out) out(x) label("x_1_during")
    x = x + h;

#pragma approx snapshot(out) out(s1) label("s1_1_during")
    s1 = s1 + 2.0 * fun(x);
  }

#pragma approx snapshot(out) out(s1) label("s1_after")
  s1 = s1 * h * pi / 3.0;

  HPACRegisterApplicationOutput(&s1, sizeof(double), "area", HDOUBLE);
  printf("ans: %.6e\n", s1);

  return 0;
}

/*
  Replace variable b           max error introduced: 0.000000e+00  count: 1           totalerr: 0.000000e+00
  Replace variable a           max error introduced: 0.000000e+00  count: 1           totalerr: 0.000000e+00
  Replace variable h           max error introduced: 4.152677e-15  count: 1           totalerr: 4.152677e-15
  Replace variable pi          max error introduced: 9.154282e-14  count: 1           totalerr: 9.569550e-14
  Replace variable result      max error introduced: 2.967209e-11  count: 2000002     totalerr: 2.976779e-11
  DO NOT replace   x           max error introduced: 2.397519e-07  count: 2000001     totalerr: 2.397817e-07
  DO NOT replace   s1          max error introduced: 8.160601e-05  count: 2000002     totalerr: 8.184579e-05
*/