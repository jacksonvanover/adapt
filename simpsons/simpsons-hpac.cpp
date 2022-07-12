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