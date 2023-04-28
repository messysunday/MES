#include "funkcje.h"
#include <cmath>

using namespace std;

double f2d(double x)
{
	double wynik = 2.0 * pow(x, 2.0) + 3.0 * x + 8;
	return wynik;
}

double f3d(double x, double y)
{
	double wynik = -2.0 * y * pow(x, 2.0) + 2.0 * x * y + 4;
	return wynik;
}

double f2d2(double x)
{
	return pow(x, 2) - 3.0 * x + 6.0;
}


