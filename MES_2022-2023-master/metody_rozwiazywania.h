#include <iostream>
#include <vector>
#include "struktury.h"

using namespace std;

#ifndef METODY_ROZWIAZYWANIA_H
#define METODY_ROZWIAZYWANIA_H

double kwadratury2d2(double (f)(double), double pocz, double kon);

double kwadratury2d3(double (f)(double));

void skaluj_przedzial_2d2(double pocz, double kon, vector <double>& pkt);

double kwadraturyGaussa(double (f)(double, double), Element4 & element4);


#endif
