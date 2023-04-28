#include "metody_rozwiazywania.h"

using namespace std;

double kwadratury2d2(double(f)(double), double pocz, double kon)
{
	vector <double> x = { (-1.0) * sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) }; //wspolrzedne x wezlow calkowania (dla 2 punktow tylko 2 wspolrzedne)
	vector <double> waga = { 1.0, 1.0 }; //wagi
	skaluj_przedzial_2d2(pocz, kon, x);
	double detJ = (kon - pocz) / 2;
	double wynik = 0.0; //zmienna na wynik calkowania
	for (int i = 0; i < 2; ++i) //dla kazdego punktu w pionie
	{
		wynik += f(x[i]) * waga[i]; //oblicz powierzchnie - wysokosc figury (wartosc funkcji f w punkcie calkowania) * waga i(czyli dl podstawy)
	}
	return detJ * wynik;
}

double kwadratury2d3(double(f)(double))
{
	vector <double> x = { (-1.0) * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) }; ////wspolrzedne x wezlow calkowania (dla 9 punktow tylko 3 rozne wspolrzedne)
	vector <double> waga = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 }; //wagi
	double wynik = 0.0;
	//tutaj pola podstaw bryl sa rozne - trzeba uwazac ktore pole z jaka wysokascia
	for (int i = 0; i < 3; ++i)
	{
		wynik += f(x[i]) * waga[i];//oblicz powierzchnie - wysokosc figury (wartosc funkcji f w punkcie calkowania) * waga i(czyli dl podstawy)
	}
	return wynik;
}

/*
double kwadraturyGaussa(double(f)(double, double), Element4 & element4)
{
	double wynik = 0.0;
	for (int i = 0; i < element4.rzad; ++i)
	{
		for (int j = 0; j < element4.rzad; ++j)
		{
			wynik += f(element4.ksi[i], element4.eta[j]) * element4.waga[i] * element4.waga[j];//oblicz objetosc - wysokosc figury (wartosc funkcji f w punkcie calkowania) * waga i * waga j (czyli pole podstawy)
		}
	}
	return wynik;
}
*/

void skaluj_przedzial_2d2(double pocz, double kon, vector <double>& pkt) //do schematu 2 pkt
{
	pkt[0] = (1.0 - pkt[0]) / 2.0 * pocz + (pkt[0] + 1.0) / 2.0 * kon;
	pkt[1] = (1.0 - pkt[1]) / 2.0 * pocz + (pkt[1] + 1.0) / 2.0 * kon;
}
