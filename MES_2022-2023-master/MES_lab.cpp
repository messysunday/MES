#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>

#include "struktury.h"
#include "wczytywanie.h"
#include "funkcje.h"

using namespace std;

//wizualizacja w paraview!
//warstwy - od x jakiegoś do jakiegoś maja okreslone parametry, wartosci wsp dac do elementu

int main()
{
	Grid grid;
	Global_Data data;
	//read(grid, data);
	generuj_siatke(grid, data);
	cout << "done\n";
	Element4 element4(2, 2); //to tyko raz
	Equations equations(data, grid);
	wypisz_macierz(equations.temps0);
	string symulacja = "NewKacperSym";
	while (data.simulationTime > data.currentTime)
	{
		wyczysc_macierz(equations.HG);
		wyczysc_macierz(equations.CG);
		wyczysc_macierz(equations.p);
		for (auto i : grid.Elements)
		{
			i.macierzH_macierzC(element4, grid, data);
			i.macierzHBC_wektorP(element4, grid, data, equations);
			equations.agregacjaH(i);
			equations.agregacjaC(i);
			equations.agregacjaP(i);
		}
		zapisz_do_pliku(symulacja + to_string(data.currentTime) + ".vtk", grid, equations);
		equations.solve();
		data.currentTime += data.simulationStepTime;
	}
}
