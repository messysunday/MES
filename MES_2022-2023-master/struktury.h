#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <iomanip>

using namespace std;

#ifndef STRUKTURY_H
#define STRUKTURY_H

struct Global_Data;
struct Node;
struct Grid;
struct Element;
struct Element4;
struct Equations;

void wypisz_macierz(vector <vector <double>>& matrix);
void wypisz_macierz(vector <double>& matrix);
void wyczysc_macierz(vector <vector<double>>& matrix);
void generuj_siatke(Grid& grid, Global_Data& data);
void zapisz_do_pliku(string nazwa, Grid& grid, Equations & equations);


struct Global_Data
{
	int simulationTime = 0; //czas symulacji
	int simulationStepTime = 0; //czas jednego kroku symulacji
	int currentTime = 0;
	double conductivity = 0.0; //przewodniosc cieplna
	double alfa = 0.0; //wspolczynnik konwekcyjnej wymiany ciepla
	double tot = 0.0; //temperatura otoczenia
	double initialTemperature = 0.0; //startowa temperatura obiektu
	double density = 0.0; //gestosc obiektu
	double specificHeat = 0.0; //cieplo wlasciwe
	void wypisz();
};

struct Node
{
	double x = 0.0;
	double y = 0.0; //wspolrzedne
	double t = 0.0; //temperatura
	bool BC = false; //czy na brzegu
	void wypisz();
};

struct Grid
{
	int nN = 0; //ilosc wierzcholkow
	int nE = 0; //ilosc elementow
	vector <Node> Nodes; //tablica wezlow
	vector <Element> Elements; //tablica elementow
	void wypisz();
};

struct Bound
{
	vector <vector<double>> macierz_N;
	vector <double> ksi;
	vector <double> eta;
	vector <double> waga;
	int rzad = 0;
	int nrS = 0;
	Bound() {};
	Bound(int newRzad, int newNrS);
	void generuj_macierz_N();
};

struct Element4 //tego tylko 1 instancja, tu sa wartosci tabelaryczne
{
	vector <vector <double>> macierz_ksi; //macierz dNi/ksi
	vector <vector <double>> macierz_eta; //macierz dNi/eta
	vector <vector <double>> CN;
	vector <double> ksi; //vectory wspolrzednych wezlow calkowania
	vector <double> eta;
	vector <double> waga; //wagi wezlow calkowania
	vector <Bound> bounds;
	int rzad = 0; //ilu wezlowa kwadratura
	Element4(int newRzad, int newRzad2); //konstruktor
	void generuj_macierz_ksi(); //funkcje generujace macierze pochodnych
	void generuj_macierz_eta();
	void generuj_macierz_CN();
	void wypisz(); //funkcja wypisujaca wszystkie macierze i vectory
};

struct Element //tego bedzie duzo
{
	vector <int> Node_ID; //tablica wierzcholkow elementow
	vector <vector <double>> jacobi; //macierz jacobiego, odwracana w trakcie obliczen
	vector <vector <double>> H; //macierz H, pozniej H + HBC
	vector <vector <double>> C; //macierz C
	vector <vector <double>> HBC; //macierz HBC
	vector <vector <double>> P; //Wektor P
	vector <vector <double>> dNdx; //pochodne funkcji kszta³tu
	vector <vector <double>> dNdy;
	double detJ = 0.0; //wyznacznik macierzy Jacobiego
	double conductivity = 0.0; //w³asnoœci materia³u, wykorzystywane w trakcie obliczeñ
	double density = 0.0;
	double specificHeat = 0.0;
	double tot = 0.0;
	double alfa = 0.0;
	Element();
	void wypisz();
	void macierzJacobiego(Element4 & element4,int nr, Grid & grid); //funkcja obliczaj¹ca macierz jacobiego
	void jacobian(); //funkcja obliczaj¹ca wyznacznik macierzy jacobiego
	void macierzOdwrotnaJacobiego(); //funkcja odwracaj¹ca macierz jacobiego
	void macierzH_macierzC(Element4& element4, Grid& grid, Global_Data& data); //odrazu macierz h razy waga, macierzC
	void macierzHBC_wektorP(Element4& element4, Grid& grid, Global_Data& data, Equations & equations);
};

struct Equations 
{
	vector<vector <double>> HG;
	vector<vector<double>> CG;
	vector <double> temps0;
	vector <double> temps1;
	vector <vector <double>> p;
	Equations(Global_Data & global_data, Grid & grid);
	void agregacjaH(const Element & element);
	void agregacjaC(const Element& element);
	void agregacjaP(const Element& element);
	void solve();
	void gauss(vector<vector<double>>& wspolczynniki, vector<vector<double>>& wolne, vector<double>& wyniki);
};

#endif // !STRUKTURY_H

