#include "wczytywanie.h"

void read(Grid& siatka, Global_Data& dane)
{
	string zmiennaNaTekst;
	ifstream odczyt("Test1_4_4.txt");
	//ifstream odczyt("Test2_4_4_MixGrid.txt");
	//ifstream odczyt("Test3_31_31_kwadrat.txt");
	if (!odczyt.is_open())
	{
		cout << "Nie udalo sie otworzyc pliku!\n";
		return;
	}
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.simulationTime;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.simulationStepTime;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.conductivity;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.alfa;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.tot;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.initialTemperature;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.density;
	odczyt >> zmiennaNaTekst;
	odczyt >> dane.specificHeat;
	odczyt >> zmiennaNaTekst;
	odczyt >> zmiennaNaTekst;
	odczyt >> siatka.nN;
	odczyt >> zmiennaNaTekst;
	odczyt >> zmiennaNaTekst;
	odczyt >> siatka.nE;
	odczyt >> zmiennaNaTekst;
	for (int i = 0; i < siatka.nN; ++i)
	{
		int pomoc = 0.0;
		string pomoc2;
		Node newNode;
		odczyt >> pomoc >> pomoc2 >> newNode.x >> pomoc2 >> newNode.y;
		newNode.t = dane.initialTemperature;
		siatka.Nodes.push_back(newNode);
	}
	getline(odczyt, zmiennaNaTekst);
	odczyt >> zmiennaNaTekst;
	odczyt >> zmiennaNaTekst;
	for (int i = 0; i < siatka.nE; ++i)
	{
		Element newElement;
		newElement.conductivity = dane.conductivity;
		newElement.density = dane.density;
		newElement.specificHeat = dane.specificHeat;
		newElement.tot = dane.tot;
		newElement.alfa = dane.alfa;
		vector <int> nodesID(4, 0);
		int pomoc = 0.0;
		string pomoc2;
		Node newNode;
		odczyt >> pomoc >> pomoc2 >> nodesID[0] >> pomoc2 >> nodesID[1] >> pomoc2 >> nodesID[2] >> pomoc2 >> nodesID[3]; //uwaga, tutaj id ida od 1 a nie od 0
		newElement.Node_ID = nodesID;
		siatka.Elements.push_back(newElement);
	}
	odczyt >> zmiennaNaTekst;
	int id;
	odczyt >> id;
	siatka.Nodes[id - 1].BC = true;
	while (!odczyt.fail())
	{
		odczyt >> zmiennaNaTekst;
		odczyt >> id;
		siatka.Nodes[id - 1].BC = true;

	}
}