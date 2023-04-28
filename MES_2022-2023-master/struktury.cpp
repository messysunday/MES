#include "struktury.h"
#include <fstream>

void wypisz_macierz (vector <vector <double>> & matrix)
{
	for (auto& i : matrix)
	{
		for (auto j : i)
		{
			cout <<setw(10) << j << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}

void wypisz_macierz(vector <double>& matrix)
{
	for (auto i : matrix)
	{
		//cout << fixed;
		cout << setw(10) << i << " ";
	}
	cout << "\n";
}

void wyczysc_macierz(vector <vector<double>>& matrix)
{
	for (int i = 0; i < matrix.size(); ++i) 
	{
		for (int j = 0; j < matrix[0].size(); ++j) 
		{
			matrix[i][j] = 0.0;
		}
	}
}

void generuj_siatke(Grid & grid, Global_Data & data) 
{
	//dane wejsciowe do zadania - sciana domu
	//warstwa 1: mur z pustakow ceramicznych grubosc 24 cm
	//warstwa 2: welna mineralna 12 cm
	//warstwa 3: cegla 12 cm
	//grubosc sciany: 48 cm
	//najwezszy element ma 12 cm grubosci, aby miec min 3 elementy na scianke, wybrano szerokosc elementu = 4 cm
	data.simulationStepTime = 3600 * 12; //12 h
	data.simulationTime = data.simulationStepTime*2*30; //30 dni
	double ceglaC = 0.4, ceglaD = 800, ceglaSH = 2510; //deska
	double izolacjaC = 0.041 , izolacjaD = 78.0, izolacjaSH = 750; //welna
	double pustakC = 0.105, pustakD = 385, pustakSH = 840; //HH
	double tot1 = 9.0; //temp na zewnatrz - bedzie brana pod uwage w warunku brzegowym dla sciany ceglanej
	double alfa1 = 27.0; //alfa cegla-powietrze
	double tot2 = 25.0; //temp wewnatrz - bedzie brana pod uwage w warunlu brzegowym dla pustaka
	double alfa2 = 0.45; //alfa pustak-powietrze
	double h = 2.2; 
	double b = 2.2; 
	double ceglax = 0.15; 
	double pustakx = 1.75;
	double izolacjax = b - ceglax - pustakx;
	double dx = 0.05; //3 elementy w sciance z cegly
	double dy = 0.05;
	double nhH = h/dy;
	double nhB = b/dx;
	grid.nN = (nhH + 1)* (nhB + 1);
	grid.nE = nhH * nhB;
	double x = 0;
	double y = 0;
	int id = 0;
	for (int i = 1; i <= nhH + 1; ++i) 
	{
		x = 0;
		for (int j = 1; j <= nhB + 1; ++j) 
		{
			id++;
			grid.Nodes.push_back(Node());
			grid.Nodes[id - 1].x = x;
			grid.Nodes[id - 1].y = y;
			x += dx;
		}
		y += dy;
	}
	int ide = 0;
	for (int i = 1; i <= nhB; ++i) 
	{
		for (int j = 1; j <= nhH; ++j)
		{
			ide++;
			int i1 = (i - 1) * (nhH + 1) + j;
			int i2 = i* (nhH + 1) + j;
			int i3 = i * (nhH + 1) + j + 1;
			int i4 = (i - 1) * (nhH + 1) + j + 1;
			grid.Elements.push_back(Element());
			grid.Elements[ide - 1].Node_ID.push_back(i1);
			grid.Elements[ide - 1].Node_ID.push_back(i4);
			grid.Elements[ide - 1].Node_ID.push_back(i3);
			grid.Elements[ide - 1].Node_ID.push_back(i2);
			
			
		}
	}
	for (int i = 0; i < grid.Nodes.size(); ++i) {
		y = grid.Nodes[i].y;
		x = grid.Nodes[i].x;
		if (x <= 0.00001 || x >= b-0.00001)
		{
			grid.Nodes[i].BC = true;
		}
		grid.Nodes[i].t = tot1;
	}
	for (auto& i : grid.Elements)
	{
		double xmin = b;
		double xmax = 0.0;
		double ymin = h;
		double ymax = 0.0;
		for (auto j : i.Node_ID) {
			if (grid.Nodes[j - 1].x < xmin)
			{
				xmin = grid.Nodes[j].x;
			}
			if (grid.Nodes[j - 1].x > xmax)
			{
				xmax = grid.Nodes[j].x;
			}
		}
		if ((xmin >= 0 && xmax <= pustakx))
		{
			i.conductivity = pustakC;
			i.density = pustakD;
			i.specificHeat = pustakSH;
			i.tot = tot2;
			i.alfa = alfa2;
		}
		else if ((xmin >= pustakx && xmax <= pustakx + izolacjax))
		{
			i.conductivity = izolacjaC;
			i.density = izolacjaD;
			i.specificHeat = izolacjaSH;
		}
		else
		{
			i.conductivity = ceglaC;
			i.density = ceglaD;
			i.specificHeat = ceglaSH;
			i.tot = tot1;
			i.alfa = alfa1;
		}
	}
}

void zapisz_do_pliku(string nazwa, Grid & grid, Equations & equations) 
{
	ofstream zapis;
	zapis.open(nazwa);
	zapis << "# vtk DataFile Version 2.0\n";
	zapis << "Unstructured Grid Example\n";
	zapis << "ASCII\n";
	zapis << "DATASET UNSTRUCTURED_GRID\n\n";
	zapis << "POINTS " << grid.nN << " double\n";
	for (auto i : grid.Nodes) 
	{
		zapis << i.x << " " << i.y << " " << 0 << "\n";
	}
	zapis << "\n";
	zapis << "CELLS " << grid.nE << " " << grid.nE * 5 << "\n";
	for (auto i : grid.Elements) 
	{
		zapis << 4 << " ";
		for (auto j : i.Node_ID)
		{
			zapis << j - 1 << " ";
		}
		zapis << "\n";
	}
	zapis << "\nCELL_TYPES " << grid.nE << "\n";
	for (int i = 0; i < grid.nE; ++i) 
	{
		zapis << 9 << "\n";
	}
	zapis << "\nPOINT_DATA " << grid.nN << "\n";
	zapis << "SCALARS Temp double 1\n";
	zapis << "LOOKUP_TABLE default\n";
	for (double i : equations.temps0) {
		zapis << i << "\n";
	}
}

vector<vector<double>> operator * (vector<vector<double>>& a, vector<vector<double>>& b)
{
	vector <double> pom(b[0].size(), 0.0);
	vector <vector <double>> wynik(a.size(), pom);
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b[0].size(); ++j)
		{
			for (int k = 0; k < b.size(); ++k)
			{
				wynik[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return wynik;
}

vector<vector<double>> operator * (vector<vector<double>>& a, vector<double>& b)
{
	vector <double> pom(b.size(), 0.0);
	vector <vector <double>> wynik(a.size(), pom);
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b.size(); ++j)
		{
			for (int k = 0; k < a[0].size(); ++k)
			{
				wynik[i][j] += a[i][k] * b[j];
			}
		}
	}
	return wynik;
}

vector<vector<double>> operator * (vector <double>& a, vector <vector <double>>& b) {
	vector <double> pom(a.size(), 0.0);
	vector <vector <double>> wynik(1, pom);
	for (int j = 0; j < b.size(); ++j)
	{
		for (int k = 0; k < a.size(); ++k)
		{
			wynik[0][j] += a[k] * b[k][j];
		}
	}
	return wynik;
}

vector<vector<double>> operator * (double& a, vector<vector<double>>& b)
{
	vector <vector <double>> wynik = b;
	for (int i = 0; i < b.size(); ++i)
	{
		for (int j = 0; j < b[i].size(); ++j)
		{
			wynik[i][j] = a * b[i][j];
		}
	}
	return wynik;
}

vector <vector<double>> T(vector <vector<double>>& matrix)
{
	vector <double> pom(matrix.size(), 0.0);
	vector <vector <double>> wynik (matrix[0].size(), pom);
	for (int i = 0; i < wynik.size(); ++i)
	{
		for (int j = 0; j < wynik[i].size(); ++j)
		{
			wynik[i][j] = matrix[j][i];
		}
	}
	return wynik;
}

vector <vector<double>> T(vector<double>& matrix)
{
	vector <vector <double>> wynik;
	for (int i = 0; i < matrix.size(); ++i)
	{
		wynik.push_back({ matrix[i] });
	}
	return wynik;
}

vector <vector <double>> operator + (vector<vector<double>> a, vector<vector<double>> b)
{
	vector <double> pom(b[0].size(), 0.0);
	vector <vector <double>> wynik(a.size(), pom);
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b[i].size(); ++j)
		{
			wynik[i][j] = a[i][j] + b[i][j];
		}
	}
	return wynik;
}

vector <vector <double>> operator - (vector<vector<double>> a, vector<vector<double>> b)
{
	vector <double> pom(b[0].size(), 0.0);
	vector <vector <double>> wynik(a.size(), pom);
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b[i].size(); ++j)
		{
			wynik[i][j] = a[i][j] - b[i][j];
		}
	}
	return wynik;
}

void Global_Data::wypisz()
{
	cout << simulationTime << "\n" << simulationStepTime << "\n" << conductivity << "\n" << alfa << "\n" << tot << "\n" << initialTemperature << "\n"
		<< density << "\n" << specificHeat << "\n";
}

void Node::wypisz()
{
	cout << "x: " << setw(10) << setprecision(10) << x << "\ty: " << setw(10) << setprecision(10) << y << "\ttemp: " << setw(10) << setprecision(10) << t << "\tBC: " << BC << "\n";
}

void Element::wypisz()
{
	for (auto i : Node_ID)
	{
		cout << i << " ";
	}
	cout << "\n";
}

Element::Element()
{
	jacobi.push_back({ 0.0, 0.0 });
	jacobi.push_back({ 0.0, 0.0 });
	dNdx = { {0.0, 0.0, 0.0, 0.0} };
	dNdy = dNdx;
	for (int i = 0; i < 4; ++i)
	{
		H.push_back({ 0.0, 0.0, 0.0, 0.0 });
		C.push_back({ 0.0, 0.0, 0.0, 0.0 });
		HBC.push_back({ 0.0, 0.0, 0.0, 0.0 });
		P.push_back({ 0.0 });
	}
}

void Element::macierzJacobiego(Element4& element4, int nr,  Grid & grid)
{
	for (int i = 0; i < 4; ++i)
	{
		jacobi[0][0] += element4.macierz_ksi[nr][i] * grid.Nodes[Node_ID[i] - 1].x; //dni/deta * xi dla punktu od 0 do 3
		jacobi[0][1] += element4.macierz_ksi[nr][i] * grid.Nodes[Node_ID[i] - 1].y;
		jacobi[1][0] += element4.macierz_eta[nr][i] * grid.Nodes[Node_ID[i] - 1].x;
		jacobi[1][1] += element4.macierz_eta[nr][i] * grid.Nodes[Node_ID[i] - 1].y;
	}
}

void Element::jacobian()
{
	detJ = jacobi[0][0] * jacobi[1][1] - jacobi[1][0] * jacobi[0][1];
}

void Element::macierzOdwrotnaJacobiego()
{
	for (int i = 0; i < jacobi.size(); ++i)
	{
		for (int j = 0; j < jacobi[i].size(); ++j)
		{
			jacobi[i][j] /= detJ;
		}
	}
	double pom = jacobi[1][1];
	jacobi[1][1] = jacobi[0][0];
	jacobi[0][0] = pom;
	jacobi[1][0] *= (-1.0);
	jacobi[0][1] *= (-1.0);
}

void Element::macierzH_macierzC(Element4& element4, Grid& grid, Global_Data & data)
{
	for (int j = 0; j < element4.rzad * element4.rzad; ++j)
	{
		wyczysc_macierz(jacobi);
		macierzJacobiego(element4, j, grid);
		jacobian();
		macierzOdwrotnaJacobiego();
		for (int i = 0; i < 4; ++i)
		{
			dNdx[0][i] = element4.macierz_ksi[j][i] * jacobi[0][0] + element4.macierz_eta[j][i] * jacobi[0][1];
			dNdy[0][i] = element4.macierz_ksi[j][i] * jacobi[1][0] + element4.macierz_eta[j][i] * jacobi[1][1];
		}
		vector <vector <double>> dNdxT = T(dNdx);
		vector <vector <double>> dNdyT = T(dNdy);
		vector <vector <double>> mnozeniex = dNdxT * dNdx;
		vector <vector <double>> mnozeniey = dNdyT * dNdy;
		vector <vector <double>> sumaMacierzy = mnozeniex + mnozeniey;
		double iloczyn = conductivity * detJ;
		double waga = element4.waga[j / element4.rzad] * element4.waga[j % element4.rzad];
		double mnoznik = waga * iloczyn;
		vector <vector <double>> Hj = mnoznik * sumaMacierzy;
		H = H + Hj;


		vector <vector <double>> N = { element4.CN[j] };
		vector <vector <double>> NT = T(N);
		vector <vector <double>> CJ = NT * N;
		mnoznik = waga * detJ;
		CJ = mnoznik * CJ;
		C = C + CJ;
	}
	double stala = (density * specificHeat) / (double)data.simulationStepTime;
	C = stala  * C;
}

void Element::macierzHBC_wektorP(Element4& element4, Grid& grid, Global_Data& data, Equations & equations) 
{	
	int bound1 = 0, bound2 = 0;
	if (grid.Nodes[Node_ID[0] - 1].BC && grid.Nodes[Node_ID[1] - 1].BC) //sciana 1
	{
		bound1 = 0;
		bound2 = 1;
		double detj = sqrt(pow((grid.Nodes[Node_ID[bound1] - 1].x - grid.Nodes[Node_ID[bound2] - 1].x), 2) + pow((grid.Nodes[Node_ID[bound1] - 1].y - grid.Nodes[Node_ID[bound2] - 1].y), 2))/2.0;
		for (int i = 0; i < element4.bounds[bound1].rzad; ++i)
		{
			vector <double> N = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> NT = T(N);
			vector <vector <double>> iloczynMacierzy = NT * N;
			iloczynMacierzy = element4.bounds[bound1].waga[i] * iloczynMacierzy;
			iloczynMacierzy = detj * iloczynMacierzy;
			iloczynMacierzy = alfa * iloczynMacierzy;
			HBC = HBC + iloczynMacierzy;
			vector <double> Ppom = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> PT = T(Ppom);
			//wypisz_macierz(PT);
			PT = alfa * PT;
			PT = detj * PT;
			PT = element4.bounds[bound1].waga[i] * PT;
			PT = tot * PT;
			P = P + PT;
		}
	}
	if (grid.Nodes[Node_ID[1] - 1].BC && grid.Nodes[Node_ID[2] - 1].BC)
	{
		bound1 = 1;
		bound2 = 2;
		double detj = sqrt(pow((grid.Nodes[Node_ID[bound1] - 1].x - grid.Nodes[Node_ID[bound2] - 1].x), 2) + pow((grid.Nodes[Node_ID[bound1] - 1].y - grid.Nodes[Node_ID[bound2] - 1].y), 2)) / 2.0;
		for (int i = 0; i < element4.bounds[bound1].rzad; ++i)
		{
			vector <double> N = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> NT = T(N);
			vector <vector <double>> iloczynMacierzy = NT * N;
			iloczynMacierzy = element4.bounds[bound1].waga[i] * iloczynMacierzy;
			iloczynMacierzy = detj * iloczynMacierzy;
			iloczynMacierzy = alfa * iloczynMacierzy;
			HBC = HBC + iloczynMacierzy;
			vector <double> Ppom = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> PT = T(Ppom);
			//wypisz_macierz(PT);
			PT = alfa * PT;
			PT = detj * PT;
			PT = element4.bounds[bound1].waga[i] * PT;
			PT = tot * PT;
			P = P + PT;
		}
	}
	if (grid.Nodes[Node_ID[2] - 1].BC && grid.Nodes[Node_ID[3] - 1].BC)
	{
		bound1 = 2;
		bound2 = 3;
		double detj = sqrt(pow((grid.Nodes[Node_ID[bound1] - 1].x - grid.Nodes[Node_ID[bound2] - 1].x), 2) + pow((grid.Nodes[Node_ID[bound1] - 1].y - grid.Nodes[Node_ID[bound2] - 1].y), 2)) / 2.0;
		for (int i = 0; i < element4.bounds[bound1].rzad; ++i)
		{
			vector <double> N = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> NT = T(N);
			vector <vector <double>> iloczynMacierzy = NT * N;
			iloczynMacierzy = element4.bounds[bound1].waga[i] * iloczynMacierzy;
			iloczynMacierzy = detj * iloczynMacierzy;
			iloczynMacierzy = alfa * iloczynMacierzy;
			HBC = HBC + iloczynMacierzy;

			vector <double> Ppom = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> PT = T(Ppom);
			//wypisz_macierz(PT);
			PT = alfa * PT;
			PT = detj * PT;
			PT = element4.bounds[bound1].waga[i] * PT;
			PT = tot * PT;
			P = P + PT;
		}
	}
	if (grid.Nodes[Node_ID[3] - 1].BC && grid.Nodes[Node_ID[0] - 1].BC)
	{
		bound1 = 3;
		bound2 = 0;
		double detj = sqrt(pow((grid.Nodes[Node_ID[bound1] - 1].x - grid.Nodes[Node_ID[bound2] - 1].x), 2) + pow((grid.Nodes[Node_ID[bound1] - 1].y - grid.Nodes[Node_ID[bound2] - 1].y), 2)) / 2.0;
		for (int i = 0; i < element4.bounds[bound1].rzad; ++i)
		{
			vector <double> N = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> NT = T(N);
			vector <vector <double>> iloczynMacierzy = NT * N;
			iloczynMacierzy = element4.bounds[bound1].waga[i] * iloczynMacierzy;
			iloczynMacierzy = detj * iloczynMacierzy;
			iloczynMacierzy = alfa * iloczynMacierzy;
			HBC = HBC + iloczynMacierzy;
			vector <double> Ppom = element4.bounds[bound1].macierz_N[i];
			vector <vector <double>> PT = T(Ppom);
			//wypisz_macierz(PT);
			PT = alfa * PT;
			PT = detj * PT;
			PT = element4.bounds[bound1].waga[i] * PT;
			PT = tot * PT;
			P = P + PT;
		}
	}
	H = HBC + H;
}

void Grid::wypisz()
{
	cout << nN << "\n";
	for (auto i : Nodes)
	{
		i.wypisz();
	}
	cout << nE << "\n";
	for (auto i : Elements)
	{
		i.wypisz();
	}
}

Element4::Element4(int newRzad, int newRzad2)
{
	bounds = vector <Bound>(4, Bound());
	rzad = newRzad;		//ilu wezlowa kwadratura
	if (rzad == 2) //definiowanie wspolrzednych i wag wezlow kwadratur w zaleznosci od ilosci wezlow
	{
		eta = { (-1.0) * sqrt(1.0 / 3.0),sqrt(1.0 / 3.0) };
		ksi = { (-1.0) * sqrt(1.0 / 3.0),sqrt(1.0 / 3.0) };
		waga = { 1.0, 1.0 };
	}
	else if (rzad == 3)
	{
		ksi = { (-1.0) * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
		eta = { (-1.0) * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
		waga = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	}
	else
	{
		rzad = 4;
		eta = { (-1.0) * sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))), (-1.0)*sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))),sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))) };
		ksi = { (-1.0) * sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))), (-1.0) * sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))),sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))) };
		waga = { (18.0-sqrt(30.0))/36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
	}
	vector <double> pom(4, 0.0); //wektor pomocniczy do tworzenia macierzy pochodnych po ksi i eat
	for (int i = 0; i < rzad * rzad; ++i) //generowanie pustych macierzy na pochodne ksi i eta (1 vector dla 1 pkt, mamy rzad^2 pkt)
	{
		macierz_ksi.push_back(pom);
		macierz_eta.push_back(pom);
		CN.push_back(pom);
	}
	generuj_macierz_ksi(); //generowanie macierzy pochodnych
	generuj_macierz_eta();
	generuj_macierz_CN();
	for (int i = 0; i < 4; ++i)
	{
		bounds[i] = Bound(newRzad2, i);
	}
}

void Element4::generuj_macierz_ksi() //pochodne Ni po ksi
{
	for (int i = 0; i < eta.size() * eta.size(); ++i) //petla po punktach calkowania (dla 2 wezlowej 4 pkt, dla 3 wezlowej 9 pkt)
	{
		macierz_ksi[i][0] = (-0.25) * (1 - eta[i / eta.size()]);	//dN1/dksi
		macierz_ksi[i][1] = (0.25) * (1 - eta[i / eta.size()]);		//dN2/dksi
		macierz_ksi[i][2] = (0.25) * (1 + eta[i / eta.size()]);		//dN3/dksi
		macierz_ksi[i][3] = (-0.25) * (1 + eta[i / eta.size()]);	//dN4/dksi
	}
	
}

void Element4::generuj_macierz_eta() //pochodne Ni po eta
{
	for (int i = 0; i < ksi.size() * ksi.size(); ++i) //petla po punktach calkowania (dla 2 wezlowej 4 pkt, dla 3 wezlowej 9 pkt)
	{
		macierz_eta[i][0] = (-0.25) * (1 - ksi[i % ksi.size()]);	//dN1/deta
		macierz_eta[i][1] = (-0.25) * (1 + ksi[i % ksi.size()]);	//dN2/deta
		macierz_eta[i][2] = (0.25) * (1 + ksi[i % ksi.size()]);		//dN3/deta
		macierz_eta[i][3] = (0.25) * (1 - ksi[i % ksi.size()]);		//dN4/deta
	}
}

void Element4::generuj_macierz_CN()
{
	for (int i = 0; i < CN.size(); ++i)
	{
		CN[i][0] = 0.25 * (1 - ksi[i % ksi.size()]) * (1 - eta[i / eta.size()]);
		CN[i][1] = 0.25 * (1 + ksi[i % ksi.size()]) * (1 - eta[i / eta.size()]);
		CN[i][2] = 0.25 * (1 + ksi[i % ksi.size()]) * (1 + eta[i / eta.size()]);
		CN[i][3] = 0.25 * (1 - ksi[i % ksi.size()]) * (1 + eta[i / eta.size()]);
	}
}

void Element4::wypisz()
{
	cout << "macierz ksi:\n";
	for (int i = 0; i < eta.size() * eta.size(); ++i)
	{
		cout << eta[i / eta.size()] << ": ";
		for (auto j : macierz_ksi[i])
		{
			cout << j << " ";
		}
		cout << "\n";
	}
	cout << "macierz eta: \n";
	for (int i = 0; i < ksi.size() * ksi.size(); ++i)
	{
		cout << ksi[i % ksi.size()] << ": ";
		for (auto j : macierz_eta[i])
		{
			cout << j << " ";
		}
		cout << "\n";
	}

}

Equations::Equations(Global_Data & global_data, Grid & grid)
{
	vector <double> pom (grid.nN, 0.0);
	temps0 = pom;
	for (int i = 0; i < grid.nN; ++i)
		{
		temps0[i] = grid.Nodes[i].t;
		}
	temps1 = pom;
	for (int i = 0; i < grid.nN; ++i)
	{
		HG.push_back(pom);
		CG.push_back(pom);
		p.push_back({ 0.0 });
	}
}

void Equations::agregacjaP(const Element& element)
{
	for (int i = 0; i < element.Node_ID.size(); ++i)
	{
		int ID = element.Node_ID[i];
		p[ID - 1][0] += element.P[i][0];
	}
}

void Equations::gauss(vector<vector<double>>& wspolczynniki, vector<vector<double>>& wolne, vector<double>& wyniki)
{
	int n = wspolczynniki.size();
	double mnoznik = 0.0;
	for (int i = 0; i < n - 1; ++i) //petla po wierszach
	{
		for (int j = i + 1; j < n; ++j) //petla po elementach wiersza != 0
		{
			mnoznik = wspolczynniki[j][i] / wspolczynniki[i][i];
			for (int k = i; k < n; ++k)
			{
				wspolczynniki[j][k] -= mnoznik * wspolczynniki[i][k];
			}
			wolne[0][j] -= wolne[0][i] * mnoznik;
		}
	}

	//czesc 2: postepowanie odwrotne
	for (int i = n - 1; i >= 0; --i)
	{
		double skladnik = 0;
		for (int k = i; k < n; k++)
		{
			skladnik += wyniki[k] * wspolczynniki[i][k];
		}
		wyniki[i] = (wolne[0][i] - skladnik) / wspolczynniki[i][i];
	}
}

void Equations::agregacjaH(const Element & element)
{
	for (int i = 0; i < element.Node_ID.size(); ++i)
	{
		for (int j = 0; j < element.Node_ID.size(); ++j)
		{
			int ID1 = element.Node_ID[i];
			int ID2 = element.Node_ID[j];
			HG[ID1 - 1][ID2 - 1] += element.H[i][j];
		}
	}
}

void Equations::solve()
{
	vector <double> pom(temps0.size(), 0.0);
	vector <vector <double>> wspolczynniki = HG + CG;
	vector <vector <double>> PT = T(p);
	vector <vector <double>> iloczyn = (temps0 * CG);
	vector <vector <double>> wolne = PT + iloczyn;
	gauss(wspolczynniki, wolne, temps1);
	wypisz_macierz(temps1);
	temps0 = temps1;
	temps1 = pom;
}

void Equations::agregacjaC(const Element& element)
{
	for (int i = 0; i < element.Node_ID.size(); ++i)
	{
		for (int j = 0; j < element.Node_ID.size(); ++j)
		{
			int ID1 = element.Node_ID[i];
			int ID2 = element.Node_ID[j];
			CG[ID1 - 1][ID2 - 1] += element.C[i][j];
		}
	}
}

Bound::Bound(int newRzad, int newNrS)
{
	rzad = newRzad;		//ilu wezlowa kwadratura
	nrS = newNrS; //ktora sciana
	vector <double> pom1;
	vector <double> pom2(rzad, -1);
	vector <double> pom3(rzad, 1);
	vector <double> pom4;
	if (rzad == 2) //definiowanie wspolrzednych i wag wezlow kwadratur w zaleznosci od ilosci wezlow
	{
		pom1 = { (-1.0) * sqrt(1.0 / 3.0),sqrt(1.0 / 3.0) };
		pom4 = { sqrt(1.0 / 3.0), (-1.0) * sqrt(1.0 / 3.0) };
		waga = { 1.0, 1.0 };
	}
	else if (rzad == 3)
	{
		pom1 = { (-1.0) * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
		pom4 = { sqrt(3.0 / 5.0), 0.0, (-1.0) * sqrt(3.0 / 5.0) };
		waga = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	}
	else
	{
		rzad = 4;
		pom1 = { (-1.0) * sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))), (-1.0) * sqrt((3.0 / 7.0) - (2.0 / 7.0) * (sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + (2.0 / 7.0) * (sqrt(6.0 / 5.0))) };
		pom4 = pom1;
		reverse(pom4.begin(), pom4.end());
		waga = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
	}
	switch (nrS) {
	case 0:
	{
		ksi = pom1;
		eta = pom2;
		break;
	}
	case 1:
	{
		ksi = pom3;
		eta = pom1;
		break;
	}
	case 2:
	{
		ksi = pom4;
		eta = pom3;
		break;
	}
	case 3:
	{
		ksi = pom2;
		eta = pom4;
		break;
	}
	}
	vector <double> pom(4, 0.0); //wektor pomocniczy do tworzenia macierzy pochodnych po ksi i eat
	for (int i = 0; i < rzad; ++i) //generowanie pustych macierzy na pochodne ksi i eta (1 vector dla 1 pkt, mamy rzad^2 pkt)
	{
		macierz_N.push_back(pom);
	}
	generuj_macierz_N(); //generowanie macierzy pochodnych
}

void Bound::generuj_macierz_N()
{
	for (int i = 0; i < rzad; ++i)
	{
		macierz_N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
		macierz_N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
		macierz_N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
		macierz_N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
	}
}
