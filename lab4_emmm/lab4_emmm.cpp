//Вариант 7. Задача складирования.

#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>
using namespace std;

//-----ПО ВАРИАНТУ----------------------
const int MONTH = 8;										//кол-во месяцев
const vector<int> ALPHA = {16,13,13,8,9,17,7,6};			//затраты на пополнение продукции
const vector<int> BETA = {6,5,10,8,10,8,15,17};				//доход от расхода продукции со склада
vector<int> X(MONTH), Y(MONTH), P(MONTH);					//глобальное объявление переменных х, у и к
const int C = 36;											//ёмкость склада в ед. продукции
const int K1 = 26;											//кол-во единиц продукции на складе в начале 1го месяца
//--------------------------------------

double P_function(int n, int k) {
	vector<double> p_dop;	//вектор, хранящий Р для каждой из 5-и точек

	if (n == MONTH - 1) {
		p_dop.push_back(0);	//О(0,0)
		p_dop.push_back(BETA[n] * k); //А(0,ki)
		p_dop.push_back(BETA[n] * C - ALPHA[n] * (C - k)); //B1(c-ki,c)
		p_dop.push_back(BETA[n] * k - ALPHA[n] * C); //B2(c,ki)
		p_dop.push_back(-ALPHA[n] * (C - k)); //C(c-ki,0)

		int p_max = p_dop[0];
		int x, y;
		if (p_max <= p_dop[0]) {
			p_max = p_dop[0];
			x = 0;
			y = 0;
		}
		if (p_max <= p_dop[1]) {
			p_max = p_dop[1];
			x = 0;
			y = k;
		}
		if (p_max <= p_dop[2]) {
			p_max = p_dop[2];
			x = C - k;
			y = C;
		}
		if (p_max <= p_dop[3]) {
			p_max = p_dop[3];
			x = C;
			y = k;
		}
		if (p_max <= p_dop[4]) {
			p_max = p_dop[4];
			x = C - k;
			y = 0;
		}
		X[n] = x;
		Y[n] = y;
		P[n] = p_max;
		p_dop.clear();
		return p_max;
		
	} else {
		p_dop.push_back(P_function(n + 1, k));	//О(0,0)
		p_dop.push_back(BETA[n] * k + P_function(n + 1, 0)); //А(0,ki)
		p_dop.push_back(ALPHA[n] * k + (BETA[n] - ALPHA[n]) * C + P_function(n + 1, 0)); //B1(c-ki,c)
		p_dop.push_back(BETA[n] * k - ALPHA[n] * C + P_function(n + 1, C)); //B2(c,ki)
		p_dop.push_back(ALPHA[n] * k - ALPHA[n] * C + P_function(n + 1, C)); //C(c-ki,0)

		int p_max = p_dop[0];
		int x, y;
		if (p_max <= p_dop[0]) {
			p_max = p_dop[0];
			x = 0;
			y = 0;
		}
		if (p_max <= p_dop[1]) {
			p_max = p_dop[1];
			x = 0;
			y = k;
		}
		if (p_max <= p_dop[2]) {
			p_max = p_dop[2];
			x = C - k;
			y = C;
		}
		if (p_max <= p_dop[3]) {
			p_max = p_dop[3];
			x = C;
			y = k;
		}
		if (p_max <= p_dop[4]) {
			p_max = p_dop[4];
			x = C - k;
			y = 0;
		}
		X[n] = x;
		Y[n] = y;
		P[n] = p_max;
		p_dop.clear();
		return p_max;
	}
}

//Вызов основных функций:
int Calculating() {
	
	//О(0,0)
	//А(0,ki)
	//B1(c-ki,c)
	//B2(c,ki)
	//C(c-ki,0)


	//вывод в файл
	ofstream file;
	file.open("output.txt");

	//вид таблицы:
	//i ki xi yi betai*yi-alphai*xi
	//будет 8 строк,тк 8 месяцев
	
	file << "i\tk[i]\tx[i]\ty[i]\tbeta[i]*y[i]-alpha[i]*x[i]" << endl;
	int k = K1;
	for (int i = 0; i < MONTH; i++) {
		P[i] = P_function(i, k);
		file << (i + 1) << '\t' << k << '\t' << X[i] << '\t' << Y[i] << '\t' << (BETA[i] * Y[i] - ALPHA[i] * X[i]) << '\n';
		k += X[i] - Y[i];
	}
	file << "Максимальная прибыль P1(" << K1 << ") = " << P[0];
	file.close();

	return 5;
}

int main() {
	setlocale(LC_ALL, "Russian");

	cout << "Происходит вычисление..." << endl;

	int logic = Calculating();

	if (logic == 5)
		cout << "Вычисления произведены!" << endl;

	return 0;
}