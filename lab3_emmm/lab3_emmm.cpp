//Вариант 7.
//Решить бесконечношаговую задачу о замене оборудования.

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>

using namespace std;

//------ПО ВАРИАНТУ----------------
const int P = 252;
const double ALPHA = 0.8;

double r_t(double t) {
	return (190 - (146 / (t + 1)));
}

double phi_t(double t) {
	return (127 / (t + 1));
}
//---------------------------------
//Вывод в файл под названием output
void output_file(double P0, int month, vector<double> p0_i) {
	ofstream file;
	file.open("output.txt");
	
	file << "Значения Р(0) в каждом месяце до N = " << p0_i.size() << endl;
	file << "\tN\tP(0)" << endl;
	for (int i = 0; i < p0_i.size(); i++) {
		file << "\t" << (i + 1) << "\t" << p0_i[i] << endl;
	}

	file << endl << "Ответ: P(0) = " << P0 << ", при N = " << month;

	file.close();
}

//Вычисление P0
double P_0(int N) {
	double sum_p0 = 0;

	for (int i = 0; i < N; i++) {
		sum_p0 += (1 / (1 - pow(ALPHA, N))) * pow(ALPHA, i) * r_t(i);
	}

	sum_p0 += (1 / (1 - pow(ALPHA, N))) * pow(ALPHA, N) * (P - phi_t(N));
	return sum_p0;
}

//Вызов основных методов:
int Calculating() {
	vector<double> p0_i;							//таблица всех вычесленных Р(0)
	int N = 1;										//количество месяцев, ч/з которые делается замена оборудования

	p0_i.push_back(P_0(N));
	while (true) {
		N++;
		p0_i.push_back(P_0(N));

		if (p0_i[N - 1] > p0_i[N - 2]) break;		//если значение Р(0) начало возрастать, то завершить цикл
	}

	double min_p0 = p0_i[0];
	int month = 1;

	for (int i = 0; i < p0_i.size(); i++) {			//поиск минимального Р(0)
		if (p0_i[i] < min_p0) {
			min_p0 = p0_i[i];
			month = i + 1;
		}
	}

	output_file(min_p0, month, p0_i);

	return 7;
}
//----------------------------------

int main() {
	setlocale(LC_ALL, "Russian");

	cout << "Происходит вычисление..." << endl;

	int logic = Calculating();

	if (logic == 7)
		cout << "Вычисления произведены!" << endl;

	return 0;
}