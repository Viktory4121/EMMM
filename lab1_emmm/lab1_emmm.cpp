//Вариант 7. n=5 m=5 лямбда=60 мю=6 j=500 k=4

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>

#define NUMBER_OF_SERIES 100	//количество серий
#define NUMBER_OF_TESTS 1000	//количество испытаний
#define NUMBER_OF_TABLE 98		//таблица, которая будет выведена в файл
#define ERR -1000.0				//говорит о том, что заявка была не принята на обслуживание

#define N 5			//количество каналов
#define M 5			//количество мест в очереди
#define LAMBDA 60	//интенсивность приёма заявок
#define MU 6		//интенсивность работы канала
#define J 500		//заявка, после которой каналов становится К
#define K 4			//количество каналов после J-й заявки

using namespace std;

vector<double> Generating_numbers() {	//генерация случайных чисел
	vector<double> vv;
	for (int i = 0; i < NUMBER_OF_TESTS; ++i) {
		vv.push_back((rand() / float(32767)) + 0.000001);
	}
	return vv;
}

vector<double> Delta_t(vector<double> rr) {
	vector<double> vv;
	vv.push_back(ERR);
	for (int i = 1; i < NUMBER_OF_TESTS; ++i) vv.push_back(-(log(rr[i]) / LAMBDA));
	return vv;
}

vector<double> T_served(vector<double> qq) {	//время обслуживания заявки
	vector<double> vv;
	for (int i = 0; i < NUMBER_OF_TESTS; ++i) vv.push_back(-(log(qq[i]) / MU));
	return vv;
}

vector<double> T_arrival(vector<double> d_t) {	//время прибытия заявки
	vector<double> vv;
	vv.push_back(0);
	for (int i = 1; i < NUMBER_OF_TESTS; ++i) vv.push_back(d_t[i] + vv[i - 1]);
	return vv;
}

double T_end_of_dervice(double t_ss, double t_s) {	//время конца обслуживания одной заявки(принимает на вход 2 числа, а не массивы)
	return (t_ss + t_s);
}

void Output_(vector<double> CKO, vector<double> MO) {	//вывод в файл output.txt
	ofstream file;
	file.open("output.txt");

	file << "Число обслуженных заявок:" << endl;
	file << "СКО = " << CKO[6] << endl;
	file << "М(х) = " << MO[6] << endl << endl;

	file << "Относительная пропускная способность" << endl;
	file << "СКО = " << CKO[0] << endl;
	file << "М(х) = " << MO[0] << endl << endl;

	file << "Абсолютная пропускная способность:" << endl;
	file << "СКО = " << CKO[1] << endl;
	file << "М(х) = " << MO[1] << endl << endl;

	file << "Среднее время обслуживания:" << endl;
	file << "СКО = " << CKO[2] << endl;
	file << "М(х) = " << MO[2] << endl << endl;

	file << "Среднее время нахождения заявки в очереди: " << endl;
	file << "СКО = " << CKO[3] << endl;
	file << "М(х) = " << MO[3] << endl << endl;

	file << "Среднее число занятых каналов:" << endl;
	file << "СКО = " << CKO[4] << endl;
	file << "М(х) = " << MO[4] << endl << endl;

	file << "Среднее число занятых мест в очереди:" << endl;
	file << "СКО = " << CKO[5] << endl;
	file << "М(х) = " << MO[5] << endl << endl;

	file.close();
}

void Output_table(vector<double> r, vector<double> q, vector<double> d_t, vector<double> t_a, vector<double> t_s,
	vector<double> n_e, vector<double> m_e, vector<double> ch_num, vector<double> t_ss, vector<double> t_es) {	//вывод в файл table.txt
	//'-'
	ofstream file;
	file.open("table.txt");

	file << " i\tr\tq\tdt\ttпр\tn\tm\tNкан.\ttн.о\ttобс\ttк.о" << endl;
	for (int i = 0; i < NUMBER_OF_TESTS; ++i) {

		file << fixed << setprecision(3) << " " << (i + 1) << "\t" << (r[i]) << "\t" << (q[i]) << "\t";

		if (d_t[i] == ERR) file << "--" << "\t";
		else file << fixed << setprecision(3) << (d_t[i]) << "\t";

		file << fixed << setprecision(3) << (t_a[i]) << "\t" << int(n_e[i]) << "\t" << int(m_e[i]) << "\t";

		if (ch_num[i] == ERR) file << "--" << "\t";
		else file << int(ch_num[i]) << "\t";

		if (t_ss[i] == ERR) file << "--" << "\t";
		else file << fixed << setprecision(3) << (t_ss[i]) << "\t";

		file << fixed << setprecision(3) << (t_s[i]) << "\t";

		if (t_es[i] == ERR) file << "--" << endl;
		else file << fixed << setprecision(3) << (t_es[i]) << endl;
	}

	file.close();

}

int Calculating_tables() {	//функция для заполнения таблицы
	vector<vector<double>> series_res(NUMBER_OF_SERIES, vector<double>(7, 0));
	vector<double> CKO(7), MO(7);
	int sum_served_requests = 1;
	int sum_served_t_a;
	double sum_t_s, sum_t_m;

	for (int lll = 0; lll < NUMBER_OF_SERIES; ++lll) {

		vector<double> r, q, d_t, t_a, t_s, n_e, m_e, ch_num, t_ss, t_es;
		r = Generating_numbers();	//случайное число
		q = Generating_numbers();	//случайное число
		d_t = Delta_t(r);			//дельта t
		t_s = T_served(q);			//время обслуживания
		t_a = T_arrival(d_t);		//время прибытия

		//заполнение первой строки:
		n_e.push_back(0);
		m_e.push_back(0);
		t_ss.push_back(0);
		t_es.push_back(T_end_of_dervice(t_ss[0], t_s[0]));
		ch_num.push_back(1);
		//------------------------

		int ch_off; //номер канала, который перестанет работать после J-ого канала

		//заполнение n_e, m_e, ch_num, t_ss, t_es:
		for (int i = 1; i < NUMBER_OF_TESTS; i++) {

			if (i < J) { //условие, что до J-ой заявки количество каналов = N

				int log_int = 1;
				vector<int> log_vector(N); //содержит информацию о каждом канале, что он уже работал (или работает)

				for (int k = 1; k <= N; k++) { //проверка на то, что до текущей заявки работали все каналы
					for (int j = 0; j < i; j++) {
						if (ch_num[j] == k) {
							log_vector[k - 1] = 1;
							break;
						}
						else {
							log_vector[k - 1] = 0;
							log_int = 0;
						}
					}
				}

				double max_n_t_es, min_n_t_es = 1000;
				int sum_z = 0, number_ch, index_n_min, index_n_max, logic = 0;

				for (int k = 1; k <= N; k++) { //прогонка по каналам
					max_n_t_es = -1.0;

					for (int j = 0; j < i; j++) { //прогонка по заявкам, которые были до текущей

						if (k == ch_num[j] && max_n_t_es <= t_es[j]) { //поиск максимального времени к.о. для текущего канала
							max_n_t_es = t_es[j];
							index_n_max = j;
							logic = 1;
						}

						if (t_a[i] < t_es[j] && k == ch_num[j]) { //расчёт занятых каналов и мест в очереди в период прихода текущей заявки
							sum_z++;
						}

						if (log_int == 0 && !log_vector[k - 1]) { //если есть хотя бы 1 канал, который до этого не использовался
							number_ch = k;
							log_int = 7;
						}
					}
					if (min_n_t_es > max_n_t_es && logic) {
						min_n_t_es = max_n_t_es;
						index_n_min = index_n_max;
					}
					logic = 0;
				}

				if (sum_z >= N) {
					n_e.push_back(N);
					m_e.push_back(sum_z - n_e[i]);
				}
				else {
					n_e.push_back(sum_z);
					m_e.push_back(0);
				}

				if (m_e[i] != M) { //текущая заявка обслужится

					if (t_a[i] >= t_es[index_n_min]) {
						//если заявка пришла позже, чем освободится канал
						ch_num.push_back(ch_num[index_n_min]);
						t_ss.push_back(t_a[i]);
						t_es.push_back(T_end_of_dervice(t_ss[i], t_s[i]));
						sum_served_requests++;
						continue;

					}
					else if (t_a[i] < t_es[index_n_min] && log_int == 7) {
						//если заявка прибыла раньше, чем ближайщий по времени к.о. канал завершит обслуживание
						//если какой-то канал до этой заявки не использовался
						ch_num.push_back(number_ch);
						t_ss.push_back(t_a[i]);
						t_es.push_back(T_end_of_dervice(t_ss[i], t_s[i]));
						sum_served_requests++;
						continue;

					}
					else if (t_a[i] < t_es[index_n_min]) {
						//если заявка пришла раньше, чем освободится канал
						ch_num.push_back(ch_num[index_n_min]);
						t_ss.push_back(t_es[index_n_min]);
						t_es.push_back(T_end_of_dervice(t_ss[i], t_s[i]));
						sum_served_requests++;
						continue;
					}
				}
				else if (m_e[i] == M) {
					//заявка не обслуживается
					ch_num.push_back(ERR);
					t_ss.push_back(ERR);
					t_es.push_back(ERR);
				}

				//-------------------------------------------------------------------------------------------------------------------------------------------------------------

			}
			else { //условие, что с J-ой заявки количество каналов = К (уходит заявка, которая раньше закончит обслуживание)

				if (i == J) {
					//определение, какая заявка уходит

					double max_n_t_es, min_n_t_es = 1000;
					int sum_z = 0, index_n_min, index_n_max, logic = 0;

					for (int k = 1; k <= N; k++) { //прогонка по каналам
						max_n_t_es = -1.0;

						for (int j = 0; j < i; j++) { //прогонка по заявкам, которые были до текущей

							if (k == ch_num[j] && max_n_t_es <= t_es[j]) { //поиск максимального времени к.о. для текущего канала
								max_n_t_es = t_es[j];
								index_n_max = j;
								logic = 1;
							}

							if (t_a[i] < t_es[j] && k == ch_num[j]) { //расчёт занятых каналов и мест в очереди в период прихода текущей заявки
								sum_z++;
							}
						}
						if (min_n_t_es > max_n_t_es && logic) {
							min_n_t_es = max_n_t_es;
							index_n_min = index_n_max;
						}
						logic = 0;
					}
					ch_off = ch_num[index_n_min];
				}

				//продолжение расчётов

				double max_n_t_es, min_n_t_es = 1000;
				int sum_z = 0, index_n_min, index_n_max, logic = 0;

				for (int k = 1; k <= N; k++) { //прогонка по каналам
					max_n_t_es = -1.0;

					if (k == ch_off) continue; //пропуск исключенного канала

					for (int j = 0; j < i; j++) { //прогонка по заявкам, которые были до текущей

						if (k == ch_num[j] && max_n_t_es <= t_es[j]) { //поиск максимального времени к.о. для текущего канала
							max_n_t_es = t_es[j];
							index_n_max = j;
							logic = 1;
						}

						if (t_a[i] < t_es[j] && k == ch_num[j]) { //расчёт занятых каналов и мест в очереди в период прихода текущей заявки
							sum_z++;
						}
					}
					if (min_n_t_es > max_n_t_es && logic) {
						min_n_t_es = max_n_t_es;
						index_n_min = index_n_max;
					}
					logic = 0;
				}

				if (sum_z >= K) {
					n_e.push_back(K);
					m_e.push_back(sum_z - n_e[i]);
				}
				else {
					n_e.push_back(sum_z);
					m_e.push_back(0);
				}


				if (m_e[i] != M) {//текущая заявка обслужится

					if (t_a[i] >= t_es[index_n_min]) {
						//если заявка пришла позже, чем освободится канал
						ch_num.push_back(ch_num[index_n_min]);
						t_ss.push_back(t_a[i]);
						t_es.push_back(T_end_of_dervice(t_ss[i], t_s[i]));
						sum_served_requests++;
						continue;

					}
					else if (t_a[i] < t_es[index_n_min]) {
						//если заявка пришла раньше, чем освободится канал
						ch_num.push_back(ch_num[index_n_min]);
						t_ss.push_back(t_es[index_n_min]);
						t_es.push_back(T_end_of_dervice(t_ss[i], t_s[i]));
						sum_served_requests++;
						continue;
					}
				}
				else if (m_e[i] == M) {
					//заявка не обслуживается
					ch_num.push_back(ERR);
					t_ss.push_back(ERR);
					t_es.push_back(ERR);
				}
			}

		}

		//расчёт характеристик:

		sum_served_t_a = 0;
		sum_t_s = 0;
		sum_t_m = 0;
		for (int h = 0; h < NUMBER_OF_TESTS; h++) {
			if (t_es[h] < t_a[NUMBER_OF_TESTS - 1]) {
				if (ch_num[h] != ERR) {
					sum_t_s += t_s[h];
					sum_t_m += t_ss[h] - t_a[h];
					sum_served_t_a++;
				}
			}
		}

		series_res[lll][6] = sum_served_t_a;
		series_res[lll][0] = double(sum_served_t_a) / NUMBER_OF_TESTS;
		series_res[lll][1] = sum_served_t_a / t_a[NUMBER_OF_TESTS - 1];
		series_res[lll][2] = sum_t_s / sum_served_t_a;
		series_res[lll][3] = sum_t_m / sum_served_t_a;
		series_res[lll][4] = sum_t_s / t_a[NUMBER_OF_TESTS - 1];
		series_res[lll][5] = sum_t_m / t_a[NUMBER_OF_TESTS - 1];

		//--------------------

		if (lll == NUMBER_OF_TABLE - 1) { //вывод текущей таблицы
			Output_table(r, q, d_t, t_a, t_s, n_e, m_e, ch_num, t_ss, t_es);
		}
	}

	for (int h = 0; h < 7; h++) {
		double sum = 0;
		for (int g = 0; g < NUMBER_OF_SERIES; g++) {
			sum += series_res[g][h];
		}
		MO[h] = sum / NUMBER_OF_SERIES;
		sum = 0;
		for (int g = 0; g < NUMBER_OF_SERIES; g++) {
			sum += (series_res[g][h] - MO[h]) * (series_res[g][h] - MO[h]);
		}
		CKO[h] = sqrt(sum / NUMBER_OF_SERIES);
	}

	Output_(CKO, MO); //вывод характеристик

	return 5;
}

int main() {
	setlocale(LC_ALL, "Russian");

	cout << "Происходит вычисление..." << endl;

	int logic = Calculating_tables();

	if (logic == 5)
		cout << "Вычисления произведены!" << endl;

	return 0;
}