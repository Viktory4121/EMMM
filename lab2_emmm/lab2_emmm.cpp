//Вариант 7. Динамическая задача распределения ресурсов.

#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>

using namespace std;

const int MONTH = 4;   //Кол-во месяцев

//---------------------------------РЕАЛИЗАЦИЯ ОСНОВНЫХ АЛГОРИТМОВ----------------------------------
//чтение 4 таблиц из файла
void read_file(vector<double> *x, vector<double> *f, vector<double> *g, 
               vector<double> *phi, vector<double> *psi) {

    string line;
    ifstream file("Вариант 7.txt"); // окрываем файл для чтения

    if (file.is_open()){

        if(getline(file, line) && line == "   x      f(x)"){ //считывание первой таблицы
                
            int i = 0;
            double xx, ff;

            //putchar('\n'); //перевод курсора на следующую строку

            while (getline(file, line) && line != "") {
                istringstream iss(line);
                iss >> xx >> ff;

                (*x).push_back(xx);
                (*f).push_back(ff);
                //putchar('\n');
                i++;
            }
        }

       // putchar('\n'); //переходим к следующей таблице

        if (getline(file, line) && line == "   x      g(x)") { //считывание первой таблицы

            int i = 0;
            double xx, gg;

            //putchar('\n'); //перевод курсора на следующую строку

            while (getline(file, line) && line != "") {
                istringstream iss(line);
                iss >> xx >> gg;

                (*g).push_back(gg);
                //putchar('\n');
                i++;
            }
        }

        //putchar('\n'); //переходим к следующей таблице

        if (getline(file, line) && line == "   x    phi(x)") { //считывание первой таблицы

            int i = 0;
            double xx, phii;

            //putchar('\n'); //перевод курсора на следующую строку

            while (getline(file, line) && line != "") {
                istringstream iss(line);
                iss >> xx >> phii;

                (*phi).push_back(phii);

                //putchar('\n');
                i++;
            }
        }

        //putchar('\n'); //переходим к следующей таблице

        if (getline(file, line) && line == "   x    psi(x)") { //считывание первой таблицы

            int i = 0;
            double xx, psii;

            //putchar('\n'); //перевод курсора на следующую строку

            while (getline(file, line) && line != "") {
                istringstream iss(line);
                iss >> xx >> psii;

                (*psi).push_back(psii);

                //putchar('\n');
                i++;
            }
        }
    }

    file.close();     // закрываем файл
}

//вычисление аппроксимации
double approximation(vector<double> x, vector<double> fun, double x0) {

    double max = 0.0, min = 0.0;
    int index1 = -1, index2 = -1;

    for (int i = 0; i < x.size(); i++) {
        if (x[i] == x0) {                   //если такое значение х0 есть в таблице общих х
            return fun[i];
        } else if (x[i] < x0) {             //поиск левой границы для аппроксимации
            min = x[i];
            index1 = i;
        } else {                            //поиск правой границы
            max = x[i];
            index2 = i;
            break;
        }
    }

    return (((fun[index1] * (max - x0)) + (fun[index2] * (x0 - min))) / (max - min));
}

//вычисление границ диапазонов для всех с К2 до Кn - НЕ НАДО
/*void Search_K(vector<int> *min, vector<int> *max, vector<double> x, vector<double> phi, vector<double> psi) {
    //ceil() округляет в сторону большего, а floor() - в сторону меньшего

    double buf = x[x.size() - 1];               //перебор со значения К1

    for (int k = MONTH - 2; k >= 0; k--) {

        double mmax = 0, mmin = 0;
        for (int i = x.size() - 1; i >= 0; i--) {

            if (x[i] == buf) {
                buf = double(floor(phi[i]));
                (*min)[k] = floor(phi[i]);
                break;

            } else {

                if (x[i] < buf) { 
                    mmin = x[i];
                    break;
                } else {
                    mmax = x[i];
                }
            }
        }

        if ((*min)[k] == 0) {
            (*min)[k] = floor(approximation(x, phi, buf));          //Вычисление аппроксимации
            buf = (*min)[k];
        }
    }

    //тот же самый цикл для пси
    double duf = x[x.size() - 1];                                   //перебор со значения К1

    for (int k = MONTH - 2; k >= 0; k--) {

        double mmax = 0, mmin = 0;
        for (int i = x.size() - 1; i >= 0; i--) {

            if (x[i] == duf) {
                duf = double(floor(psi[i]));
                (*max)[k] = floor(psi[i]);
                break;

            } else {

                if (x[i] < duf) {
                    mmin = x[i];
                    break;
                } else {
                    mmax = x[i];
                }
            }
        }

        if ((*max)[k] == 0) {
            (*max)[k] = floor(approximation(x, psi, duf));          //Вычисление аппроксимации
            duf = (*max)[k];
        }
    }
}*/

//вывод таблицы обратного планирования - ПЕРЕДЕЛАТЬ
void Output_inverse(vector<double> k, vector<double> x, vector<double> p, vector<int> sum, vector<int> index) {

    ofstream file;
    file.open("inverse.txt");

    int m = 0;
    for (int i = MONTH - 1; i >= 0; i--) {
        file << "Месяц " << i + 1 << ":" << endl;
        file << "\tk\tx\tP" << endl;
        for (int t = index[m]; t < (index[m] + sum[m]); t++) {
            file << fixed << setprecision(3) << "\t" << int(k[t]) << "\t" << int(x[t]) << "\t" << p[t] << endl;
        }
        file << endl << endl;
        m++;
    }
    file.close();
}

//вывод таблицы прямого планирования
void Output_direct(vector<vector<double>> v) {
    ofstream file;
    file.open("direct.txt");
    double sum = 0;

    file << "Прямое планирование:" << endl << endl;

    file << "i  |k\t\t|x\t\t|y\t\t|f(x)\t\t|g(y)\t\t|phi(x)\t\t|psi(y)" << endl;

    for (int i = 0; i < MONTH; i++) {
        file << (i + 1) << "  ";
        for (int j = 0; j < 7; j++) {
            if(v[i][j] != -1.0) file << fixed << setprecision(3) << "|" << v[i][j] << "  \t";
            else file << "|--\t\t";
        }
        file << endl;
        sum += v[i][3] + v[i][4];
    }

    file << endl << "Ответ: P1(k1) = P1(" << int(v[0][0]) << ") = " << sum;

    file.close();
}

//прямое планирование:
void Direct_planning(vector<double> k, vector<double> x, vector<double> p, vector<int> sum, vector<int> ind_sum,            //новые таблицы
    vector<double> x0, vector<double> f, vector<double> g, vector<double> phi, vector<double> psi) {       //таблицы по варианту

    //индексы двумерного вектора:
    //0  1  2   3    4     5      6
    //k  x  y  f(x) g(y) phi(x) psi(y)
    vector<vector<double>> table(MONTH, vector<double>(7));

    for (int i = 0; i < MONTH; i++) {
        if (i == 0) {                                               //первый месяц---------------------------

            table[i][0] = k[k.size() - 1];                          //k1
            table[i][1] = x[x.size() - 1];                          //x1
            table[i][2] = table[i][0] - table[i][1];                //y1

            //далее - аппроксимация:
            table[i][3] = approximation(x0, f, table[i][1]);        //f(x1)
            table[i][4] = approximation(x0, g, table[i][2]);        //g(y1)
            table[i][5] = approximation(x0, phi, table[i][1]);      //phi(x1)
            table[i][6] = approximation(x0, psi, table[i][2]);      //psi(y1)

        } else if (i != MONTH - 1) {                                //i-й месяц, кроме первого и крайнего----

            table[i][0] = table[i - 1][5] + table[i - 1][6];        //ki = phi_i-1 + psi_i-1

            //==создание доп.векторов для аппроксимации==
            vector <double> dop_x, dop_k;
            for (int v = ind_sum[MONTH - 1 - i]; v <= ind_sum[MONTH - i] - 1; v++) {
                dop_x.push_back(x[v]);
                dop_k.push_back(k[v]);
            }
            //вычисление аппроксимации для xi
            table[i][1] = approximation(dop_k, dop_x, table[i][0]); //xi
            dop_x.clear();
            dop_k.clear();
            //===========================================
            
            table[i][2] = table[i][0] - table[i][1];                //yi
            table[i][3] = approximation(x0, f, table[i][1]);        //f(xi)
            table[i][4] = approximation(x0, g, table[i][2]);        //g(yi)
            table[i][5] = approximation(x0, phi, table[i][1]);      //phi(xi)
            table[i][6] = approximation(x0, psi, table[i][2]);      //psi(yi)

        } else {                                                    //крайний месяц-----------------------------

            table[i][0] = table[i - 1][5] + table[i - 1][6];        //ki = phi_i-1 + psi_i-1

            //==создание доп.векторов для аппроксимации==
            vector <double> dop_x, dop_k;
            for (int v = ind_sum[MONTH - 1 - i]; v <= ind_sum[MONTH - i] - 1; v++) {
                dop_x.push_back(x[v]);
                dop_k.push_back(k[v]);
            }
            //вычисление аппроксимации для xi
            table[i][1] = approximation(dop_k, dop_x, table[i][0]); //xi
            dop_x.clear();
            dop_k.clear();
            //===========================================

            table[i][2] = table[i][0] - table[i][1];                //yi
            table[i][3] = approximation(x0, f, table[i][1]);        //f(xi)
            table[i][4] = approximation(x0, g, table[i][2]);        //g(yi)
            table[i][5] = -1.0;                                     //phi(xi)
            table[i][6] = -1.0;                                     //psi(yi)
        }
    }

    Output_direct(table);
}

//вычисление обратного планирования
void Inverse_planning(vector<double> x, vector<double> f, vector<double> phi, vector<double> g,
                      vector<double> psi, int K1) {

    vector<double> kk_, xx_, pp_;   //здесь хранятся все таблицы по всем месяцам (они просто будут добавляться одни за другими)
    vector<int> sum, index_sum;     //количество строк у каждой таблицы каждого месяца и индекс начала этой таблице в общем векторе

    int com_sum = 0;
    int step = x[x.size() - 1] - x[x.size() - 2];   //вычисление шага у х в таблицах по варианту

    vector<double> dop_x, dop_p, dop_k;             //дополнительные вектора для исключения путаницы

    index_sum.push_back(0);                         //первая таблица начинается с 0 индекса
    for (int i = 0; i < MONTH; i++) {
        sum.push_back(0);

        if (i == 0) {               //для последнего месяца

            for (int k = x[0]; k <= x[x.size() - 1]; k += step) {            //перебор по ki значениям соотв. диапазонов

                kk_.push_back(k);

                double maximum = -10.0, pn, g_fun;
                double x_;
                for (int t = 0; t < x.size(); t++) {                                          //поиск максимального P

                    if (x[t] > k) continue;                                //проверка на то, чтобы х > к

                    g_fun = approximation(x, g, (k - x[t]));

                    pn = f[t] + g_fun;
                    if (maximum <= pn) {
                        maximum = pn;
                        x_ = x[t];
                    }
                }
                xx_.push_back(x_);
                pp_.push_back(maximum);

                sum[i]++;
                com_sum++;
            }

            index_sum.push_back(com_sum);

        } else if (i != MONTH - 1) { //для всех месяцев, кроме первого и последнего
                
            for (int k = x[0]; k <= x[x.size() - 1]; k += step) {                     //перебор по ki значениям соотв. диапазонов

                kk_.push_back(k);

                double maximum = -10.0, pn, g_fun, psi_fun, p_fun;
                double x_;

                for (int tt = 0; tt < x.size(); tt++) {                            //поиск максимального P

                    if (x[tt] > k) continue;                                       //проверка на то, чтобы х > к

                    g_fun = approximation(x, g, (k - x[tt]));
                    psi_fun = approximation(x, psi, (k - x[tt]));

                    //--поиск P--
                    for (int v = index_sum[i - 1]; v <= index_sum[i] - 1; v++) {    //заполнение доп. векторов таблицами за будущий месяц(т.е. прошлая таблица)
                        dop_x.push_back(xx_[v]);
                        dop_p.push_back(pp_[v]);
                        dop_k.push_back(kk_[v]);
                    }
                    p_fun = approximation(dop_k, dop_p, (phi[tt] + psi_fun));

                    dop_x.clear();
                    dop_p.clear();
                    dop_k.clear();
                    //-----------

                    pn = f[tt] + g_fun + p_fun;
                    if (maximum <= pn) {
                        maximum = pn;
                        x_ = x[tt];
                    }
                }
                //после того, как длины векторов xx_, pp_ = 173, продожлать по шагам протыкивать
                xx_.push_back(x_);
                pp_.push_back(maximum);

                sum[i]++;
                com_sum++;
            }
            index_sum.push_back(com_sum);

        } else { //если дошли до первого месяца, то
            
            kk_.push_back(K1);
            
            double maximum = -10.0, pn, g_fun, psi_fun, p_fun;
            double x_;

            for (int tr = 0; tr < x.size(); tr++) { //поиск максимального P
                
                g_fun = approximation(x, g, (K1 - x[tr]));
                psi_fun = approximation(x, psi, (K1 - x[tr]));

                //--поиск Р-----
                for (int v = index_sum[i - 1]; v <= index_sum[i] - 1; v++) {
                    dop_x.push_back(xx_[v]);
                    dop_p.push_back(pp_[v]);
                    dop_k.push_back(kk_[v]);
                }
                p_fun = approximation(dop_k, dop_p, (phi[tr] + psi_fun));

                dop_x.clear();
                dop_p.clear();
                dop_k.clear();
                //--------------

                pn = f[tr] + g_fun + p_fun;
                if (maximum <= pn) {
                    maximum = pn;
                    x_ = x[tr];
                }
            }
            xx_.push_back(x_);
            pp_.push_back(maximum);

            sum[i]++;
        }
    }

    Output_inverse(kk_, xx_, pp_, sum, index_sum);      //Вывод таблиц обратного планирования

    Direct_planning(kk_, xx_, pp_, sum, index_sum, x, f, g, phi, psi);     //Вызов функции для составления прямого планирования
}





//функция для прямой оптимизация - 1 таблица и вывод в файл: функция для расчетов, функция для вывода

//-------------------------------------------------------------------------------------------------
//Вызов функций
int Calculating() {
    vector<double> x;                                  //Вектор х для всех функций одинаковый
    vector<double> f, phi;                             //первое производство
    vector<double> g, psi;                             //второе производство
    read_file(&x, &f, &g, &phi, &psi);

    int K1 = x[x.size() - 1];                          //Общее количество ресурсов
    
    //Обратное планирование:
    Inverse_planning(x, f, phi, g, psi, K1);
    //Вызов функции для прямого планирования производится из функции для обратного планирования

    //Очищение векторов:
    x.clear();
    f.clear();
    g.clear();
    phi.clear();
    psi.clear();
	return 5;
}
//-------------------------------------------------------------------------------------------------

int main() {
	setlocale(LC_ALL, "Russian");

	cout << "Происходит вычисление..." << endl;

	int logic = Calculating();

	if (logic == 5)
		cout << "Вычисления произведены!" << endl;

	return 0;
}