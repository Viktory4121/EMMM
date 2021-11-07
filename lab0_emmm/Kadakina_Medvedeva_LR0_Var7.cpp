//Kadakina_Medvedeva_LR0_Var7.cpp
//[6,11]  f(x)=abs(cos(6*x+4)/x)

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random>
#include <ctime>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>

#define NUMBER_OF_TABLE 57      //таблица, которая будет выведена в файл
#define NUMBER_OF_TESTS 100     //количество испытаний
#define NUMBER_OF_SAMPLES 1000  //количество выборок
#define A 6.0     //начало промежутка
#define B 11.0    //конец промежутка
#define Max_y 1/A //функция на данном промежутке будет находится ниже этого значения (точка не касается функцию, но приближена))

using namespace std;

//функция по варианту
vector<double> Function(vector<double> xx, vector<double> ff) { // функция по варианту
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) {
        ff.push_back(abs(cos(6 * xx[i] + 4) / xx[i]));
    }
    return ff;
}

//генерация случайных чисел
vector<double> GenerateRandomNumbers(vector<double> vv) {
    default_random_engine generator(time(0));
    uniform_real<double> dist(0.0, 1.0);
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) vv.push_back(dist(generator));
    return vv;
}

//вычисление х
vector<double> Generate_x(vector<double> rr, vector<double> xx) {
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) {
        xx.push_back(A + (B - A) * rr[i]);
    }
    return xx;
}

//вычисление у
vector<double> Generate_y(vector<double> qq, vector<double> yy) {
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) {
        yy.push_back(Max_y * qq[i]);
    }
    return yy;
}

//проверка на принадлежность точки к области
vector<int> Checking_Numbers(vector<double> yy, vector<double> ff, vector<int> mm) {
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) {
        if (yy[i] <= ff[i]) {
            mm.push_back(1);
        } else {
            mm.push_back(0);
        }
    }
    return mm;
}

//вывод таблицы в файл
void Output_table(vector<double> rr, vector<double> qq, vector<double> xx, vector<double> yy, vector<double> ff, vector<int> mm) {

    //оформление вывода (1 таблица из 1000 строчек и 7 столбцов, мат ожидание и СКО в конце)
    //i, r, q, x, y, f(x), m  - столбцы

    ofstream file;
    file.open("output.txt");

    file << " i\tr\tq\tx\ty\tf(x)\tm" << endl;
    for (int i = 0; i < NUMBER_OF_SAMPLES; ++i) {
        file << fixed << setprecision(3) << " " << (i + 1) << "\t" << (rr[i]) << "\t" << (qq[i]) << "\t" << (xx[i]) << "\t" << (yy[i]) << "\t" << (ff[i]) << "\t" << (mm[i]) << endl;
    }

    file.close();
}

//вывод значений СКО и М(х) в файл
void Output_CKO_MO(double CKO, double MO) {
    ofstream file;
    file.open("output.txt", ios_base::app);

    file << endl << "CKO = " << CKO;
    file << endl << "M(x) = " << MO;

    file.close();
}

//пересчёт таблицы
int Recalculating_the_Table(vector<double> rr, vector<double> qq, vector<double> xx, vector<double> yy, vector<double> ff, vector<int> mm) {
    double CKO = 0, MO = 0;
    vector<double> integ; //приближённое значение интеграла
    double sum = 0.0; //количество точек, принадлежащих области

    for (int i = 0; i < NUMBER_OF_TESTS; ++i) {
        rr = GenerateRandomNumbers(rr);
        qq = GenerateRandomNumbers(qq);
        xx = Generate_x(rr, xx);
        yy = Generate_y(qq, yy);
        ff = Function(xx, ff);
        mm = Checking_Numbers(yy, ff, mm);

        if (i == NUMBER_OF_TABLE - 1) {
            Output_table(rr, qq, xx, yy, ff, mm);
        }

        sum = 0;
        for (int j = 0; j < NUMBER_OF_SAMPLES; ++j) {
            if (mm[j] == 1) ++sum;
        }

        integ.push_back((sum / NUMBER_OF_SAMPLES) * (B - A) * Max_y);
        MO += integ[i] / NUMBER_OF_TESTS;

        rr.clear();
        qq.clear();
        xx.clear();
        yy.clear();
        ff.clear();
        mm.clear();
    }

    for (int i = 0; i < NUMBER_OF_TESTS; ++i) {
        CKO += (integ[i] - MO) * (integ[i] - MO);
    }

    Output_CKO_MO(CKO, MO);

    return 7;
}


int main() {
    setlocale(LC_ALL, "Russian");
    vector<double> r, q, x, y, func;
    vector<int> m;
    int logic = Recalculating_the_Table(r, q, x, y, func, m);

    if (logic == 7)
        cout << "Вычисления произведены!" << endl;

    return 0;
}
