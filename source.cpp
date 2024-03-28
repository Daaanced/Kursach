#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "source.h"

using namespace std;

int left_con, right_con, N = 0, nodes;
double sum1 = 0, sum2 = 0;

vector<int> ia;
vector<double> di, al; // матрица А
vector<double> M_di, M_al; // матрица M
vector<double> G_di, G_al; //матрица G
vector<double> grid; // пространственная сетка
vector<double> grid_time; // пространственная сетка
vector<double> q; // вектора решения
vector<double> d; // вектор правой части изменённый
vector<double> qM; // q * M * 1 / deltat вектор
double lambda = 1.0;
double gamma = 1.0;
double beta = 2.0;

double f_function(double r, double t) { return  exp(t); } // значение функции правой части в заданной точке

double u_function(double r, double t) { return exp(t); } // аналитическое значение функции

double du_function(double r, double t) { return   0; } // производная аналитического значения функции

double ubeta_function(double r, double t) { return r * t + pow(r, 2) * t; } // функция ubeta для третьих краевых

void left_con1(double t) // левые краевые условия 1-го рода
{
	di[0] = 1;
	d[0] = 1 * u_function(grid[0], t);
	d[1] = d[1] - (al[0] * d[0]);
	d[2] = d[2] - (al[1] * d[0]);
	al[0] = 0;
	al[1] = 0;
}

void right_con1(double t) // правые краевые условия 1-го рода
{
	di.back() = 1;
	d.back() = 1 * u_function(grid.back(), t);
	d[d.size() - 2] = d[d.size() - 2] - (al[al.size() - 1] * d.back());
	d[d.size() - 3] = d[d.size() - 3] - (al[al.size() - 2] * d.back());
	al[al.size() - 1] = 0;
	al[al.size() - 2] = 0;
}

void left_con2(double t) // левые краевые условия 2-го рода
{
	d[0] += -lambda * du_function(grid[0], t);
}

void right_con2(double t) // правые краевые условия 2-го рода
{
	d.back() += lambda * du_function(grid.back(), t);
}

//void right_con3(double t) // правые краевые условия 3-го рода
//{
//	di.back() += beta;
//	d.back() += beta * ubeta_function(grid.back(), t);
//}


void build_matrix() // портрет матрицы и считывание сетки из файла
{
	int i;
	ifstream cons, grid_r, grid_t;
	double buf;
	grid_r.open("nodes.txt"); //открытие файла
	if (!(grid_r.is_open()))
	{
		printf_s("Ошибка открытия файла nodes.txt\n");
		exit;
	}
	if (!grid_r.eof()) grid_r >> buf;
	else
	{
		printf_s("Некорректное количество узлов");
		exit;
	}
	while (!grid_r.eof()) //считывание сетки из файла
	{
		grid.push_back(buf);
		grid_r >> buf;
		grid.push_back((grid.back() + buf) / 2.0);
		N++;
	}
	grid.push_back(buf);
	nodes = grid.size();
	grid_r.close();

	grid_t.open("time.txt"); //открытие файла
	if (!(grid_t.is_open()))
	{
		printf_s("Ошибка открытия файла time.txt\n");
		exit;
	}
	if (!grid_t.eof()) grid_t >> buf;
	else
	{
		printf_s("Некорректное количество узлов ПО ВРЕМЕНИ");
		exit;
	}
	while (!grid_t.eof()) //считывание сетки из файла
	{
		grid_time.push_back(buf);
		grid_t >> buf;
	}
	grid_time.push_back(buf);
	int timesize = grid_time.size();
	grid_t.close();

	for (int i = 0; i < grid_time.size(); i++) cout << grid_time[i] << endl;

	if (timesize < 3)
	{
		printf_s("Некорректное количество узлов ПО ВРЕМЕНИ");
		exit;
	}

	if (nodes < 3)
	{
		printf_s("Некорректное количество узлов");
		exit;
	}
	cons.open("data.txt");
	if (!cons.is_open())
	{
		printf_s("Ошибка открытия файла data.txt\n");
		exit;
	}
	if (cons.eof())
	{
		printf_s("Отсутствуют данные в файле data.txt\n");
		exit;
	}
	cons >> left_con >> right_con;
	if (left_con < 1 || left_con > 3)
	{
		cout << "Некорректное левое краевое условие";
		exit;
	}
	if (right_con < 1 || right_con > 3)
	{
		cout << "Некорректное левое правое условие";
		exit;
	}
	cons.close();
	ia.resize(nodes + 1);
	ia[0] = 1;
	ia[1] = 1;
	for (i = 0; i < nodes - 1; i++) ia[i + 2] = ia[i + 1] + 1 + i % 2;

	di.resize(nodes);
	al.resize(ia.back() - 1);
	M_di.resize(nodes);
	M_al.resize(ia.back() - 1);
	G_di.resize(nodes);
	G_al.resize(ia.back() - 1);
	q.resize(nodes);
	d.resize(nodes);
	grid_r.close();
}

void build_A(vector<double>& di, vector<double>& al, vector<double>& M_di, vector<double>& M_al, vector<int>& ia)
{
	for (int i = 0; i < nodes; i++) di[i] += M_di[i];
	for (int i = 0; i < al.size(); i++) al[i] += M_al[i];
}

void matrix_M(double delta) //добавляет в локальную матрицу матрицу масс
{
	int i;
	double h;
	for (i = 0; i < N; i++)
	{
		h = (grid[2 * i + 2] - grid[2 * i]); //вычисляем шаг
		M_di[2 * i] += (gamma * h / 30.0) * ((grid[2 * i] * 4.0) + (h / 2.0));
		M_di[2 * i + 1] += (gamma * h / 30.0) * ((grid[2 * i] * 16.0) + (h / 2.0) * 16.0);
		M_di[2 * i + 2] += (gamma * h / 30.0) * ((grid[2 * i] * 4.0) + (h / 2.0) * 7.0);
		M_al[ia[2 * i + 2] - 2] += (gamma * h / 30.0) * (grid[2 * i] * 2.0);
		M_al[ia[2 * i + 3] - 3] += (gamma * h / 30.0) * ((grid[2 * i] * -1.0) + (h / 2.0) * -1.0);
		M_al[ia[2 * i + 3] - 2] += (gamma * h / 30.0) * ((grid[2 * i] * 2.0) + (h / 2.0) * 4.0);
	}

	for (int i = 0; i < M_di.size(); i++) M_di[i] = M_di[i] / delta;
	for (int i = 0; i < M_al.size(); i++) M_al[i] = M_al[i] / delta;
}

void matrix_G()	//добавляет в локальную матрицу матрицу жесткости
{
	double h;
	for (int i = 0; i < N; i++)
	{
		h = (grid[2 * i + 2] - grid[2 * i]); //вычисляем шаг
		di[2 * i] += lambda / (3 * h) * ((grid[2 * i] * 7.0) + (h / 2.0) * 3.0);
		di[2 * i + 1] += lambda / (3 * h) * ((grid[2 * i] * 16.0) + (h / 2.0) * 16.0);
		di[2 * i + 2] += lambda / (3 * h) * ((grid[2 * i] * 7.0) + (h / 2.0) * 11.0);
		al[ia[2 * i + 2] - 2] += lambda / (3 * h) * ((grid[2 * i] * -8.0) + (h / 2.0) * -4.0);
		al[ia[2 * i + 3] - 3] += lambda / (3 * h) * (grid[2 * i] + (h / 2.0));
		al[ia[2 * i + 3] - 2] += lambda / (3 * h) * ((grid[2 * i] * -8.0) + (h / 2.0) * -12.0);
	}
}

void build_right(vector<double>& M_di, vector<double>& d, double t) //генерирует вектор правой части
{
	int i;
	vector<double> sum;
	sum.resize(nodes);
	double h, a, b, c, c11, c22, c33, c21, c31, c32;
	for (i = 0; i < N; i++)
	{
		h = (grid[2 * i + 2] - grid[2 * i]); //вычисляем шаг
		a = f_function(grid[2 * i], t);
		b = f_function(grid[2 * i + 1], t);
		c = f_function(grid[2 * i + 2], t);
		c11 = h / 30.0 * ((grid[2 * i] * 4.0) + (h / 2.0));
		c22 = h / 30.0 * ((grid[2 * i] * 16.0) + (h / 2.0) * 16.0);
		c33 = h / 30.0 * ((grid[2 * i] * 4.0) + (h / 2.0) * 7.0);
		c21 = h / 30.0 * (grid[2 * i] * 2.0);
		c31 = h / 30.0 * ((grid[2 * i] * -1.0) + (h / 2.0) * -1.0);
		c32 = h / 30.0 * ((grid[2 * i] * 2.0) + (h / 2.0) * 4.0);
		d[2 * i] += c11 * a + c21 * b + c31 * c;
		d[2 * i + 1] += c21 * a + c22 * b + c32 * c;
		d[2 * i + 2] += c31 * a + c32 * b + c33 * c;
	}

	for (int i = 0; i < d.size(); i++)
	{
		d[i] += qM[i];
	}
}

void conditions(double t)
{
	if (left_con == 1) left_con1(t);
	if (left_con == 2) left_con2(t);
	//if (left_con == 3) left_con3(t);
	if (right_con == 1) right_con1(t);
	if (right_con == 2) right_con2(t);
	//if (right_con == 3) right_con3(t);
}

int LLt()
{
	int i, j, k;
	double per;
	int a;
	int b;
	di[0] = sqrt(di[0]); // l11 = sqrt(a11)
	for (i = 1; i < nodes; i++) // пробегаем по строкам матрицы
	{
		a = i - ia[i + 1] + ia[i]; // начало i-й строки в абсолютной нумерации
		for (j = 0; j < ia[i + 1] - ia[i]; j++) // пробегаем по столбцам до диагонали
		{
			b = a + j - ia[a + j + 1] + ia[a + j]; // начало j-й строки в абсолютной нумерации
			per = al[ia[i] + j - 1];
			if (a < b) for (k = ia[a + j + 1] - ia[a + j] - 1; k >= 0; k--) per -= al[ia[a + j] + k - 1] * al[ia[i] + b - a + k - 1]; // если i-я строка началась раньше j-й строки
			else for (k = a - b; k < ia[a + j + 1] - ia[a + j]; k++) per -= al[ia[i] + k - 1 - (a - b)] * al[ia[a + j] + k - 1]; // если j-я строка началась раньше i-й строки
			al[ia[i] + j - 1] = per / di[a + j];
		}
		per = di[i]; // диагональные элементы
		for (k = 0; k < ia[i + 1] - ia[i]; k++) per -= al[ia[i] + k - 1] * al[ia[i] + k - 1];
		di[i] = sqrt(per);
		if (di[i] == 0) return 1;
	}
	return 0;
}

void gauss()
{
	int i, j;
	q[0] = d[0] / di[0]; // прямой ход метода Гаусса
	for (i = 0; i < nodes; i++)
	{
		q[i] = d[i];
		for (j = 0; j < ia[i + 1] - ia[i]; j++)
			q[i] -= al[ia[i] + j - 1] * q[i - ia[i + 1] + ia[i] + j];
		q[i] = q[i] / di[i];
	}
	for (i = 0; i < nodes; i++)
		d[i] = q[i]; // заносим в вектор правой части решение обратного хода метода Гаусса
	for (i = nodes - 1; i >= 0; i--) // обратный ход метода Гаусса
	{
		q[i] = d[i] / di[i];
		for (j = 0; j < ia[i + 1] - ia[i]; j++) d[i - ia[i + 1] + ia[i] + j] -= q[i] * al[ia[i] + j - 1];
	}
}



void mult_matrix_vector(vector<double> dd, vector<double> all, vector<double> qq, vector<double>& result, int elems) // умножение матрицы M на вектор q(i - 1)
{
	vector<double> res(elems);
	int length = res.size();
	for (int i = 0; i < N; i++)
	{
		res[2 * i] += qq[2 * i] * dd[2 * i] + qq[2 * i + 1] * all[3 * i] + qq[2 * i + 2] * all[3 * i + 1];
		res[2 * i + 1] += qq[2 * i] * all[3 * i] + qq[2 * i + 1] * dd[2 * i + 1] + qq[2 * i + 2] * all[3 * i + 2];
		res[2 * i + 2] += qq[2 * i] * all[3 * i + 1] + qq[2 * i + 1] * all[3 * i + 2];
	}

	res[length - 1] += qq[length - 1] * dd[length - 1];
	result = res;
}



void clearM()
{
	d.assign(d.size(), 0);
	di.assign(di.size(), 0);
	M_di.assign(M_di.size(), 0);
	G_di.assign(G_di.size(), 0);
	al.assign(al.size(), 0);
	M_al.assign(M_al.size(), 0);
	G_al.assign(G_al.size(), 0);
}

void output()
{
	int j;
	double qz = 0;
	FILE* fout;
	if ((fopen_s(&fout, "output.txt", "w")) != 0)
	{
		cout << "Не удаётся открыть файл output.txt";
		return;
	}

	for (int i = 0; i < q.size(); i++) q[i] = u_function(grid[i], grid_time[0]);
	fprintf_s(fout, "r (узлы)\t\tq = u(r) \t\t\tq* \t\t\t\t\t||q-q*|| (погрешность)\n");
	for (int j = 0; j < nodes; j++)
	{
		qz = u_function(grid[j], grid_time[0]);
		sum1 += pow(q[j] - qz, 2);
		sum2 += pow(qz, 2);
		fprintf_s(fout, "%.2lf\t\t\t%2.10lf\t\t%2.10lf\t\t%2.10lf\n", grid[j], q[j], qz, abs(q[j] - qz));
	}
	for (int i = 1; i < grid_time.size(); i++)
	{
		double delta = grid_time[i] - grid_time[i - 1];
		matrix_M(delta);
		matrix_G();
		mult_matrix_vector(M_di, M_al, q, qM, q.size());
		q.assign(d.size(), 0);
		build_right(M_di, d, grid_time[i]);
		build_A(di, al, M_di, M_al, ia);
		conditions(grid_time[i]);
		LLt();
		gauss();
		fprintf_s(fout, "r (узлы)\t\tq = u(r) \t\t\tq* \t\t\t\t\t||q-q*|| (погрешность)\n");
		for (int j = 0; j < nodes; j++)
		{
			qz = u_function(grid[j], grid_time[i]);
			sum1 += pow(q[j] - qz, 2);
			sum2 += pow(qz, 2);
			fprintf_s(fout, "%.2lf\t\t\t%2.10lf\t\t%2.10lf\t\t%2.10lf\n", grid[j], q[j], qz, abs(q[j] - qz));
		}
		fprintf_s(fout, "Относительная погрешность ||q*-q|| / ||q*|| = %.2le\n", sqrt(sum1) / sqrt(sum2));
		clearM();
	}

	fclose(fout);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	build_matrix();
	output();

}

