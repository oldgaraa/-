// Реализация метода оптимального управления (метод Понтрягина) при выведении РН на орбиту
#include "stdafx.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "math.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "vecSost.h"
const char FILENAME[] = { "Result.txt" };
std::ofstream txtfile(FILENAME);
std::ofstream txtfile2("2.txt");

void Integr(double p_m0, double p_y0, double Tetta0, double *yk, double *tk, double *mk, bool flagPrint = true);
void inversionMatrix(double **A, int N);
double sgn(double A);

int main()
{
	char *s = new char[5];
	std::cout << s;
	//Мой и Оличкин
	double yk1 = 245E3, yk2 = 420E3;
	//double yk1 = 234E3, yk2 = 405E3;
	std::vector <double> X, Y, Yfinal;
	std::vector <double> dX, delta;
	dX.resize(3);
	delta.resize(3);
	dX[0] = 1E-5;
	dX[1] = 1E-6;
	dX[2] = 1E-3;
	//НУ						
	X.resize(3);
	Y.resize(3);
	Yfinal.resize(3);
	X[2] = 500;	//715.88; //1195.152708;						// 2 p_m0
	X[1] = 0.0026;//0.00358997;//0.001910665 ;					// 1 p_y0
	X[0] = 58.81442961 * M_PI / 180;//60.488* M_PI / 180;//56.14268853 * M_PI / 180;			// 0 Tetta0
	Yfinal[2] = 283.925;						// 2 Vyk (tk) //277.967706
	Yfinal[1] = 0;					// 1 mk (ro) 0.13563
	Yfinal[0] = yk1;						// 0 yk

	std::vector <double> dFi;
	dFi.resize(3); //невязки
	double Xprir[3][3] = { 0 }, Yprir[3][3] = { 0 };
	int N = 3;
	double **Jak = new double *[N];
	for (int i = 0; i < N; i++) Jak[i] = new double[N];

	//варьируем p_y0, Tetta0, p_m0;
	// получаем yk, Vyk, mk

	//варьируем p_y0, Tetta0
	// получаем yk, mk
	for (;;) {
		txtfile.close();
		txtfile.open(FILENAME, std::ofstream::out);
		Integr(X[2], X[1], X[0], &Y[0], &Y[2], &Y[1], 0);

		//считаем приращения
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				Xprir[i][j] = (i == j) ? X[j] + dX[j] : X[j];

		for (int i = 0; i < N; i++) {
			Integr(Xprir[i][2], Xprir[i][1], Xprir[i][0], &Yprir[i][0], &Yprir[i][2], &Yprir[i][1], 0);
		}

		//считаем Якобиан
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
 				Jak[j][i] = (Yprir[i][j] - Y[j]) / (dX[i]);
			}
		}

		/*for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				std::cout << std::setw(20) << Jak[i][j];

			std::cout << std::endl ;
		}
		std::cout<< std::endl;*/


		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				if (abs(Jak[i][j]) <= 1E-8)
					std::cout << "JAK<0!!!!!!!!!" << std::endl;
		}
		inversionMatrix(Jak, N);

		/*for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				std::cout << std::setw(20) << Jak[i][j];

			std::cout << std::endl;
		}
		*/

		//невязки
		for (int i = 0; i < N; i++) dFi[i] = (Y[i] - Yfinal[i]);

		if ((abs(dFi[0]) < 1E-3) && (abs(dFi[1]) < 1E-8) && (abs(dFi[2]) < 1E-3)) break;

		for (int i = 0; i < N; i++) {
			delta[i] = 0;
			for (int j = 0; j < N; j++) 
				delta[i] += dFi[j] * Jak[i][j];
		}

		double stepmax = 1;
		//if (abs(dFi[0]) < 10000) stepmax = 0.1;
		//if (abs(dFi[0]) < 1000) stepmax = 0.01;
		//if (abs(dFi[0]) < 100) stepmax = 0.001;
		//if (abs(dFi[0]) < 10) stepmax = 0.00001;
		//if (abs(dFi[0]) < 1) stepmax = 0.000001;
		//if (abs(dFi[0]) < 0.1) stepmax = 0.000001;
		if (abs(delta[0]) > stepmax * M_PI / 180.0) delta[0] = stepmax * M_PI / 180.0*sgn(delta[0]);
		if (abs(delta[1]) > stepmax*1E-4) delta[1] = stepmax*1E-4*sgn(delta[1]);
		if (abs(delta[2]) > 10*stepmax) delta[2] = 10*stepmax*sgn(delta[2]);
		for (int i = 0; i < N; i++) X[i] -= delta[i];
		//if (X[0] < 10 * M_PI / 180) X[0] = 10 * M_PI / 180;
		//if (X[0] > 80 * M_PI / 180) X[0] = 80 * M_PI / 180;
		//if (X[1] < 0) X[1] = 0;
		//if (X[1] > 0.0023) X[1] = 0.0023;
		//if (X[2] < 500) X[2] = 500;
		//if (X[2] > 3000) X[2] = 3000;
		std::cout << std::endl;
		std::cout << std::setw(20) << X[0] * 180 / M_PI << std::setw(20) << X[1] << std::setw(20) << X[2] << std::endl;
		std::cout << std::setw(20) << dFi[0] * 180 / M_PI << std::setw(20) << dFi[1] << std::setw(20) << dFi[2] << std::endl;
		std::cout << std::endl;
		txtfile2 << std::setw(20) << dFi[0]/Yfinal[0] * 180 / M_PI << std::setw(20) << dFi[1]/Yfinal[1] << std::setw(20) << dFi[2] / Yfinal[2] << std::endl;

	}
	std::cout << std::endl;
	std::cout << std::setw(20) << X[0] * 180 / M_PI << std::setw(20) << X[1] << std::setw(20) << X[2] << std::endl;
	std::cout << std::setw(20) << dFi[0] * 180 / M_PI << std::setw(20) << dFi[1] << std::setw(20) << dFi[2] << std::endl;
	std::cout << std::endl;
	txtfile2 << std::setw(20) << dFi[0] / Yfinal[0] * 180 / M_PI << std::setw(20) << dFi[1] / Yfinal[1] << std::setw(20) << dFi[2] / Yfinal[2] << std::endl;

	txtfile.close();
	txtfile.open(FILENAME, std::ofstream::out);
	Integr(X[2], X[1], X[0], &Y[0], &Y[2], &Y[1]);

	txtfile.close();
	for (int i = 0; i < N; i++)
		delete[] Jak[i];
	delete[] Jak;
	system("pause");
	return 0;
}

//Функция интегрирования
void Integr(double p_m0, double p_y0, double Tetta0, double *yk, double *tk, double *mk, bool flagPrint) {
	vecSost Rock(p_m0, p_y0, Tetta0), buf(0, 0, 0);
	//double tk = 430, 
		double dt0 = 1, dt = dt0;
	double ro_min = 1E100;
	bool flag_ro,  //знак функции ro на предыдущем шаге
		flag_ronew,	//знак функции ro на этом шаге
		flag_end = 0;   //второй переход
	int step = 0,
		st = 25;

	for (;;) {
		
		//дробление шага на выход
		if (Rock.Parametrs[5] < 0) {
		//if (Rock.Parametrs[0] > tk) {
			Rock = buf;
			dt /= 2;
		}
		else if (flagPrint == 1) {
			if (step%st == 0) Rock.print(txtfile);
			step++;
		}
		if (abs(Rock.Parametrs[5] - 0) <= 1E-10) {
		//if (abs(Rock.Parametrs[0] - tk) <= 1E-10) {
			if (flagPrint == 1) Rock.print(txtfile);
			break;
		}

		flag_ro = Rock.ro > 0 ? 1 : 0;
		buf = Rock;
		Rock.kutaStep(dt);
		flag_ronew = Rock.ro > 0 ? 1 : 0;

		//дробление шага на ro
		if (flag_ro != flag_ronew) {
			if (abs(Rock.ro) < 1E-13) {
				dt = dt0;
				flag_ronew = (flag_ro == 1) ? 1 : 0;
			}
			else {
				Rock = buf;
				dt /= 2;
			}
		}
		//std::cout << Rock.Parametrs[0] << std::endl;
		if (Rock.ro < ro_min) ro_min = Rock.ro;
	}
	*yk = Rock.Parametrs[3];
	*tk = Rock.Parametrs[0];
	//*mk = Rock.Parametrs[1];
	*mk = ro_min;
}

//Функция инверсии матрицы
void inversionMatrix(double **A, int N)
{
	double temp;

	double **E = new double *[N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}

double sgn(double A) {
	return A >= 0 ? 1 : -1;
}