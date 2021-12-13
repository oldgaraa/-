#include "stdafx.h"
#include "vecSost.h"
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
#include "math.h"

vecSost::vecSost(double p_m0, double p_y0, double Tang0)
{
	NU(p_m0, p_y0, Tang0);
}


vecSost::~vecSost()
{
}

void vecSost::kutaStep(double step)
{
	std::vector <double> par(Parametrs.size());
	std::vector <double> K1(d_dt.size());
	par = Parametrs;
	vecSost::otherCalculations();
	K1 = vecSost::df_dt(par);
	for (int i = 0; i < K1.size(); i++) K1[i] *= step;

	std::vector<double> K2(d_dt.size());
	par = vecSost::Sum(Parametrs, K1, 0.5);
	vecSost::otherCalculations();
	K2 = vecSost::df_dt(par);
	for (int i = 0; i < K2.size(); i++) K2[i] *= step;

	std::vector<double> K3(d_dt.size());
	par = vecSost::Sum(Parametrs, K2, 0.5);
	vecSost::otherCalculations();
	K3 = vecSost::df_dt(par);
	for (int i = 0; i < K3.size(); i++) K3[i] *= step;

	std::vector<double> K4(d_dt.size());
	par = vecSost::Sum(Parametrs, K3, 1);
	vecSost::otherCalculations();
	K4 = vecSost::df_dt(par);
	for (int i = 0; i < K4.size(); i++) K4[i] *= step;

	Parametrs = vecSost::Sum(Parametrs, K1, 1.0 / 6.0);
	Parametrs = vecSost::Sum(Parametrs, K2, 2.0 / 6.0);
	Parametrs = vecSost::Sum(Parametrs, K3, 2.0 / 6.0);
	Parametrs = vecSost::Sum(Parametrs, K4, 1.0 / 6.0);
	vecSost::otherCalculations();
}

//Правые части
std::vector <double> vecSost::df_dt(std::vector <double> par)
{

	std::vector <double> d_dt(7);
	//t
	d_dt[0] = 1;

	//m
	d_dt[1] = -betta;

	//x
	d_dt[2] = par[4];

	//y
	d_dt[3] = par[5];

	//Vx
	d_dt[4] = betta*W / par[1] * p_Vx0 / sqrt(p_Vx0*p_Vx0 + pow(-p_y0*par[0] + p_Vy0, 2));

	//Vy
	d_dt[5] = betta*W / par[1] * (-p_y0*par[0] + p_Vy0) / sqrt(p_Vx0*p_Vx0 + pow(-p_y0*par[0] + p_Vy0, 2)) - gc;

	//pm
	d_dt[6] = betta*W / par[1]/ par[1] * sqrt(p_Vx0*p_Vx0 + pow(-p_y0*par[0] + p_Vy0, 2));

	return d_dt;
}

//Прочие расчеты
void vecSost::otherCalculations()
{
	p[4] = -p[2] * Parametrs[0] + p_Vx0;
	p[5]= -p[3] * Parametrs[0] + p_Vy0;
	Tang = atan2(p[5], p[4]);
	ro = (p[4] * cos(Tang) + p[5] * sin(Tang) - Parametrs[6] * Parametrs[1] / W);
	betta = ro > 0 ? bettamax : 0;
	//betta = Parametrs[1] <= 0 ? 0 : betta;
	H = (p[4] * cos(Tang) + p[5] * sin(Tang))*betta*W / Parametrs[1] - p[5] * gc + p[2] * Parametrs[4] + p[3] * Parametrs[5] - Parametrs[6] * betta;
	double ppmm = ((p[4] * cos(Tang) + p[5] * sin(Tang))*betta*W / Parametrs[1] - p[5] * gc + p[2] * Parametrs[4] + p[3] * Parametrs[5]) / betta;
}

//Оператор сложения
std::vector <double> vecSost::Sum(const std::vector <double> par, const std::vector <double> K, const double C) {
	
	std::vector <double> res(par.size());
	for (int i = 0; i < par.size(); i++) res[i] = par[i] + K[i] * C;
	return res;
}

// вывод в файл
void vecSost::print(std::ofstream& txtfile)
{
	for (int i = 0; i < Parametrs.size(); i++) txtfile << std::setw(25) << Parametrs[i];
	txtfile << std::setprecision(10) << std::setw(25) << Tang*180/M_PI << std::setw(25) << ro << std::setw(25) << p[3] << std::setw(25) << p[4] << std::setw(25) << p[5] << std::setw(25) << H << std::setw(25) << betta ;
	txtfile << std::endl;
}

//Начальные условия
void vecSost::NU(double p_m0, double p_y0, double Tang0)
{
	Parametrs.resize(7);
	d_dt.resize(7);
	p.resize(7);

	//t
	Parametrs[0] = 0;
	//m
	Parametrs[1] = 1;
	//x
	Parametrs[2] = 0;
	//y
	Parametrs[3] = 44500;
	//Vx
	Parametrs[4] = 1650;
	//Vy
	Parametrs[5] = 900;
	//pm
	Parametrs[6] = p_m0;

	p_Vx0 = cos(Tang0);
	p_Vy0 = sin(Tang0);

	//px
	p[2] = 0;
	//py
	p[3] = p_y0;
	//pVx
	p[4] = p_Vx0;
	//pVy
	p[5] = p_Vy0;


	bettamax = n0 * Parametrs[1] * gc / W;

	//this->p_m0 = ((p_Vx0 * cos(Tang0) + p_Vy0 * sin(Tang0))*bettamax*W / Parametrs[1] - p[5] * gc + p[2] * Parametrs[4] + p[3] * Parametrs[5]) / bettamax;
	this->p_m0 = p_m0;
	Parametrs[6] = this->p_m0;
	this->p_y0 = p_y0;
	this->Tang0 = Tang0;
	otherCalculations();
}

vecSost& vecSost::operator= (const vecSost& A) {
	//проверка на самоприсваивание
	if (this == &A) {
		return *this;
	}
	Parametrs = A.Parametrs;
	d_dt = A.d_dt;
	bettamax = A.bettamax;
	betta = A.betta;
	p_Vx0 = A.p_Vx0;
	p_Vy0 = A.p_Vy0;
	H = A.H;
	Tang = A.Tang;
	ro = A.ro;
	p_m0 = A.p_m0;
	Tang0 = A.Tang0;
	p_y0 = A.p_y0;
	p = A.p;
	return *this;
}
