#pragma once
//����� ������� ���������
//#include "GOST4401-81.h"
#include <vector>

class vecSost
{
public:
	vecSost(double p_m0, double p_y0, double Tang0);
	~vecSost();
	
	//��� � ������� ������
	double n0 = 0.99, W = 3189.05, gc = 9.80665, bettamax, betta;
	//double n0 = 1.01, W = 3185.20, gc = 9.80665, bettamax, betta;
	//double n0 = 0.95, W = 3174.90, gc = 9.80665, bettamax, betta, mk= 0.13563;
	double p_Vx0, p_Vy0, H, Tang, ro;
	double p_m0, p_y0, Tang0;

	std::vector <double> Parametrs; //������ ���������
	std::vector <double> d_dt;     //������ ����� ��
	std::vector <double> p;

	void kutaStep(double step);    //���� ��� ��������������
	std::vector <double> df_dt(std::vector <double> par);  // ������ ������ ������
	void otherCalculations();		// ������ ��������������� ����������	
	std::vector <double> Sum (const std::vector <double> par, const std::vector <double> K, const double C); //�������� ��������
	void print(std::ofstream& txtfile); //�����
	void NU(double p_m0, double p_y0, double Tang0);									//��������� �������
	vecSost& operator= (const vecSost& a);
	//0 t
	//1 m
	//2 x
	//3 y
	//4 Vx
	//5 Vy
	//6 pm
	//
};

