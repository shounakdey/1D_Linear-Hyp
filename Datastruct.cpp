/*
This program solves the following simple linear hyperbolic equation - 

du/dt + (c)du/dx = 0

with c = 1 and initial condition - exp(-100(x-0.5)^2), using the following methods

1. Upwind Scheme
2. Lax-Wendroff
3. MacCormack

The different methods can be selected and the parameters can be changed by changing
the values in "Read_File.txt". The values are printed in DAT format in fixed intervals 
which can be modified in the "print_file" subroutine in "Functions.cpp".

Developed by - Shounak Dey, PhD Candidate, JNCASR, 3/4/2026
*/

/**********************************************************************************************************************************************************
											 DEFINITION OF DATA TYPES AND DECLARATION OF DATA STRUCTURES
***********************************************************************************************************************************************************/

/*Header Files*/
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<climits>
#include<fstream>
#include "Matrix.h"
typedef double dataType;

using namespace std;


/*Grid variables*/
struct gridparam
{
	AMatrix<dataType> x;				//1D array containing the grid-points

	gridparam(unsigned N): x(N){}
};


/*Solver parameters*/
struct simuparam
{
	unsigned int N;						//Number of points in the domain
	double x_ini;						//Initial point of domain
	double x_fin;						//Final point of domain
	double c;							//Velocity
	double delt;
	double delx;
	double t_final;
	double a;							//a = c*delt/delx
	int selector;
};


/*Flow variables*/
struct flowvar
{
	AMatrix<dataType> u;				//Array to hold values in n-th time step
	AMatrix<dataType> u1;				//Array to hold values in (n+1)-th time step

	flowvar(unsigned N):u(N),u1(N){}
};