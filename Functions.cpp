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

/***********************************************************************************************************************************************************
Subroutines
/***********************************************************************************************************************************************************

/*Functions Declaration*/

void Read_Data(simuparam&);
void Static_Print(simuparam&);
void Initialize_and_Grid_formation(simuparam&, gridparam&, flowvar&);
void Upwind_Scheme(simuparam&, flowvar&);
void LW_Scheme(simuparam&, flowvar&);
void MacCormack_Scheme(simuparam&, flowvar&);
void print_file(int,simuparam&, flowvar&, gridparam&);


/*Function to read simulation parameters from file*/

void Read_Data(simuparam &sim)
{
	ifstream in("Read_File.txt");
	in>>sim.N;
	in>>sim.x_ini;
	in>>sim.x_fin;
	in>>sim.c;
	in>>sim.delt;
	in>>sim.t_final;
	in>>sim.selector;

	sim.delx = (sim.x_fin - sim.x_ini)/(sim.N-1);
	sim.a = (sim.c*sim.delt)/sim.delx;
}


/*Printing name of method and simulation parameters on the terminal*/

void Static_Print(simuparam &sim)
{
	if(sim.selector == 0)
		cout<<"\n******Upwind Scheme******\n\n";

	if(sim.selector == 1)
		cout<<"\n******Lax-Wendroff Scheme******\n\n";

	if(sim.selector == 2)
		cout<<"\n******MacCormack Scheme******\n\n";

	cout<<"Number of points - "<<sim.N<<endl
		<<"Time step - "<<sim.delt<<endl
		<<"Step size - "<<sim.delx<<endl
		<<"Constant Velocity - "<<sim.c<<"\n\n";
}


/*Formation of 1D grid and initialization of values*/

void Initialize_and_Grid_formation(simuparam &sim, gridparam &grid, flowvar &flow)
{
	for(int i=0; i<sim.N; i++)
	{
		grid.x(i) = double(i)/double(sim.N-1);
		flow.u(i) = exp(-100*pow((grid.x(i)-0.5),2));
	}

	ofstream file;
	file.open("Initial.dat");
	for(int i=0; i<sim.N; i++)
		file<<grid.x(i)<<"\t"<<flow.u(i)<<endl;
	file.close();

	flow.u1(0) = flow.u(0);
	flow.u1(sim.N-1) = flow.u(sim.N-1);
}


/*Function for the Upwind Scheme*/

void Upwind_Scheme(simuparam &sim, flowvar &flow)
{
	for(int i=1; i<sim.N; i++)		
		flow.u1(i) = flow.u(i) - sim.a*(flow.u(i) - flow.u(i-1));

	flow.u1(0) = flow.u(sim.N-1);
}


/*Function for the Lax-Wendroff Scheme*/

void LW_Scheme(simuparam &sim, flowvar &flow)
{
	for(int i=1; i<sim.N-1; i++)		
		flow.u1(i) = flow.u(i) - 0.5*sim.a*(flow.u(i+1) - flow.u(i-1)) + 0.5*pow(sim.a,2)*(flow.u(i+1) - 2*flow.u(i) + flow.u(i-1));			
	
	flow.u1(0) = flow.u(sim.N-1);
	flow.u1(sim.N-1) = flow.u(sim.N-2);
}


/*Function for the MacCormack Scheme*/

void MacCormack_Scheme(simuparam &sim, flowvar &flow)
{
	double f1 = 0;
	double f2 = 0;

	for(int i = 1; i<sim.N-1; i++)
	{
		f1 = flow.u(i-1,0) - sim.a*(flow.u(i,0) - flow.u(i-1,0));
		f2 = flow.u(i,0) - sim.a*(flow.u(i+1,0) - flow.u(i,0));

		flow.u1(i,0) = 0.5*(flow.u(i,0) + f2 - sim.a*(f2 - f1));
	}

	flow.u1(0,0) = flow.u(sim.N-1,0);
	flow.u1(sim.N-1,0) = flow.u1(sim.N-2,0);
}


/*Function to print u-values at regular intervals during the solving process*/

void print_file(int ctr, simuparam &sim, flowvar &flow, gridparam &grid)
{
	ostringstream fname;
	fname << "data_" << setw(4) << setfill('0') << ctr/1000 << ".dat";

	ofstream file0(fname.str());
    // if (!file0) return 1;

    for(int i=0; i<sim.N; i++)
		file0<<grid.x(i)<<"\t"<<flow.u1(i)<<endl;

	file0.close();
}