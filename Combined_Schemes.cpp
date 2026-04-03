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


/*Source Files*/
#include "Datastruct.cpp"
#include "Functions.cpp"

/***********************************************************************************************************************************************************
The Master Function
************************************************************************************************************************************************************/

int main()
{
	simuparam sim;
	double t = 0;
	int ctr = 0;

	Read_Data(sim);

	gridparam grid(sim.N);
	flowvar flow(sim.N);

	Static_Print(sim);

	Initialize_and_Grid_formation(sim,grid,flow);

	while(t<sim.t_final)
	{
		if(sim.selector == 0)		
			Upwind_Scheme(sim,flow);
		

		if(sim.selector == 1)		
			LW_Scheme(sim,flow);
		

		if(sim.selector == 2)		
			MacCormack_Scheme(sim,flow);


		for(int i=0; i<sim.N; i++)
			flow.u(i) = flow.u1(i);

		if(ctr%1000 == 0)
		{
			print_file(ctr,sim,flow,grid);
			cout<<"Iteration Number - "<<ctr<<"\t\tSimulation Time - "<<t<<endl;
		}

		t = t+sim.delt;
		ctr++;
	}

	return 0;
}