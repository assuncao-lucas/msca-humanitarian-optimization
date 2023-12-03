#include <iostream>
#include <algorithm>
#include <queue>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <dirent.h>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <boost/dynamic_bitset.hpp>
#include "graph.h"
#include "instance.h"
#include "CommodityFormulationsForm.h"
#include "timer.h"
#include "matrix.hpp"
#include "general.h"
#include "graph_algorithms.h"
#include "solution.hpp"
#include "heuristic_solution.h"
#include "feasibility_pump.h"
#include "ALNS.h"

std::pair<char,std::pair<int,int>> DetermineVariable(Instance& instance, std::string index, bool bianchessi)
{
	std::pair<char,std::pair<int,int>> variable;
	Graph * graph = instance.graph();
	int num_vertices = graph->num_vertices();
	int num_arcs = graph->num_arcs();
	int index_num = std::stoi(index,NULL);

	if(index_num <= num_vertices)
	{
		variable.first = 'y';
		variable.second.first = index_num -1;
		variable.second.second = -1;
	}else if(index_num <= num_vertices + num_arcs)
	{
		variable.first = 'x';
		variable.second = (graph->arcs_map_)[index_num - num_vertices -1];
	}else
	{
		if(bianchessi) variable.first = 'z';
		else variable.first = 'f';
		variable.second = (graph->arcs_map_)[index_num - num_vertices - num_arcs -1];
	}

	return variable;
}

// y_i = i (os primeiros num_vertices elementos)
// x_ij = num_vertices + pos(i,j)
// f_ij = num_vertices + num_arcs + pos(i,j)
void ParsePORTAFile(Instance& instance, std::string file_name)
{
	std::fstream input, output;
	input.open("./PORTA/src/example.poi.ieq", std::fstream::in);

	if(!input.is_open()){
		std::cout << "Could not open file!" << std::endl;
		throw 1;
	}

	output.open(file_name.c_str(), std::fstream::out);

	if(!output.is_open()){
		std::cout << "Could not open file!" << std::endl;
		throw 1;
	}

	std::string line;
	while(std::getline(input, line))
	{
		std::size_t prev = 0, pos;
		while ((pos = line.find_first_of("x", prev)) != std::string::npos)
		{
			output<< line.substr(prev,pos-prev);
			// acha índice aqui
			prev = pos + 1;
			pos = line.find_first_of("><=+- ", prev);
			if(pos != std::string::npos)
			{
				std::string index = line.substr(prev,pos-prev+1);
				std::pair<char,std::pair<int,int>> variable = DetermineVariable(instance,index,false);

				switch (variable.first)
				{
					case 'y': output << "y(" << variable.second.first << ")"; break;
					case 'x': output << "x(" << variable.second.first << ")(" << variable.second.second << ")"; break;
					case 'f': output << "f(" << variable.second.first << ")(" << variable.second.second << ")"; break;
					case 'z': output << "z(" << variable.second.first << ")(" << variable.second.second << ")"; break;
					default: break;
				}
				char separator = line[pos];
				switch (separator)
				{
					case '>': output << ">="; prev = pos + 2; break;
					case '<': output << "<="; prev = pos + 2; break;
					case '=': output << "=="; prev = pos + 2; break;
					case '+': output << "+"; prev = pos + 1; break;
					case '-': output << "-"; prev = pos + 1; break;
					case ' ': output << " "; prev = pos + 1; break;
					default: break;
				}
			}else
			{
				throw 2;
				std::cout << "Invalid PORTA file." << std::endl;
			}
			prev = pos+1;
		}
		if (prev < line.length())
		output << line.substr(prev, std::string::npos);
		output << std::endl;
	}

	input.close();
	output.close();
}

void WritePORTAIneqFile(Instance& instance, double * R0, double* Rn, std::string porta_file)
{
	int num_vehicles = instance.num_vehicles();
	int num_vertices = (instance.graph())->num_vertices();
	int num_mandatory = (instance.graph())->num_mandatory();
	GArc * curr_arc = NULL;
	Graph* graph = instance.graph();
	int num_arcs = graph->num_arcs();

	std::fstream file;
	file.open(porta_file.c_str(),std::fstream::out);

	file << "DIM = " << num_vertices + 2*num_arcs << std::endl;

	file << "LOWER_BOUNDS" << std::endl;

	file << "0";
	for(int i = 1; i < num_vertices; i++)
	{
		file << " 0";
	}

	for(int i = 0; i < 2*num_arcs; i++)
	{
		file << " 0";
	}

	file << std::endl;

	file << "UPPER_BOUNDS" << std::endl;

	file << "1";
	for(int i = 1; i < num_vertices; i++)
	{
		file << " 1";
	}

	for(int i = 0; i < num_arcs; i++)
	{
		file << " 1";
	}

	for(int i = 0; i < num_arcs; i++)
	{
		file << " " << instance.limit();
	}

	file << std::endl;

	file << "INEQUALITIES_SECTION" << std::endl;

	int num_adj_arc = 0;
	file << "x1 == 1" << std::endl;
	file << "x" << num_vertices << " == 1" << std::endl;

	for(int i = 1; i < num_vertices - 1; i++)
	{
		num_adj_arc = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[i][j] != NULL)
			{
				file << "+x"<< 1 + num_vertices + graph->pos(i,j);
				num_adj_arc++;
			}
		}

		if (i <= num_mandatory)
		{
			if(num_adj_arc == 0) file << "x" << 1 +i << " == 0" << std::endl;
			else file << " == 1" << std::endl;
			file << "x" << 1 +i << " == 1" << std::endl;
		}
		else
		{
			if(num_adj_arc == 0) file << "x" << 1 +i << " == 0" << std::endl;
			else  file << "-x" << 1 +i << " == 0" << std::endl;
		}
	}


	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[0][j] != NULL) file << "+x" << 1 +num_vertices + graph->pos(0,j);
	}

	file << " == " << num_vehicles << std::endl;

	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[j][num_vertices-1] != NULL) file << "+x" << 1 +num_vertices + graph->pos(j,num_vertices-1);
	}

	file << " == " << num_vehicles << std::endl;

	int cont = 0;
	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[j][0] != NULL)
		{
			file << "+x" << 1 +num_vertices + graph->pos(j,0);
			cont++;
		}
	}

	if(cont != 0) file << " == 0" << std::endl;
	cont = 0;

	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[num_vertices-1][j] != NULL)
		{
			file << "+x" << 1 +num_vertices + graph->pos(num_vertices-1,j);
			cont++;
		}
	}

	if(cont != 0) file << " == 0" << std::endl;


	for(int i = 1; i < num_vertices-1; i++)
	{
		cont = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[i][j] != NULL)
			{
				file << "+x" << 1 +num_vertices + graph->pos(i,j);
				cont++;
			}

			if((*(instance.graph()))[j][i] != NULL)
			{
				file << "-x" << 1 +num_vertices + graph->pos(j,i);
				cont++;
			}
		}
		if(cont != 0) file << " == 0" << std::endl;
	}


	for(int i = 1; i < num_vertices; i++)
	{
		curr_arc = (*(instance.graph()))[0][i];

		if(curr_arc != NULL)
		{
			int coef = (instance.limit() - curr_arc->dist())*K_PRECISION;
			file << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(0,i) << "-x" << 1 +num_vertices + num_arcs + graph->pos(0,i) << " == 0" << std::endl;

		}
	}

	for(int i = 1; i < num_vertices-1; i++)
	{
		cont = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[j][i] != NULL)
			{
				file << "+x" << 1 +num_vertices + num_arcs + graph->pos(j,i);
				cont++;
			}

			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				file << "-x" << 1 +num_vertices + num_arcs + graph->pos(i,j);
				int coef = curr_arc->dist()*K_PRECISION;
				file << "-" << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(i,j);
				cont++;
			}
		}

		if(cont != 0) file << " == 0" << std::endl;
	}


	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				int coef = 0;
				coef = - (instance.limit() - R0[i] - curr_arc->dist())*K_PRECISION;
				file << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(i,j) << "+x" << 1 +num_vertices +  num_arcs + graph->pos(i,j) << " <= 0" << std::endl;
			}
		}
	}

	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				int coef = 0;
				coef = - (Rn[j])*K_PRECISION;
				file << coef << "/" << K_PRECISION << "x" << 1 + num_vertices + graph->pos(i,j) << "+x" << 1 +num_vertices +  num_arcs + graph->pos(i,j) << " >= 0" << std::endl;
			}
		}
	}


	file << "END" << std::endl;
	file.close();
}

void WritePORTAIneqFileBianchessi(Instance& instance, double * R0, double* Rn, std::string file_name)
{
	int num_vehicles = instance.num_vehicles();
	int num_vertices = (instance.graph())->num_vertices();
	int num_mandatory = (instance.graph())->num_mandatory();
	GArc * curr_arc = NULL;
	Graph* graph = instance.graph();
	int num_arcs = graph->num_arcs();

	std::fstream file;
	file.open(file_name.c_str(),std::fstream::out);

	file << "DIM = " << num_vertices + 2*num_arcs << std::endl;

	file << "LOWER_BOUNDS" << std::endl;

	file << "0";
	for(int i = 1; i < num_vertices; i++)
	{
		file << " 0";
	}

	for(int i = 0; i < 2*num_arcs; i++)
	{
		file << " 0";
	}

	file << std::endl;

	file << "UPPER_BOUNDS" << std::endl;

	file << "1";
	for(int i = 1; i < num_vertices; i++)
	{
		file << " 1";
	}

	for(int i = 0; i < num_arcs; i++)
	{
		file << " 1";
	}

	for(int i = 0; i < num_arcs; i++)
	{
		file << " " << instance.limit();
	}

	file << std::endl;

	file << "INEQUALITIES_SECTION" << std::endl;

	int num_adj_arc = 0;
	file << "x1 == 1" << std::endl;
	file << "x" << num_vertices << " == 1" << std::endl;

	for(int i = 1; i < num_vertices - 1; i++)
	{
		num_adj_arc = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[i][j] != NULL)
			{
				file << "+x"<< 1 + num_vertices + graph->pos(i,j);
				num_adj_arc++;
			}
		}

		if (i <= num_mandatory)
		{
			if(num_adj_arc == 0) file << "x" << 1 +i << " == 0" << std::endl;
			else file << " == 1" << std::endl;
			file << "x" << 1 +i << " == 1" << std::endl;
		}
		else
		{
			if(num_adj_arc == 0) file << "x" << 1 +i << " == 0" << std::endl;
			else  file << "-x" << 1 +i << " == 0" << std::endl;
		}
	}


	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[0][j] != NULL) file << "+x" << 1 +num_vertices + graph->pos(0,j);
	}

	file << " == " << num_vehicles << std::endl;

	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[j][num_vertices-1] != NULL) file << "+x" << 1 +num_vertices + graph->pos(j,num_vertices-1);
	}

	file << " == " << num_vehicles << std::endl;

	int cont = 0;
	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[j][0] != NULL)
		{
			file << "+x" << 1 +num_vertices + graph->pos(j,0);
			cont++;
		}
	}

	if(cont != 0) file << " == 0" << std::endl;
	cont = 0;

	for(int j = 0; j < num_vertices; j++)
	{
		if((*(instance.graph()))[num_vertices-1][j] != NULL)
		{
			file << "+x" << 1 +num_vertices + graph->pos(num_vertices-1,j);
			cont++;
		}
	}

	if(cont != 0) file << " == 0" << std::endl;


	for(int i = 1; i < num_vertices-1; i++)
	{
		cont = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[i][j] != NULL)
			{
				file << "+x" << 1 +num_vertices + graph->pos(i,j);
				cont++;
			}

			if((*(instance.graph()))[j][i] != NULL)
			{
				file << "-x" << 1 +num_vertices + graph->pos(j,i);
				cont++;
			}
		}
		if(cont != 0) file << " == 0" << std::endl;
	}


	for(int i = 1; i < num_vertices; i++)
	{
		curr_arc = (*(instance.graph()))[0][i];

		if(curr_arc != NULL)
		{
			int coef = (curr_arc->dist())*K_PRECISION;
			file << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(0,i) << "-x" << 1 +num_vertices + num_arcs + graph->pos(0,i) << " == 0" << std::endl;

		}
	}

	for(int i = 1; i < num_vertices-1; i++)
	{
		cont = 0;
		for(int j = 0; j < num_vertices; j++)
		{
			if((*(instance.graph()))[j][i] != NULL)
			{
				file << "-x" << 1 +num_vertices + num_arcs + graph->pos(j,i);
				cont++;
			}

			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				file << "+x" << 1 +num_vertices + num_arcs + graph->pos(i,j);
				int coef = curr_arc->dist()*K_PRECISION;
				file << "-" << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(i,j);
				cont++;
			}
		}

		if(cont != 0) file << " == 0" << std::endl;
	}


	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				int coef = 0;
				coef = - (instance.limit()-Rn[j])*K_PRECISION;
				file << coef << "/" << K_PRECISION << "x" << 1 +num_vertices + graph->pos(i,j) << "+x" << 1 +num_vertices +  num_arcs + graph->pos(i,j) << " <= 0" << std::endl;
			}
		}
	}

	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			curr_arc = (*(instance.graph()))[i][j];
			if(curr_arc != NULL)
			{
				int coef = 0;
				coef = - (R0[i]+curr_arc->dist())*K_PRECISION;
				file << coef << "/" << K_PRECISION << "x" << 1 + num_vertices + graph->pos(i,j) << "+x" << 1 +num_vertices +  num_arcs + graph->pos(i,j) << " >= 0" << std::endl;
			}
		}
	}

	file << "END" << std::endl;
	file.close();
}

void CapacitatedMultiCommodity(Instance & inst, double * R, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vehicles);
	NumVar3Matrix f(env, num_vehicles);

	CapacitatedMultiCommodityForm(env,model,x,f,inst,R,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeCapacitatedMultiCommodity(env,model,x,f,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);
	env.end();
	//return sol;
}

void CapacitatedMultiTwoCommodity(Instance & inst, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();
	//int num_mandatory = graph->num_mandatory();

	//Solution<int>  *sol = new Solution<int>(num_vertices-1,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVarMatrix f_x0_source(env, num_vehicles);

	NumVar3Matrix x(env, num_vehicles);
	NumVar3Matrix f(env, num_vehicles);

	//IloNumVarArray y(env, num_vertices - num_mandatory-1, 0.0, 1.0, ILOINT);

	CapacitatedMultiTwoCommodityForm(env,model,x,f,f_x0_source, inst,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeCapacitatedMultiTwoCommodity(env,model,x,f,f_x0_source,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();

	//return sol;
}

void CapacitatedTwoCommodity(Instance & inst, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();

	//Solution<int> * sol = new Solution<int>(num_vertices-1,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	IloNumVarArray x0_source;
	IloNumVarArray x0_splited;

	if(solve_relax)
	{
		x0_source = IloNumVarArray(env, num_vertices, 0, 1, ILOFLOAT);
		x0_splited = IloNumVarArray(env, num_vertices, 0, 1, ILOFLOAT);
	}else
	{
		x0_source = IloNumVarArray(env, num_vertices, 0, 1, ILOINT);
		x0_splited = IloNumVarArray(env, num_vertices, 0, 1, ILOINT);
	}

	NumVarMatrix f_x0_source(env, 2);
	NumVarMatrix f_x0_splited(env, 2);

	//IloNumVarArray y(env, num_vertices - num_mandatory-1, 0.0, 1.0, ILOINT);
	NumVarMatrix x(env, num_vertices);
	NumVarMatrix f(env, num_vertices);

	CapacitatedTwoCommodityForm(env,model,x0_source,x0_splited,x,f,f_x0_source,f_x0_splited,inst,solve_relax);

	Matrix<int>* sol_matrix = new Matrix<int>(num_vertices,num_vertices,0);

	optimizeCapacitatedTwoCommodity(env,model,x0_source,x0_splited,x,f,f_x0_source,f_x0_splited,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;

	//sol_matrix->Print();
	delete sol_matrix;
	sol_matrix = NULL;
	//(sol->solution())->Print();

	env.end();

	//return sol;
}

void KulkarniBhave(Instance& inst, double * R0, double * Rn, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVarMatrix x(env, num_vertices);
	IloNumVarArray u(env, num_vertices, 0.0, IloInfinity, ILOFLOAT);

	KulkarniBhaveForm(env,model,x,u,inst,R0,Rn,solve_relax);

	Matrix<int>* sol_matrix = new Matrix<int>(num_vertices,num_vertices,0);

	optimizeKulkarniBhave(env,model,x,u,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;

	//sol_matrix->Print();
	//(sol->solution())->Print();

	delete sol_matrix;
	sol_matrix = NULL;

	env.end();

	//return sol;
}

std::list<UserCutGeneral*> * Bianchessi(Instance& inst, double* R0, double * Rn, double time_limit, Solution<int>* sol, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts = NULL)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_arcs = graph->num_arcs();
	int num_routes = inst.num_vehicles();

	std::list<UserCutGeneral*> * cuts = NULL;
	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	IloNumVarArray x(env);
	IloNumVarArray y(env);
	IloNumVarArray z(env, num_arcs, 0, IloInfinity, ILOFLOAT);
	IloNumVar slack(env, 0, num_routes, ILOFLOAT);

	if(solve_relax)
	{
		x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
		y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOFLOAT);
	}
	else
	{
		x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOINT);
		y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOINT);
	}

	BianchessiForm(env,model,slack,y,x,z,inst,R0,Rn,solve_relax);

	Matrix<int>* sol_matrix = new Matrix<int>(num_vertices,num_vertices,0);

	cuts = optimizeBianchessi(env,model,y,x,z,inst,*sol,sol_matrix,time_limit,solve_relax,callback,find_root_cuts, R0, Rn, initial_cuts);

	//std::cout << cost << std::endl;

	//sol_matrix->Print();
	//(sol->solution())->Print();

	delete sol_matrix;
	sol_matrix = NULL;

	env.end();
	return cuts;
}

void BianchessiIndexedInVehicles(Instance& inst, double * R0, double * Rn, double time_limit, Solution<int>* sol, bool solve_relax, bool callback = false)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vertices);
	NumVarMatrix z(env, num_vertices);
	//IloNumVarArray y;

	//if(solve_relax) y = IloNumVarArray(env, num_vertices, 0, 1, ILOFLOAT);
	//else y = IloNumVarArray(env, num_vertices, 0, 1, ILOINT);

	BianchessiIndexedInVehiclesForm(env,model,x,z,inst,R0,Rn,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeBianchessiIndexedInVehicles(env,model,x,z,inst,*sol,sol_matrix,time_limit,solve_relax,callback);

	//std::cout << cost << std::endl;

	//sol_matrix->Print();
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();

	//return sol;
}

void CapacitatedSingleCommodityIndexedInVehicles(Instance& inst, double * R0, double * Rn, double time_limit, Solution<int>* sol, bool solve_relax, bool callback = false)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vertices);
	NumVarMatrix f(env, num_vertices);

	CapacitatedSingleCommodityIndexedInVehiclesForm(env,model,x,f,inst,R0, Rn,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeCapacitatedSingleCommodityIndexedInVehicles(env,model,x,f,inst,*sol,sol_matrix,time_limit,solve_relax, callback);

	//std::cout << cost << std::endl;

	//sol_matrix->Print();
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();

	//return sol;
}

void PreprocessInstance(Instance& inst, double * R0, double * Rn, double time_limit, Solution<int>* sol, double primal_bound, std::list<UserCutGeneral*> * initial_cuts = NULL)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_arcs = graph->num_arcs();
	int num_mandatory = graph->num_mandatory();
	int num_routes = inst.num_vehicles();
	bool find_root_cuts = true;
	int num_vertices_removed = 0, num_arcs_removed = 0;

	std::list<UserCutGeneral*> * cuts = NULL;

	Solution<int> * sol_iter = new Solution<int>(num_vertices);

	IloEnv env;
	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());
	IloModel model(env);

	IloNumVar slack(env, 0, num_routes, ILOFLOAT);
	IloNumVarArray x(env);
	IloNumVarArray y(env);
	IloNumVarArray f(env, num_arcs, 0, IloInfinity, ILOFLOAT);

	x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
	y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOFLOAT);

	CapacitatedSingleCommodityForm(env,model,slack,y,x,f,inst,R0, Rn,true);

	cuts = optimizeCapacitatedSingleCommodity(cplex,env,model,y,x,f,inst,*sol_iter,NULL,time_limit,true,false,find_root_cuts, R0, Rn, NULL,NULL);

	if(sol_iter->is_feasible_)
	{
		IloNumArray x_red_costs(env);
		IloNumArray y_red_costs(env);

		IloNumArray x_values(env);
		IloNumArray y_values(env);

		IloCplex::BasisStatusArray y_statuses(env);
		IloCplex::BasisStatusArray x_statuses(env);

		cplex.getReducedCosts(y_red_costs,y);
		cplex.getValues(y_values,y);
		cplex.getBasisStatuses(y_statuses,y);

		for(int i = 0; i < num_vertices; i++)
		{
			//std::cout << "y[" << i << "]: " << y_values[i] << " " << y_red_costs[i] << " status: "  << cplex.getBasisStatus(y[i]) << std::endl;

			if( ((y_statuses)[i] == IloCplex::BasisStatus::AtLower) && double_equals(1.0*(y_values[i]),0.0) && double_less(sol_iter->lp_ - 1.0*(y_red_costs[i]), primal_bound)) num_vertices_removed++;
			//getchar();getchar();
		}

		cplex.getReducedCosts(x_red_costs,x);
		cplex.getValues(x_values,x);
		cplex.getBasisStatuses(x_statuses,x);

		for(int i = 0; i < num_arcs; i++)
		{
			//std::cout << "x[" << i << "]: " << x_values[i] << " " << x_red_costs[i] << " status: "  << cplex.getBasisStatus(x[i]) << std::endl;

			/*if(double_equals(1.0*(x_values[i]),0.0))
			{
			std::cout << sol_iter->lp_ - 1.0*(x_red_costs[i]) << " < " << primal_bound << std::endl;
		}*/

		if( ((x_statuses)[i] == IloCplex::BasisStatus::AtLower) && double_equals(1.0*(x_values[i]),0.0) && double_less(sol_iter->lp_ - 1.0*(x_red_costs[i]), primal_bound)) num_arcs_removed++;
		//getchar();getchar();
	}

	if(num_vertices_removed + num_arcs_removed > 0)
	{
		std::cout << "* num_vertices_removed: " << num_vertices_removed << std::endl;
		std::cout << "* num_arcs_removed: " << num_arcs_removed << std::endl;
		getchar();getchar();
	}

	x_red_costs.end();
	y_red_costs.end();

	x_values.end();
	y_values.end();

	x_statuses.end();
	y_statuses.end();
}

if(cuts != NULL)
{
	DeleteCuts(cuts);
	delete cuts;
	cuts  = NULL;
}

delete sol_iter;
sol_iter = NULL;
cplex.end();
env.end();
}

std::list<UserCutGeneral*> * CapacitatedSingleCommodityExtended(Instance& inst, double * R0, double * Rn, double time_limit, Solution<int>* sol, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts = NULL)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_arcs = graph->num_arcs();
	int num_vehicles = inst.num_vehicles();

	std::list<UserCutGeneral*> * cuts = NULL;

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	IloNumVarArray x(env);
	IloNumVarArray y(env);
	IloNumVarArray f(env, num_arcs, 0, IloInfinity, ILOFLOAT);
	if(solve_relax)
	{
		x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
		y = IloNumVarArray(env, num_vehicles*num_vertices, 0.0, 1.0, ILOFLOAT);
	}
	else
	{
		x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOINT);
		y = IloNumVarArray(env, num_vehicles*num_vertices, 0.0, 1.0, ILOINT);
	}

	IloNumVarArray f_art(env, num_vehicles*num_vertices, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray x_art(env, num_vehicles*num_vertices, 0, 1, ILOINT);
	IloNumVarArray y_art(env, num_vehicles*num_vehicles, 0, 1, ILOINT);

	CapacitatedSingleCommodityExtendedForm(env,model,y,x,f,y_art,x_art,f_art,inst,R0, Rn,solve_relax);

	Matrix<int>* sol_matrix = new Matrix<int>(num_vertices,num_vertices,0);
	cuts = optimizeCapacitatedSingleCommodityExtended(env,model,y,x,f,y_art,x_art,f_art,inst,*sol,sol_matrix,time_limit,solve_relax,callback,find_root_cuts, R0, Rn, initial_cuts);
	//std::cout << cost << std::endl;
	//sol_matrix->Print();
	//(sol->solution())->Print();

	delete sol_matrix;
	sol_matrix = NULL;

	env.end();

	return cuts;

	//return sol;
}

void ButtRyan94MultiTwoCommodity(Instance& inst, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices-1,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vehicles);
	NumVar3Matrix f(env, num_vehicles);

	ButtRyan94MultiTwoCommodityForm(env,model,x,f,inst,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeButtRyan94MultiTwoCommodity(env,model,x,f,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();
	//return sol;
}

void ButtRyan94TwoCommodity(Instance& inst, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices-1,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vehicles);
	NumVarMatrix f(env, num_vertices);

	ButtRyan94TwoCommodityForm(env,model,x,f,inst,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeButtRyan94TwoCommodity(env,model,x,f,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();
	//return sol;
}

void ButtRyan94MultiCommodity(Instance& inst, double * r, int type, double time_limit, Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vehicles);
	NumVar3Matrix f(env, num_vehicles);

	ButtRyan94MultiCommodityForm(env,model,x,f,inst,r,type,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeButtRyan94MultiCommodity(env,model,x,f,inst,*sol,sol_matrix,time_limit, solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();
	//return sol;
}

void ButtRyan94SingleCommodity(Instance& inst, double * r, int type, double time_limit,Solution<int>* sol, bool solve_relax)
{
	Graph * graph = inst.graph();
	int num_vertices = graph->num_vertices();
	int num_vehicles = inst.num_vehicles();

	//Solution<int> * sol = new Solution<int>(num_vertices,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

	IloEnv env;
	IloModel model(env);

	NumVar3Matrix x(env, num_vehicles);

	NumVarMatrix f(env, num_vertices);

	ButtRyan94SingleCommodityForm(env,model,x,f,inst,r,type,solve_relax);

	Matrix<int>** sol_matrix = Allocate3Matrix<int>(num_vehicles, num_vertices, num_vertices);

	optimizeButtRyan94SingleCommodity(env,model,x,f,inst,*sol,sol_matrix,time_limit,solve_relax);

	//std::cout << cost << std::endl;
	//Print3Matrix<int>(sol_matrix, num_vehicles);
	//(sol->solution())->Print();

	Delete3Matrix<int>(sol_matrix, num_vehicles);

	env.end();
}

Instance* GenerateFilesFromOPinstancesIter(std::string dir_path, std::string file_name, double mandatory_percentage)
{
	std::string curr_file = dir_path;
	curr_file.append(file_name);
	Graph * graph = NULL;
	Instance * new_instance = NULL;

	std::fstream file;
	file.open(curr_file.c_str(), std::fstream::in);
	if(!file.is_open()){
		std::cout << "Could not open file" << std::endl;
		throw 1;
		return NULL;
	}

	int num_vertices = 0, num_vehicles = 0, num_mandatory = 0;
	double limit = 0.0;

	/*std::cout << file_name << std::endl;
	std::cout << num_vertices << std::endl;
	std::cout << num_mandatory << std::endl;
	std::cout << num_vehicles << std::endl;
	std::cout << limit << std::endl;
	std::cout << "*******************" << std::endl;*/

	std::vector<std::pair<double,double>> coordinates;
	std::pair<double,double> curr_coordinate;
	int curr_profit;

	std::vector<int> * profits_vec = new std::vector<int>();

	file >> limit >> num_vehicles;

	// skips position 1 in coordinates vector: will be filled with last node (will be a mandatory node)
	while(file >> curr_coordinate.first)
	{
		num_vertices++;
		file >> curr_coordinate.second >> curr_profit;
		coordinates.push_back(curr_coordinate);
		profits_vec->push_back(curr_profit);
	}

	file.close();

	// relocates destination to the end of the vector
	curr_coordinate.first = coordinates[1].first;
	curr_coordinate.second = coordinates[1].second;
	curr_profit = (*profits_vec)[1];

	coordinates[1].first = coordinates[num_vertices-1].first;
	coordinates[1].second = coordinates[num_vertices-1].second;
	(*profits_vec)[1] = (*profits_vec)[num_vertices-1];

	coordinates[num_vertices-1].first = curr_coordinate.first;
	coordinates[num_vertices-1].second = curr_coordinate.second;
	(*profits_vec)[num_vertices-1] = curr_profit;

	num_mandatory = (int)ceil(mandatory_percentage * num_vertices);

	graph = new Graph(num_vertices,num_mandatory,(&((*profits_vec)[0])));

	/*for(int i = 0; i < num_vertices; i++)
	{
	std::cout << i << " : (" << coordinates[i].first << ", " << coordinates[i].second << ")";
	if(i <= num_mandatory) std::cout << " | profit: 0.0" << std::endl;
	else std::cout << " | profit: " << profits[i - num_mandatory-1] << std::endl;
	getchar();
}*/

for(int i = 0; i < num_vertices; i++)
{
	for(int j = i+1; j < num_vertices; j ++)
	{
		graph->AddEdge(i,j,euclidian_distance(coordinates[i], coordinates[j]));
	}
}

new_instance = new Instance(graph, num_vehicles,limit);

std::vector<bool> selected_vertices(num_vertices, false);
int iter_mandatory = 0, cont_mandatory = 0;
while(cont_mandatory < num_mandatory)
{
	iter_mandatory = rand()%(num_vertices-2) + 1;
	if(!selected_vertices[iter_mandatory])
	{
		selected_vertices[iter_mandatory] = true;
		cont_mandatory++;
		new_instance->mandatory_list_.push_back(iter_mandatory);
	}
}

return new_instance;
}

Instance* GenerateFilesFromTOPinstancesIter(std::string dir_path, std::string file_name, double mandatory_percentage)
{
	std::string curr_file = dir_path;
	curr_file.append(file_name);
	Graph * graph = NULL;
	Instance * new_instance = NULL;

	std::cout << file_name << std::endl;

	std::fstream file;
	file.open(curr_file.c_str(), std::fstream::in);
	if(!file.is_open()){
		std::cout << "Could not open file" << std::endl;
		throw 1;
		return NULL;
	}

	int num_vertices = 0, num_vehicles = 0, num_mandatory = 0;
	double limit = 0.0;

	std::string line;
	getline(file,line);

	std::size_t pos_1 = line.find_first_of(" ");

	std::stringstream parameter;
	parameter << line.substr(pos_1 + 1);
	parameter >> num_vertices;

	num_mandatory = (int)ceil(mandatory_percentage * num_vertices);

	getline(file,line);
	pos_1 = line.find_first_of(" ");

	parameter.clear();
	parameter << line.substr(pos_1 + 1);
	parameter >> num_vehicles;

	getline(file,line);
	pos_1 = line.find_first_of(" ");

	parameter.clear();
	parameter << line.substr(pos_1 + 1);
	parameter >> limit;

	/*std::cout << file_name << std::endl;
	std::cout << num_vertices << std::endl;
	std::cout << num_mandatory << std::endl;
	std::cout << num_vehicles << std::endl;
	std::cout << limit << std::endl;
	std::cout << "*******************" << std::endl;*/

	std::vector<std::pair<double,double>> coordinates(num_vertices);
	int * profits = new int[num_vertices];

	// skips position 1 in coordinates vector: will be filled with last node (will be a mandatory node)
	for(int i = 0; i < num_vertices; i++)
	{
		file >> (coordinates[i]).first >> (coordinates[i]).second >> profits[i];
		//std::cout << (coordinates[i]).first << " " << (coordinates[i]).second << " " << profits[i];
		//getchar();getchar();
	}

	file.close();

	graph = new Graph(num_vertices,num_mandatory,profits);

	/*for(int i = 0; i < num_vertices; i++)
	{
	std::cout << i << " : (" << coordinates[i].first << ", " << coordinates[i].second << ")";
	if(i <= num_mandatory) std::cout << " | profit: 0.0" << std::endl;
	else std::cout << " | profit: " << profits[i - num_mandatory-1] << std::endl;
	getchar();
}*/

for(int i = 0; i < num_vertices; i++)
{
	for(int j = i+1; j < num_vertices; j ++)
	{
		graph->AddEdge(i,j,euclidian_distance(coordinates[i], coordinates[j]));
	}
}

new_instance = new Instance(graph, num_vehicles,limit);

std::vector<bool> selected_vertices(num_vertices, false);
int iter_mandatory = 0, cont_mandatory = 0;
while(cont_mandatory < num_mandatory)
{
	iter_mandatory = rand()%(num_vertices-2) + 1;
	if(!selected_vertices[iter_mandatory])
	{
		selected_vertices[iter_mandatory] = true;
		cont_mandatory++;
		new_instance->mandatory_list_.push_back(iter_mandatory);
	}
}

return new_instance;
}

int GenerateFilesFromOPinstances(std::string dir_path, double mandatory_percentage)
{
	Instance * inst = NULL;
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (dir_path.c_str())) != NULL)
	{
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL)
		{
			std::string curr_file(ent->d_name);
			if((curr_file != ".")&&(curr_file != "..")&&( (curr_file.size() > 0)&& (curr_file[curr_file.size() -1] != '~') )) inst = GenerateFilesFromOPinstancesIter(dir_path, curr_file, mandatory_percentage);
			if(inst != NULL)
			{
				std::string folder;
				std::stringstream s_percentage;
				std::string percentage;
				std::fstream output;

				std::size_t pos = dir_path.find_first_of("OP");

				folder = dir_path.substr(0,pos);
				folder.append("ST");
				folder.append(dir_path.substr(pos));
				folder.append("Set_");

				s_percentage << mandatory_percentage;
				percentage = s_percentage.str();

				folder.append(percentage);
				folder.append("//");
				folder.append(curr_file);
				std::cout << folder << std::endl;
				inst->WriteToFile(folder,curr_file,mandatory_percentage);
				delete inst;
				inst = NULL;
			}
		}
		closedir (dir);
	}else
	{
		return -1;
	}
	return 0;
}

int GenerateFilesFromTOPinstances(std::string dir_path, double mandatory_percentage, bool repair = false)
{
	Instance * inst = NULL;
	DIR *dir;

	struct dirent *ent;
	if ((dir = opendir (dir_path.c_str())) != NULL)
	{
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL)
		{
			std::string curr_file(ent->d_name);
			if((curr_file != ".")&&(curr_file != "..")&&( (curr_file.size() > 0)&& (curr_file[curr_file.size() -1] != '~') )) inst = GenerateFilesFromTOPinstancesIter(dir_path, curr_file, mandatory_percentage);
			if(inst != NULL)
			{
				std::string folder;
				std::stringstream s_percentage;
				std::string percentage;
				std::fstream output;

				std::size_t pos = dir_path.find_first_of("TOP");

				folder = dir_path.substr(0,pos);
				folder.append("S");
				folder.append(dir_path.substr(pos));
				folder.append("Set_");

				s_percentage << mandatory_percentage;
				percentage = s_percentage.str();

				folder.append(percentage);
				folder.append("//");
				folder.append(curr_file);

				std::string folder2 = folder;
				folder2.replace(14,4,"STOP_backup");

				//std::cout << folder << std::endl;
				//std::cout << folder2 << std::endl;

				if(repair)
				{
					Instance inst_old(folder2);

					inst->mandatory_list_ = inst_old.mandatory_list_;
				}

				inst->WriteToFile(folder,curr_file,mandatory_percentage);
				//getchar(); getchar();
				delete inst;
				inst = NULL;
			}
		}
		closedir (dir);
	}else
	{
		return -1;
	}
	return 0;
}

void SolveCutAndBranch(std::vector<std::string>& instances)
{
	double time_limit = -1, original_time_limit = 0.0;
	int option = -1;
	Graph * graph = NULL;
	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

	std::string folder;
	std::string file_name;

	std::list<UserCutGeneral*> * cuts = NULL;

	std::cout << "Cut and Branch:" << std::endl;

	do{
		std::cout << "*******************************************" << std::endl
		<< " * Time limit in seconds (-1 to be unlimited): " << std::endl;
		std::cin >> time_limit;
	}while( double_less(time_limit,0.0) && !(double_equals(time_limit,-1)));

	original_time_limit = time_limit;
	do{
		std::cout << "*******************************************" << std::endl
		<< "Select cuts:" << std::endl
		<< " 0 - ALL" << std::endl
		<< " 1 - Generalized Connectivity Constraints (Grazia)" << std::endl
		<< " 2 - Flow bound cuts" << std::endl
		<< " 3 - Conflict cuts" << std::endl
		<< " 4 - LCI cuts" << std::endl
		<< " 5 - Clique Conflict cuts" << std::endl
		<< " 6 - Clique Conflict cuts Extended" << std::endl
		<< " 7 - LCI cuts Extended" << std::endl
		<< "Option: " ;
		std::cin >> option;
		switch(option)
		{
			case 0:{
				//for(int i = 0; i < K_NUM_TYPES_CALLBACKS; i++) (*CALLBACKS_SELECTION)[i] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
				(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
				(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
				(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
				(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;
				break;
			}
			case 1:{
				(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
				break;
			}
			case 2:{
				(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
				break;
			}
			case 3:{
				(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
				break;
			}
			case 4:{
				(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
				break;
			}
			case 5:{
				(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
				break;
			}
			case 6:{
				(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT_EXT] = true;
				break;
			}
			case 7:{
				(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT_EXT] = true;
				break;
			}
			default: option = -1; std::cout << "Invalid option!" << std::endl; break;
		}
	}while(option == -1);

	do{
		std::cout << "*******************************************" << std::endl
		<< " 0 - Run ALL" << std::endl
		<< " 1 - Capacitated single commodity" << std::endl
		<< " 2 - Run Bianchessi et al." << std::endl
		<< " 3 - Capacitated single commodity Extended" << std::endl
		<< "Option: " ;
		std::cin >> option;
		switch(option)
		{
			case 0:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;
					Instance inst(instances[i]);
					graph = inst.graph();
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);

					cuts = CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,true,NULL,NULL);

					if(time_limit != -1) time_limit = std::max(0.0, original_time_limit - sol->root_time_);
					sol->milp_time_ = sol->root_time_;

					CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false,cuts,NULL);

					sol->write_to_file("cb_csc3",folder,file_name);

					DeleteCuts(cuts);
					delete cuts;
					cuts = NULL;

					sol->reset();

					cuts = Bianchessi(inst,R,Rn,time_limit,sol,true,false,true);

					if(time_limit != -1) time_limit = std::max(0.0, original_time_limit - sol->root_time_);
					sol->milp_time_ = sol->root_time_;

					Bianchessi(inst,R,Rn,time_limit,sol,false,false,false,cuts);

					sol->write_to_file("cb_b1",folder,file_name);

					DeleteCuts(cuts);
					delete cuts;
					cuts = NULL;

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
			}
			case 1:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;

					Instance inst(instances[i]);
					graph = inst.graph();
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					/*inst.ComputeConflictGraph();
					sol->pre_processing_time_ = inst.time_spent_in_preprocessing_;
					sol->num_maximal_cliques_ = (inst.conflicts_list_).size();
					sol->write_to_file("find_all_maximal_cliques_tomita",folder,file_name);
					delete sol;
					sol = NULL;
					continue;*/

					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);

					cuts = CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,true,NULL,NULL);

					//double primal_bound = 0.0;
					//PreprocessInstance(inst,R,Rn,time_limit,sol,primal_bound,cuts);

					if(time_limit != -1) time_limit = std::max(0.0, original_time_limit - sol->root_time_);
					sol->milp_time_ = sol->root_time_;

					CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false,cuts,NULL);

					sol->write_to_file("teste",folder,file_name);

					DeleteCuts(cuts);
					delete cuts;
					cuts = NULL;

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
			}
			case 2:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;
					Instance inst(instances[i]);
					graph = inst.graph();
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);

					cuts = Bianchessi(inst,R,Rn,time_limit,sol,true,false,true);

					if(time_limit != -1) time_limit = std::max(0.0, original_time_limit - sol->root_time_);
					sol->milp_time_ = sol->root_time_;

					Bianchessi(inst,R,Rn,time_limit,sol,false,false,false,cuts);

					sol->write_to_file("stop_cb_b1",folder,file_name);

					DeleteCuts(cuts);
					delete cuts;
					cuts = NULL;

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
				case 3:{
					for (size_t i = 0; i < instances.size(); i++)
					{
						split_file_path(instances[i],folder,file_name);
						std::cout << "* " << file_name << std::endl;

						Instance inst(instances[i]);
						graph = inst.graph();
						Solution<int> * sol = new Solution<int>(graph->num_vertices());

						/*inst.ComputeConflictGraph();
						sol->pre_processing_time_ = inst.time_spent_in_preprocessing_;
						sol->num_maximal_cliques_ = (inst.conflicts_list_).size();
						sol->write_to_file("find_all_maximal_cliques_tomita",folder,file_name);
						delete sol;
						sol = NULL;
						continue;*/

						double * R = Dijkstra(graph,false,false);
						double * Rn = Dijkstra(graph,true,false);

						cuts = CapacitatedSingleCommodityExtended(inst,R,Rn,time_limit,sol,true,false,true);

						//double primal_bound = 0.0;
						//PreprocessInstance(inst,R,Rn,time_limit,sol,primal_bound,cuts);

						if(time_limit != -1) time_limit = std::max(0.0, original_time_limit - sol->root_time_);
						sol->milp_time_ = sol->root_time_;

						//CapacitatedSingleCommodityExtended(inst,R,Rn,time_limit,sol,false,false,false,cuts);

						sol->write_to_file("relax_clique_ccs_lcis_cb_csc3e",folder,file_name);

						DeleteCuts(cuts);
						delete cuts;
						cuts = NULL;

						delete sol;
						sol = NULL;

						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;
					}
					break;
				}
			}
			default: option = -1; std::cout << "Invalid option!" << std::endl; break;
		}
	}while(option == -1);

	if(cuts != NULL)
	{
		DeleteCuts(cuts);
		delete cuts;
		cuts = NULL;
	}
}

void SolveBranchAndCutCallback(std::vector<std::string>& instances)
{
	double time_limit = -1;
	int option = -1;
	Graph * graph = NULL;
	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

	std::string folder;
	std::string file_name;
	bool callback = true;

	std::cout << "Branch and Cut via Callback:" << std::endl;

	do{
		std::cout << "*******************************************" << std::endl
		<< " * Time limit in seconds (-1 to be unlimited): " << std::endl;
		std::cin >> time_limit;
	}while( double_less(time_limit,0.0) && !(double_equals(time_limit,-1)));

	do{
		std::cout << "*******************************************" << std::endl
		<< "Select callback:" << std::endl
		<< " 0 - ALL" << std::endl
		<< " 1 - Generalized Connectivity Constraints (Grazia)" << std::endl
		<< " 2 - Flow bound cuts" << std::endl
		<< " 3 - Flow bound cuts (by initial enumeration)" << std::endl
		<< " 4 - Conflict cuts" << std::endl
		<< " 5 - Cover bound Cuts" << std::endl
		<< " 6 - Path bound Cuts" << std::endl
		<< "Option: " ;
		std::cin >> option;
		switch(option)
		{
			case 0:{
				//for(int i = 0; i < K_NUM_TYPES_CALLBACKS; i++) (*CALLBACKS_SELECTION)[i] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
				(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
				//(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
				break;
			}
			case 1:{
				(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
				break;
			}
			case 2:{
				(*CALLBACKS_SELECTION)[K_TYPE_FLOW_BOUNDS_CUT] = true;
				break;
			}
			case 3:{
				(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
				callback = false;
				break;
			}
			case 4:{
				(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
				break;
			}
			case 5:{
				(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
				break;
			}
			case 6:{
				(*CALLBACKS_SELECTION)[K_TYPE_PATH_BOUND_CUT] = true;
				break;
			}
			default: option = -1; std::cout << "Invalid option!" << std::endl; break;
		}
	}while(option == -1);

	do{
		std::cout << "*******************************************" << std::endl
		<< " 0 - Run ALL" << std::endl
		<< " 1 - Capacitated single commodity" << std::endl
		<< " 2 - Run Bianchessi et al." << std::endl
		<< "Option: " ;
		std::cin >> option;
		switch(option)
		{
			case 0:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;
					Instance inst(instances[i]);
					graph = inst.graph();
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);

					CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,callback,false,NULL,NULL);
					//CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,callback,false);
					sol->write_to_file("bc_r1_csc3",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
			}
			case 1:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;
					Instance inst(instances[i]);
					graph = inst.graph();
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);

					CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,true,false,NULL,NULL);
					sol->write_to_file("bc_new_csc3",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				std::cout << "terminou!" << std::endl;
				break;
			}
			case 2:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);

					std::cout << "* " << file_name << std::endl;

					Instance inst(instances[i]);

					graph = inst.graph();
					double * R = Dijkstra(graph,false,false);
					double * Rn = Dijkstra(graph,true,false);
					Solution<int> * sol = new Solution<int>(graph->num_vertices());

					Bianchessi(inst,R,Rn,time_limit,sol,false,true,false);
					sol->write_to_file("bc_b1",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
			}
			default: option = -1; std::cout << "Invalid option!" << std::endl; break;
		}
	}while(option == -1);
}

void SolveCompactFormulations(std::vector<std::string>& instances)
{
	double time_limit = -1;
	int option = -1;
	Graph * graph = NULL;
	Solution<int> * sol = NULL;
	double * r = NULL, *R = NULL, *Rn = NULL;

	std::string folder;
	std::string file_name;

	std::cout << "COMPACT FORMULATIONS:" << std::endl;

	do{
		std::cout << "*******************************************" << std::endl
		<< " * Time limit in seconds (-1 to be unlimited): " << std::endl;
		std::cin >> time_limit;
	}while( double_less(time_limit,0.0) && !(double_equals(time_limit,-1)));

	do{
		std::cout << "*******************************************" << std::endl
		<< " 0 - Run ALL" << std::endl
		<< " 1 - Run Butt&Ryan single commodity" << std::endl
		<< " 2 - Run Butt&Ryan multi commodity" << std::endl
		<< " 3 - Run Butt&Ryan two commodity" << std::endl
		<< " 4 - Run Butt&Ryan multi two commodity" << std::endl
		<< " 5 - Run Capacitated single commodity" << std::endl
		<< " 6 - Run Capacitated multi commodity" << std::endl
		<< " 7 - Run Capacitated two commodity" << std::endl
		<< " 8 - Run Capacitated multi two commodity" << std::endl
		<< " 9 - Run Kulkarni and Bhave" << std::endl
		<< " 10 - Run Bianchessi et al." << std::endl
		<< " 11 - Run Capacitated single commodity Extended" << std::endl
		<< "Option: " ;
		std::cin >> option;
		switch(option){
			case 0:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					try{
						split_file_path(instances[i],folder,file_name);
						std::cout << "* " << file_name << std::endl;
						Instance inst(instances[i]);

						graph = inst.graph();

						sol = new Solution<int>(graph->num_vertices());

						r = Dijkstra(graph,false,true);
						R = Dijkstra(graph,false,false);
						Rn = Dijkstra(graph,true,false);

						//ButtRyan94SingleCommodity(inst,NULL,0,time_limit,sol,true);
						//ButtRyan94SingleCommodity(inst,NULL,0,time_limit,sol,false);
						//sol->write_to_file("bcsc1",folder,file_name);
						//sol->reset();
						//ButtRyan94SingleCommodity(inst,r,0,time_limit,sol,true);
						//ButtRyan94SingleCommodity(inst,r,0,time_limit,sol,false);
						//sol->write_to_file("bcsc2",folder,file_name);
						//sol->reset();
						//ButtRyan94SingleCommodity(inst,NULL,1,time_limit,sol,true);
						//ButtRyan94SingleCommodity(inst,NULL,1,time_limit,sol,false);
						//sol->write_to_file("bcsc3",folder,file_name);
						//sol->reset();
						/*ButtRyan94SingleCommodity(inst,r,1,time_limit,sol,true);
						ButtRyan94SingleCommodity(inst,r,1,time_limit,sol,false);
						sol->write_to_file("bcsc4",folder,file_name);
						sol->reset();*/

						//ButtRyan94MultiCommodity(inst,NULL,0,time_limit,sol,true);
						//ButtRyan94MultiCommodity(inst,NULL,0,time_limit,sol,false);
						//sol->write_to_file("bcmc1",folder,file_name);
						//sol->reset();
						//ButtRyan94MultiCommodity(inst,r,0,time_limit,sol,true);
						//ButtRyan94MultiCommodity(inst,r,0,time_limit,sol,false);
						//sol->write_to_file("bcmc2",folder,file_name);
						//sol->reset();
						//ButtRyan94MultiCommodity(inst,NULL,1,time_limit,sol,true);
						//ButtRyan94MultiCommodity(inst,NULL,1,time_limit,sol,false);
						//sol->write_to_file("bcmc3",folder,file_name);
						//sol->reset();
						/*ButtRyan94MultiCommodity(inst,r,1,time_limit,sol,true);
						ButtRyan94MultiCommodity(inst,r,1,time_limit,sol,false);
						sol->write_to_file("bcmc4",folder,file_name);
						sol->reset();

						ButtRyan94TwoCommodity(inst,time_limit,sol,true);
						ButtRyan94TwoCommodity(inst,time_limit,sol,false);
						sol->write_to_file("bctc",folder,file_name);
						sol->reset();

						ButtRyan94MultiTwoCommodity(inst,time_limit,sol,true);
						ButtRyan94MultiTwoCommodity(inst,time_limit,sol,false);
						sol->write_to_file("bcmtc",folder,file_name);
						sol->reset();*/

						//CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,true);
						//CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,false);
						//sol->write_to_file("csc1",folder,file_name);
						//sol->reset();

						//CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,true);
						//CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,false);
						//sol->write_to_file("csc2",folder,file_name);
						//sol->reset();

						CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,false,NULL,NULL);
						CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false,NULL,NULL);
						sol->write_to_file("csc3",folder,file_name);
						sol->reset();

						/*CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,true);
						CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("csc4",folder,file_name);
						sol->reset();*/

						/*CapacitatedMultiCommodity(inst,NULL,time_limit,sol,true);
						CapacitatedMultiCommodity(inst,NULL,time_limit,sol,false);
						sol->write_to_file("cmc1",folder,file_name);
						sol->reset();*/

						/*CapacitatedMultiCommodity(inst,R,time_limit,sol,true);
						CapacitatedMultiCommodity(inst,R,time_limit,sol,false);
						sol->write_to_file("cmc2",folder,file_name);
						sol->reset();

						CapacitatedTwoCommodity(inst,time_limit,sol,true);
						CapacitatedTwoCommodity(inst,time_limit,sol,false);
						sol->write_to_file("ctc1",folder,file_name);
						sol->reset();

						CapacitatedMultiTwoCommodity(inst,time_limit,sol,true);
						CapacitatedMultiTwoCommodity(inst,time_limit,sol,false);
						sol->write_to_file("cmtc1",folder,file_name);
						sol->reset();

						KulkarniBhave(inst,R,Rn,time_limit,sol,true);
						KulkarniBhave(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("kb",folder,file_name);
						sol->reset();*/

						/*Bianchessi(inst,R,Rn,time_limit,sol,true,false,false);
						//Bianchessi(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("b1_relax",folder,file_name);
						sol->reset();*/

						/*BianchessiIndexedInVehicles(inst,R,Rn,time_limit,sol,true);
						BianchessiIndexedInVehicles(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("b2",folder,file_name);*/

						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;

						delete [] r;
						r = NULL;

						delete sol;
						sol = NULL;
					}
					catch (const std::bad_alloc& ex)
					{
						sol->out_of_memory_ = true;
						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;

						delete [] r;
						r = NULL;

						delete sol;
						sol = NULL;
						std::cout << "Out of memory exception" << std::endl;
						continue;
					}
				}
				break;
			}
			case 1:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					std::cout << " * " << file_name << std::endl;
					graph = inst.graph();
					r = Dijkstra(graph,false,true);
					sol = new Solution<int>(graph->num_vertices());

					/*ButtRyan94SingleCommodity(inst,NULL,0,time_limit,sol,true);
					ButtRyan94SingleCommodity(inst,NULL,0,time_limit,sol,false);
					sol->write_to_file("bcsc1",folder,file_name);
					sol->reset();

					ButtRyan94SingleCommodity(inst,r,0,time_limit,sol,true);
					ButtRyan94SingleCommodity(inst,r,0,time_limit,sol,false);
					sol->write_to_file("bcsc2",folder,file_name);
					sol->reset();

					ButtRyan94SingleCommodity(inst,NULL,1,time_limit,sol,true);
					ButtRyan94SingleCommodity(inst,NULL,1,time_limit,sol,false);
					sol->write_to_file("bcsc3",folder,file_name);
					sol->reset();*/

					ButtRyan94SingleCommodity(inst,r,1,time_limit,sol,true);
					ButtRyan94SingleCommodity(inst,r,1,time_limit,sol,false);
					sol->write_to_file("stop_bcsc4",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] r;
					r = NULL;
				}
				break;
			}
			case 2:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					graph = inst.graph();
					r = Dijkstra(graph,false,true);

					sol = new Solution<int>(graph->num_vertices());

					ButtRyan94MultiCommodity(inst,NULL,0,time_limit,sol,true);
					ButtRyan94MultiCommodity(inst,NULL,0,time_limit,sol,false);
					sol->write_to_file("bcmc1",folder,file_name);
					sol->reset();

					ButtRyan94MultiCommodity(inst,r,0,time_limit,sol,true);
					ButtRyan94MultiCommodity(inst,r,0,time_limit,sol,false);
					sol->write_to_file("bcmc2",folder,file_name);
					sol->reset();

					ButtRyan94MultiCommodity(inst,NULL,1,time_limit,sol,true);
					ButtRyan94MultiCommodity(inst,NULL,1,time_limit,sol,false);
					sol->write_to_file("bcmc3",folder,file_name);
					sol->reset();

					ButtRyan94MultiCommodity(inst,r,1,time_limit,sol,true);
					ButtRyan94MultiCommodity(inst,r,1,time_limit,sol,false);
					sol->write_to_file("bcmc4",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] r;
					r = NULL;
				}
				break;
			}
			case 3:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					ButtRyan94TwoCommodity(inst,time_limit,sol,true);
					ButtRyan94TwoCommodity(inst,time_limit,sol,false);
					sol->write_to_file("bctc",folder,file_name);

					delete sol;
					sol = NULL;
				}
				break;
			}
			case 4:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					ButtRyan94MultiTwoCommodity(inst,time_limit,sol,true);
					ButtRyan94MultiTwoCommodity(inst,time_limit,sol,false);
					sol->write_to_file("bcmtc",folder,file_name);

					delete sol;
					sol = NULL;
				}
				break;
			}
			case 5:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					std::cout << "* " << file_name << std::endl;
					Instance inst(instances[i]);
					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					R = Dijkstra(graph,false,false);
					Rn = Dijkstra(graph,true,false);
					/*CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,true);
					CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,false);
					sol->write_to_file("csc1",folder,file_name);
					sol->reset();

					CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,true);
					CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,false);
					sol->write_to_file("csc2",folder,file_name);
					sol->reset();*/

					CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,false,NULL,NULL);
					//CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false);
					sol->write_to_file("relax_csc3_reinforced2",folder,file_name);
					sol->reset();

					/*CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,true);
					CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,false);
					sol->write_to_file("csc4",folder,file_name);*/

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;

					delete [] Rn;
					Rn = NULL;
				}
				break;
			}
			case 6:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					R = Dijkstra(graph,false,false);
					CapacitatedMultiCommodity(inst,NULL,time_limit,sol,true);
					CapacitatedMultiCommodity(inst,NULL,time_limit,sol,false);
					sol->write_to_file("cmc1",folder,file_name);
					sol->reset();

					CapacitatedMultiCommodity(inst,R,time_limit,sol,true);
					CapacitatedMultiCommodity(inst,R,time_limit,sol,false);
					sol->write_to_file("cmc2",folder,file_name);

					delete sol;
					sol = NULL;

					delete [] R;
					R = NULL;
				}
				break;
			}
			case 7:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);
					Instance inst(instances[i]);
					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					CapacitatedTwoCommodity(inst,time_limit,sol,true);
					CapacitatedTwoCommodity(inst,time_limit,sol,false);
					sol->write_to_file("ctc1",folder,file_name);

					delete sol;
					sol = NULL;
				}
				break;
			}
			case 8:{
				for (size_t i = 0; i < instances.size(); i++)
				{
					split_file_path(instances[i],folder,file_name);

					//std::cout << folder << std::endl;
					//std::cout << file_name << std::endl;

					Instance inst(instances[i]);

					graph = inst.graph();
					sol = new Solution<int>(graph->num_vertices());

					CapacitatedMultiTwoCommodity(inst,time_limit,sol,true);
					CapacitatedMultiTwoCommodity(inst,time_limit,sol,false);
					sol->write_to_file("cmtc1",folder,file_name);

					delete sol;
					sol = NULL;
				}
				break;
				case 9:
				{
					for (size_t i = 0; i < instances.size(); i++)
					{
						split_file_path(instances[i],folder,file_name);

						//std::cout << folder << std::endl;
						std::cout << file_name << std::endl;

						Instance inst(instances[i]);

						graph = inst.graph();
						sol = new Solution<int>(graph->num_vertices());
						R = Dijkstra(graph,false,false);
						Rn = Dijkstra(graph,true,false);

						KulkarniBhave(inst,R,Rn,time_limit,sol,true);
						KulkarniBhave(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("stop_kb",folder,file_name);

						delete sol;
						sol = NULL;

						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;
					}
				}
				break;
				case 10:
				{
					for (size_t i = 0; i < instances.size(); i++)
					{
						split_file_path(instances[i],folder,file_name);

						std::cout << "* " << file_name << std::endl;

						Instance inst(instances[i]);

						graph = inst.graph();
						R = Dijkstra(graph,false,false);
						Rn = Dijkstra(graph,true,false);
						sol = new Solution<int>(graph->num_vertices());

						Bianchessi(inst,R,Rn,time_limit,sol,true,false,false);
						//Bianchessi(inst,R,Rn,time_limit,sol,false,false,false);
						sol->write_to_file("stop_relax_b1",folder,file_name);
						//sol->reset();

						/*BianchessiIndexedInVehicles(inst,R,Rn,time_limit,sol,true);
						BianchessiIndexedInVehicles(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("b2",folder,file_name);*/

						delete sol;
						sol = NULL;

						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;
					}
				}
				break;
				case 11:{
					for (size_t i = 0; i < instances.size(); i++)
					{
						split_file_path(instances[i],folder,file_name);
						std::cout << "* " << file_name << std::endl;
						Instance inst(instances[i]);
						graph = inst.graph();
						sol = new Solution<int>(graph->num_vertices());

						R = Dijkstra(graph,false,false);
						Rn = Dijkstra(graph,true,false);
						/*CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,true);
						CapacitatedSingleCommodity(inst,NULL,NULL,time_limit,sol,false);
						sol->write_to_file("csc1",folder,file_name);
						sol->reset();

						CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,true);
						CapacitatedSingleCommodity(inst,R,NULL,time_limit,sol,false);
						sol->write_to_file("csc2",folder,file_name);
						sol->reset();*/

						CapacitatedSingleCommodityExtended(inst,R,Rn,time_limit,sol,true,false,false);
						//CapacitatedSingleCommodityExtended(inst,R,Rn,time_limit,sol,false,false,false);
						sol->write_to_file("relax_csc3e",folder,file_name);
						sol->reset();

						/*CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,true);
						CapacitatedSingleCommodityIndexedInVehicles(inst,R,Rn,time_limit,sol,false);
						sol->write_to_file("csc4",folder,file_name);*/

						delete sol;
						sol = NULL;

						delete [] R;
						R = NULL;

						delete [] Rn;
						Rn = NULL;
					}
					break;
				}
			}
			default: option = -1; std::cout << "Invalid option!" << std::endl; break;
		}
	}while(option == -1);
}

void GenrateRelaxCSVTable(std::vector<std::string> dirs, bool stop, double total_time_limit)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_none");
	/*algorithms.push_back("relax_all_gccs");
	//algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_all_ccs");
	//algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_all_gccs_ccs");
	//algorithms.push_back("relax_gccs_ccs_cb_csc3");
	//algorithms.push_back("relax_all_gccs_ccs_lcis");
	//algorithms.push_back("relax_clique_ccs_active_vertices");
	algorithms.push_back("relax_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_clique_ccs_total_active_vertices_unilateral");
	algorithms.push_back("relax_maximum_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_lcis_clique_ccs_active_vertices_unilateral");

	algorithms.push_back("relax_angle_0.1_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.05_clique_ccs_active_vertices_unilateral");*/
	algorithms.push_back("cb_csc3");
	algorithms.push_back("relax_angle_0.03_clique_ccs_active_vertices_unilateral");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//CSV//stop_relax_table.csc";
	else output_name = ".//tables//CSV//relax_table.csv";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_cuts_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(),0.0);
	std::vector<double> total_avg_cuts(algorithms.size(),0.0);
	int total_num_instances = 0;

	output << "set";
	for(size_t j = 0; j < algorithms.size(); j++)
	{
		output << ";" << algorithms[j] << ";" << algorithms[j];
	}
	output << std::endl;
	for(size_t j = 0; j < algorithms.size(); j++)
	{
		output << ";tempo(s);#cuts";
	}
	output << std::endl;
	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> time_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> cuts_per_algo(algorithms.size(),std::vector<double>());
		std::vector<double> avg_time(algorithms.size(),0.0);
		std::vector<double> avg_cuts(algorithms.size(),0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

		total_num_instances+=(instances.size());
		for(size_t i = 0; i < instances.size(); i++)
		{
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_time;
				std::string status;
				double time = 0.0;
				int num_cuts = 0;
				std::string line;

				for(int i = 0; i < 5; i++) getline(input,line);

				size_t pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;
				//std::cout << line << "    " << time << std::endl;
				//getchar(); getchar();

				for(int i = 0; i < 4; i++) getline(input,line);

				for(int i = 0; i < K_NUM_TYPES_CALLBACKS; i++)
				{
					std::stringstream s_cuts;
					int cuts_iter = 0;
					getline(input,line);
					pos = line.find_first_of(":");
					size_t pos2 = line.find_first_of("/");
					s_cuts << line.substr(pos + 2, pos2-pos -2);
					s_cuts >> cuts_iter;
					num_cuts+=cuts_iter;

					//std::cout << line << "    " << cuts_iter << std::endl;
					//getchar(); getchar();
				}
				time_per_algo[j].push_back(time);
				total_time_per_algo[j].push_back(time);
				total_avg_time[j] += time;
				avg_time[j] += time;

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				cuts_per_algo[j].push_back(num_cuts);
				total_cuts_per_algo[j].push_back(num_cuts);

				total_avg_cuts[j] += num_cuts;

				avg_cuts[j] += num_cuts;
				input.close();
			}

			//getchar();getchar();
		}

		output << j+1;

		for(size_t j = 0; j < algorithms.size(); j++)
		{
			if((time_per_algo[j]).size() > 0) avg_time[j]/=(1.0*((time_per_algo[j]).size()));
			else avg_time[j] = -1;
			avg_cuts[j]/=(1.0*((cuts_per_algo[j]).size()));
			output << ";" << avg_time[j] << ";" << avg_cuts[j];
		}
		output << std::endl;
	}

	output << "Total";
	for(size_t j = 0; j < algorithms.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_time[j]/=(1.0*((total_time_per_algo[j]).size()));
		total_avg_cuts[j]/=(1.0*((total_cuts_per_algo[j]).size()));
		output << ";" << total_avg_time[j] << ";" << total_avg_cuts[j];
	}

	output << std::endl;

	output.close();
}

void GenerateAlgorithmsLatexTable(std::vector<std::string> dirs, bool stop, double total_time_limit)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	//algorithms.push_back("bc_b1");
	//algorithms.push_back("csc3");
	//algorithms.push_back("csc3_2");
	//algorithms.push_back("cb_csc3");
	//algorithms.push_back("cb2_csc3");

	/*std::string algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);*/

	std::string algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_algorithms_new.txt";
	else output_name = ".//tables//latex//table_algorithms_new.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(),0.0);
	std::vector<double> total_avg_gap(algorithms.size(),0.0);
	std::vector<int> total_num_optimal(algorithms.size(),0);
	int total_num_instances = 0;

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> time_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo(algorithms.size(),std::vector<double>());
		std::vector<double> avg_time(algorithms.size(),0.0);
		std::vector<double> avg_gap(algorithms.size(),0.0);
		std::vector<double> st_dev(algorithms.size(),0.0);

		std::vector<int> num_optimal(algorithms.size(),0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

		total_num_instances+=(instances.size());
		for(size_t i = 0; i < instances.size(); i++)
		{
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb, s_ub, s_time;
				std::string status;
				double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if(s_lb.str() == "-inf") lb = -1;
				else s_lb >> lb;

				getline(input,line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if(s_ub.str() == "inf") ub = -1;
				else s_ub >> ub;

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if((status != "OPTIMAL") && (status != "INFEASIBLE"))
				{
					if((!(double_equals(lb,-1))) && (!(double_equals(ub,-1))))
					{
						if(double_less(time,total_time_limit))
						{
							(num_optimal[j])++;
							(total_num_optimal[j])++;
							time_per_algo[j].push_back(time);
							total_time_per_algo[j].push_back(time);

							total_avg_time[j] += time;
							avg_time[j] += time;
						}else gap = (100.0*(ub-lb))/ub;
					}else gap = 100.0;
				}else
				{
					(num_optimal[j])++;
					(total_num_optimal[j])++;
					time_per_algo[j].push_back(time);
					total_time_per_algo[j].push_back(time);

					total_avg_time[j] += time;
					avg_time[j] += time;
				}

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				if(!double_equals(gap,0.0))
				{
					gap_per_algo[j].push_back(gap);
					total_gap_per_algo[j].push_back(gap);

					total_avg_gap[j] += gap;

					avg_gap[j] += gap;
				}
				input.close();
			}

			//getchar();getchar();
		}

		output << j+1;

		if(stop) output << "\\_5\\%";

		for(size_t j = 0; j < algorithms.size(); j++)
		{
			if((time_per_algo[j]).size() > 0) avg_time[j]/=(1.0*((time_per_algo[j]).size()));
			else avg_time[j] = -1;
			if((gap_per_algo[j]).size() > 0)
			{
				avg_gap[j]/=(1.0*((gap_per_algo[j]).size()));
				st_dev[j] = StDev(gap_per_algo[j],avg_gap[j]);
			}else avg_gap[j] = st_dev[j] = -1;
			output << " & & " << num_optimal[j] << "/" << instances.size() << " & " << avg_time[j] << " & " << avg_gap[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t j = 0; j < algorithms.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		if((total_time_per_algo[j]).size() > 0) total_avg_time[j]/=(1.0*((total_time_per_algo[j]).size()));
		else total_avg_time[j] = -1;
		if((total_gap_per_algo[j]).size() > 0) total_avg_gap[j]/=(1.0*((total_gap_per_algo[j]).size()));
		else total_avg_gap[j] = -1;
		output << "& & " << total_num_optimal[j] << "/" << total_num_instances << " & " << total_avg_time[j] << " & " << total_avg_gap[j] << " & " << StDev(total_gap_per_algo[j],total_avg_gap[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateValidIneqsActiveLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_b1");
	algorithms.push_back("relax_none_cb_csc3");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_valid_ineq_active.txt";
	else output_name = ".//tables//latex//table_valid_ineq_active.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_percentage_active_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_percentage_active(algorithms.size(),0.0);

	std::vector<int> total_num_active_per_algo(algorithms.size(), 0);
	std::vector<int> total_count(algorithms.size(),0);

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> percentage_active_per_algo(algorithms.size(),std::vector<double>());
		std::vector<double> avg_percentage_active(algorithms.size(),0.0);
		std::vector<double> st_dev(algorithms.size(),0.0);

		std::vector<int> num_active_per_algo(algorithms.size(), 0);
		std::vector<int> count(algorithms.size(),0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);
		for(size_t i = 0; i < instances.size(); i++)
		{
			//std::cout << "* " << instances[i] << std::endl;
			int num_active = 0, num_total = 0;
			double percentage_active = 0.0;
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_active, s_total;
				int curr_active = 0, curr_total = 0;
				std::string line;

				for(int k = 0; k < 12; k++) getline(input,line);
				size_t pos = line.find_first_of(":"), pos2 = line.find_first_of("/");

				s_active << line.substr(pos + 2,pos2 - pos -2);
				s_active >> num_active;

				s_total << line.substr(pos2 + 1);
				s_total >> num_total;

				//std::cout << num_active << "/" << num_total << std::endl;
				//getchar();getchar();

				if(num_total != 0)
				{
					percentage_active = (num_active*1.0)/num_total;
					percentage_active_per_algo[j].push_back(percentage_active);
					total_percentage_active_per_algo[j].push_back(percentage_active);
					total_avg_percentage_active[j] += percentage_active;
					avg_percentage_active[j] += percentage_active;

					(num_active_per_algo[j])+=num_active;
					(total_num_active_per_algo[j])+=num_active;
					(count[j])++;
					(total_count[j])++;
				}

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;
				input.close();
			}
		}

		output << j+1;

		for(size_t j = 0; j < algorithms.size(); j++)
		{
			avg_percentage_active[j]/=(1.0*((percentage_active_per_algo[j]).size()));
			st_dev[j] = StDev(percentage_active_per_algo[j],avg_percentage_active[j]);
			output << " & & " << (1.0*num_active_per_algo[j])/(count[j]) << " & " << avg_percentage_active[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";

	for(size_t j = 0; j < algorithms.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_percentage_active[j]/=(1.0*((total_percentage_active_per_algo[j]).size()));
		output << "& & " << (1.0*total_num_active_per_algo[j])/(total_count[j]) << " & " << total_avg_percentage_active[j] << " & " << StDev(total_percentage_active_per_algo[j],total_avg_percentage_active[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateValidIneqsImprovementsLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;
	std::vector<std::string> algorithms2;

	algorithms.push_back("relax_cb_csc3");
	algorithms.push_back("relax_FBs_cb_csc3");

	algorithms2.push_back("relax_b1_sem");
	algorithms2.push_back("relax_b1");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_valid_ineq_improvements.txt";
	else output_name = ".//tables//latex//table_valid_ineq_improvements.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(),std::vector<double>()), total_improvement_per_algo2(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(),0.0), total_avg_improvement2(algorithms.size(),0.0);

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> improvement_per_algo(algorithms.size(),std::vector<double>()), improvement_per_algo2(algorithms.size(),std::vector<double>());
		std::vector<double> avg_improvement(algorithms.size(),0.0),avg_improvement2(algorithms.size(),0.0);
		std::vector<double> st_dev(algorithms.size(),0.0), st_dev2(algorithms.size(),0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);
		for(size_t i = 0; i < instances.size(); i++)
		{
			double original_lp = 0.0, original_lp2 = 0.0;;
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);
				s_lp >> lp;

				if(j == 0)
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE"))
					{
						original_lp = -1.0;
					}else original_lp = lp;
					curr_improvement = 0.0;
				}else
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE") || (double_equals(original_lp,-1)))
					{
						curr_improvement = 0.0;
					}else
					{
						if(double_equals(original_lp,0.0)) curr_improvement = 0;
						else curr_improvement = (100*(original_lp-lp))/original_lp;
					}
				}

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				improvement_per_algo[j].push_back(curr_improvement);
				total_improvement_per_algo[j].push_back(curr_improvement);
				total_avg_improvement[j] += curr_improvement;
				avg_improvement[j] += curr_improvement;
				input.close();
			}

			for(size_t j = 0; j < algorithms2.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms2[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);
				s_lp >> lp;

				if(j == 0)
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE"))
					{
						original_lp2 = -1.0;
					}else original_lp2 = lp;
					curr_improvement = 0.0;
				}else
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE") || (double_equals(original_lp2,-1)))
					{
						curr_improvement = 0.0;
					}else
					{
						if(double_equals(original_lp2,0.0)) curr_improvement = 0;
						else curr_improvement = (100*(original_lp2-lp))/original_lp2;
					}


				}

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				improvement_per_algo2[j].push_back(curr_improvement);
				total_improvement_per_algo2[j].push_back(curr_improvement);
				total_avg_improvement2[j] += curr_improvement;
				avg_improvement2[j] += curr_improvement;
				input.close();
			}

			//getchar();getchar();
		}

		output << j+1;

		for(size_t j = 1; j < algorithms2.size(); j++)
		{
			avg_improvement2[j]/=(1.0*(instances.size()));
			st_dev2[j] = StDev(improvement_per_algo2[j],avg_improvement2[j]);
			output << " & & "<< avg_improvement2[j] << " & " << st_dev2[j];
		}
		for(size_t j = 1; j < algorithms.size(); j++)
		{
			avg_improvement[j]/=(1.0*(instances.size()));
			st_dev[j] = StDev(improvement_per_algo[j],avg_improvement[j]);
			output << " & & "<< avg_improvement[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t j = 1; j < algorithms2.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement2[j]/=(1.0*((total_improvement_per_algo2[j]).size()));
		output << "& & "<< total_avg_improvement2[j] << " & " << StDev(total_improvement_per_algo2[j],total_avg_improvement2[j]);
	}

	for(size_t j = 1; j < algorithms.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement[j]/=(1.0*((total_improvement_per_algo[j]).size()));
		output << "& & "<< total_avg_improvement[j] << " & " << StDev(total_improvement_per_algo[j],total_avg_improvement[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateGapsNodesCutsLatexTable(std::vector<std::string> dirs, bool stop, double time_limit)
{
	//std::cout << "3" << std::endl; getchar(); getchar();
	std::vector<std::string> algorithms;

	algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	algorithms.push_back("cb2_csc3");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_gaps_nodes_cuts.txt";
	else output_name = ".//tables//latex//table_gaps_nodes_cuts.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	for(size_t k = 0; k < algorithms.size(); k++)
	{
		std::vector<double> total_gaps;
		double total_avg_gaps = 0.0;
		std::vector<double> total_nodes;
		double total_avg_nodes = 0.0;
		std::vector<std::vector<double>> total_cuts(K_NUM_TYPES_CALLBACKS,std::vector<double>());
		std::vector<double> total_avg_cuts(K_NUM_TYPES_CALLBACKS,0.0);

		output << algorithms[k] << std::endl;

		for(size_t j = 0; j < dirs.size(); j++)
		{
			std::vector<double> gaps;
			double avg_gaps = 0.0;
			std::vector<double> nodes;
			double avg_nodes = 0.0;
			std::vector<std::vector<double>> cuts(K_NUM_TYPES_CALLBACKS,std::vector<double>());
			std::vector<double> avg_cuts(K_NUM_TYPES_CALLBACKS,0.0);

			std::vector<std::string> instances;
			std::string folder = dirs[j].substr(20);
			AddInstancesFromDirectory(dirs[j],instances,false);

			for(size_t i = 0; i < instances.size(); i++)
			{
				double lb = 0.0, ub = 0.0;
				double curr_gap = 0.0, curr_nodes = 0.0, curr_cuts = 0.0, time = 0.0;

				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[k]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;

				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_nodes, s_lb, s_ub, s_time;
				std::string status;

				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if(s_lb.str() == "-inf") lb = -1;
				else s_lb >> lb;

				getline(input,line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if(s_ub.str() == "inf") ub = -1;
				else s_ub >> ub;

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if((status != "OPTIMAL") && (status != "INFEASIBLE"))
				{
					if((!(double_equals(lb,-1))) && (!(double_equals(ub,-1))))
					{
						if(double_less(time,time_limit))
						{
							curr_gap = 0.0;
						}else curr_gap = (100.0*(ub-lb))/ub;
					}else curr_gap = 100.0;
				}else
				{
					curr_gap = 0.0;
				}

				getline(input,line);
				getline(input,line);
				getline(input,line);

				for(int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l) getline(input,line);

				getline(input,line);
				getline(input,line);

				// read cuts here
				for(int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
				{
					std::stringstream s_cuts;
					getline(input,line);
					pos = line.find_first_of(":");
					size_t pos2 = line.find_first_of("/");
					s_cuts << line.substr(pos + 2, pos2 - pos -2);
					curr_cuts = 0.0;
					s_cuts >> curr_cuts;

					cuts[l].push_back(curr_cuts);
					total_cuts[l].push_back(curr_cuts);
					avg_cuts[l] += curr_cuts;
					total_avg_cuts[l] += curr_cuts;
					//std::cout << line << " " << curr_cuts << std::endl;
					//getchar();getchar();
				}

				getline(input,line);
				pos = line.find_first_of(":");
				s_nodes << line.substr(pos + 2);
				s_nodes >> curr_nodes;

				gaps.push_back(curr_gap);
				total_gaps.push_back(curr_gap);
				avg_gaps += curr_gap;
				total_avg_gaps += curr_gap;

				nodes.push_back(curr_nodes);
				total_nodes.push_back(curr_nodes);
				avg_nodes += curr_nodes;
				total_avg_nodes += curr_nodes;

				input.close();
			}
			//getchar();getchar();

			avg_gaps /= (gaps.size());
			avg_nodes /= (nodes.size());
			output << j+1;
			if(stop) output << "\\_5\\%";

			output << " & & " << avg_gaps << " & " << StDev(gaps,avg_gaps) << " & & " << avg_nodes << " &";

			for(int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
			{
				if(algorithms[k] == "bc_b1")
				{
					if(l == K_TYPE_GCC_CUT)
					{
						avg_cuts[l]/= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}

				if(algorithms[k] == "cb_csc3")
				{
					if((l == K_TYPE_GCC_CUT) || (l == K_TYPE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
					{
						avg_cuts[l]/= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}

				if(algorithms[k] == "cb2_csc3")
				{
					if((l == K_TYPE_CLIQUE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
					{
						avg_cuts[l]/= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}

			}

			output <<  "\\\\" << std::endl;
		}

		total_avg_gaps /= (total_gaps.size());
		total_avg_nodes /= (total_nodes.size());
		output << " & & " << total_avg_gaps << " & " << StDev(total_gaps,total_avg_gaps) << " & & " << total_avg_nodes << " &";

		for(int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
		{
			if(algorithms[k] == "bc_b1")
			{
				if(l == K_TYPE_GCC_CUT)
				{
					total_avg_cuts[l]/= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}

			if(algorithms[k] == "cb_csc3")
			{
				if((l == K_TYPE_GCC_CUT) || (l == K_TYPE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
				{
					total_avg_cuts[l]/= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}

			if(algorithms[k] == "cb2_csc3")
			{
				if((l == K_TYPE_CLIQUE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
				{
					total_avg_cuts[l]/= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}
		}

		output <<  "\\\\" << std::endl;
	}
	output.close();
}

void GenerateRootPrimalAndDualGapsLatexTable(std::vector<std::string> dirs, bool stop, double time_limit, bool user_cuts)
{
	//std::cout << "3" << std::endl; getchar(); getchar();
	std::vector<std::pair<std::string,std::string>> algorithms;

	//algorithms.push_back(std::pair<std::string,std::string>("bc_b1_initial_primal_bound","bc_b1"));
	algorithms.push_back(std::pair<std::string,std::string>("cb_csc3_initial_primal_bound","cb_csc3"));
	algorithms.push_back(std::pair<std::string,std::string>("cb2_csc3_initial_primal_bound","cb2_csc3"));

	std::vector<std::pair<std::string,std::string>> algorithms2;
	//algorithms2.push_back(std::pair<std::string,std::string>("bc_b1_initial_primal_bound","bc_b1"));
	//algorithms2.push_back(std::pair<std::string,std::string>("cb_csc3_initial_primal_bound","cb_csc3"));
	algorithms2.push_back(std::pair<std::string,std::string>("cb2_csc3_initial_primal_bound","cb2_csc3"));;

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_primal_dual_gaps";
	else output_name = ".//tables//latex//table_primal_dual_gaps";

    if(!user_cuts) output_name += "_no_user_cuts";

    output_name += ".txt";

	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_initial_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_initial_gap_per_algo_cplex_cuts(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_primal_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_dual_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_dual_gap_per_algo_no_cplex_cuts(algorithms.size(),std::vector<double>());

	std::vector<double> total_avg_initial_gap(algorithms.size(),0.0);
	std::vector<double> total_avg_initial_gap_cplex_cuts(algorithms.size(),0.0);
	std::vector<double> total_avg_gap(algorithms.size(),0.0);
	std::vector<double> total_avg_primal_gap(algorithms.size(),0.0);
	std::vector<double> total_avg_dual_gap(algorithms.size(),0.0);
	std::vector<double> total_avg_dual_gap_no_cplex_cuts(algorithms.size(),0.0);

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> initial_gap_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> initial_gap_per_algo_cplex_cuts(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> primal_gap_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> dual_gap_per_algo(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> dual_gap_per_algo_no_cplex_cuts(algorithms.size(),std::vector<double>());

		std::vector<double> avg_initial_gap(algorithms.size(),0.0);
		std::vector<double> avg_initial_gap_cplex_cuts(algorithms.size(),0.0);
		std::vector<double> avg_gap(algorithms.size(),0.0);
		std::vector<double> avg_primal_gap(algorithms.size(),0.0);
		std::vector<double> avg_dual_gap(algorithms.size(),0.0);
		std::vector<double> avg_dual_gap_no_cplex_cuts(algorithms.size(),0.0);

		std::vector<double> st_dev_initial_gap(algorithms.size(),0.0);
		std::vector<double> st_dev_initial_gap_cplex_cuts(algorithms.size(),0.0);
		std::vector<double> st_dev_gap(algorithms.size(),0.0);
		std::vector<double> st_dev_primal_gap(algorithms.size(),0.0);
		std::vector<double> st_dev_dual_gap(algorithms.size(),0.0);
		std::vector<double> st_dev_dual_gap_no_cplex_cuts(algorithms.size(),0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

		for(size_t i = 0; i < instances.size(); i++)
		{
			double root_lb = 0.0, root_ub2 = 0.0, root_ub = 0.0, best_lb = -1.0;
			double lb = 0.0, ub = 0.0;

			// find best primal bound
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append((algorithms[j]).second);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb1;
				std::string status;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				//std::cout << status << std::endl;
				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_lb1 << line.substr(pos + 2);

				if(s_lb1.str() == "-inf") lb = -1.0;
				else s_lb1 >> lb;

				if((!double_equals(lb,-1.0)) && (status != "INFEASIBLE") && (double_greater(lb,best_lb))) best_lb = lb;

				//std::cout << "best_lb: " << best_lb << std::endl;
				input.close();
			}

			for(size_t j = 0; j < algorithms2.size(); j++)
			{
				double curr_initial_gap_cplex_cuts = 0.0, curr_initial_gap = 0.0, curr_gap = 0.0, curr_primal_gap = 0.0, curr_dual_gap = 0.0, curr_dual_gap_no_cplex_cuts = 0.0, time = 0.0;
				ub = 0.0;
				root_ub = 0.0;
				root_ub2 = 0.0;
				lb = 0.0;
				root_lb = 0.0;

				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append((algorithms2[j]).first);

                if(!user_cuts) curr_file += "_no_user_cuts";
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}


				std::string status, status2;
				std::string line;
				std::stringstream s_root_lb, s_root_ub2;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				getline(input,line);

				pos = line.find_first_of(":");
				s_root_lb << line.substr(pos + 2);
				if(s_root_lb.str() == "-inf") root_lb = -1.0;
				else s_root_lb >> root_lb;

				getline(input,line);

				pos = line.find_first_of(":");
				s_root_ub2 << line.substr(pos + 2);
				if(s_root_ub2.str() == "inf") root_ub2 = -1.0;
				else s_root_ub2 >> root_ub2;

				input.close();

				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append((algorithms2[j]).second);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_root_ub, s_lb, s_ub, s_time;
				getline(input,line);
				pos = line.find_first_of(":");
				status2 = line.substr(pos + 2);

				std::string::iterator end_pos2 = std::remove(status2.begin(),status2.end(),' ');
				status2.erase(end_pos2, status2.end());

				getline(input,line);

				pos = line.find_first_of(":");
				s_root_ub << line.substr(pos + 2);
				if(s_root_ub.str() == "inf") root_ub = -1.0;
				else s_root_ub >> root_ub;

				getline(input,line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if(s_lb.str() == "-inf") lb = -1.0;
				else s_lb >> lb;

				getline(input,line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if(s_ub.str() == "inf") ub = -1.0;
				else s_ub >> ub;


				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if((status2 != "OPTIMAL") && (status2 != "INFEASIBLE"))
				{
					if((double_equals(ub,-1)) || (double_equals(lb,-1))) curr_gap = 100.0;
					else if(double_equals(ub,0.0) || double_less(time,time_limit)) curr_gap = 0.0;
					else curr_gap = (100*(ub - lb))/ub;
				}else
				{
					curr_gap = 0.0;
				}

				if(status != "INFEASIBLE")
				{
					if((double_equals(root_ub,-1)) || (double_equals(root_lb,-1))) curr_initial_gap = 100.0;
					else if(double_equals(root_ub,0.0)) curr_initial_gap = 0.0;
					else curr_initial_gap = (100*(root_ub - root_lb))/root_ub;

					if((double_equals(root_ub2,-1)) || (double_equals(root_lb,-1))) curr_initial_gap_cplex_cuts = 100.0;
					else if(double_equals(root_ub2,0.0)) curr_initial_gap_cplex_cuts = 0.0;
					else curr_initial_gap_cplex_cuts = (100*(root_ub2 - root_lb))/root_ub2;
				}else curr_initial_gap = curr_initial_gap_cplex_cuts = 0.0;

				if(!double_equals(best_lb,-1.0))
				{
					if((double_equals(root_lb,-1)) || (double_equals(root_ub2,-1))) curr_primal_gap = 100.0;
					else if(double_equals(root_ub2,0.0)) curr_primal_gap = 0.0;
					else curr_primal_gap = (100*(best_lb - root_lb))/root_ub2;

					if(double_equals(root_ub2,-1)) curr_dual_gap = 100.0;
					else if(double_equals(root_ub2,0.0)) curr_dual_gap = 0.0;
					else curr_dual_gap = (100*(root_ub2 - best_lb))/root_ub2;

					if(double_equals(root_ub,-1)) curr_dual_gap_no_cplex_cuts = 100.0;
					else if(double_equals(root_ub,0.0)) curr_dual_gap_no_cplex_cuts = 0.0;
					else curr_dual_gap_no_cplex_cuts = (100*(root_ub - best_lb))/root_ub;
				}else curr_primal_gap = curr_dual_gap = curr_dual_gap_no_cplex_cuts = 0.0;

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				initial_gap_per_algo[j].push_back(curr_initial_gap);
				initial_gap_per_algo_cplex_cuts[j].push_back(curr_initial_gap_cplex_cuts);
				gap_per_algo[j].push_back(curr_gap);
				primal_gap_per_algo[j].push_back(curr_primal_gap);
				dual_gap_per_algo[j].push_back(curr_dual_gap);
				dual_gap_per_algo_no_cplex_cuts[j].push_back(curr_dual_gap_no_cplex_cuts);

				total_initial_gap_per_algo[j].push_back(curr_initial_gap);
				total_initial_gap_per_algo_cplex_cuts[j].push_back(curr_initial_gap_cplex_cuts);
				total_gap_per_algo[j].push_back(curr_gap);
				total_primal_gap_per_algo[j].push_back(curr_primal_gap);
				total_dual_gap_per_algo[j].push_back(curr_dual_gap);
				total_dual_gap_per_algo_no_cplex_cuts[j].push_back(curr_dual_gap_no_cplex_cuts);

				total_avg_initial_gap[j] += curr_initial_gap;
				total_avg_initial_gap_cplex_cuts[j] += curr_initial_gap_cplex_cuts;
				total_avg_gap[j] += curr_gap;
				total_avg_primal_gap[j] += curr_primal_gap;
				total_avg_dual_gap[j] += curr_dual_gap;
				total_avg_dual_gap_no_cplex_cuts[j] += curr_dual_gap_no_cplex_cuts;

				avg_initial_gap[j] += curr_initial_gap;
				avg_initial_gap_cplex_cuts[j] += curr_initial_gap_cplex_cuts;
				avg_gap[j] += curr_gap;
				avg_primal_gap[j] += curr_primal_gap;
				avg_dual_gap[j] += curr_dual_gap;
				avg_dual_gap_no_cplex_cuts[j] += curr_dual_gap_no_cplex_cuts;

				input.close();
			}
			//getchar();getchar();
		}

		output << j+1;
		if(stop) output << "\\_5\\%";

		for(size_t k = 0; k < algorithms2.size(); k++)
		{
			avg_initial_gap[k]/=(1.0*(instances.size()));
			st_dev_initial_gap[k] = StDev(initial_gap_per_algo[k],avg_initial_gap[k]);
			avg_initial_gap_cplex_cuts[k]/=(1.0*(instances.size()));
			st_dev_initial_gap_cplex_cuts[k] = StDev(initial_gap_per_algo_cplex_cuts[k],avg_initial_gap_cplex_cuts[k]);
			avg_gap[k]/=(1.0*(instances.size()));
			st_dev_gap[k] = StDev(gap_per_algo[k],avg_gap[k]);
			avg_primal_gap[k]/=(1.0*(instances.size()));
			st_dev_primal_gap[k] = StDev(primal_gap_per_algo[k],avg_primal_gap[k]);
			avg_dual_gap[k]/=(1.0*(instances.size()));
			st_dev_dual_gap[k] = StDev(dual_gap_per_algo[k],avg_dual_gap[k]);
			avg_dual_gap_no_cplex_cuts[k]/=(1.0*(instances.size()));
			st_dev_dual_gap_no_cplex_cuts[k] = StDev(dual_gap_per_algo_no_cplex_cuts[k],avg_dual_gap_no_cplex_cuts[k]);
			output << " & & " << avg_initial_gap_cplex_cuts[k] << " & " << st_dev_initial_gap_cplex_cuts[k] << " & & " << avg_primal_gap[k] << " & " << st_dev_primal_gap[k] << " & & " << avg_dual_gap[k] << " & " << st_dev_dual_gap[k];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t j = 0; j < algorithms2.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_initial_gap[j]/=(1.0*((total_initial_gap_per_algo[j]).size()));
		total_avg_initial_gap_cplex_cuts[j]/=(1.0*((total_initial_gap_per_algo_cplex_cuts[j]).size()));
		total_avg_gap[j]/=(1.0*((total_gap_per_algo[j]).size()));
		total_avg_primal_gap[j]/=(1.0*((total_primal_gap_per_algo[j]).size()));
		total_avg_dual_gap[j]/=(1.0*((total_dual_gap_per_algo[j]).size()));
		total_avg_dual_gap_no_cplex_cuts[j]/=(1.0*((total_dual_gap_per_algo_no_cplex_cuts[j]).size()));
		output << "& & "<< total_avg_initial_gap_cplex_cuts[j] << " & " << StDev(total_initial_gap_per_algo_cplex_cuts[j],total_avg_initial_gap_cplex_cuts[j]) << "& & " << total_avg_primal_gap[j] << " & " << StDev(total_primal_gap_per_algo[j],total_avg_primal_gap[j]) << " & & " << total_avg_dual_gap[j] << " & " << StDev(total_dual_gap_per_algo[j],total_avg_dual_gap[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateLPImprovementsLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_FBs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_cb_csc3");
	algorithms.push_back("relax_FBs_CCs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_cb_csc3");
	algorithms.push_back("relax_FBs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_AVICs_cb_csc3");

	algorithms.push_back("relax_FBs_GCCs_CCs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_AVICs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_LCIs_AVICs_cb_csc3");*/
	//algorithms.push_back("relax_5_cb_csc3");
	/*algorithms.push_back("relax_6_cb_csc3");
	algorithms.push_back("relax_7_cb_csc3");
	algorithms.push_back("relax_8_cb_csc3");
	algorithms.push_back("relax_9_cb_csc3");*/


	//algorithms.push_back("cb_csc3");
	//algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");
	/*algorithms.push_back("relax_all_gccs");
	//algorithms.push_back("relax_all_gccs2");
	//algorithms.push_back("relax_all_gccs3");
	//algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_all_ccs");
	//algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_all_gccs_ccs");
	//algorithms.push_back("relax_gccs_ccs_cb_csc3");
	//algorithms.push_back("relax_all_gccs_ccs_lcis");
	//algorithms.push_back("relax_clique_ccs");
	//algorithms.push_back("relax_clique_ccs_active_vertices");
	algorithms.push_back("relax_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_clique_ccs_total_active_vertices_unilateral");
	algorithms.push_back("relax_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.1_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.05_clique_ccs_active_vertices_unilateral");*/
	//algorithms.push_back("relax_csc3e");
	//algorithms.push_back("relax_cb_csc3e"); // é o relax_clique_ccs_cb_csc3e
	//algorithms.push_back("relax_clique_ccs_lcis_cb_csc3e");
	//algorithms.push_back("relax_angle_0.03_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_lcis_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_clique_ccs");
	//algorithms.push_back("relax_all_cb_csc3");
	//algorithms.push_back("relax_all_maximal_cliques_ccs_cb_csc3");
	//algorithms.push_back("relax_all_clique_ccs_cb_csc3");
	//algorithms.push_back("relax_all_clique_ccs_cb_csc3_new");
	//algorithms.push_back("relax_clique_ccs_cb_csc3");
	//algorithms.push_back("relax_clique_ccs_lcis_cb_csc3");
	//algorithms.push_back("relax_clique_ccs_cb_csc3_new_unilateral");
	//algorithms.push_back("relax_clique_ccs_cb_csc3");
	/*algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_lcis_cb_csc3");
	algorithms.push_back("relax_gccs_ccs_cb_csc3");
	algorithms.push_back("relax_all_cb_csc3");*/
	//algorithms.push_back("relax_cccs_cb_csc3");
	//algorithms.push_back("relax_all_maximal_cliques_ccs_cb_csc3");
	//algorithms.push_back("relax_gccs_all_maximal_cliques_ccs_cb_csc3");

	//algorithms.push_back("relax_3gccs_cb_csc3");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//latex//stop_table_LP_improvements.txt";
	else output_name = ".//tables//latex//table_LP_improvements.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(),0.0);

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> improvement_per_algo(algorithms.size(),std::vector<double>());
		std::vector<double> avg_improvement(algorithms.size(),0.0);
		std::vector<double> st_dev(algorithms.size(),0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);
		for(size_t i = 0; i < instances.size(); i++)
		{
			std::cout << instances[i] << std::endl;
			double original_lp = 0.0;
			for(size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);

				if(s_lp.str() == "inf") lp = -1;
				else s_lp >> lp;

				if(j == 0)
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE"))
					{
						original_lp = -1.0;
					}else original_lp = lp;
					curr_improvement = 0.0;
				}else
				{
					if((double_equals(lp,-1)) || (status == "INFEASIBLE") || (double_equals(original_lp,-1)))
					{
						curr_improvement = -1;
					}else
					{
						if(double_equals(original_lp,0.0)) curr_improvement = 0.0;
						else curr_improvement = (100*(original_lp-lp))/original_lp;

						if(double_less(curr_improvement,0)) std::cout << original_lp << " - " << lp << std::endl;
					}
				}

				//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				if(!double_equals(curr_improvement,-1))
				{
					improvement_per_algo[j].push_back(curr_improvement);
					total_improvement_per_algo[j].push_back(curr_improvement);
					total_avg_improvement[j] += curr_improvement;
					avg_improvement[j] += curr_improvement;
				}
				input.close();
			}

			//getchar();getchar();
		}

		output << j+1;

		for(size_t j = 1; j < algorithms.size(); j++)
		{
			avg_improvement[j]/=(1.0*((improvement_per_algo[j]).size()));
			st_dev[j] = StDev(improvement_per_algo[j],avg_improvement[j]);
			output << " & & "<< avg_improvement[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t j = 1; j < algorithms.size(); j++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement[j]/=(1.0*((total_improvement_per_algo[j]).size()));
		output << "& & "<< total_avg_improvement[j] << " & " << StDev(total_improvement_per_algo[j],total_avg_improvement[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateCutsConfigurationsPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_FBs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_cb_csc3");
	algorithms.push_back("relax_FBs_CCs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_cb_csc3");
	algorithms.push_back("relax_FBs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_AVICs_cb_csc3");

	algorithms.push_back("relax_FBs_GCCs_CCs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_AVICs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_LCIs_cb_csc3");
	//algorithms.push_back("relax_FBs_CCCs_LCIs_AVICs_cb_csc3");

	std::fstream output;
	std::string output_name;
	if(stop) output_name = ".//tables//CSV//stop_table_LP_bounds.csv";
	else output_name = ".//tables//CSV//table_LP_bounds.csv";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

  output << "Configuration;";

  for(size_t i = 0;  i < algorithms.size(); ++i)
  {
    output << i << ";";
  }

  output << std::endl;

  output << "Instance;";

  for(size_t i = 0;  i < algorithms.size(); ++i)
    output << "upper bound;";

  output << std::endl;
	output << std::setprecision(2) << std::fixed;

	for(size_t j = 0; j < dirs.size(); j++)
	{

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

    std::sort(instances.begin(),instances.end());

		for(size_t i = 0; i < instances.size(); i++)
		{
			std::cout << instances[i] << std::endl;


      if(stop) output << instances[i].substr(0,instances[i].size() - 4) << "_5%;";
			else output << instances[i].substr(0,instances[i].size() - 4) << ";";

			for(size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(),std::fstream::in);

				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				std::string line;

				getline(input,line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);

        if(status == "INFEASIBLE") output << " - ;";
        else output << s_lp.str() << ";";

        //output << s_lp.str() << ";";
				input.close();
			}

			//getchar();getchar();
      output << std::endl;
		}

	}

	output.close();
}

void GenerateExactAlgorithmsPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop, double time_limit)
{
	std::string curr_file;
	double time = 0.0;
	double ub = 0.0, lb = 0.0;

    std::vector<std::string> algorithms;
    algorithms.push_back("bc_b1");
    algorithms.push_back("cb_csc3");
    /*algorithms.push_back("cb2_csc3");*/

		/*std::string algo = "cb2_csc3_";
		if(stop) algo += "stop_";
		algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
		algorithms.push_back(algo);

		algo = "cb2_csc3_";
		if(stop) algo += "stop_";
		algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
		algorithms.push_back(algo);*/

	std::fstream output;
	std::string output_name;

	//std::cout << output_name << std::endl;

	if(stop) output_name = ".//tables//CSV//stop_exact_algorithms_per_instance.csv";
	else output_name = ".//tables//CSV//exact_algorithms_per_instance.csv";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	output << std::setprecision(2) << std::fixed;

	output << "Algorithm;";
	for(size_t i = 0; i < algorithms.size(); ++i)
	{
		for(size_t j = 0; j < 3; ++j)
		{
			output << algorithms[i] << ";";
		}
	}

	output << std::endl;

	output << "Instance;";
	for(size_t i = 0; i < algorithms.size(); ++i)
	{
			output << "lb;ub;time(s);";
	}

	output << std::endl;

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

		std::sort(instances.begin(),instances.end());

		//output << folder << std::endl;

		for(size_t i = 0; i < instances.size(); i++)
		{
			if(stop) output << instances[i].substr(0,instances[i].size() - 4) << "_5%;";
			else output << instances[i].substr(0,instances[i].size() - 4) << ";";

            for(size_t k = 0; k < algorithms.size(); ++k)
            {
						    curr_file = ".//solutions//";
						    curr_file.append(folder);
						    curr_file.append("s_");
						    if(stop) curr_file.append("stop_");
			                curr_file.append(algorithms[k]);
						    curr_file.append("_");
						    curr_file.append(instances[i]);
						    //std::cout << curr_file << std::endl;

						    std::fstream input;
						    input.open(curr_file.c_str(),std::fstream::in);

						    if(!input.is_open())
						    {
							    std::cout << "Could not open file " << curr_file << std::endl;
							    throw 4;
						    }

						    std::stringstream s_lb, s_ub, s_time;
						    std::string status;
						    std::string line;

						    getline(input,line);
						    size_t pos = line.find_first_of(":");
						    status = line.substr(pos + 2);

						    std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
						    status.erase(end_pos, status.end());

						    getline(input,line);
						    getline(input,line);
						    pos = line.find_first_of(":");
						    s_lb << line.substr(pos + 2);
						    getline(input,line);
						    pos = line.find_first_of(":");
						    s_ub << line.substr(pos + 2);

						    getline(input,line);
						    getline(input,line);
						    pos = line.find_first_of(":");
						    s_time << line.substr(pos + 2);

						    s_time >> time;
						    s_lb >> lb;
						    s_ub >> ub;
						    //gap1 = 0.0;

						    //if((status != "OPTIMAL") && (status != "INFEASIBLE"))
						    //{
							    //if(!double_less(time,time_limit)) gap1 = (100.0*(ub-lb))/ub;
						    //}

			                if( ( (status == "OPTIMAL") || (double_less(time,7200))) && (!double_equals(lb,ub)) ) std::cout << algorithms[k] << "_" << instances[i] << " " << status << " " << time << "s " << lb << " != " << ub << std::endl;

						    if(double_greater(time,7200)) time = 7200;
						    if(status == "INFEASIBLE") output << " - ; - ; " << time << ";";
						    else output << s_lb.str() << " ; " << s_ub.str() << " ; " << time << ";";
						    input.close();
            }
            output << std::endl;
		}
	}
	output.close();
}

void GenerateAppendixLatexTable(std::vector<std::string> dirs, bool stop, double time_limit)
{
	std::string curr_file;
	double time = 0.0;
	double ub = 0.0, lb = 0.0;

    std::vector<std::string> algorithms;
    /*algorithms.push_back("bc_b1");
    algorithms.push_back("cb_csc3");
    algorithms.push_back("cb2_csc3");*/

		std::string algo = "cb2_csc3_";
		if(stop) algo += "stop_";
		algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
		algorithms.push_back(algo);

		algo = "cb2_csc3_";
		if(stop) algo += "stop_";
		algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
		algorithms.push_back(algo);

	std::fstream output;
	std::string output_name;

	//std::cout << output_name << std::endl;

	if(stop) output_name = ".//tables//latex//stop_appendix.txt";
	else output_name = ".//tables//latex//appendix.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	output << std::setprecision(2) << std::fixed;

	for(size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j],instances,false);

		std::sort(instances.begin(),instances.end());

		//output << folder << std::endl;

		for(size_t i = 0; i < instances.size(); i++)
		{
			if(stop) output << instances[i].substr(0,instances[i].size() - 4) << "\\_5\\%";
			else output << instances[i].substr(0,instances[i].size() - 4);

            for(size_t k = 0; k < algorithms.size(); ++k)
            {
			    curr_file = ".//solutions//";
			    curr_file.append(folder);
			    curr_file.append("s_");
			    if(stop) curr_file.append("stop_");
                curr_file.append(algorithms[k]);
			    curr_file.append("_");
			    curr_file.append(instances[i]);
			    //std::cout << curr_file << std::endl;

			    std::fstream input;
			    input.open(curr_file.c_str(),std::fstream::in);

			    if(!input.is_open())
			    {
				    std::cout << "Could not open file " << curr_file << std::endl;
				    throw 4;
			    }

			    std::stringstream s_lb, s_ub, s_time;
			    std::string status;
			    std::string line;

			    getline(input,line);
			    size_t pos = line.find_first_of(":");
			    status = line.substr(pos + 2);

			    std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
			    status.erase(end_pos, status.end());

			    getline(input,line);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_lb << line.substr(pos + 2);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_ub << line.substr(pos + 2);

			    getline(input,line);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_time << line.substr(pos + 2);

			    s_time >> time;
			    s_lb >> lb;
			    s_ub >> ub;
			    //gap1 = 0.0;

			    //if((status != "OPTIMAL") && (status != "INFEASIBLE"))
			    //{
				    //if(!double_less(time,time_limit)) gap1 = (100.0*(ub-lb))/ub;
			    //}

                if( ( (status == "OPTIMAL") || (double_less(time,7200))) && (!double_equals(lb,ub)) ) std::cout << algorithms[k] << "_" << instances[i] << " " << status << " " << time << "s " << lb << " != " << ub << std::endl;


                if(k != 0) output << " & ";
			    if(double_greater(time,7200)) time = 7200;
			    if(status == "INFEASIBLE") output << " & -- & -- & " << time;
			    else output << " & " << s_lb.str() << " & " << s_ub.str() << " & " << time;
			    input.close();

			    /*curr_file = ".//solutions//";
			    curr_file.append(folder);
			    curr_file.append("s_");
			    if(stop) curr_file.append("stop_");
			    curr_file.append("relax_none_cb_csc3_");
			    curr_file.append(instances[i]);
			    //std::cout << curr_file << std::endl;

			    input.open(curr_file.c_str(),std::fstream::in);
			    if(!input.is_open())
			    {
				    std::cout << "Could not open file " << curr_file << std::endl;
				    throw 4;
			    }

			    std::stringstream s_lp, s_lp_cuts;

			    getline(input,line);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_lp << line.substr(pos + 2);

			    input.close();

			    curr_file = ".//solutions//";
			    curr_file.append(folder);
			    curr_file.append("s_");
			    if(stop) curr_file.append("stop_");
			    curr_file.append("cb_csc3_");
			    curr_file.append(instances[i]);
			    //std::cout << curr_file << std::endl;

			    input.open(curr_file.c_str(),std::fstream::in);
			    if(!input.is_open())
			    {
				    std::cout << "Could not open file " << curr_file << std::endl;
				    throw 4;
			    }

			    std::stringstream s_lb2, s_ub2, s_time2;

			    getline(input,line);
			    pos = line.find_first_of(":");
			    status = line.substr(pos + 2);

			    end_pos = std::remove(status.begin(),status.end(),' ');
			    status.erase(end_pos, status.end());

			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_lp_cuts << line.substr(pos + 2);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_lb2 << line.substr(pos + 2);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_ub2 << line.substr(pos + 2);

			    getline(input,line);
			    getline(input,line);
			    pos = line.find_first_of(":");
			    s_time2 << line.substr(pos + 2);

			    s_time2 >> time;

			    s_lb2 >> lb;
			    s_ub2 >> ub;
			    gap2 = 0.0;

			    if((status != "OPTIMAL") && (status != "INFEASIBLE"))
			    {
				    if(!double_less(time,time_limit)) gap2 = (100.0*(ub-lb))/ub;
			    }

			    if(double_greater(time,7200)) time = 7200;
			    if(status == "INFEASIBLE") output << "-- & -- & " << time << " & -- & --";
			    else output << s_lb2.str() << " & " << s_ub2.str() << " & " << time << " & " << s_lp.str() << " & " << s_lp_cuts.str();*/
			    //for(int i = 0; i < 10; i++) getline(input,line);
            }
            output << "\\\\" << std::endl;
		}

		if(j != dirs.size() - 1) output << "\\\\" << std::endl;
	}
	output.close();
}

void GenerateCSVFile(std::string dir, double time_limit, bool stop)
{
	std::string folder = dir.substr(20);

	std::vector<std::string> instances;
	AddInstancesFromDirectory(dir,instances,false);
	std::sort(instances.begin(),instances.end());
	std::string curr_file;

	std::vector<std::string> algorithms;
	/*algorithms.push_back("bcsc1");
	algorithms.push_back("bcsc2");
	algorithms.push_back("bcsc3");*/
	//algorithms.push_back("bcsc4");

	/*algorithms.push_back("bcmc1");
	algorithms.push_back("bcmc2");
	algorithms.push_back("bcmc3");*/
	//algorithms.push_back("bcmc4");

	//algorithms.push_back("bctc");
	//algorithms.push_back("bcmtc");

	//algorithms.push_back("csc1");
	//algorithms.push_back("csc3");
	//algorithms.push_back("csc4");

	//algorithms.push_back("cmc1");
	//algorithms.push_back("cmc2");

	//algorithms.push_back("ctc1");
	//algorithms.push_back("ctc2");

	//algorithms.push_back("cmtc1");
	//algorithms.push_back("cmtc2");

	//algorithms.push_back("kb");

	//algorithms.push_back("b1");
	//algorithms.push_back("b2");

	/*algorithms.push_back("bcsc4");
	algorithms.push_back("bcmc4");
	algorithms.push_back("bctc");
	algorithms.push_back("bcmtc");*/
	//algorithms.push_back("csc3_relax1");

	/*algorithms.push_back("csc3_raw");
	algorithms.push_back("csc3");
	algorithms.push_back("csc3_avics");
	algorithms.push_back("avics_csc3");*/

	//algorithms.push_back("bc_csc3_heuristic_all_cuts");
	//algorithms.push_back("csc3");
	//algorithms.push_back("csc3_heuristic_advance_2");
	//algorithms.push_back("cb_csc3_heuristic_avics_hard_coded_advance_2");
	//algorithms.push_back("cb_csc3_heuristic_avics_cuts_advance_2");

	//algorithms.push_back("cb_csc3");
	//algorithms.push_back("cb_csc3_heuristic");
	//algorithms.push_back("cb_csc3_heuristic_advance_2");
	//algorithms.push_back("cb_csc3_new_heuristic_advance_2");

	/*algorithms.push_back("cb_avics_gccs_ccs_lcis_csc3");
	algorithms.push_back("cb_avics_clique_ccs_lcis_csc3");
	algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");*/
	//algorithms.push_back("cb_avics_clique_ccs_lcis_csc3");
	//algorithms.push_back("cb_heuristic");
	//algorithms.push_back("cb_heuristic_advance_2");
	//algorithms.push_back("cb_csc3_heuristic_avics_cover_advance_2");
	algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	//algorithms.push_back("cb2_csc3");
	//algorithms.push_back("bc_b1");

	std::string algo = "cb_csc3_";
	/*if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	/*algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);*/

	algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	//algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");
	/*algorithms.push_back("bc_branchcallback_csc3");
	algorithms.push_back("bc_cutcallback_branchcallback_csc3");
	algorithms.push_back("bc_cutcallback_csc3");
	algorithms.push_back("bc_initialflowbounds_cutcallback_branchcallback_csc3");*/
	//algorithms.push_back("st_bc_b1");
	//algorithms.push_back("st_cb_csc3");
	//algorithms.push_back("csc3_relax3");
	//algorithms.push_back("cb_r3");
	//algorithms.push_back("cb_r4");
	/*algorithms.push_back("bc1_0.05_PREPROCESS_WITH_Y_b1");
	algorithms.push_back("bc1_0.05_PREPROCESS_WITH_Y_csc3");
	algorithms.push_back("bc4_0.6_0.1_PREPROCESS_WITH_Y_b1");
	algorithms.push_back("bc4_0.6_0.1_PREPROCESS_WITH_Y_csc3");*/
	//algorithms.push_back("csc4");
	/*algorithms.push_back("cmc2");
	algorithms.push_back("ctc1");
	algorithms.push_back("cmtc1");
	algorithms.push_back("kb");*/
	//algorithms.push_back("b2");

	std::fstream output;
	std::string output_name = ".//tables//CSV//";
	output_name.append(folder);
	if(stop) output_name.append("correct_table_stop.csv");
	else output_name.append("correct_table_top.csv");
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;
	output << "Formulation;;;;";

	for(size_t j = 0; j < algorithms.size(); j++)
	{
		if(stop) algorithms[j] = "stop_" + algorithms[j];
		output << algorithms[j] <<";" << algorithms[j] <<";" << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";;";
	}

	output << std::endl << "Instance;T;m;;";
	for(size_t j = 0; j < algorithms.size(); j++)
	{
		output << "status;lp;lb;ub;gap;time(s);;";
	}

	output << std::endl;

	for(size_t i = 0; i < instances.size(); i++)
	{
		curr_file = dir + instances[i];
		double limit = 0.0;
		int num_vehicles = 0;

		std::fstream file;
		file.open(curr_file.c_str(),std::fstream::in);

		if(file.is_open())
		{
			std::stringstream parameter;
			std::string line;
			getline(file,line);

			getline(file,line);

			getline(file,line);
			std::size_t pos_1 = line.find_first_of(" ");

			parameter.clear();
			parameter << line.substr(pos_1 + 1);
			parameter >> num_vehicles;

			getline(file,line);
			pos_1 = line.find_first_of(" ");

			parameter.clear();
			parameter << line.substr(pos_1 + 1);
			parameter >> limit;

			file.close();
		}

		output << instances[i] << ";" << limit << ";" << num_vehicles << ";;";
		for(size_t j = 0; j < algorithms.size(); j++)
		{
			curr_file = ".//solutions//";
			curr_file.append(folder);
			curr_file.append("s_");
			curr_file.append(algorithms[j]);
			curr_file.append("_");
			curr_file.append(instances[i]);
			//std::cout << curr_file << std::endl;

			std::fstream input;
			std::cout << curr_file << std::endl;
			input.open(curr_file.c_str(),std::fstream::in);

			if(!input.is_open())
			{
				output << ";" << ";" << ";" << ";" << ";" << ";;";
				continue;
			}

			std::stringstream s_lp, s_lb, s_ub, s_time;
			std::string status;
			double lp, lb, ub, time, gap;
			std::string line;

			getline(input,line);
			size_t pos = line.find_first_of(":");
			status = line.substr(pos + 2);

			getline(input,line);
			pos = line.find_first_of(":");
			s_lp << line.substr(pos + 2);

			getline(input,line);
			pos = line.find_first_of(":");
			s_lb << line.substr(pos + 2);

			getline(input,line);
			pos = line.find_first_of(":");
			s_ub << line.substr(pos + 2);

			getline(input,line);
			getline(input,line);
			pos = line.find_first_of(":");
			s_time << line.substr(pos + 2);

			s_lp >> lp;
			s_lb >> lb;
			s_ub >> ub;
			s_time >> time;

			std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
			status.erase(end_pos, status.end());

			if(double_equals(lp,-1)) output << status << ";" << "-" << ";";
			else output << status << ";" << lp << ";";
			if(double_less(time,time_limit)) status = "OPTIMAL";

			if(status != "INFEASIBLE")
			{
				gap = (100.0*(ub - lb))/ub;

				output << lb << ";" << ub << ";" <<  gap << ";" << time << ";;";
			}else output << "-" << ";" << "-" << ";" <<  "-" << ";" << time << ";;";

		}
		output << std::endl;
	}

	output.close();
}

double RetrieveTimeSpentInFP(Instance& inst, std::string algo, std::string folder, std::string file_name)
{
		double t1 = 0.0, t2 = 0.0;
		std::stringstream s_t1, s_t2;

    std::fstream input;
    std::string path = ".//solutions";
    path.append(folder);
    //struct stat sb;
    //if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
    path.append("s_");
    path.append(algo);
    path.append("_");
    path.append(file_name);
	//std::cout << path << std::endl;

    input.open(path.c_str(),std::fstream::in);

		if(!(input.is_open())) throw "Error opening file 2";
		std::string line;
		for(int i = 1; i <= 7; ++i) getline(input,line);

		size_t pos = line.find_first_of(":");
		s_t1 << line.substr(pos + 2);
		s_t1 >> t1;

		for(int i = 1; i <= 5; ++i) getline(input,line);

		pos = line.find_first_of(":");
		s_t2 << line.substr(pos + 2);
		s_t2 >> t2;

		input.close();

		//std::cout << path << std::endl;
		//std::cout << t1 << std::endl;
		//std::cout << t2 << std::endl;
		return t1+t2;
}

double RetrieveTimeSpentInLNS(Instance& inst, std::string algo, std::string folder, std::string file_name)
{
		double time = 0.0;
		std::stringstream s_time;

    std::fstream input;
    std::string path = ".//solutions";
    path.append(folder);
    //struct stat sb;
    //if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
    path.append("s_");
    path.append(algo);
    path.append("_");
    path.append(file_name);
	//std::cout << path << std::endl;

    input.open(path.c_str(),std::fstream::in);

		if(!(input.is_open())) throw "Error opening file 2";
		std::string line;
		for(int i = 1; i <= 5; ++i) getline(input,line);

		size_t pos = line.find_first_of(":");
		s_time << line.substr(pos + 2);
		s_time >> time;
		input.close();

		//std::cout << path << std::endl;
		//std::cout << time << std::endl;
		return time;
}

void AddInstances(std::vector<std::string> & instances)
{
	//AddInstancesFromDirectory(".//instances//STOP//unsolved//",instances);
	AddInstancesFromDirectory(".//instances//STOP//current//",instances);
	/*AddInstancesFromDirectory(".//instances//STOP//Set_21_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_32_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_33_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_64_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_66_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_100_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_102_234//Set_0.05//",instances);*/
}

std::string GenerateAlgorithmName()
{
	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
	std::string algo;

	if((*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT]) algo += "FBs_";
	if((*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT]) algo += "GCCs_";
	if((*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT]) algo += "CCs_";
	if((*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT]) algo += "CCCs_";
	if((*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT]) algo += "LCIs_";
	if((*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT]) algo += "AVICs_";

	return algo;
}

static const struct option longOpts[] = {
	{ "compact", no_argument, NULL, 'a' },
	{ "cutting-plane", no_argument, NULL, 'b' },
	{ "branch-and-cut", no_argument, NULL, 'c' },
	{ "bianchessi", no_argument, NULL, 'd' },
	{ "capacity-based", no_argument, NULL, 'e' },
	{ "time-limit", required_argument, NULL, 'f' },
	{ "instance", required_argument, NULL, 'g' },
	{ "generate-convex-hull", no_argument, NULL, 'h' },
	{ "FBs", no_argument, NULL, 'i' },
	{ "GCCs", no_argument, NULL, 'j' },
	{ "CCs", no_argument, NULL, 'k' },
	{ "CCCs", no_argument, NULL, 'l' },
	{ "LCIs", no_argument, NULL, 'm' },
	{ "AVICs", no_argument, NULL, 'n' },
	{ NULL, no_argument, NULL, 0 }
};

void ParseArgumentsAndRun(int argc, char* argv[] )
{
	std::string instance, folder, file_name;
	int c;
	bool solve_compact = false, solve_cb = false, solve_bc = false, bianchessi = false, capacity_based = false, generate_convex_hull = false;
	double time_limit = 0.0;

	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

	while( ( c = getopt_long (argc, argv, "012345:6:7", longOpts, NULL) ) != -1 )
	{
		switch(c)
		{
			case 'a':
			solve_compact = true;
			break;
			case 'b':
			solve_cb = true;
			break;
			case 'c':
			solve_bc = true;
			break;
			case 'd':
			bianchessi = true;
			break;
			case 'e':
			capacity_based = true;
			break;
			case 'f':
			if(optarg) time_limit = std::atoi(optarg);
			break;
			case 'g':
			if(optarg) instance = std::string(optarg);
			break;
			case 'h':
			generate_convex_hull = true;
			break;
			case 'i':
			(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
			break;
			case 'j':
			(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
			break;
			case 'k':
			(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
			break;
			case 'l':
			(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
			break;
			case 'm':
			(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
			break;
			case 'n':
			(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;
			break;
		}
	}

	if((solve_compact && solve_bc) || (solve_compact && solve_cb) || (solve_bc && solve_cb) || (generate_convex_hull && (solve_compact || solve_bc || solve_cb))) throw 2;
	if(bianchessi && capacity_based) throw 3;

	split_file_path(instance,folder,file_name);
	std::cout << "* " << file_name << std::endl;

	Instance inst(instance);
	Graph* graph = inst.graph();

	Solution<int> * sol = new Solution<int>(graph->num_vertices());

	double * r = Dijkstra(graph,false,true);
	double * R = Dijkstra(graph,false,false);
	double * Rn = Dijkstra(graph,true,false);

	if(solve_compact)
	{
		if(bianchessi)
		{
			Bianchessi(inst,R,Rn,time_limit,sol,true,false,false);
			//Bianchessi(inst,R,Rn,time_limit,sol,false,false,false);
			std::string algo;
			if(K_STOP) algo += "stop_";
			algo += "relax_b1";
			std::cout << algo << std::endl;
			sol->write_to_file(algo,folder,file_name);
		}

		if(capacity_based)
		{
			//std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

			//(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
			//(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;

			//CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,false);
			CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false,NULL,NULL);
			sol->write_to_file("csc3_initial_primal_bound",folder,file_name);
		}
	}

	if(solve_cb)
	{
		//std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
		std::list<UserCutGeneral*> * cuts = NULL;

		(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
		(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
		(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
		(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
		//(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
		//(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;

		if(bianchessi)
		{
			cuts = Bianchessi(inst,R,Rn,time_limit,sol,true,false,true);

			if(time_limit != -1) time_limit = std::max(0.0, time_limit - sol->root_time_);
			sol->milp_time_ = sol->root_time_;

			Bianchessi(inst,R,Rn,time_limit,sol,false,false,false,cuts);

			sol->write_to_file("stop_cb_b1",folder,file_name);
		}

		if(capacity_based)
		{
			HeuristicSolution * initial_sol = NULL;

			if(K_PREPROCESS_REDUCED_COSTS)
			{
				int seed = 1;
				std::string algo = ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed);
				initial_sol = new HeuristicSolution();
				initial_sol->ReadFromFile(inst,algo,folder,file_name);

				PreprocessInstance(inst,R,Rn,time_limit,sol,1.0*(initial_sol->profits_sum_));

				delete initial_sol;
				initial_sol = NULL;
			}

			if(K_ADD_INITIAL_HEURISTIC_SOLUTION)
			{
				srand(time(NULL));
				int seed = rand()%10 + 1;
				std::string algo1 = FPHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed);
				std::string algo2 = ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed);

				double time_spent = RetrieveTimeSpentInFP(inst,algo1,folder,file_name);
				time_spent += RetrieveTimeSpentInLNS(inst,algo2,folder,file_name);

				initial_sol = new HeuristicSolution();
				initial_sol->ReadFromFile(inst,algo2,folder,file_name);

				if(!double_equals(time_limit,-1)) time_limit = std::max(0.0, time_limit - time_spent);
				(sol->milp_time_) += time_spent;
				//std::cout << time_limit << std::endl;
			}

			cuts = CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,false,true,NULL,NULL);

			if(!double_equals(time_limit,-1)) time_limit = std::max(0.0, time_limit - sol->root_time_);
			(sol->milp_time_) += (sol->root_time_);

			CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,false,false,cuts,initial_sol);

			std::string algo;
			if(K_STOP) algo = "stop_";
			algo += "cb_csc3";
			if(K_ADD_INITIAL_HEURISTIC_SOLUTION) algo += ("_" + ALNSHeuristicSolution::GenerateFileName());
			std::cout << algo << std::endl;
			//getchar(); getchar();
			sol->write_to_file(algo,folder,file_name);
			//sol->write_to_file("cb_angle_0.03_clique_ccs_lcis_csc3",folder,file_name);
			if(initial_sol != NULL)
			{
				delete initial_sol;
				initial_sol = NULL;
			}
		}

		DeleteCuts(cuts);
		if(cuts != NULL)
		{
			delete cuts;
			cuts = NULL;
		}
	}

	if(solve_bc)
	{
		//std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

		//(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
		//(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
		//(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
		//(*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
		//			(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
		//			(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;

		if(bianchessi)
		{
			//Bianchessi(inst,R,Rn,time_limit,sol,true,true,false);
			Bianchessi(inst,R,Rn,time_limit,sol,false,true,false);

			std::string algo;
			if(K_STOP) algo = "stop_";
			algo += "bc_b1_initial_primal_bound_no_user_cuts";
			std::cout << algo << std::endl;

			sol->write_to_file(algo,folder,file_name);
		}

		if(capacity_based)
		{
			HeuristicSolution * initial_sol = NULL;
			if(K_ADD_INITIAL_HEURISTIC_SOLUTION)
			{
				int seed = 1;
				std::string algo = ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed);
				initial_sol = new HeuristicSolution();
				initial_sol->ReadFromFile(inst,algo,folder,file_name);
			}
			//CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,true,true,false);
			CapacitatedSingleCommodity(inst,R,Rn,time_limit,sol,false,true,false,NULL,initial_sol);
			sol->write_to_file("bc_csc3_heuristic_all_cuts",folder,file_name);
			std::cout << "bc_csc3_heuristic_all_cuts" << std::endl;

			if(initial_sol != NULL)
			{
				delete initial_sol;
				initial_sol = NULL;
			}
		}
	}

	if(generate_convex_hull)
	{
		if(bianchessi)
		{
			std::string porta_file = "example.ieq";
			WritePORTAIneqFileBianchessi(inst,R,Rn,"./PORTA/src/"+porta_file);
			std::string command = "cd PORTA/src/; ./vint "+porta_file;
			int res = std::system(command.c_str());

			porta_file = "example.poi";
			command = "cd PORTA/src/; ./traf "+porta_file;
			res = std::system(command.c_str());

			ParsePORTAFile(inst,".//PORTA//porta_b1_"+file_name);
		}
		if(capacity_based)
		{
			std::string porta_file = "example.ieq";
			WritePORTAIneqFile(inst,R,Rn,"./PORTA/src/"+porta_file);
			std::string command = "cd PORTA/src/; ./vint "+porta_file;
			int res = std::system(command.c_str());

			porta_file = "example.poi";
			command = "cd PORTA/src/; ./traf "+porta_file;
			res = std::system(command.c_str());

			ParsePORTAFile(inst,".//PORTA//porta_csc3_"+file_name);
		}
	}

	delete [] R;
	R = NULL;

	delete [] Rn;
	Rn = NULL;

	delete [] r;
	r = NULL;

	delete sol;
	sol = NULL;
}

int main(int argc, char* argv[])
{
	/*std::list<std::list<int>> conflicts_list;
	std::vector<std::list<int>> conflict_graph(3,std::list<int>());
	conflict_graph[0].push_back(1);
	conflict_graph[1].push_back(0);

	FindAllMaximalConflictCliques2(&conflict_graph,conflicts_list);
	return 0;*/

	double time_limit = 7200;
	try
	{

		if(K_GETOPT)
		{
			ParseArgumentsAndRun(argc, argv);
            DeleteTimer();
            DeleteCallbackSelection();
			return 0;
		}
		std::vector<std::string> instances;

		int option = -1;
		double mandatory_percentage = 0.05;

		std::cout << "#    INSTANCES SELECTED HARD CODED." << std::endl;
		do{
			std::cout << "*******************************************" << std::endl
			<< " 0 - Generate STOP instances from TOP/OP ones" << std::endl
			<< " 1 - Generate result tables" << std::endl
			<< "*******************************************" << std::endl
			<< " 2 - Perform ALL EXPERIMENTS" << std::endl
			<< " 3 - Compact Formulations" << std::endl
			<< " 4 - Branch-and-Cut via Callback" << std::endl
			<< " 5 - Cut-and-Branch (only add cuts at root node)" << std::endl
			<< "Option: ";
			std::cin >> option;
			switch(option){
				case 0:{
					break;
					//GenerateFilesFromTOPinstances(".//instances//TOP//current//",mandatory_percentage);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_21_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_32_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_33_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_64_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_66_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_100_234//",mandatory_percentage,true);
					GenerateFilesFromTOPinstances(".//instances//TOP//Set_102_234//",mandatory_percentage,true);

					//GenerateFilesFromTOPinstances(".//instances//TOP//other_new_instances//",mandatory_percentage);

					/*GenerateFilesFromOPinstances(".//instances//OP//set_64_1//",mandatory_percentage);
					GenerateFilesFromOPinstances(".//instances//OP//set_66_1//",mandatory_percentage);
					GenerateFilesFromOPinstances(".//instances//OP//Tsiligirides 1//",mandatory_percentage);
					GenerateFilesFromOPinstances(".//instances//OP//Tsiligirides 2//",mandatory_percentage);
					GenerateFilesFromOPinstances(".//instances//OP//Tsiligirides 3//",mandatory_percentage);*/

					break;
				}
				case 1:{
					//GenerateCSVFile(".//instances//STOP//unsolved//",time_limit,false);
					/*GenerateCSVFile(".//instances//STOP//Set_21_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_32_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_33_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_64_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_66_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_100_234//Set_0.05//",time_limit,true);
					GenerateCSVFile(".//instances//STOP//Set_102_234//Set_0.05//",time_limit,true);

					GenerateCSVFile(".//instances//STOP//Set_21_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_32_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_33_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_64_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_66_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_100_234//Set_0.05//",time_limit,false);
					GenerateCSVFile(".//instances//STOP//Set_102_234//Set_0.05//",time_limit,false);*/

					/*GenerateCSVFile(".//instances//STOP//Set_21_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_32_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_33_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_64_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_66_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_100_234//Set_0.1//",time_limit);
					GenerateCSVFile(".//instances//STOP//Set_102_234//Set_0.1//",time_limit);*/

					std::vector<std::string> dirs;
					//dirs.push_back(".//instances//STOP//current//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_32_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_21_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_33_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_100_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_66_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_64_234//Set_0.1//");
					//dirs.push_back(".//instances//STOP//Set_102_234//Set_0.1//");
					//GenerateLPImprovementsLatexTable(dirs,false);
					//GenerateAlgorithmsLatexTable(dirs,false,time_limit);
					//GenerateAppendixLatexTable(dirs,false, time_limit);

					dirs.clear();
					dirs.push_back(".//instances//STOP//Set_32_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_21_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_33_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_100_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_66_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_64_234//Set_0.05//");
					dirs.push_back(".//instances//STOP//Set_102_234//Set_0.05//");
					//GenerateRootPrimalAndDualGapsLatexTable(dirs,false,time_limit,true);
					//GenerateRootPrimalAndDualGapsLatexTable(dirs,true,time_limit,true);
					//GenerateGapsNodesCutsLatexTable(dirs,false,time_limit);
					//GenerateLPImprovementsLatexTable(dirs,false);
					//GenerateLPImprovementsLatexTable(dirs,true);
					//GenerateRelaxCSVTable(dirs,false,time_limit);
					//GenerateValidIneqsImprovementsLatexTable(dirs,false);
					//GenerateLPImprovementsLatexTable(dirs,true);
					//GenerateValidIneqsImprovementsLatexTable(dirs,true);
					//GenerateAlgorithmsLatexTable(dirs,false,time_limit);
					//GenerateAlgorithmsLatexTable(dirs,true,time_limit);
					//GenerateAlgorithmsLatexTable(dirs,true,time_limit);
					//GenerateAppendixLatexTable(dirs,false,time_limit);
					//GenerateAppendixLatexTable(dirs,true,time_limit);

					//GenerateExactAlgorithmsPerInstanceResultsCSV(dirs,true,time_limit);
					//GenerateExactAlgorithmsPerInstanceResultsCSV(dirs,false,time_limit);

					GenerateCutsConfigurationsPerInstanceResultsCSV(dirs,true);
					GenerateCutsConfigurationsPerInstanceResultsCSV(dirs,false);

					//dirs.clear();
					//dirs.push_back(".//instances//STOP//current//");
					//GenerateValidIneqsActiveLatexTable(dirs,false);
					//GenerateValidIneqsActiveLatexTable(dirs,true);
					break;
				}
				case 2:{
					AddInstances(instances);
					SolveCompactFormulations(instances);
					SolveBranchAndCutCallback(instances);
					break;
				}
				case 3:{
					AddInstances(instances);
					SolveCompactFormulations(instances);
					break;
				}
				case 4:{
					AddInstances(instances);
					SolveBranchAndCutCallback(instances);
					break;
				}
				case 5:{
					AddInstances(instances);
					SolveCutAndBranch(instances);
					break;
				}
				default: option = -1; std::cout << "Invalid option!" << std::endl; break;
			}
		}while(option == -1);

		DeleteTimer();
		DeleteCallbackSelection();
	}
	catch(const std::runtime_error& re)
	{
		std::cout << "Runtime error: " << re.what() << std::endl;
	}
	catch(const std::exception& ex)
	{
		std::cout << "Error occurred: " << ex.what() << std::endl;
	}
	catch(const int& error)
	{
		std::cout << "Error occurred: " << error << std::endl;
	}
	catch (IloException& e)
	{
		std::cout << "Concert Exception: " << e << std::endl;
	}
	catch (const char * e)
	{
		std::cout << e << std::endl;
	}
	catch (const std::string& e)
	{
		std::cout << e << std::endl;
	}
	catch(...)
	{
		std::cout << "Unknown failure occurred. Possible memory corruption" << std::endl;
	}

	return 0;
}
