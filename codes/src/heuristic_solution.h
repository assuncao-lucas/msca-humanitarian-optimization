#ifndef HEURISTIC_SOLUTION_H_
#define HEURISTIC_SOLUTION_H_

#include <list>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "route.h"
#include "graph.h"
#include "instance.h"

class VertexStatus
{
public:
	VertexStatus();
	~VertexStatus();
	bool selected_;
	int route_;
	std::list<int>::iterator pos_;
};

class HeuristicSolution
{
public:
	HeuristicSolution();
	HeuristicSolution(int dimension, int dimension2, int num_routes);
	virtual ~HeuristicSolution();
	virtual void Reset(int dimension, int dimension2, int num_routes);
	bool is_infeasible_;
	bool is_feasible_;
	bool is_optimal_;
	int profits_sum_;
	int dimension_;
	int dimension2_;
	int num_routes_;
	std::vector<Route> routes_vec_;
	std::vector<VertexStatus> vertex_status_vec_;
	std::list<int> unvisited_vertices_;
	bool PreviewAddVertex(Instance &, int vertex, int route, int pos, std::list<int>::iterator& it, int& profit_variation, double& time_variation);
	void AddVertex(int vertex, int route, std::list<int>::iterator pos_it, int profit_variation, double time_variation);
	bool PreviewRemoveVertex(Instance &, int vertex, int& profit_variation, double& time_variation, bool allow_infeasible_routes = false);
	void RemoveVertex(int vertex, int profit_variation, double time_variation);
	bool PreviewInterRouteMoveVertex(Instance &, int vertex, int r2, int pos, std::list<int>::iterator& it, int& profit_variation1, double& time_variation1, int& profit_variation2, double& time_variation2);
	void InterRouteMoveVertex(int vertex, int r2, std::list<int>::iterator it, int profit_variation1, double time_variation1, int profit_variation2, double time_variation2);
	bool PreviewInterRouteSwap(Instance &, int r1, int pos_i1, int pos_f1, int r2, int pos_i2, int pos_f2, std::list<int>::iterator& it_i1, std::list<int>::iterator& it_f1,
		int& profit_variation1, double& time_variation1, std::list<int>::iterator& it_i2, std::list<int>::iterator& it_f2, int& profit_variation2, double& time_variation2);
	void InterRouteSwap(int r1, int r2, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1,
			int profit_variation1, double time_variation1, std::list<int>::iterator it_i2, std::list<int>::iterator it_f2, int profit_variation2, double time_variation2);

	bool PreviewInterRouteSwapUnrouted(Instance& inst, int r1, int pos_i1, int pos_f1, std::list<int>::iterator& it_i1, std::list<int>::iterator& it_f1,int& profit_variation1, double& time_variation1, std::list<int>::iterator it_i2, int& profit_variation2);
	void InterRouteSwapUnrouted(int r1, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1, int profit_variation1, double time_variation1, std::list<int>::iterator it_i2);

	bool PreviewAddVertexToRouteWithinMinimumDistanceIncrease(Instance &, int vertex, int route, std::list<int>::iterator& it, int& profit_variation, double& time_variation, bool allow_infeasible_routes = false);
	bool PreviewAddVertexWithinMinimumDistanceIncrease(Instance &, int vertex, int& route, std::list<int>::iterator& it, int& profit_variation, double& time_variation);

	bool Do2OptImprovement(const Graph* graph, const int& route);
	bool Do3OptImprovement(const Graph* graph, const int& route);
	bool operator==(HeuristicSolution&);

	friend std::ostream& operator<< (std::ostream &out, HeuristicSolution &sol);
	bool CheckCorrectness(Instance& instance);
	virtual void WriteToFile(std::string algo, std::string folder, std::string file_name);
	virtual void ReadFromFile(Instance & inst, std::string algo, std::string folder, std::string file_name);

	boost::dynamic_bitset<> bitset_arcs_;
    boost::dynamic_bitset<> bitset_vertices_;

	std::vector<double> x_values_;
	std::vector<double> y_values_;
	double dual_bound_ = 0.0;
	void BuildBitset(Instance &);
};

class FPHeuristicSolution: public HeuristicSolution
{
public:
	FPHeuristicSolution();
	FPHeuristicSolution(int dimension, int dimension2, int num_routes);
	~FPHeuristicSolution();
	int num_iterations_stage1_;
	int num_iterations_stage2_;
	int num_perturbations_stage1_;
	int num_perturbations_stage2_;
	int num_restarts_stage1_;
	int num_restarts_stage2_;
	bool found_y_integer_;
	bool found_x_integer_;
	double time_stage1_;
	double time_stage2_;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	static std::string GenerateFileName();
	void WriteToFile(std::string algo, std::string folder, std::string file_name);
};

class LBHeuristicSolution: public HeuristicSolution
{
public:
	LBHeuristicSolution();
	LBHeuristicSolution(int dimension, int dimension2, int num_routes);
	~LBHeuristicSolution();
	int num_iterations_;
	double time_;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	static std::string GenerateFileName();
	void WriteToFile(std::string algo, std::string folder, std::string file_name);
};

class ALNSHeuristicSolution: public HeuristicSolution
{
public:
	ALNSHeuristicSolution();
	ALNSHeuristicSolution(ALNSHeuristicSolution*);
	ALNSHeuristicSolution(int dimension, int dimension2, int num_routes);
	~ALNSHeuristicSolution();
	int num_iterations_;
	int initial_solution_profits_sum_;
	int last_improve_iteration_;
	double elapsed_time_;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	void WriteToFile(std::string algo, std::string folder, std::string file_name);
	void ReadFromFile(Instance & inst, std::string algo, std::string folder, std::string file_name);
	static std::string GenerateFileName();
};

#endif
