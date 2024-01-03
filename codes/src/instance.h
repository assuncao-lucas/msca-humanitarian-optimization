#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include "src/matrix.hpp"
#include "src/graph.h"

class Instance{
private:
	Graph* graph_ = nullptr;
	int num_vehicles_ = 0;
	int uncertainty_budget_ = 0;
	double service_time_deviation_ = 0.0;
	double limit_ = 0.0;
	std::string raw_file_name_;
	int num_mandatory_ = 0;
	bool found_maximal_cliques_ = false;
	double time_spent_in_preprocessing_ = 0.0;
	std::vector<std::list<int>> * conflict_graph_ = nullptr;
	Matrix<bool> conflict_matrix_;
	std::list<int> mandatory_list_;
	std::vector<std::list<int>> conflicts_list_; // list of cliques

	std::vector<std::list<int>> map_vertices_to_cliques_;
	boost::dynamic_bitset<> active_conflict_cliques_;
	std::vector<int> map_reordered_vertices_to_original_positions_;

	void GeneratePreProcessedGraph();
	void ReorderMandatoryVertices();
public:

	int num_opt_cuts_ = 0;
	int num_feas_cuts_ = 0;
	explicit Instance() = default;
	explicit Instance(std::string dir_path, std::string file_name, int num_vehicles, double service_time_deviation, int uncertainty_budged, bool pre_process_graph);
	virtual ~Instance();
	std::string GetInstanceName() const;
	void FillInstanceFromFile(std::string dir_path, std::string file_name, double service_time_deviation);
	void AddMandatoryVertices(double mandatory_percentage);
	void WriteToFile(std::string dir_path, std::string curr_file);
	void BuildVerticesToCliquesMapping();
	void FindAllMaximalConflictCliquesBronKerosch();
	void FindAllMaximalConflictCliquesTomita();
	void SelectMaximumCliquesPerVertex();
	void ComputeConflictGraph();
	void ResetConflictsCliques();
	int num_vehicles() const {return num_vehicles_;}
	int num_mandatory() const {return num_mandatory_;}
	bool found_maximal_cliques(){return found_maximal_cliques_;}
	bool FoundConflictGraph(){return conflict_graph_!= nullptr;}
	double limit() const {return limit_;}
	double time_spent_in_preprocessing() const {return time_spent_in_preprocessing_;}
	int uncertainty_budget() const {return uncertainty_budget_;}
	double service_time_deviation() const{return service_time_deviation_;}
	boost::dynamic_bitset<> active_conflict_cliques() { return active_conflict_cliques_;}
	const std::vector<std::list<int>>& map_vertices_to_cliques() const {return map_vertices_to_cliques_;}
	const std::vector<std::list<int>>& conflicts_list() const{ return conflicts_list_;}
	void set_graph(Graph* graph);
	const Graph* graph() const {return graph_;}
	friend std::ostream& operator<< (std::ostream &out, const Instance &instance);
};
