#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include "heuristic_solution.h"
#include "general.h"

VertexStatus::VertexStatus()
{
	this->selected_ = false;
	this->route_ = -1;
}

VertexStatus::~VertexStatus()
{
}

HeuristicSolution::HeuristicSolution()
{
	this->profits_sum_ = 0;
	this->dimension_ = 0;
	this->dimension2_ = 0;
	this->is_infeasible_ = false;
	this->is_feasible_ = false;
	this->is_optimal_ = false;
}

HeuristicSolution::~HeuristicSolution()
{
}

HeuristicSolution::HeuristicSolution(int dimension, int dimension2, int num_routes)
{
	HeuristicSolution::Reset(dimension, dimension2, num_routes);
}

void HeuristicSolution::Reset(int dimension, int dimension2, int num_routes)
{
	this->profits_sum_ = 0;
	this->dimension_ = dimension;
	this->dimension2_ = dimension2;
	this->num_routes_ = num_routes;
	this->routes_vec_ = std::vector<Route>(num_routes);
	this->vertex_status_vec_ = std::vector<VertexStatus>(dimension);
	this->is_infeasible_ = false;
	this->is_feasible_ = false;
	this->is_optimal_ = false;
	(this->unvisited_vertices_).clear();

	this->bitset_arcs_ = boost::dynamic_bitset<>(dimension2, 0);
	this->bitset_vertices_ = boost::dynamic_bitset<>(dimension, 0);

	for (int i = 1; i < dimension; ++i)
	{
		(this->unvisited_vertices_).push_front(i);
		VertexStatus *status = &((this->vertex_status_vec_)[i]);
		status->selected_ = false;
		status->route_ = -1;
		status->pos_ = (this->unvisited_vertices_).begin();
	}
}

bool HeuristicSolution::Do2OptImprovement(const Graph *graph, const int &route)
{
	Route *curr_route = &((this->routes_vec_)[route]);
	int size_route = (int)((curr_route->vertices_).size());

	if (size_route < 2)
		return false;

	std::list<int>::iterator vi_it;
	std::list<int>::iterator vk_it;
	std::list<int>::iterator pre_vi_it, pos_vk_it, it_init, it_end;
	int vi = 0, vk = 0, pre_vi = 0, pos_vk = 0, v1 = 0, v2 = 0;
	int num_swaps = 0, temp = 0, num_vertices = graph->num_vertices();
	GArc *pre_arc = nullptr, *pos_arc = nullptr, *new_pre_arc = nullptr, *new_pos_arc = nullptr, *curr_arc = nullptr;

	double time_variation = 0.0;

	vi_it = (curr_route->vertices_).begin();
	for (int i = 0; i < size_route - 1; i++)
	{
		vk_it = vi_it;
		vi = (*vi_it);

		if (i > 0)
		{
			pre_vi_it = vi_it;
			--pre_vi_it;
			pre_vi = (*pre_vi_it);
		}
		else
		{
			pre_vi = 0;
		}

		for (int k = i + 1; k < size_route; k++)
		{
			++vk_it;
			vk = (*vk_it);

			if (k < size_route - 1)
			{
				pos_vk_it = vk_it;
				++pos_vk_it;
				pos_vk = (*pos_vk_it);
			}
			else
			{
				pos_vk = num_vertices - 1;
			}

			// compute time variation properly!!
			new_pre_arc = (*graph)[pre_vi][vk];
			new_pos_arc = (*graph)[vi][pos_vk];

			if ((new_pre_arc == nullptr) || (new_pos_arc == nullptr))
				continue;

			pre_arc = (*graph)[pre_vi][vi];
			pos_arc = (*graph)[vk][pos_vk];

			// std::cout << "pre_vi: " << pre_vi << " vi: " << vi << " vk:" << vk << " pos_vk: " << pos_vk << std::endl;

			// compute time variation after 2-opt
			time_variation = (new_pre_arc->distance() + new_pos_arc->distance() - pre_arc->distance() - pos_arc->distance());

			it_init = vi_it;
			v2 = *it_init;
			it_end = vk_it;
			bool invalid_swap = false;

			do
			{
				v1 = v2;
				++it_init;
				v2 = *it_init;

				curr_arc = (*graph)[v1][v2];
				time_variation -= curr_arc->distance();
				curr_arc = (*graph)[v2][v1];
				if (curr_arc == nullptr)
				{
					invalid_swap = true;
					break;
				}
				else
					time_variation += curr_arc->distance();
			} while (it_init != it_end);
			// se há melhora no tempo da rota, faz o swap
			if ((!invalid_swap) && double_less(time_variation, 0.0))
			{
				// std::cout << *curr_route << std::endl;
				// std::cout << "após reverse Swap de " << vi << " a " << vk << std::endl;
				(curr_route->time_) += (time_variation);
				num_swaps = (k - i + 1) / 2;
				for (int cont = 0; cont < num_swaps; cont++)
				{
					temp = *vi_it;

					*vi_it = *vk_it;
					*vk_it = temp;

					// atualiza os iterators no PoleStatus
					((this->vertex_status_vec_)[*vi_it]).pos_ = vi_it;
					((this->vertex_status_vec_)[*vk_it]).pos_ = vk_it;

					++vi_it;
					--vk_it;
				}

				// std::cout << *curr_route << std::endl;
				// getchar(); getchar();

				return true;
			}
		}
		++vi_it;
	}
	return false;
}

bool HeuristicSolution::Do3OptImprovement(const Graph *graph, const int &route)
{
	Route *curr_route = &((this->routes_vec_)[route]);
	size_t size_route = ((curr_route->vertices_).size());

	if (size_route < 4)
		return false;

	std::list<int>::iterator vi_it;
	std::list<int>::iterator vj_it;
	std::list<int>::iterator vk_it;
	std::list<int>::iterator pre_vi_it, pos_vk_it, pos_vj_it, it_init, it_end;
	int vi = 0, vj = 0, vk = 0, pre_vi = 0, pos_vj = 0, pos_vk = 0, v1 = 0, v2 = 0;
	int num_vertices = graph->num_vertices();
	GArc *new_arc1 = nullptr, *new_arc2 = nullptr, *new_arc3 = nullptr, *curr_arc = nullptr;
	bool invalid_swap = false;

	double time_variation = 0.0, minus_variation = 0.0;

	std::list<int> temp;
	std::list<int>::iterator new_route_it, temp_it;

	vi_it = (curr_route->vertices_).begin();
	for (int i = 0; i <= size_route - 4; i++)
	{
		vi = (*vi_it);

		if (i > 0)
		{
			pre_vi_it = vi_it;
			--pre_vi_it;
			pre_vi = (*pre_vi_it);
		}
		else
		{
			pre_vi = 0;
		}

		vj_it = vi_it;
		++vj_it; // because next loop begins from i + 1
		for (int j = i + 1; j <= size_route - 3; j++)
		{
			vj = (*vj_it);

			pos_vj_it = vj_it;
			++pos_vj_it;
			pos_vj = (*pos_vj_it);

			vk_it = vj_it;
			std::advance(vk_it, 2); // because next loop begins from j + 2
			for (int k = j + 2; k <= size_route - 1; k++)
			{
				vk = (*vk_it);

				if (k < size_route - 1)
				{
					pos_vk_it = vk_it;
					pos_vk_it++;
					pos_vk = (*pos_vk_it);
				}
				else
				{
					pos_vk = num_vertices - 1;
				}

				minus_variation = ((*graph)[pre_vi][vi])->distance() + ((*graph)[vj][pos_vj])->distance() + ((*graph)[vk][pos_vk])->distance();

				// primeira possibilidade de reconexão
				invalid_swap = false;
				new_arc1 = (*graph)[pre_vi][pos_vj];
				new_arc2 = (*graph)[vk][vi];
				new_arc3 = (*graph)[vj][pos_vk];

				if ((new_arc1 == nullptr) || (new_arc2 == nullptr) || (new_arc3 == nullptr))
					invalid_swap = true;
				// time_variation = ((*min_dists)[actual_pre_vi][actual_pos_vj] + (*min_dists)[actual_vk][actual_vi] + (*min_dists)[actual_vj][actual_pos_vk] - minus_variation)/((this->curr_instance_)->speed_vehicle());

				if (!invalid_swap)
				{
					// in this case, does not need to reverse any part of the route
					time_variation = new_arc1->distance() + new_arc2->distance() + new_arc3->distance() - minus_variation;
				}

				// se há melhora no tempo da rota, faz o swap
				if ((!invalid_swap) && double_less(time_variation, 0.0))
				{
					// std::cout << *curr_route << std::endl;
					// std::cout << "após 3-opt com:" << std::endl;
					// std::cout << "vi: " << vi << ", pre_vi:" << pre_vi <<", vj: " << vj << ", pos_vj:" << pos_vj << ", vk: " << vk << ",pos_vk: " << pos_vk << std::endl;
					(curr_route->time_) += (time_variation);
					temp = curr_route->vertices_;

					new_route_it = vi_it;

					temp_it = temp.begin();
					std::advance(temp_it, j + 1);

					for (int cont = j + 1; cont <= k; cont++)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it++;
					}

					temp_it = temp.begin();
					std::advance(temp_it, i);

					for (int cont = i; cont <= j; cont++)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it++;
					}

					// std::cout << "case 1" << std::endl;
					// std::cout << "i: " << vi << " j: " << vj << " k: " << vk << std::endl;
					// std::cout << *curr_route << std::endl;
					// getchar();getchar();

					return true;
				}

				// segunda possibilidade de reconexão]
				invalid_swap = false;
				new_arc1 = (*graph)[pre_vi][vj];
				new_arc2 = (*graph)[vi][vk];
				new_arc3 = (*graph)[pos_vj][pos_vk];

				if ((new_arc1 == nullptr) || (new_arc2 == nullptr) || (new_arc3 == nullptr))
					invalid_swap = true;
				// time_variation = ((*min_dists)[actual_pre_vi][actual_vj] + (*min_dists)[actual_vi][actual_vk] + (*min_dists)[actual_pos_vj][actual_pos_vk] - minus_variation)/((this->curr_instance_)->speed_vehicle());

				if (!invalid_swap)
				{
					// in this case, does not need to reverse any part of the route
					time_variation = new_arc1->distance() + new_arc2->distance() + new_arc3->distance() - minus_variation;

					it_init = vi_it;
					v2 = *it_init;
					it_end = vj_it;
					// bool invalid_swap = false;

					while (it_init != it_end)
					{
						v1 = v2;
						++it_init;
						v2 = *it_init;

						curr_arc = (*graph)[v1][v2];
						time_variation -= curr_arc->distance();
						curr_arc = (*graph)[v2][v1];
						if (curr_arc == nullptr)
						{
							invalid_swap = true;
							break;
						}
						else
							time_variation += curr_arc->distance();
					}

					if (!invalid_swap)
					{
						it_init = pos_vj_it;
						v2 = *it_init;
						it_end = vk_it;
						// bool invalid_swap = false;

						while (it_init != it_end)
						{
							v1 = v2;
							++it_init;
							v2 = *it_init;

							curr_arc = (*graph)[v1][v2];
							time_variation -= curr_arc->distance();
							curr_arc = (*graph)[v2][v1];
							if (curr_arc == nullptr)
							{
								invalid_swap = true;
								break;
							}
							else
								time_variation += curr_arc->distance();
						}
					}
				}

				// se há melhora no tempo da rota, faz o swap
				if ((!invalid_swap) && double_less(time_variation, 0.0))
				{
					// std::cout << *curr_route << std::endl;
					// std::cout << "após 3-opt com:" << std::endl;
					// std::cout << "vi: " << vi << ", pre_vi:" << pre_vi <<", vj: " << vj << ", pos_vj:" << pos_vj << ", vk: " << vk << ",pos_vk: " << pos_vk << std::endl;
					(curr_route->time_) += (time_variation);
					temp = curr_route->vertices_;

					new_route_it = vi_it;

					temp_it = temp.begin();
					std::advance(temp_it, j);

					for (int cont = j; cont >= i; cont--)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it--;
					}

					temp_it = temp.begin();
					std::advance(temp_it, k);

					for (int cont = k; cont >= j + 1; cont--)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it--;
					}

					// std::cout << *curr_route << std::endl;
					// getchar();getchar();
					// std::cout << "case 2" << std::endl;
					// std::cout << "i: " << vi << " j: " << vj << " k: " << vk << std::endl;

					return true;
				}

				// terceira possibilidade de reconexão
				invalid_swap = false;
				new_arc1 = (*graph)[pre_vi][vk];
				new_arc2 = (*graph)[pos_vj][vi];
				new_arc3 = (*graph)[vj][pos_vk];

				if ((new_arc1 == nullptr) || (new_arc2 == nullptr) || (new_arc3 == nullptr))
					invalid_swap = true;
				// time_variation = ((*min_dists)[actual_pre_vi][actual_vk] + (*min_dists)[actual_pos_vj][actual_vi] + (*min_dists)[actual_vj][actual_pos_vk] - minus_variation)/((this->curr_instance_)->speed_vehicle());

				if (!invalid_swap)
				{
					// in this case, does not need to reverse any part of the route
					time_variation = new_arc1->distance() + new_arc2->distance() + new_arc3->distance() - minus_variation;

					it_init = pos_vj_it;
					v2 = *it_init;
					it_end = vk_it;
					// bool invalid_swap = false;

					while (it_init != it_end)
					{
						v1 = v2;
						++it_init;
						v2 = *it_init;

						curr_arc = (*graph)[v1][v2];
						time_variation -= curr_arc->distance();
						curr_arc = (*graph)[v2][v1];
						if (curr_arc == nullptr)
						{
							invalid_swap = true;
							break;
						}
						else
							time_variation += curr_arc->distance();
					}
				}

				// se há melhora no tempo da rota, faz o swap
				if ((!invalid_swap) && double_less(time_variation, 0.0))
				{
					// std::cout << *curr_route << std::endl;
					// std::cout << "após 3-opt com:" << std::endl;
					// std::cout << "vi: " << vi << ", pre_vi:" << pre_vi <<", vj: " << vj << ", pos_vj:" << pos_vj << ", vk: " << vk << ",pos_vk: " << pos_vk << std::endl;
					(curr_route->time_) += (time_variation);
					temp = curr_route->vertices_;

					new_route_it = vi_it;

					temp_it = temp.begin();
					std::advance(temp_it, k);

					for (int cont = k; cont >= j + 1; cont--)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it--;
					}

					temp_it = temp.begin();
					std::advance(temp_it, i);

					for (int cont = i; cont <= j; cont++)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it++;
					}

					// std::cout << *curr_route << std::endl;
					// getchar();getchar();
					// std::cout << "case 3" << std::endl;
					// std::cout << "i: " << vi << " j: " << vj << " k: " << vk << std::endl;

					return true;
				}

				// quarta possibilidade de reconexão
				invalid_swap = false;
				new_arc1 = (*graph)[pre_vi][pos_vj];
				new_arc2 = (*graph)[vk][vj];
				new_arc3 = (*graph)[vi][pos_vk];

				if ((new_arc1 == nullptr) || (new_arc2 == nullptr) || (new_arc3 == nullptr))
					invalid_swap = true;
				// time_variation = ((*min_dists)[actual_pre_vi][actual_pos_vj] + (*min_dists)[actual_vk][actual_vj] + (*min_dists)[actual_vi][actual_pos_vk] - minus_variation)/((this->curr_instance_)->speed_vehicle());

				if (!invalid_swap)
				{
					// in this case, does not need to reverse any part of the route
					time_variation = new_arc1->distance() + new_arc2->distance() + new_arc3->distance() - minus_variation;

					it_init = vi_it;
					v2 = *it_init;
					it_end = vj_it;
					// bool invalid_swap = false;

					while (it_init != it_end)
					{
						v1 = v2;
						++it_init;
						v2 = *it_init;

						curr_arc = (*graph)[v1][v2];
						time_variation -= curr_arc->distance();
						curr_arc = (*graph)[v2][v1];
						if (curr_arc == nullptr)
						{
							invalid_swap = true;
							break;
						}
						else
							time_variation += curr_arc->distance();
					}
				}

				// se há melhora no tempo da rota, faz o swap
				if ((!invalid_swap) && double_less(time_variation, 0.0))
				{
					// std::cout << *curr_route << std::endl;
					// std::cout << "após 3-opt com:" << std::endl;
					// std::cout << "vi: " << vi << ", pre_vi:" << pre_vi <<", vj: " << vj << ", pos_vj:" << pos_vj << ", vk: " << vk << ",pos_vk: " << pos_vk << std::endl;
					(curr_route->time_) += (time_variation);
					temp = curr_route->vertices_;

					new_route_it = vi_it;

					temp_it = temp.begin();
					std::advance(temp_it, j + 1);

					for (int cont = j + 1; cont <= k; cont++)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it++;
					}

					temp_it = temp.begin();
					std::advance(temp_it, j);

					for (int cont = j; cont >= i; cont--)
					{
						*new_route_it = *temp_it;

						((this->vertex_status_vec_)[*new_route_it]).pos_ = new_route_it;

						new_route_it++;
						temp_it--;
					}

					// std::cout << *curr_route << std::endl;
					// getchar();getchar();
					// std::cout << "case 4" << std::endl;
					// std::cout << "i: " << vi << " j: " << vj << " k: " << vk << std::endl;

					return true;
				}

				++vk_it;
			}
			vj_it++;
		}
		vi_it++;
	}
	return false;
}

void HeuristicSolution::AddVertex(int vertex, int route, std::list<int>::iterator it, int profit_variation, double time_variation)
{
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);

	// remove from list of unvisited_vertices
	this->unvisited_vertices_.erase(status->pos_);

	// add to route
	Route *curr_route = &((this->routes_vec_)[route]);
	int max_pos = (int)((curr_route->vertices_).size());
	// std::list<int>::iterator it = (curr_route->vertices_).begin();
	// if((pos == -1)||(pos > max_pos)) pos = max_pos;
	// if(pos > 0) std::advance(it,pos);

	status->selected_ = true;
	status->route_ = route;
	status->pos_ = (curr_route->vertices_).insert(it, vertex);

	(curr_route->time_) += time_variation;
	(curr_route->sum_profits_) += profit_variation;

	(this->profits_sum_) += profit_variation;
}

void HeuristicSolution::InterRouteSwap(int r1, int r2, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1,
									   int profit_variation1, double time_variation1, std::list<int>::iterator it_i2, std::list<int>::iterator it_f2, int profit_variation2, double time_variation2)
{
	Route *route1 = &((this->routes_vec_)[r1]), *route2 = &((this->routes_vec_)[r2]);

	// change the routes of the status of the vertices moved!!!!!!!!!!
	std::list<int>::iterator it = it_i1;
	do
	{
		((this->vertex_status_vec_)[*it]).route_ = r2;
		++it;
	} while (it != it_f1);

	it = it_i2;
	do
	{
		((this->vertex_status_vec_)[*it]).route_ = r1;
		++it;
	} while (it != it_f2);

	(route1->vertices_).splice(it_i1, route2->vertices_, it_i2, it_f2);
	(route2->vertices_).splice(it_f2, route1->vertices_, it_i1, it_f1);

	(route1->time_) += time_variation1;
	(route1->sum_profits_) += profit_variation1;

	(route2->time_) += time_variation2;
	(route2->sum_profits_) += profit_variation2;
}

void HeuristicSolution::InterRouteSwapUnrouted(int r1, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1, int profit_variation1, double time_variation1, std::list<int>::iterator it_i2)
{
	Route *route1 = &((this->routes_vec_)[r1]);

	// change the routes of the status of the vertices moved!!!!!!!!!!
	std::list<int>::iterator it = it_i1;
	do
	{
		((this->vertex_status_vec_)[*it]).route_ = -1;
		((this->vertex_status_vec_)[*it]).selected_ = false;
		//(this->unvisited_vertices_).push_front(*it);
		//((this->vertex_status_vec_)[*it]).pos_ = (this->unvisited_vertices_).begin();
		++it;
	} while (it != it_f1);

	((this->vertex_status_vec_)[*it_i2]).route_ = r1;
	((this->vertex_status_vec_)[*it_i2]).selected_ = true;

	(route1->vertices_).splice(it_i1, this->unvisited_vertices_, it_i2);
	(this->unvisited_vertices_).splice((this->unvisited_vertices_).begin(), route1->vertices_, it_i1, it_f1);

	(route1->time_) += time_variation1;
	(route1->sum_profits_) += profit_variation1;

	(this->profits_sum_) += profit_variation1;

	//(route2->time_) += time_variation2;
	//(route2->sum_profits_) += profit_variation2;
}

bool HeuristicSolution::PreviewInterRouteSwap(Instance &inst, int r1, int pos_i1, int pos_f1, int r2, int pos_i2, int pos_f2, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1,
											  int &profit_variation1, double &time_variation1, std::list<int>::iterator &it_i2, std::list<int>::iterator &it_f2, int &profit_variation2, double &time_variation2)
{
	const Graph *graph = inst.graph();
	profit_variation1 = profit_variation2 = 0;
	time_variation1 = time_variation2 = 0.0;
	int pre_vertex1 = 0, pos_vertex1 = 0, pre_vertex2 = 0, pos_vertex2 = 0;
	Route *route1 = &((this->routes_vec_)[r1]), *route2 = &((this->routes_vec_)[r2]);
	int max_pos1 = (int)((route1->vertices_).size()) - 1, max_pos2 = (int)((route2->vertices_).size()) - 1;
	GArc *new_pre_arc = nullptr, *new_pos_arc = nullptr;
	int begin_segment_vertex1, end_segment_vertex1, begin_segment_vertex2, end_segment_vertex2;

	int num_vertices = graph->num_vertices();

	int intra_route_profit_loss1 = 0, intra_route_profit_loss2 = 0;
	double intra_route_time_variation1 = 0.0, intra_route_time_variation2 = 0.0;

	if ((pos_i1 > pos_f1) || (pos_i1 > max_pos1) || (pos_f1 > max_pos1))
		return false;
	if ((pos_i2 > pos_f2) || (pos_i2 > max_pos2) || (pos_f2 > max_pos2))
		return false;

	// ROUTE 1
	// compute pre_vertex
	it_i1 = (route1->vertices_).begin();
	if (pos_i1 == 0)
		pre_vertex1 = 0;
	else
	{
		if (pos_i1 > 1)
			std::advance(it_i1, pos_i1 - 1);
		pre_vertex1 = *it_i1;
		++it_i1;
	}

	begin_segment_vertex1 = *it_i1;
	// std::cout << "pre_vertex1: " << pre_vertex1 << std::endl;
	// std::cout << "begin_seg_vertex1: " << begin_segment_vertex1 << std::endl;

	it_f1 = it_i1;
	int v1 = pre_vertex1, v2 = 0;

	for (int i = 0; i <= pos_f1 - pos_i1; i++)
	{
		v2 = *it_f1;

		intra_route_profit_loss1 += (graph->vertices_info())[v2].profit_;
		time_variation1 -= ((*graph)[v1][v2])->distance();
		if (i > 0)
			intra_route_time_variation1 += ((*graph)[v1][v2])->distance();
		if (i == pos_f1 - pos_i1)
			end_segment_vertex1 = v2;
		++it_f1;
		v1 = v2;
	}

	// compute pos_vertex
	if (pos_f1 == max_pos1)
		pos_vertex1 = num_vertices - 1; // at this point, it_f1 can be .end()!
	else
		pos_vertex1 = *it_f1;

	// std::cout << "pos_vertex1: " << pos_vertex1 << std::endl;
	// std::cout << "end_seg_vertex1: " << end_segment_vertex1 << std::endl;

	time_variation1 -= ((*graph)[v1][pos_vertex1])->distance();

	// ROUTE 2
	// compute pre_vertex
	it_i2 = (route2->vertices_).begin();
	if (pos_i2 == 0)
		pre_vertex2 = 0;
	else
	{
		if (pos_i2 > 1)
			std::advance(it_i2, pos_i2 - 1);
		pre_vertex2 = *it_i2;
		++it_i2;
	}

	begin_segment_vertex2 = *it_i2;
	// std::cout << "pre_vertex2: " << pre_vertex2 << std::endl;
	// std::cout << "begin_seg_vertex2: " << begin_segment_vertex2 << std::endl;
	it_f2 = it_i2;
	v1 = pre_vertex2;
	v2 = 0;

	for (int i = 0; i <= pos_f2 - pos_i2; i++)
	{
		v2 = *it_f2;

		intra_route_profit_loss2 += (graph->vertices_info())[v2].profit_;
		time_variation2 -= ((*graph)[v1][v2])->distance();
		if (i > 0)
			intra_route_time_variation2 += ((*graph)[v1][v2])->distance();
		if (i == pos_f2 - pos_i2)
			end_segment_vertex2 = v2;
		++it_f2;
		v1 = v2;
	}

	// compute pos_vertex
	if (pos_f2 == max_pos2)
		pos_vertex2 = num_vertices - 1; // at this point, it_f2 can be .end()!
	else
		pos_vertex2 = *it_f2;

	// std::cout << "pos_vertex2: " << pos_vertex2 << std::endl;
	// std::cout << "end_seg_vertex2: " << end_segment_vertex2 << std::endl;

	time_variation2 -= ((*graph)[v1][pos_vertex2])->distance();

	// route 1
	profit_variation1 = intra_route_profit_loss2 - intra_route_profit_loss1;

	new_pre_arc = ((*graph)[pre_vertex1][begin_segment_vertex2]);
	new_pos_arc = ((*graph)[end_segment_vertex2][pos_vertex1]);

	if ((new_pre_arc == nullptr) || (new_pos_arc == nullptr))
		return false;
	time_variation1 += (new_pre_arc->distance() + intra_route_time_variation2 + new_pos_arc->distance());

	if (double_greater(route1->time_ + time_variation1, inst.limit()))
		return false;

	// route 2
	profit_variation2 = intra_route_profit_loss1 - intra_route_profit_loss2;

	new_pre_arc = ((*graph)[pre_vertex2][begin_segment_vertex1]);
	new_pos_arc = ((*graph)[end_segment_vertex1][pos_vertex2]);

	if ((new_pre_arc == nullptr) || (new_pos_arc == nullptr))
		return false;
	time_variation2 += (new_pre_arc->distance() + intra_route_time_variation1 + new_pos_arc->distance());

	if (double_greater(route2->time_ + time_variation2, inst.limit()))
		return false;

	return true;
}

bool HeuristicSolution::PreviewInterRouteSwapUnrouted(Instance &inst, int r1, int pos_i1, int pos_f1, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1,
													  int &profit_variation1, double &time_variation1, std::list<int>::iterator it_i2, int &profit_variation2)
{
	const Graph *graph = inst.graph();
	profit_variation1 = profit_variation2 = 0;
	time_variation1 = 0.0;
	int pre_vertex1 = 0, pos_vertex1 = 0;
	Route *route1 = &((this->routes_vec_)[r1]);
	int max_pos1 = (int)((route1->vertices_).size()) - 1;
	GArc *new_pre_arc = nullptr, *new_pos_arc = nullptr;
	int begin_segment_vertex1, end_segment_vertex1;

	int num_vertices = graph->num_vertices();

	int intra_route_profit_loss1 = 0, intra_route_profit_loss2 = 0;
	double intra_route_time_variation1 = 0.0;

	if ((pos_i1 > pos_f1) || (pos_i1 > max_pos1) || (pos_f1 > max_pos1))
		return false;

	// ROUTE 1
	// compute pre_vertex
	it_i1 = (route1->vertices_).begin();
	if (pos_i1 == 0)
		pre_vertex1 = 0;
	else
	{
		if (pos_i1 > 1)
			std::advance(it_i1, pos_i1 - 1);
		pre_vertex1 = *it_i1;
		++it_i1;
	}

	begin_segment_vertex1 = *it_i1;
	// std::cout << "pre_vertex1: " << pre_vertex1 << std::endl;
	// std::cout << "begin_seg_vertex1: " << begin_segment_vertex1 << std::endl;

	it_f1 = it_i1;
	int v1 = pre_vertex1, v2 = 0;

	for (int i = 0; i <= pos_f1 - pos_i1; i++)
	{
		v2 = *it_f1;

		intra_route_profit_loss1 += (graph->vertices_info())[v2].profit_;
		time_variation1 -= ((*graph)[v1][v2])->distance();
		if (i > 0)
			intra_route_time_variation1 += ((*graph)[v1][v2])->distance();
		if (i == pos_f1 - pos_i1)
			end_segment_vertex1 = v2;
		++it_f1;
		v1 = v2;
	}

	// compute pos_vertex
	if (pos_f1 == max_pos1)
		pos_vertex1 = num_vertices - 1; // at this point, it_f1 can be .end()!
	else
		pos_vertex1 = *it_f1;

	// std::cout << "pos_vertex1: " << pos_vertex1 << std::endl;
	// std::cout << "end_seg_vertex1: " << end_segment_vertex1 << std::endl;

	time_variation1 -= ((*graph)[v1][pos_vertex1])->distance();

	// ROUTE 2
	intra_route_profit_loss2 = (graph->vertices_info())[*it_i2].profit_;

	// route 1
	profit_variation1 = intra_route_profit_loss2 - intra_route_profit_loss1;

	new_pre_arc = ((*graph)[pre_vertex1][*it_i2]);
	new_pos_arc = ((*graph)[*it_i2][pos_vertex1]);

	if ((new_pre_arc == nullptr) || (new_pos_arc == nullptr))
		return false;
	time_variation1 += (new_pre_arc->distance() + new_pos_arc->distance());

	if (double_greater(route1->time_ + time_variation1, inst.limit()))
		return false;

	// route 2
	profit_variation2 = intra_route_profit_loss1 - intra_route_profit_loss2;

	return true;
}

bool HeuristicSolution::PreviewAddVertexWithinMinimumDistanceIncrease(Instance &inst, int vertex, int &route, std::list<int>::iterator &it, int &profit_variation, double &time_variation)
{
	profit_variation = 0;
	time_variation = std::numeric_limits<double>::infinity();
	int curr_profit_variation = 0;
	double curr_time_variation = 0.0;
	std::list<int>::iterator curr_it;
	bool can_add = false;

	for (int curr_route = 0; curr_route < this->num_routes_; ++curr_route)
	{
		if (PreviewAddVertexToRouteWithinMinimumDistanceIncrease(inst, vertex, curr_route, curr_it, curr_profit_variation, curr_time_variation))
		{
			if (double_less(curr_time_variation, time_variation))
			{
				can_add = true;
				route = curr_route;
				it = curr_it;
				profit_variation = curr_profit_variation;
				time_variation = curr_time_variation;
			}
		}
	}

	return can_add;
}

bool HeuristicSolution::PreviewAddVertexToRouteWithinMinimumDistanceIncrease(Instance &inst, int vertex, int route, std::list<int>::iterator &it, int &profit_variation, double &time_variation, bool allow_infeasible_routes)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = std::numeric_limits<double>::infinity();
	double curr_time_variation = 0.0;
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((this->routes_vec_)[route]);
	size_t max_pos = (int)((curr_route->vertices_).size());
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	int num_vertices = graph->num_vertices();
	bool can_add = false;

	if (status->selected_)
		return false;

	std::list<int>::iterator curr_it = (curr_route->vertices_).begin();
	for (int pos = 0; pos <= max_pos; ++pos)
	{
		// compute pos_vertex
		if (pos == max_pos)
			pos_vertex = num_vertices - 1;
		else
			pos_vertex = *curr_it;

		GArc *pre_arc = (*graph)[pre_vertex][vertex], *pos_arc = (*graph)[vertex][pos_vertex], *curr_arc = (*graph)[pre_vertex][pos_vertex];

		if ((pre_arc != nullptr) && (pos_arc != nullptr) && (curr_arc != nullptr))
		{
			curr_time_variation = (pre_arc->distance() + pos_arc->distance() - curr_arc->distance());
			if (!(double_greater(curr_route->time_ + curr_time_variation, inst.limit())) || allow_infeasible_routes)
			{
				// simulate addition to route
				if (double_less(curr_time_variation, time_variation))
				{
					can_add = true;
					time_variation = curr_time_variation;
					it = curr_it;
				}
			}
		}
		++curr_it;
		pre_vertex = pos_vertex;
	}

	if (can_add)
	{
		profit_variation = (graph->vertices_info())[vertex].profit_;
		return true;
	}

	return false;
}

bool HeuristicSolution::PreviewAddVertex(Instance &inst, int vertex, int route, int pos, std::list<int>::iterator &it, int &profit_variation, double &time_variation)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = 0.0;
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((this->routes_vec_)[route]);
	int max_pos = (int)((curr_route->vertices_).size());
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	int num_vertices = graph->num_vertices();

	if ((pos == -1) || (pos > max_pos))
		pos = max_pos;

	// compute pre_vertex
	it = (curr_route->vertices_).begin();
	if (pos == 0)
		pre_vertex = 0;
	else
	{
		if (pos > 1)
			std::advance(it, pos - 1);
		pre_vertex = *it;
		++it;
	}

	// compute pos_vertex
	if (pos == max_pos)
		pos_vertex = num_vertices - 1;
	else
		pos_vertex = *it;

	/*if(pos == 0) pre_vertex = 0;
	else
	{
		std::list<int>::iterator pre_vertex_it = (curr_route->vertices_).begin();
		if(pos > 1) std::advance(pre_vertex_it,pos-1);
		pre_vertex = *pre_vertex_it;
	}

	//std::cout << "pre_vertex: " << pre_vertex << std::endl;

	if(pos == max_pos) pos_vertex = num_vertices - 1;
	else
	{
		std::list<int>::iterator pos_vertex_it = (curr_route->vertices_).begin();
		if(pos > 0) std::advance(pos_vertex_it,pos);
		pos_vertex = *pos_vertex_it;
	}*/

	// std::cout << "pos_vertex: " << pos_vertex << std::endl;

	GArc *pre_arc = (*graph)[pre_vertex][vertex], *pos_arc = (*graph)[vertex][pos_vertex], *curr_arc = (*graph)[pre_vertex][pos_vertex];

	if ((!(status->selected_)) && (pre_arc != nullptr) && (pos_arc != nullptr) && (curr_arc != nullptr))
	{
		// simulate addition to route
		profit_variation = (graph->vertices_info())[vertex].profit_;

		// if route is currently empty
		// if(max_pos == 0) time_variation = pre_arc->distance() + pos_arc->distance();
		// else
		time_variation = (pre_arc->distance() + pos_arc->distance() - curr_arc->distance());

		if (!double_greater(curr_route->time_ + time_variation, inst.limit()))
			return true;
	}
	return false;
}

void HeuristicSolution::RemoveVertex(int vertex, int profit_variation, double time_variation)
{
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	// std::cout << "removing from route " << status->route_ << std::endl;
	//  remove from route
	Route *curr_route = &((this->routes_vec_)[status->route_]);
	// std::cout << *curr_route << std::endl;
	(curr_route->vertices_).erase(status->pos_);
	(curr_route->time_) += time_variation;
	(curr_route->sum_profits_) += profit_variation;

	// add to list of unvisited_vertices_
	(this->unvisited_vertices_).push_front(vertex);
	status->selected_ = false;
	status->route_ = -1;
	status->pos_ = (this->unvisited_vertices_).begin();

	(this->profits_sum_) += profit_variation;
}

bool HeuristicSolution::PreviewRemoveVertex(Instance &inst, int vertex, int &profit_variation, double &time_variation, bool allow_infeasible_routes)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = 0.0;
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	if (!(status->selected_))
		return false;

	int num_vertices = graph->num_vertices();
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((this->routes_vec_)[status->route_]);
	std::list<int>::iterator pos = status->pos_;
	int max_pos = (int)((curr_route->vertices_).size());

	if (max_pos == 0)
		return false;

	if (pos == (curr_route->vertices_).begin())
		pre_vertex = 0;
	else
	{
		std::list<int>::iterator pre_vertex_it = pos;
		--pre_vertex_it;
		pre_vertex = *pre_vertex_it;
	}
	// std::cout << "pre_vertex: " << pre_vertex << std::endl;

	std::list<int>::iterator pos_vertex_it = pos; // at this point, pos is never the .end()!!!
	++pos_vertex_it;
	if (pos_vertex_it == (curr_route->vertices_).end())
		pos_vertex = num_vertices - 1;
	else
		pos_vertex = *pos_vertex_it;

	// std::cout << "pos_vertex: " << pos_vertex << std::endl;

	GArc *pre_arc = (*graph)[pre_vertex][vertex], *pos_arc = (*graph)[vertex][pos_vertex], *new_arc = (*graph)[pre_vertex][pos_vertex];

	// if((new_arc != nullptr)||(max_pos == 1))
	if (new_arc != nullptr)
	{
		// simulate removal from route
		profit_variation = -(graph->vertices_info())[vertex].profit_;

		// if the route becomes empty
		// if(max_pos == 1) time_variation = -(pre_arc->distance() + pos_arc->distance());
		// else
		time_variation = -(pre_arc->distance() + pos_arc->distance() - new_arc->distance());

		if (!double_greater(curr_route->time_ + time_variation, inst.limit()) || allow_infeasible_routes)
			return true;
	}
	return false;
}

bool HeuristicSolution::PreviewInterRouteMoveVertex(Instance &inst, int vertex, int r2, int pos, std::list<int>::iterator &it, int &profit_variation1, double &time_variation1, int &profit_variation2, double &time_variation2)
{
	const Graph *graph = inst.graph();
	bool can_move = false;
	profit_variation1 = profit_variation2 = 0;
	time_variation1 = time_variation2 = 0.0;
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	if (status->route_ != r2)
	{
		if (PreviewRemoveVertex(inst, vertex, profit_variation1, time_variation1))
		{
			// to do the preview of addition!!!
			status->selected_ = false;
			if (PreviewAddVertex(inst, vertex, r2, pos, it, profit_variation2, time_variation2))
			{
				can_move = true;
			}
			status->selected_ = true;
		}
	}
	return can_move;
}

void HeuristicSolution::InterRouteMoveVertex(int vertex, int r2, std::list<int>::iterator it, int profit_variation1, double time_variation1, int profit_variation2, double time_variation2)
{
	VertexStatus *status = &((this->vertex_status_vec_)[vertex]);
	// remove from route r1
	Route *curr_route = &((this->routes_vec_)[status->route_]);
	(curr_route->vertices_).erase(status->pos_);
	(curr_route->time_) += time_variation1;
	(curr_route->sum_profits_) += profit_variation1;

	// add to route r2
	curr_route = &((this->routes_vec_)[r2]);
	int max_pos = (int)((curr_route->vertices_).size());
	// std::list<int>::iterator it = (curr_route->vertices_).begin();
	// if((pos == -1)||(pos > max_pos)) pos = max_pos;
	// if(pos > 0) std::advance(it,pos);

	status->route_ = r2;
	status->pos_ = (curr_route->vertices_).insert(it, vertex);

	(curr_route->time_) += time_variation2;
	(curr_route->sum_profits_) += profit_variation2;
}

bool HeuristicSolution::CheckCorrectness(Instance &instance)
{
	const Graph *graph = instance.graph();
	int num_vertices = graph->num_vertices(), num_mandatory = instance.num_mandatory(), num_vehicles = instance.num_vehicles(), num_arcs = graph->num_arcs();
	int total_profits = 0, curr_profits = 0;
	double curr_time = 0.0;
	Route *curr_route = nullptr;
	VertexStatus *curr_status = nullptr;
	GArc *curr_arc = nullptr;
	int v1 = 0, v2 = 0;
	int count_mandatory = 0;
	boost::dynamic_bitset<> visited_vertices(num_vertices, 0);
	boost::dynamic_bitset<> visited_arcs(num_arcs, 0);

	visited_vertices[0] = visited_vertices[num_vertices - 1] = 1;

	if (this->num_routes_ != num_vehicles)
		throw "Invalid solution 1";

	if (this->is_infeasible_)
		return true;

	for (int i = 0; i < num_vehicles; i++)
	{
		curr_route = &((this->routes_vec_)[i]);
		curr_profits = v1 = 0;
		curr_time = 0.0;
		if (!((curr_route->vertices_).empty()))
		{
			for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			{
				v2 = *it;
				if (visited_vertices[v2])
					throw "Invalid solution 2";
				visited_vertices[v2] = 1;

				if ((v2 > 0) && (v2 <= num_mandatory))
					++count_mandatory;
				curr_status = &((this->vertex_status_vec_)[v2]);
				if ((!(curr_status->selected_)) || (it != curr_status->pos_) || (curr_status->route_ != i))
					throw "Invalid solution 3";
				if (v2 > num_mandatory)
				{
					total_profits += (graph->vertices_info())[v2].profit_;
					curr_profits += (graph->vertices_info())[v2].profit_;
				}

				curr_arc = (*graph)[v1][v2];
				if (curr_arc == nullptr)
					throw "Invalid solution 4";
				curr_time += (curr_arc->distance());
				visited_arcs[graph->pos(v1, v2)] = 1;
				v1 = v2;
			}

			v2 = num_vertices - 1;
			curr_arc = (*graph)[v1][v2];
			if (curr_arc == nullptr)
				throw "Invalid solution 5";
			curr_time += (curr_arc->distance());
			visited_arcs[graph->pos(v1, v2)] = 1;
		}
		else
		{
			curr_arc = (*graph)[0][num_vertices - 1];
			if (curr_arc == nullptr)
				throw "Invalid solution 6";
			curr_time = curr_arc->distance();
			visited_arcs[graph->pos(0, num_vertices - 1)] = 1;
		}

		if (!double_equals(curr_time, curr_route->time_))
		{
			std::cout << curr_time << " " << curr_route->time_ << std::endl;
			throw "Invalid solution 7";
		}
		if (double_greater(curr_route->time_, instance.limit()))
			throw "Invalid solution 7.2";
		if (curr_profits != curr_route->sum_profits_)
			throw "Invalid solution 8";
	}

	if (total_profits != this->profits_sum_)
		throw "Invalid solution 9";

	if (visited_vertices.count() != ((size_t)num_vertices - (this->unvisited_vertices_).size()))
		throw "Invalid solution 10";

	for (auto it = (this->unvisited_vertices_).begin(); it != (this->unvisited_vertices_).end(); ++it)
	{
		curr_status = &((this->vertex_status_vec_)[*it]);
		if ((curr_status->selected_) || (curr_status->pos_ != it) || (curr_status->route_ != -1) || visited_vertices[*it])
			throw "Invalid solution 11";
	}

	if (count_mandatory != num_mandatory)
		throw "invalid solution 12";

	if (visited_arcs != this->bitset_arcs_)
		throw "Invalid solution 13";
	if (visited_vertices != this->bitset_vertices_)
		throw "Invalid solution 14";
	return true;
}

void HeuristicSolution::BuildBitset(const Instance &instance)
{
	this->bitset_arcs_ = boost::dynamic_bitset<>(this->dimension2_, 0);
	this->bitset_vertices_ = boost::dynamic_bitset<>(this->dimension_, 0);
	const Graph *graph = instance.graph();
	int num_vertices = graph->num_vertices(), num_vehicles = instance.num_vehicles();
	Route *curr_route = nullptr;
	int v1 = 0, v2 = 0;

	if (this->is_infeasible_)
		return;

	(this->bitset_vertices_)[0] = 1;

	for (int i = 0; i < num_vehicles; i++)
	{
		v1 = 0;
		curr_route = &((this->routes_vec_)[i]);
		if (!((curr_route->vertices_).empty()))
		{
			for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			{
				v2 = *it;
				(this->bitset_vertices_)[v2] = 1;

				(this->bitset_arcs_)[graph->pos(v1, v2)] = 1;
				v1 = v2;
			}

			v2 = 0;
			(this->bitset_arcs_)[graph->pos(v1, v2)] = 1;
		}
	}
}

bool HeuristicSolution::operator==(HeuristicSolution &other)
{
	return (this->bitset_arcs_ == other.bitset_arcs_);
}

std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol)
{
	if (sol.is_infeasible_)
		out << "[INFEASIBLE]";
	if (sol.is_optimal_)
		out << "[OPTMAL]";
	out << "[profits_sum: " << sol.profits_sum_ << "]" << std::endl;
	for (int i = 0; i < sol.num_routes_; ++i)
		out << ((sol.routes_vec_)[i]);
	out << "unvisited: ";
	for (auto it = (sol.unvisited_vertices_).begin(); it != (sol.unvisited_vertices_).end(); ++it)
		out << *it << " ";
	out << std::endl;
	return out;
}

void HeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name)
{
	Route *curr_route = nullptr;
	std::fstream file;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out | std::fstream::app);

	file << std::setprecision(5) << std::fixed;

	file << K_FILE_DELIMITER << std::endl
		 << "Total profit sum: " << this->profits_sum_ << std::endl
		 << "Num routes: ";
	if ((this->is_infeasible_) || (!(this->is_feasible_)))
		file << "0" << std::endl;
	else
	{
		file << this->num_routes_ << std::endl;

		for (int i = 0; i < this->num_routes_; i++)
		{
			curr_route = &((this->routes_vec_)[i]);
			file << curr_route->sum_profits_ << " " << curr_route->time_ << " ";

			for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			{
				file << instance.getOriginalVertexPosition(*it) << " ";
			}

			// file << this->dimension_ - 1;
			file << std::endl;
		}
	}

	file.close();
}

void HeuristicSolution::ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	const Graph *graph = inst.graph();
	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();

	this->Reset(num_vertices, num_arcs, inst.num_vehicles());
	Route *curr_route = nullptr;
	int num_routes = 0;
	std::fstream input;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);
	std::string line;

	if (!(input.is_open()))
		throw "Error opening file 1";
	while (!(input.eof()))
	{
		getline(input, line);
		if (line == K_FILE_DELIMITER)
		{
			std::stringstream s_profits_sum, s_num_routes;

			getline(input, line);
			size_t pos = line.find_first_of(":");
			s_profits_sum << line.substr(pos + 2);
			s_profits_sum >> this->profits_sum_;
			// std::cout << this->profits_sum_ << std::endl;

			getline(input, line);
			pos = line.find_first_of(":");
			s_num_routes << line.substr(pos + 2);
			s_num_routes >> num_routes;

			if ((num_routes == 0) || (num_routes != this->num_routes_))
				this->is_feasible_ = false;
			else
			{
				this->is_feasible_ = true;
				int pre_vertex = 0, curr_vertex = 0;
				VertexStatus *status = nullptr;
				for (int i = 0; i < num_routes; i++)
				{
					curr_route = &((this->routes_vec_)[i]);
					getline(input, line);
					std::stringstream s_route(line);
					// std::cout << s_route.str() << std::endl;
					s_route >> curr_route->sum_profits_ >> curr_route->time_;
					// std::cout << curr_route->sum_profits_ << " " << curr_route->time_ << std::endl;
					curr_vertex = pre_vertex = 0;
					while (s_route >> curr_vertex)
					{
						status = &((this->vertex_status_vec_)[curr_vertex]);

						// remove from list of unvisited_vertices
						(this->unvisited_vertices_).erase(status->pos_);

						// adds vertex to route
						status->selected_ = true;
						status->route_ = i;
						status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), curr_vertex);

						//(this->bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
						pre_vertex = curr_vertex;
					}

					curr_vertex = num_vertices - 1;
					//(this->bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
				}
			}
			break;
		}
	}

	if (!(this->is_infeasible_) && (this->is_feasible_) && (this->unvisited_vertices_).empty())
		this->is_optimal_ = true;

	input.close();
}

ALNSHeuristicSolution::ALNSHeuristicSolution() : HeuristicSolution()
{
	this->initial_solution_profits_sum_ = 0;
	this->num_iterations_ = 0;
	this->elapsed_time_ = 0.0;
}

ALNSHeuristicSolution::ALNSHeuristicSolution(int dimension, int dimension2, int num_routes) : HeuristicSolution(dimension, dimension2, num_routes)
{
	this->initial_solution_profits_sum_ = 0;
	this->num_iterations_ = 0;
	this->last_improve_iteration_ = 0;
	this->elapsed_time_ = 0.0;
}

ALNSHeuristicSolution::ALNSHeuristicSolution(ALNSHeuristicSolution *sol)
{
	(*this) = (*sol);

	VertexStatus *status = nullptr;
	Route *curr_route = nullptr;

	for (int i = 0; i < this->num_routes_; ++i)
	{
		curr_route = &((this->routes_vec_)[i]);
		for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			((this->vertex_status_vec_)[*it]).pos_ = it;
	}

	for (auto it = (this->unvisited_vertices_).begin(); it != (this->unvisited_vertices_).end(); ++it)
		((this->vertex_status_vec_)[*it]).pos_ = it;
}

ALNSHeuristicSolution::~ALNSHeuristicSolution()
{
}

void ALNSHeuristicSolution::Reset(int dimension, int dimension2, int num_routes)
{
	HeuristicSolution::Reset(dimension, dimension2, num_routes);

	this->initial_solution_profits_sum_ = 0;
	this->num_iterations_ = 0;
	this->last_improve_iteration_ = 0;
	this->elapsed_time_ = 0.0;
}

void ALNSHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name)
{
	std::fstream file;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(5) << std::fixed;

	if (this->is_infeasible_)
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (this->is_optimal_)
		file << "STATUS: OPTIMAL" << std::endl;
	else if (this->is_feasible_)
		file << "STATUS: FEASIBLE" << std::endl;
	else
		file << "STATUS: POSSIBLY FEASIBLE" << std::endl;
	file << "profits sum: " << this->profits_sum_ << std::endl;
	file << "initial profits sum: " << this->initial_solution_profits_sum_ << std::endl;
	file << "# iterations: " << this->num_iterations_ << std::endl;
	file << "time(s): " << this->elapsed_time_ << std::endl;
	file << "last iteration of improvement: " << this->last_improve_iteration_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}

void ALNSHeuristicSolution::ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	HeuristicSolution::ReadFromFile(inst, algo, folder, file_name);

	std::fstream input;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);

	if (!(input.is_open()))
		throw "Error opening file 2";
	std::string line;
	getline(input, line);

	if (line == "STATUS: FOUND ARC INTEGER FEASIBLE")
		this->is_feasible_ = true;
	if (line == "STATUS: INFEASIBLE")
	{
		this->is_infeasible_ = true;
		this->is_feasible_ = false;
	}

	input.close();
}

std::string ALNSHeuristicSolution::GenerateFileName()
{

	std::string file_name = FPHeuristicSolution::GenerateFileName();
	if (K_ALNS_MULTI_THREAD)
		file_name += "_MT";
	if ((K_PATH_RELINKING) && (!K_STOP_AT_FIRST_PR_DETECTION))
		file_name += "_PR";
	file_name += ("_ALNS_");
	if (K_STOP_AT_FIRST_PR_DETECTION)
		file_name += "CONV";
	else
		file_name += std::to_string(K_NUM_ITERATIONS_ALNS);

	return file_name;
}

FPHeuristicSolution::FPHeuristicSolution(int dimension, int dimension2, int num_routes) : HeuristicSolution(dimension, dimension2, num_routes)
{
	this->num_iterations_stage1_ = 0;
	this->num_iterations_stage2_ = 0;
	this->num_perturbations_stage1_ = 0;
	this->num_perturbations_stage2_ = 0;
	this->num_restarts_stage1_ = 0;
	this->num_restarts_stage2_ = 0;
	this->found_y_integer_ = false;
	this->found_x_integer_ = false;
	this->time_stage1_ = 0.0;
	this->time_stage2_ = 0.0;
}

FPHeuristicSolution::FPHeuristicSolution() : HeuristicSolution()
{
	this->num_iterations_stage1_ = 0;
	this->num_iterations_stage2_ = 0;
	this->num_perturbations_stage1_ = 0;
	this->num_perturbations_stage2_ = 0;
	this->num_restarts_stage1_ = 0;
	this->num_restarts_stage2_ = 0;
	this->found_y_integer_ = false;
	this->found_x_integer_ = false;
	this->time_stage1_ = 0.0;
	this->time_stage2_ = 0.0;
}

FPHeuristicSolution::~FPHeuristicSolution()
{
}

void FPHeuristicSolution::Reset(int dimension, int dimension2, int num_routes)
{
	HeuristicSolution::Reset(dimension, dimension2, num_routes);

	this->num_iterations_stage1_ = 0;
	this->num_iterations_stage2_ = 0;
	this->num_perturbations_stage1_ = 0;
	this->num_perturbations_stage2_ = 0;
	this->num_restarts_stage1_ = 0;
	this->num_restarts_stage2_ = 0;
	this->found_y_integer_ = false;
	this->found_x_integer_ = false;
	this->time_stage1_ = 0.0;
	this->time_stage2_ = 0.0;
}

std::string FPHeuristicSolution::GenerateFileName()
{
	std::string file_name;
	if (K_STOP)
		file_name += "stop_";
	file_name += "fp";
	if (K_FEASIBILITY_PUMP_ADD_CUTS)
		file_name += "_cuts";
	if (K_SOLVE_STAGE_1)
		file_name += "_alpha1_" + std::to_string(int(100 * initial_alpha_stage1));
	file_name += ("_alpha2_" + std::to_string(int(100 * initial_alpha_stage2)) + "_p_");
	if (FIXED_FLIP_BASIS)
		file_name += ("fixed_" + std::to_string(K_FLIP_BASIS));
	else
		file_name += std::to_string(int(100 * perturbation_flip_percentage));
	return file_name;
}

void FPHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name)
{
	std::fstream file;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(5) << std::fixed;

	if (this->is_infeasible_)
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (this->found_x_integer_)
		file << "STATUS: FOUND ARC INTEGER FEASIBLE" << std::endl;
	else if (this->found_y_integer_)
		file << "STATUS: FOUND VERTEX INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;
	file << "profit sum: " << this->profits_sum_ << std::endl;
	file << "STAGE 1: " << std::endl;
	file << "# iterations: " << this->num_iterations_stage1_ << std::endl;
	file << "# perturbations: " << this->num_perturbations_stage1_ << std::endl;
	file << "# restarts: " << this->num_restarts_stage1_ << std::endl;
	file << "time(s): " << this->time_stage1_ << std::endl;
	file << "STAGE 2: " << std::endl;
	file << "# iterations: " << this->num_iterations_stage2_ << std::endl;
	file << "# perturbations: " << this->num_perturbations_stage2_ << std::endl;
	file << "# restarts: " << this->num_restarts_stage2_ << std::endl;
	file << "time(s): " << this->time_stage2_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}

LBHeuristicSolution::LBHeuristicSolution(int dimension, int dimension2, int num_routes) : HeuristicSolution(dimension, dimension2, num_routes)
{
	this->num_iterations_ = 0;
	this->time_ = 0.0;
}

LBHeuristicSolution::LBHeuristicSolution() : HeuristicSolution()
{
	this->num_iterations_ = 0;
	this->time_ = 0.0;
}

LBHeuristicSolution::~LBHeuristicSolution()
{
}

void LBHeuristicSolution::Reset(int dimension, int dimension2, int num_routes)
{
	HeuristicSolution::Reset(dimension, dimension2, num_routes);

	this->num_iterations_ = 0;
	this->time_ = 0.0;
}

std::string LBHeuristicSolution::GenerateFileName()
{
	/*
		std::string file_name;
		if(K_STOP) file_name += "stop_";
		file_name += "fp";
		if(K_FEASIBILITY_PUMP_ADD_CUTS) file_name += "_cuts";
		if(K_SOLVE_STAGE_1) file_name += "_alpha1_" + std::to_string(int(100*initial_alpha_stage1));
		file_name += ("_alpha2_" + std::to_string(int(100*initial_alpha_stage2)) + "_p_");
		if(FIXED_FLIP_BASIS) file_name += ("fixed_" + std::to_string(K_FLIP_BASIS));
		else file_name += std::to_string(int(100*perturbation_flip_percentage));
		return file_name;*/
	return "";
}

void LBHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name)
{ /*
	 std::fstream file;
	 std::string path = ".//solutions";
	 path.append(folder);
	 //struct stat sb;
	 //if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	 path.append("s_");
	 path.append(algo);
	 path.append("_");
	 path.append(file_name);
	 //std::cout << path << std::endl;

	 file.open(path.c_str(),std::fstream::out);

	 file << std::setprecision(5) << std::fixed;

		 if(this->is_infeasible_) file << "STATUS: INFEASIBLE" << std::endl;
	 else if(this->found_x_integer_) file << "STATUS: FOUND ARC INTEGER FEASIBLE" << std::endl;
		 else if(this->found_y_integer_) file << "STATUS: FOUND VERTEX INTEGER FEASIBLE" << std::endl;
		 else file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;
		 file << "profit sum: " << this->profits_sum_ << std::endl;
	 file << "STAGE 1: " << std::endl;
	 file << "# iterations: " << this->num_iterations_stage1_ << std::endl;
	 file << "# perturbations: " << this->num_perturbations_stage1_ << std::endl;
		 file << "# restarts: " << this->num_restarts_stage1_ << std::endl;
	 file << "time(s): " << this->time_stage1_ << std::endl;
		 file << "STAGE 2: " << std::endl;
	 file << "# iterations: " << this->num_iterations_stage2_ << std::endl;
	 file << "# perturbations: " << this->num_perturbations_stage2_ << std::endl;
		 file << "# restarts: " << this->num_restarts_stage2_ << std::endl;
	 file << "time(s): " << this->time_stage2_ << std::endl;

		 file.close();
 */
	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}
