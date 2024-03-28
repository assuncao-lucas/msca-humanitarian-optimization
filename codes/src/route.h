#pragma once

#include <iostream>
#include <list>

class Route
{
public:
	Route();
	Route(const Route &route);
	~Route();
	int sum_profits_;
	double time_;
	std::list<int> vertices_;
	friend std::ostream &operator<<(std::ostream &out, Route &route);
};
