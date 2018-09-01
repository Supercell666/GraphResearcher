#pragma once
#ifndef _CLUSTERISATOR_H_
#define _CLUSTERISATOR_H_

#include "graph.h"
namespace webgr {
	struct cluster
	{
		int e_count, v_count;
		unsigned id;
		long long di_sum;

		cluster(unsigned _id) :e_count(0), v_count(0),
			di_sum(0), id(_id) {}
		bool operator==(const cluster& cl) const
		{
			return (this->id == cl.id);
		}

		bool operator<(const cluster& cl) const
		{
			return (this->id < cl.id);
		}

		virtual ~cluster(){}
	};

	class clusterisator
	{
	protected:
		std::vector<cluster*> v_to_cl;
		default_graph temp_graph;
		
		virtual void _Init() = 0;

	public:
		clusterisator() {}
		void init(const default_graph& gr)
		{
			temp_graph = gr;
			v_to_cl.resize(temp_graph.capacity());
			_Init();
		}
		void init(default_graph&& gr)
		{
			temp_graph = std::move(gr);
			v_to_cl.resize(temp_graph.capacity());
			_Init();
		}
		virtual ~clusterisator()=0;
	};
}
#endif // _CLUSTERISATOR_H_