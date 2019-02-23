#pragma once
#ifndef _BOMCE_H_
#define _BOMCE_H_

#include <set>
#include <map>
#include <vector>
#include <list>
#include <memory>
#include <math.h>
#include "clusterisator.h"
#include "graphutility.h"

namespace webgr
{
	class bomce_cluster
	{
	public:
		int e_count, v_count, loops;
		unsigned id;
		long long di_sum;

		bomce_cluster(unsigned _id = 0) :id(_id), e_count(0),
			v_count(0), di_sum(0), loops(0) {}
		bomce_cluster(const bomce_cluster& cl) :id(cl.id), e_count(cl.e_count),
			v_count(cl.v_count), di_sum(cl.di_sum), loops(cl.loops) {}
		bomce_cluster(bomce_cluster&& cl) :id(cl.id), e_count(cl.e_count),
			v_count(cl.v_count), di_sum(cl.di_sum), loops(cl.loops)
		{
			cl.e_count = 0;
			cl.v_count = 0;
			cl.di_sum = 0;
			cl.loops = 0;
		}

		bool operator==(const bomce_cluster& cl) const
		{
			return (this->id == cl.id);
		}

		bool operator<(const bomce_cluster& cl) const
		{
			return (this->id < cl.id);
		}

		bomce_cluster& operator=(bomce_cluster&& cl)
		{
			e_count = cl.e_count;
			v_count = cl.v_count;
			id = cl.id;
			di_sum = cl.di_sum;
			loops = cl.loops;
			cl.e_count = 0;
			cl.v_count = 0;
			cl.di_sum = 0;
			cl.loops = 0;
			return *this;
		}

		bomce_cluster& operator=(const bomce_cluster& cl)
		{
			e_count = cl.e_count;
			v_count = cl.v_count;
			id = cl.id;
			di_sum = cl.di_sum;
			loops = cl.loops;
			return *this;
		}
	};

	struct cl_cmp
	{
		bool operator()(const std::list<bomce_cluster>::iterator& cl1,
			const std::list<bomce_cluster>::iterator& cl2) const
		{
			return (cl1->id < cl2->id);
		}
	};

	typedef std::vector<std::list<std::pair<unsigned,
		default_graph::vertex_name_type>>> bomce_tmpresult_type;

	class v_distribution
	{
	private:
		default_graph const& temp_graph;
		bomce_tmpresult_type const& temp_result;
		std::queue<unsigned> del_index;
		std::list<bomce_cluster> clusters;
		std::vector<std::list<bomce_cluster>::iterator> v_to_cl;
	public:
		v_distribution(const default_graph& gr, const bomce_tmpresult_type& tmp_result) :
			temp_graph(gr), temp_result(tmp_result) {}

		void add_to_cl(const std::list<bomce_cluster>::iterator& cl,
			unsigned gr_index, unsigned aij_sum_delt)
		{
			v_to_cl[gr_index] = cl;
			cl->v_count += temp_result[gr_index].size();
			cl->e_count += aij_sum_delt + temp_graph[gr_index].loop;
			cl->di_sum += temp_graph[gr_index].deg();
		}

		// Для тех случаев, когда у вершины нет смежных рёбер
		// Возвращает номер созданного кластера.
		void make_cluster(unsigned gr_index)
		{
			unsigned index;
			if (del_index.empty())
				index = clusters.size();
			else
			{
				index = del_index.front();
				del_index.pop();
			}
			clusters.push_front(bomce_cluster(index));
			add_to_cl(clusters.begin(), gr_index, 0);
		}
		
		void erase_from_cl(std::list<bomce_cluster>::iterator& it,
			unsigned current, unsigned aij_temp)
		{
			if ((it->v_count -= temp_result[current].size()) == 0)
			{
				del_index.push(it->id);
				clusters.erase(it);
				v_to_cl[current] = clusters.end();
				return;
			}
			it->e_count -= aij_temp + temp_graph[current].loop;
			it->di_sum -= temp_graph[current].deg();
			v_to_cl[current] = clusters.end();
		}

		//double move_vertex(unsigned);
		//double move_vertex(unsigned, const std::list<bomce_cluster>::iterator&);
			
	};

	/*double v_distribution::move_vertex(unsigned gr_index,
		const std::list<bomce_cluster>::iterator& cl)
	{
		if (v_to_cl[gr_index]->id == cl->id)
			return 0.0;
		unsigned cur_e = 0, cl_e = 0;
		temp_graph[gr_index].foreach_neighbors([&](std::pair<unsigned,
			default_graph::edge_type> pr)
		{
			if (v_to_cl[pr.first]->id == v_to_cl[gr_index]->id)
				++cur_e;
			else if (v_to_cl[pr.first]->id == cl->id)
				++cl_e;
		});
		double ret = d_modular(gr_index, cl, cl_e) -
			d_modular(v_to_cl[gr_index], cur_e, gr_index);
		erase_from_cl(v_to_cl[gr_index], gr_index, cur_e);
		add_to_cl(cl, gr_index, cl_e);
		return ret;
	}*/

	struct bomce_result
	{
		// Сруктура связей кластеров
		default_graph const& connections;
		// Параметры кластеров
		std::vector<std::list<bomce_cluster>::iterator> const& v_to_cl;
		// Состав вешин в каждом кластере
		bomce_tmpresult_type const& vertexes;
		// Список кластеров
		std::list<bomce_cluster> const& clusters;
		// число петель в исходном графе
		unsigned loop_count;

		bomce_result(default_graph const& gr, std::vector<
			std::list<bomce_cluster>::iterator> const& _v_to_cl,
			bomce_tmpresult_type const& temp_result,
			std::list<bomce_cluster> const& _clusters, unsigned loops) :
			v_to_cl(_v_to_cl), clusters(_clusters), connections(gr),
			vertexes(temp_result), loop_count(loops) {}
	};

	class bomce
	{
	private:
		unsigned edge_count, loop_count;
		std::list<bomce_cluster> clusters;
		std::vector<std::list<bomce_cluster>::iterator> v_to_cl;
		default_graph temp_graph;
		// номер кластера на список верщин в нём
		bomce_tmpresult_type temp_result;
		std::queue<unsigned> del_index;

		// Удаление вершины из кластера
		void erase_from_cl(std::list<bomce_cluster>::iterator&,
			unsigned, unsigned);

		// Добавляет вершину it в кластер с итератором cl
		// aij_sum_delt - кол. рёбер между вставляемой вершиной и
		// элементами кластера cl
		void add_to_cl(const std::list<bomce_cluster>::iterator& cl,
			unsigned gr_index, unsigned aij_sum_delt)
		{
			v_to_cl[gr_index] = cl;
			cl->v_count += temp_result[gr_index].size();
			cl->e_count += aij_sum_delt + temp_graph[gr_index].loop;
			cl->di_sum += temp_graph[gr_index].deg();
		}

		// Для тех случаев, когда у вершины нет смежных рёбер
		// Возвращает номер созданного кластера.
		void make_cluster(unsigned gr_index)
		{
			unsigned index;
			if (del_index.empty())
				index = clusters.size();
			else
			{
				index = del_index.front();
				del_index.pop();
			}
			clusters.push_front(bomce_cluster(index));
			add_to_cl(clusters.begin(), gr_index, 0);
		}

		void meta_graph();
		void _Init();
		template<class comp>
		void _Init(std::map<std::list<bomce_cluster>::iterator, unsigned, comp>&,
			const std::map<unsigned, default_graph::edge_type>&);

		// Возврат. знач. не является значением модулярность
		// Только для операций сравнения < > =
		double d_modular(const std::list<bomce_cluster>::iterator& cl,
			unsigned current_aij, unsigned gr_index)
		{
			int cur = temp_graph[gr_index].deg();
			return 2 * current_aij - (scale*(cl->di_sum - cur)*cur) / edge_count;

		}

		// Возврат. знач. не является значением модулярность
		// Только для операций сравнения < > =
		double d_modular(unsigned gr_index, const std::list<bomce_cluster>::
			iterator& it_cl, unsigned aij_temp)
		{
			return 2 * aij_temp - (scale*temp_graph[gr_index].deg()*
				it_cl->di_sum) / edge_count;
		}

		double move_vertex(unsigned);
		double move_vertex(unsigned, const std::list<bomce_cluster>::iterator&);

	public:
		double scale;

		typedef std::function<void()> distributor;
		static const distributor default_distribution;

		bomce() : edge_count(0), loop_count(0), scale(1.0) {}
		void init(default_graph&& gr)
		{
			temp_graph = std::move(gr);
			_Init();
		}

		void init(const default_graph& gr)
		{
			temp_graph = gr;
			_Init();
		}

		bool next_iteration(unsigned clusters_count = 1,
			bool need_meta_graph = true);

		bool next_iteration2(unsigned clusters_count = 1,
			bool need_meta_graph = true);

		bomce_result result() const
		{
			return bomce_result(temp_graph, v_to_cl,
				temp_result, clusters, loop_count);
		}

		virtual ~bomce() {}
	};

	double modular(const bomce_result&);

	const std::function<void()> default_distribution([]() {});

	// Удаление вершины current из кластера it
	void bomce::erase_from_cl(std::list<bomce_cluster>::iterator& it,
		unsigned current, unsigned aij_temp)
	{
		if ((it->v_count -= temp_result[current].size()) == 0)
		{
			del_index.push(it->id);
			clusters.erase(it);
			v_to_cl[current] = clusters.end();
			return;
		}
		it->e_count -= aij_temp + temp_graph[current].loop;
		it->di_sum -= temp_graph[current].deg();
		v_to_cl[current] = clusters.end();
	}

	template<class comp>
	// иннициализация cl(смежные кластеры)
	void bomce::_Init(std::map<std::list<bomce_cluster>::iterator,
		unsigned, comp>& cl, const std::map<unsigned,
		default_graph::edge_type>& graph_v)
	{
		for (auto& i : graph_v)
		{
			auto it2 = cl.find(v_to_cl[i.first]);
			if (it2 != cl.end())
				it2->second += i.second;
			else
				cl.insert({ v_to_cl[i.first], i.second });
		}
	}

	void bomce::_Init()
	{
		edge_count = temp_graph.edge_count();
		loop_count = 0;
		temp_result.resize(temp_graph.size());
		del_index = std::move(std::queue<unsigned>());
		clusters.clear();
		v_to_cl.resize(temp_graph.size());
		for (auto& i : temp_graph)
		{
			temp_result[i.first] = { std::make_pair(i.first,
				i.second.name()) };
			loop_count += i.second.loop;
			make_cluster(i.first);
		}
	}

	// Перемещает вершину it туда, куда нужно
	// Возвращает true, если перемещение произошло
	double bomce::move_vertex(unsigned gr_index)
	{
		/*
		пройтись по смежным вершинам
		получить множество смежных кластеров и смежных вершин
		определить предпочтительный для присоединения кластер
		определить предпочтительную для присоединения вершину
		выполнить либо cluster_add либо make_cluster(модулярность не должна ухудшиться),
		перед этим выполнив cluster_erase если нужно
		*/
		typedef std::map<std::list<bomce_cluster>::iterator, unsigned, cl_cmp> mp_type;
		mp_type cl; // смежный кластер на aij_temp(кол. рёбер между it и кластером)
		_Init(cl, temp_graph[gr_index].input);
		_Init(cl, temp_graph[gr_index].output);
		unsigned current_aij = 0;
		auto it = cl.find(v_to_cl[gr_index]);
		if (it != cl.end())
		{
			current_aij = it->second;
			cl.erase(v_to_cl[gr_index]);
		}
		if (cl.empty())
			return 0.0;
		std::pair<double, unsigned> max;
		// max - наилучшее внекластерное
		std::list<bomce_cluster>::iterator cl_it;
		mp_type::iterator it1 = cl.begin();
		max = { d_modular(gr_index, (cl_it = it1->first),
			it1->second), it1->second };
		for (++it1; it1 != cl.end(); ++it1)
		{
			std::list<bomce_cluster>::iterator it_tmp;
			std::pair<double, unsigned> tmp = { d_modular(gr_index,
				(it_tmp = it1->first), it1->second), it1->second };
			if (tmp.first > max.first)
			{
				max = tmp;
				cl_it = it_tmp;
			}
		}
		double p = d_modular(v_to_cl[gr_index], current_aij, gr_index);
		if (p < max.first && max.first > 0)
		{
			this->erase_from_cl(v_to_cl[gr_index], gr_index, current_aij);
			this->add_to_cl(cl_it, gr_index, max.second);
			return max.first - p;
		}
		else if (p < 0)
		{
			this->erase_from_cl(v_to_cl[gr_index], gr_index, current_aij);
			make_cluster(gr_index);
			return -p;
		}
		return 0.0;
	}

	double bomce::move_vertex(unsigned gr_index,
		const std::list<bomce_cluster>::iterator& cl)
	{
		if (v_to_cl[gr_index]->id == cl->id)
			return 0.0;
		unsigned cur_e = 0, cl_e = 0;
		temp_graph[gr_index].foreach_neighbors([&](std::pair<unsigned,
			default_graph::edge_type> pr)
		{
			if (v_to_cl[pr.first]->id == v_to_cl[gr_index]->id)
				++cur_e;
			else if (v_to_cl[pr.first]->id == cl->id)
				++cl_e;
		});
		double ret = d_modular(gr_index, cl, cl_e) -
			d_modular(v_to_cl[gr_index], cur_e, gr_index);
		erase_from_cl(v_to_cl[gr_index], gr_index, cur_e);
		add_to_cl(cl, gr_index, cl_e);
		return ret;
	}

	// Создаёт метаграф
	void bomce::meta_graph()
	{
		default_graph meta_gr;
		unsigned t = 0;
		for (auto& i : clusters)
		{
			i.id = t;
			meta_gr.new_vertex(t++);
		}
		std::vector<std::list<std::pair<unsigned, default_graph::
			vertex_name_type>>> tmp_result(clusters.size());
		for (auto& i : temp_graph)
		{
			t = v_to_cl[i.first]->id;
			tmp_result[t].splice(tmp_result[t].end(), temp_result[i.first]);
			meta_gr.connect_by_index(t, t, i.second.loop);
			for (auto& j : i.second.output)
				meta_gr.connect_by_index(t, v_to_cl[j.first]->id, j.second);
		}
		temp_result = std::move(tmp_result);
		temp_graph = std::move(meta_gr);
		v_to_cl.resize(clusters.size());
		t = 0;
		for (auto it = clusters.begin(); it != clusters.end(); ++it)
			v_to_cl[t++] = it;
		del_index = std::move(std::queue<unsigned>());
	}

	bool bomce::next_iteration(unsigned clusters_count,
		bool need_meta_graph)
	{
		if (temp_graph.size() > clusters.size() && need_meta_graph)
			meta_graph();
		double contin;
		std::vector<unsigned> prep(temp_graph.size());
		for (unsigned i = 0; i < temp_graph.size(); ++i)
			prep[i] = i;
		std::mt19937 gen(time(0));
		const double pr = 1.2;
		for (scale = 0.7; scale <= 1.3; scale += 0.2) // 0.222734
		{
			do
			{
				contin = 0;
				std::shuffle(prep.begin(), prep.end(), gen);
				for (auto& i : prep)
					contin += move_vertex(i);
			} while (contin > 1e-5);
		}
		return (temp_graph.size() > clusters.size());
	}

	bool bomce::next_iteration2(unsigned clusters_count,
		bool need_meta_graph)
	{
		if (temp_graph.size() > clusters.size() && need_meta_graph)
			meta_graph();
		double contin;
		std::vector<unsigned> prep(temp_graph.size());
		for (unsigned i = 0; i < temp_graph.size(); ++i)
			prep[i] = i;
		scale = 1;
		std::mt19937 gen(time(0));
		do
		{
			contin = 0;
			std::shuffle(prep.begin(), prep.end(), gen);
			for (auto& i : prep)
				contin += move_vertex(i);

		} while (contin > 0);
		return (temp_graph.size() > clusters.size());
	}

	double modular(const bomce_result& res)
	{
		double ret1 = 0, ret2 = 0;
		for (auto& i : res.clusters)
		{
			ret1 += i.e_count;
			ret2 += pow(i.di_sum, 2);
		}
		double t = 1.0 / (2.0 * res.connections.edge_count());
		return t*(2.0*ret1 - res.loop_count - t*ret2);
	}

	double modular2(const bomce_result& res, const default_graph& gr)
	{
		unsigned s = 0;
		for (auto& i : res.vertexes)
			s += i.size();
		std::vector<unsigned> vcl(s);
		for (unsigned i = 0; i < res.vertexes.size(); ++i)
		{
			for (auto& k : res.vertexes[i])
				vcl[k.first] = i;
		}
		double e1 = 0, d1 = 0;
		for (auto& i : gr)
		{
			for (auto& k : gr)
			{
				if (vcl[i.first] == vcl[k.first])
				{
					e1 += gr.count_by_index(i.first, k.first);
					d1 += i.second.deg() * k.second.deg();
				}
			}
		}
		return (2.0*e1 - d1 / (2.0*gr.edge_count())) / (2.0*gr.edge_count());
	}

	double ratio_cut(const bomce_result& res)
	{
		double ret = 0;
		for (auto& i : res.clusters)
			ret += (i.di_sum - 2.0*i.e_count) / res.vertexes[i.id].size();
		return ret;
	}

	double normalized_cut(const bomce_result& res)
	{
		double ret = 0;
		for (auto& i : res.clusters)
			ret += 1.0 - (2.0*i.e_count) / i.di_sum;
		return ret;
	}

	std::vector<std::set<default_graph::vertex_name_type>>
		clusters_from(const bomce_result& res)
	{
		std::vector<std::set<default_graph::vertex_name_type>> ret(res.clusters.size());
		for (unsigned i = 0; i < res.vertexes.size(); ++i)
		{
			for (auto& el : res.vertexes[i])
				ret[i].insert(el.second);
		}
		return ret;
	}
}
#endif // _bomce_H_
