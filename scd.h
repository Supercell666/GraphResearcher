#pragma once
#ifndef SCD_H_
#define SCD_H_

#include "clusterisator.h"
#include "graphutility.h"

namespace webgr
{
	struct scd_cluster : public cluster
	{
		int out_edges;

		scd_cluster(unsigned _id) : cluster(_id), out_edges(0) {}
		double density() const
		{
			if (v_count <= 1)
				return 1.0;
			return (2.0*e_count) / ((double)v_count*(v_count - 1));
		}

		double density2(unsigned aij_temp, const default_graph::
			iterator::value_type& it) const
		{
			if (v_count <= 2)
				return 1.0;
			return (2.0*(e_count - aij_temp - it.second.loop)) /
				((v_count - 1)*(v_count - 2));
		}
	};

	struct scd_result
	{
		// Состав вешин в каждом кластере
		std::map<unsigned, std::list<unsigned>> vertexes;
		// Параметры кластеров
		std::map<unsigned, scd_cluster>& clusters;
		std::vector<int>& v_to_cl;

		scd_result(std::map<unsigned, scd_cluster>& cl,
			std::vector<int>& _v_to_cl) :clusters(cl),
			v_to_cl(_v_to_cl) {}
		scd_result(const scd_result& res) :vertexes(res.vertexes),
			clusters(res.clusters), v_to_cl(res.v_to_cl) {}
		scd_result& operator=(const scd_result& res)
		{
			vertexes = res.vertexes;
			clusters = res.clusters;
			v_to_cl = res.v_to_cl;
			return *this;
		}
	};

	class scd
	{
	private:
		std::vector<int> v_to_cl;
		default_graph temp_graph;
		double glob_clc; // глобальный коэффициент класткризацци после препроцессинга
		std::map<unsigned, scd_cluster> clusters;
		std::vector<double> clc;

		void add_to_cl(const std::map<unsigned, scd_cluster>::iterator&,
			const default_graph::iterator::value_type&, unsigned);
		void erase_from_cl(const std::map<unsigned, scd_cluster>::iterator&,
			const default_graph::iterator::value_type&, unsigned, bool need_del = true);
		unsigned make_cluster(const default_graph::iterator::value_type&);
		void add_to_cl(const std::map<unsigned, scd_cluster>::iterator&,
			unsigned, unsigned);
		unsigned make_cluster(unsigned);

		int best_movement(const default_graph::iterator::value_type&);

		void _Init(std::map<unsigned, unsigned>&, const std::map<
			unsigned, default_graph::edge_type>&);
		void _Init_2(std::vector<bool>&, const std::map<unsigned,
			default_graph::edge_type>&, std::pair<std::map<unsigned,
			std::set<unsigned>>::iterator,
			std::map<unsigned, scd_cluster>::iterator>&);
		bool move_anywhere(const default_graph::iterator::value_type&);
		double wcc_i(const default_graph::iterator::value_type&,
			unsigned, unsigned, double, unsigned);
		double wcc_r(const default_graph::iterator::value_type& it,
			std::map<unsigned, scd_cluster>::iterator& cur_cl, unsigned current_aij)
		{
			unsigned out_e = cur_cl->second.out_edges -
				(int)it.second.deg() - 2 * (int)current_aij;
			double this_d = cur_cl->second.density2(current_aij, it);
			double ret = -wcc_i(it, out_e, cur_cl->second.v_count - 1,
				this_d, current_aij);
			return ret;
		}
		void del_outtriangl_edges();
	protected:
		virtual void _Init(bool);
	public:
		void init(const default_graph& gr, bool need_del_edges = false)
		{
			temp_graph = gr;
			_Init(need_del_edges);
		}
		void init(default_graph&& gr, bool need_del_edges = false)
		{
			temp_graph = std::move(gr);
			_Init(need_del_edges);
		}
		bool clusterizate();
		bool clusterizate2();
		scd_result result();
	};

	void scd::_Init(std::map<unsigned, unsigned>& cl, const std::map<
		unsigned, default_graph::edge_type>& graph_v)
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

	void scd::add_to_cl(const std::map<unsigned, scd_cluster>::iterator& cl,
		const default_graph::iterator::value_type& it, unsigned aij_sum_delt)
	{
		v_to_cl[it.first] = cl->first;
		++(cl->second.v_count);
		cl->second.e_count += (aij_sum_delt + it.second.loop);
		cl->second.out_edges += (it.second.deg() - 2 * aij_sum_delt);
	}

	void scd::erase_from_cl(const std::map<unsigned, scd_cluster>::
		iterator& cl, const default_graph::iterator::value_type& it,
		unsigned aij_temp, bool need_del)
	{
		auto& clsecond = cl->second;
		if (--(clsecond.v_count) > 0)
		{
			clsecond.e_count -= (aij_temp + it.second.loop);
			clsecond.out_edges -= it.second.deg();
			clsecond.out_edges += 2 * aij_temp;
		}
		else if (need_del)
			clusters.erase(cl);
		v_to_cl[it.first] = -1;
	}

	unsigned scd::make_cluster(const default_graph::iterator::value_type& it)
	{
		unsigned ret = (clusters.empty()) ? 0 : (clusters.rbegin()->first + 1);
		auto cl = clusters.insert({ ret, scd_cluster(0) }).first;
		add_to_cl(cl, it, 0);
		return ret;
	}

	void scd::add_to_cl(const std::map<unsigned, scd_cluster>::iterator& cl,
		unsigned v_index, unsigned aij_sum_delt)
	{
		v_to_cl[v_index] = cl->first;
		auto& clsecond = cl->second;
		++clsecond.v_count;
		clsecond.e_count += (aij_sum_delt + temp_graph[v_index].loop);
		clsecond.out_edges += temp_graph[v_index].deg();
		clsecond.out_edges -= 2 * aij_sum_delt;
	}

	unsigned scd::make_cluster(unsigned v_index)
	{
		unsigned ret = (clusters.empty()) ? 0 : (clusters.rbegin()->first + 1);
		auto cl = clusters.insert({ ret, scd_cluster(0) }).first;
		add_to_cl(cl, v_index, temp_graph[v_index].loop);
		return ret;
	}

	void scd::_Init_2(std::vector<bool>& used, const std::map<unsigned,
		default_graph::edge_type>& neighbors, std::pair<std::map<unsigned,
		std::set<unsigned>>::iterator, std::map<unsigned,
		scd_cluster>::iterator>& current_c)
	{
		for (auto& v : neighbors)
		{
			if (used[v.first])
				continue;
			used[v.first] = true;
			if (v_to_cl[v.first] == current_c.second->first)
				continue;
			unsigned aij_temp = 0;
			for (auto& k : current_c.first->second)
			{
				aij_temp += temp_graph.count_by_index(v.first, k) +
					temp_graph.count_by_index(k, v.first);
			}
			add_to_cl(current_c.second, v.first, aij_temp);
			current_c.first->second.insert(v.first);
		}
	}

	void scd::_Init(bool need_del_edges)
	{
		glob_clc = 0;
		v_to_cl.resize(temp_graph.capacity(), -1);
		for (auto& i : v_to_cl)
			i = -1;
		clusters.clear();
		clc.resize(temp_graph.capacity());
		std::vector<unsigned> prep(temp_graph.capacity());
		parallel_for_each(temp_graph.begin(), temp_graph.end(), [&prep, this](
			const default_graph::iterator::value_type& it, unsigned t_id)
		{
			clc[it.first] = cl_coef(this->temp_graph, it);
			prep[it.first] = it.first;
		});
		for (auto& t : clc)
			glob_clc += t;
		glob_clc /= temp_graph.size();
		std::function<bool(unsigned, unsigned)> comp([this](
			unsigned first, unsigned second)
		{
			if (abs(clc[first] - clc[second]) < 1e-5)
				return (this->temp_graph[first].deg() > this->temp_graph[second].deg());
			else
				return (clc[first] > clc[second]);
		});
		std::sort(prep.begin(), prep.end(), comp);
		std::vector<bool> used(temp_graph.size(), false);
		std::map<unsigned, std::set<unsigned>> mp; // содержимое кластеров, нужно только при препроцессинге
		for (auto& i : prep)
		{
			if (used[i])
				continue;
			used[i] = true;
			unsigned cl_index = make_cluster(i);
			mp.insert({ cl_index, { i } });
			auto current_c = std::make_pair(mp.find(cl_index), clusters.find(cl_index));
			_Init_2(used, temp_graph[i].output, current_c);
			_Init_2(used, temp_graph[i].input, current_c);
		}
		if (need_del_edges)
		{
			del_outtriangl_edges();
			glob_clc = 0;
			parallel_for_each(temp_graph.begin(), temp_graph.end(), [&prep, this](
				const default_graph::iterator::value_type& it, unsigned t_id)
			{
				clc[it.first] = cl_coef(this->temp_graph, it);
			});
			for (auto& t : clc)
				glob_clc += t;
			glob_clc /= temp_graph.size();
		}

		//glob_clc = webgr::cl_coef(temp_graph);
	}

	void scd::del_outtriangl_edges()
	{
		std::vector<bool> used(this->temp_graph.capacity(), false);
		for (auto& i : temp_graph)
		{
			std::vector<unsigned> t;
			merging<unsigned, std::function<void(const std::pair<unsigned,
				unsigned>&)>>(i.second.input.upper_bound(i.first), i.second.input.end(),
					i.second.output.upper_bound(i.first), i.second.output.end(), [&t, &used](
						const std::pair<unsigned, unsigned>& p)
			{
				t.push_back(p.first);
			});
			used[i.first] = true;
			auto i_it = clusters.find(v_to_cl[i.first]);
			for (auto& k : t) // соседи вершины i с индексом > i
			{
				merging<unsigned, std::function<void(const std::pair<unsigned,
					unsigned>&)>>(temp_graph[k].input.upper_bound(k), temp_graph[k].input.end(),
						temp_graph[k].output.upper_bound(k), temp_graph[k].output.end(), [
							this, &t, &i, &i_it, &k](const std::pair<unsigned, unsigned>& p)
				{
					if (!std::binary_search(t.begin(), t.end(), p.first))
					{
						unsigned sum = this->temp_graph.
							disconnect_by_index(i.first, k) +
							this->temp_graph.disconnect_by_index(
								k, i.first);
						/* = this->temp_graph.disconnect_by_index(
						p.first, i.first) + this->temp_graph.
						disconnect_by_index(i.first, k) +
						this->temp_graph.disconnect_by_index(
						k, i.first) + this->temp_graph.
						disconnect_by_index(i.first, p.first);*/

						if (sum > 0)
						{
							if (v_to_cl[p.first] == v_to_cl[i.first])
							{
								i_it->second.e_count -= sum;
							}
							else
							{
								auto p_it = clusters.find(v_to_cl[p.first]);
								i_it->second.out_edges -= sum;
								p_it->second.out_edges -= sum;
							}
						}
					}
				});
			}
		}
	}

	double scd::wcc_i(const default_graph::iterator::value_type& it,
		unsigned out_edges, unsigned v_count, double density, unsigned e_sum)
	{
		if (v_count < 2)
		{
			return 0;
		}
		double q = ((double)out_edges - e_sum) / v_count;
		double qq_1w = q*(q - 1)*glob_clc;
		double d_out = it.second.deg() - e_sum;
		double r_1r_2s3 = (v_count - 1.0)*(v_count - 2.0)*pow(density, 3);
		double t1 = ((v_count - 1)*density + 1 + q) /
			((v_count + q)*(r_1r_2s3 + (e_sum - 1)*density + qq_1w*
			(density + 1) + d_out*glob_clc));
		double t2 = -(r_1r_2s3*((v_count - 1)*density + q)) /
			((r_1r_2s3 + qq_1w + q*(v_count - 1)*density*glob_clc)*
			((v_count + q)*(v_count + q - 1)));
		double didi_1d = v_count*(v_count - 1.0)*density;
		double t3 = (didi_1d*it.second.deg()) / ((didi_1d + d_out*(d_out - 1)*
			glob_clc + d_out*e_sum*glob_clc)*(v_count + d_out));
		double ret = (e_sum*t1 + (v_count - e_sum)*t2 + t3) / temp_graph.size();
		return ret;
	}

	bool scd::move_anywhere(const default_graph::iterator::value_type& it)
	{
		std::map<unsigned, unsigned> cl;
		_Init(cl, it.second.input);
		_Init(cl, it.second.output);
		// cl - смежные кластеры(по соседним вершинам)
		unsigned index = v_to_cl[it.first];
		auto it1 = cl.find(index);
		auto cur_cl = clusters.find(index);
		unsigned current_aij = (it1 != cl.end()) ? it1->second : 0;
		double wcc_this = wcc_r(it, cur_cl, current_aij);
		if (it1 != cl.end())
			cl.erase(it1);
		if (cl.empty())
			return false;

		it1 = cl.begin();
		auto dest = clusters.find(it1->first);
		double d_it = dest->second.density();
		double wcc_max = wcc_i(it, dest->second.out_edges, dest->second.v_count,
			d_it, it1->second);
		if (cur_cl->second.v_count > 1)
			wcc_max += wcc_this;
		unsigned aij_temp = it1->second;

		for (++it1; it1 != cl.end(); ++it1)
		{
			auto t_it = clusters.find(it1->first);
			d_it = t_it->second.density();
			double wcc_t = wcc_i(it, t_it->second.out_edges, t_it->second.v_count,
				d_it, it1->second);
			if (cur_cl->second.v_count > 1)
				wcc_t += wcc_this;
			if (wcc_t > wcc_max)
			{
				wcc_max = wcc_t;
				dest = t_it;
				aij_temp = it1->second;
			}
		}
		if (wcc_max < -wcc_this)
			return false;
		if (wcc_max < wcc_this && wcc_this > 0)
		{
			//std::cout << cur_cl->second.v_count << std::ends << cur_cl->second.e_count << std::endl;
			erase_from_cl(cur_cl, it, current_aij);
			make_cluster(it);
			return true;
		}
		else if (wcc_max > 0)
		{
			//std::cout << cur_cl->second.v_count << std::ends << cur_cl->second.e_count << std::endl;
			erase_from_cl(cur_cl, it, current_aij);
			add_to_cl(dest, it, aij_temp);
			return true;
		}
		return false;
	}

	int scd::best_movement(const default_graph::iterator::value_type& it)
	{
		std::map<unsigned, unsigned> cl;
		_Init(cl, it.second.input);
		_Init(cl, it.second.output);
		// cl - смежные кластеры(по соседним вершинам)
		unsigned index = v_to_cl[it.first];
		auto it1 = cl.find(index);
		auto cur_cl = clusters.find(index);
		unsigned current_aij = (it1 != cl.end()) ? it1->second : 0;
		double wcc_this = wcc_r(it, cur_cl, current_aij);
		if (it1 != cl.end())
			cl.erase(it1);
		if (cl.empty())
			return -2;

		it1 = cl.begin();
		auto dest = clusters.find(it1->first);
		double d_it = dest->second.density();
		double wcc_max = wcc_i(it, dest->second.out_edges, dest->second.v_count,
			d_it, it1->second);
		if (cur_cl->second.v_count > 1)
			wcc_max += wcc_this;
		unsigned aij_temp = it1->second;

		for (++it1; it1 != cl.end(); ++it1)
		{
			auto t_it = clusters.find(it1->first);
			d_it = t_it->second.density();
			double wcc_t = wcc_i(it, t_it->second.out_edges, t_it->second.v_count,
				d_it, it1->second);
			if (cur_cl->second.v_count > 1)
				wcc_t += wcc_this;
			if (wcc_t > wcc_max)
			{
				wcc_max = wcc_t;
				dest = t_it;
				aij_temp = it1->second;
			}
		}
		if (wcc_max < -wcc_this)
			return -2;
		if (wcc_max < wcc_this && wcc_this > 0)
		{
			return -1;
		}
		else if (wcc_max > 0)
		{
			return dest->first;
		}
		return -2;
	}

	bool scd::clusterizate()
	{
		unsigned contin, g = 0;
		do
		{
			contin = 0;
			for (auto it = temp_graph.begin(); it != temp_graph.end(); ++it)
			{
				if (move_anywhere(*it))
					++contin;
			}
			//std::cout << '\t' << contin << std::endl;
			if (++g >= 2 * temp_graph.size())
			{
				return false;
			}
		} while (contin > 0);
		return true;
	}

	bool scd::clusterizate2()
	{
		unsigned contin;
		unsigned mv = 0;
		do
		{
			contin = 0;
			std::vector<int> moves(temp_graph.capacity());
			parallel_for_each(temp_graph.begin(), temp_graph.end(), [&moves, this](
				const default_graph::iterator::value_type& val, unsigned tid)
			{
				moves[val.first] = this->best_movement(val);
			});
			for (auto it = temp_graph.begin(); it != temp_graph.end(); ++it)
			{
				if (v_to_cl[it->first] == moves[it->first] || moves[it->first] == -2)
					continue;

				unsigned this_ = 0, dest = 0;
				merging<unsigned, std::function<void(const std::pair<unsigned,
					unsigned>&)>>(it->second.input.begin(), it->second.input.end(),
						it->second.output.begin(), it->second.output.end(),
						[this, &moves, &it, &dest, &this_](const std::pair<unsigned,
							unsigned>& p)
				{
					if (v_to_cl[p.first] == v_to_cl[it->first])
						++this_;
					else if (v_to_cl[p.first] == moves[it->first])
						++dest;
				});

				erase_from_cl(clusters.find(v_to_cl[it->first]), *it, this_, false);
				if (moves[it->first] == -1)
					make_cluster(*it);
				else
					add_to_cl(clusters.find(moves[it->first]), *it, dest);
				++contin;
			}
			//std::cout << '\t' << contin << std::endl;
			/*for (auto& k : v_to_cl)
				std::cout << k << std::ends;
			std::cout << std::endl;*/
			if (contin <= 0)
				return true;
			std::vector<unsigned> del;
			for (auto& it : clusters)
			{
				if (it.second.v_count == 0)
					del.push_back(it.first);
			}
			for (auto& i : del)
				clusters.erase(i);
		} while (++mv < temp_graph.size() / 4);
		return (mv < temp_graph.size() / 4);
	}

	scd_result scd::result()
	{
		scd_result ret(this->clusters, this->v_to_cl);
		for (unsigned i = 0; i < v_to_cl.size(); ++i)
		{
			auto it = ret.vertexes.find(v_to_cl[i]);
			if (it != ret.vertexes.end())
				it->second.push_back(i);
			else
				ret.vertexes.insert({ (unsigned)v_to_cl[i], {i} });
		}
		return ret;
	}

	double modular(const scd_result& res)
	{
		unsigned d_sum = 0, ret1 = 0;
		double ret2 = 0;
		for (auto& i : res.clusters)
		{
			ret1 += i.second.e_count;
			ret2 += pow(2 * i.second.e_count + i.second.out_edges, 2);
			d_sum += i.second.out_edges;
		}
		unsigned e_count = d_sum / 2 + ret1;
		double t = 1.0 / (2.0*e_count);
		return t*(2.0*ret1 - t*ret2);
	}
}
#endif // SCD_H_
