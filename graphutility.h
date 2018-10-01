#pragma once
#ifndef _GRAPHUTILITY_H_
#define _GRAPHUTILITY_H_

#include <fstream>
#include <random>
#include <ctime>
#include <algorithm>
#include <queue>
#include <list>
#include <thread>
#include <mutex>
#include <strchange.h>
#include <unordered_map>
#include "graph.h"
#include "bomce.h"

namespace webgr
{
	/*template<class type>
	// костыль для работы с vector<atomic<...>>
	// разрешает копирование
	class crutch_atomic : public std::atomic<type>
	{
	public:
		crutch_atomic() : std::atomic<type>() {}
		crutch_atomic(const crutch_atomic& val) :std::atomic<type>()
		{
			this->store(val.load());
		}
		crutch_atomic(const type& val) :std::atomic<type>(val) {}
		operator type() const
		{
			return this->load();
		}
	};*/

	// Модель Боллобаша-Риордана
	default_graph preferred_attachment(unsigned, unsigned);

	class top_sort
	{
	private:
		mutable std::vector<unsigned> ans;
		mutable std::vector<bool> used;
		const default_graph* gr;
		void dfs(unsigned);
	public:
		const std::vector<unsigned>& sort(const default_graph& grph)
		{
			gr = &grph;
			unsigned size = std::min((unsigned)used.size(), (unsigned)grph.size());
			used.resize(grph.size(), false);
			for (unsigned i = 0; i < size; ++i)
				used[i] = false;
			ans.clear();
			for (auto& i : grph)
			{
				if (!used[i.first])
					dfs(i.first);
			}
			std::reverse(ans.begin(), ans.end());
			return ans;
		}
	};

	void top_sort::dfs(unsigned v)
	{
		used[v] = true;
		for (auto& it : (*gr)[v].output)
		{
			if (it.first != v && !used[it.first])
				dfs(it.first);
		}
		ans.push_back(v);
	}

	struct cl_incl
	{
		std::list<unsigned> g;
		unsigned aij_sum, di_sum; // сумма рёбер и сумма степеней вершин

		cl_incl() :aij_sum(0), di_sum(0) {}
	};

	struct result_graph_generator
	{
		std::vector<cl_incl> clusters; // кластер на номер вершины
		std::vector<unsigned> vertexes; // вершина на кластер
		default_graph gr; // сам граф
		unsigned loop_count;

		result_graph_generator() : loop_count(0) {}
	};

	template<class _Iterator, class _Function>
	// Экономит память
	void parallel_for_each(const _Iterator& begin, const _Iterator& end,
		const _Function& func, unsigned thread_count =
		std::thread::hardware_concurrency())
	{
		// Требует t*logN памяти
		// Общая сложность 2*t*N (без учёта выполнения func)
		std::thread** threads = new std::thread*[thread_count];
		for (unsigned id = 0; id < thread_count; ++id)
		{
			threads[id] = new std::thread([id, &thread_count, &begin, &end, &func]() {
				_Iterator it = begin;
				for (unsigned char i = 0; i < id && it != end; ++i)
					++it;
				while (it != end)
				{
					func(*it, id);
					for (unsigned char i = 0; i < thread_count && it != end; ++i)
						++it;
				}
			});
		}
		for (unsigned i = 0; i < thread_count; ++i)
		{
			threads[i]->join();
			delete threads[i];
		}
		delete[] threads;
	}

	template<class type, class _Function>
	// Экономит память
	void parallel_for_each(const std::iterator<std::random_access_iterator_tag,
		type>& begin, const std::iterator<std::random_access_iterator_tag,
		type>& end, const _Function& func, unsigned thread_count =
		std::thread::hardware_concurrency() * 2)
	{
		// Требует t*logN памяти
		// Общая сложность 2*t*N (без учёта выполнения func)
		std::thread** threads = new std::thread*[thread_count];
		for (unsigned id = 0; id < thread_count; ++id)
		{
			threads[id] = new std::thread([id, &thread_count, &begin, &end, &func]() {
				auto it = begin;
				it = it + id;
				while (it != end)
				{
					func(*it, id);
					it = it + thread_count;
				}
			});
		}
		for (unsigned i = 0; i < thread_count; ++i)
		{
			threads[i]->join();
			delete threads[i];
		}
		delete[] threads;
	}

	template<class type, class _Function>
	void parallel_for(const type& begin, const type& end, const type& step,
		const _Function& func, unsigned thread_count =
		std::thread::hardware_concurrency() * 2)
	{
		if (step == type(0))
			throw("wrong step parameter");
		std::thread** threads = new std::thread*[thread_count];
		for (unsigned id = 0; id < thread_count; ++id)
		{
			threads[id] = new std::thread([id, &step, &thread_count,
				&begin, &end, &func, &threads]()
			{
				for (auto it = begin + step*id; it < end;
					it += step*thread_count)
				{
					func(it, id);
				}
			});
		}
		for (unsigned i = 0; i < thread_count; ++i)
		{
			threads[i]->join();
			delete threads[i];
		}
		delete[] threads;
	}

	// Модель Боллобаша-Риордана
	default_graph preferred_attachment(unsigned n, unsigned m)
	{
		/* n вершин, n*m рёбер
		граф строится до достижения n*m рёбер и n*m вершин
		затем вершины объединяются друг с другом (каждые m штук в одну), пока их не станет ровно n
		*/
		n *= m;
		std::vector<unsigned> vec; // хранит номера вершин
		// при вставке ребра между first и second вершинами
		// в vec помещаются и first и second

		default_graph gr; unsigned i;
		std::mt19937 gen(std::time(0));

		for (i = 0; i < n; ++i)
		{
			gr.new_vertex(i);
			vec.push_back(i);
			// Кидать не 1 исходящую вершину, а несколько
			std::uniform_int_distribution<unsigned> dist(0, vec.size() - 1);
			unsigned temp = vec[dist(gen)];
			gr.connect_by_index(i, temp);
			vec.push_back(temp);
		}
		if (m < 2)
			return gr;
		// объединение вершин
		for (i = 0; i < n; i += m)
		{
			for (unsigned k = i + m - 1; k > i;)
			{
				for (unsigned h = i; h < k; ++h, --k)
					gr.merge_by_index(h, k);
			}
		}
		return gr;
	}

	// Генерирует кластеризованый граф
	// v_count - количество вершин, min_cl_size - минимальный размер кластера,
	// max_cl_size- максимальный размер кластера, p_in_cl - вероятность наличия ребра в кластере,
	// p_out_cl - вероятность ребра между кластерами
	result_graph_generator clustered_graph(unsigned v_count, unsigned min_cl_size,
		unsigned max_cl_size, double p_in_cl = 0.5, double p_out_cl = 0.3)
	{
		if (min_cl_size > max_cl_size || max_cl_size > v_count ||
			min_cl_size <= 1)
		{
			throw(std::string("Wrong cl_size params"));
		}
		if (p_in_cl >= 1 || p_in_cl <= p_out_cl)
		{
			throw(std::string("Wrong p_cl params"));
		}

		std::uniform_int_distribution<int> dist(min_cl_size, max_cl_size);
		std::uniform_int_distribution<int> in_out(0, 1);
		std::uniform_real_distribution<double> q_def(0, 1);

		std::mt19937 gen(time(0));
		std::vector<unsigned> size;
		unsigned count;
		for (count = 0;;)
		{
			unsigned temp = dist(gen);
			count += temp;
			size.push_back(temp);

			if (count >= v_count - max_cl_size && count <= v_count -
				min_cl_size)
			{
				break;
			}
			if (count > v_count - min_cl_size)
				break;
		}
		size.push_back(v_count - count);
		while (size.back() < min_cl_size)
		{
			bool err = true;
			for (unsigned i = 0; i < size.size() - 1 && size.back() <
				min_cl_size; ++i)
			{
				if (size[i] > min_cl_size)
				{
					++size.back();
					--size[i];
					err = false;
				}
			}
			if (size.back() < min_cl_size && err)
				throw(std::string("Wrong cl_size params"));
		}
		// генерим кластеры размерностью size[i] каждый.
		result_graph_generator ret;
		ret.clusters.resize(size.size());
		ret.vertexes.resize(v_count);
		for (unsigned i = 0; i < v_count; ++i)
			ret.gr.new_vertex(i);
		// имена от 0 до v_count; индексы от 0 до v_count
		std::vector<std::pair<unsigned, unsigned>> clusters;
		for (unsigned b_begin = 0, cl = 0; cl < size.size(); ++cl)
		{
			unsigned b_end = b_begin + size[cl];
			clusters.push_back({ b_begin, b_end });
			for (unsigned i = b_begin; i < b_end; ++i)
			{
				ret.clusters[cl].g.push_back(i);
				ret.vertexes[i] = cl;
				for (unsigned j = i + 1; j < b_end; ++j)
				{
					if (q_def(gen) <= p_in_cl)
					{
						++ret.clusters[cl].aij_sum;
						ret.clusters[cl].di_sum += 2;
						if (in_out(gen) == 0)
							ret.gr.connect_by_index(j, i, 1);
						else
							ret.gr.connect_by_index(i, j, 1);
					}
				}
			}
			auto end = ret.gr.find(b_end);
			for (auto it = ret.gr.find(b_begin); it != end; ++it)
			{
				if (ret.gr[it->first].deg() > 0)
					continue;
				std::uniform_int_distribution<int> rand_v(b_begin, b_end - 1);
				unsigned r_v;
				do {
					r_v = rand_v(gen);
				} while (r_v == it->first);

				++ret.clusters[cl].aij_sum;
				ret.clusters[cl].di_sum += 2;
				auto it1 = ret.gr.find(r_v);
				if (in_out(gen) == 0)
					ret.gr.connect(it1, it, 1);
				else
					ret.gr.connect(it, it1, 1);
			}
			b_begin = b_end;
		}
		size.clear();
		// Зашумление
		for (unsigned i = 0; i < clusters.size(); ++i)
		{
			for (unsigned k = i + 1; k < clusters.size(); ++k)
			{
				for (unsigned it_f = clusters[i].first;
					it_f < clusters[i].second; ++it_f)
				{
					for (unsigned it_s = clusters[k].first;
						it_s < clusters[k].second; ++it_s)
					{
						if (q_def(gen) <= p_out_cl)
						{
							++ret.clusters[i].di_sum;
							++ret.clusters[k].di_sum;
							if (in_out(gen) == 0)
								ret.gr.connect_by_index(it_s, it_f, 1);
							else
								ret.gr.connect_by_index(it_f, it_s, 1);
						}
					}
				}
			}
		}
		return ret;
	}

	template<class _Iterator, class _Function>
	// Слияние 2х отсортированных контейнеров в 1 отсортированный
	void merging(const _Iterator& begin1, const _Iterator& end1,
		const _Iterator& begin2, const _Iterator& end2, const _Function& func)
	{
		auto it1 = begin1, it2 = begin2;
		while (it1 != end1 && it2 != end2)
		{
			if (it1->first < it2->first)
			{
				func(*it1);
				++it1;
			}
			else if (it1->first > it2->first)
			{
				func(*it2);
				++it2;
			}
			else
			{
				func({ it1->first, it1->second + it2->second });
				++it1; ++it2;
			}
		}
		if (it1 != end1)
		{
			do
			{
				func(*it1);
				++it1;
			} while (it1 != end1);
		}
		else if (it2 != end2)
		{
			do
			{
				func(*it2);
				++it2;
			} while (it2 != end2);
		}
	}

	template<typename type, class _Func, class compare>
	void foreach_equal(const std::iterator<std::random_access_iterator_tag, type>& begin1,
		const std::iterator<std::random_access_iterator_tag, type>& end1,
		const std::iterator<std::random_access_iterator_tag, type>& begin2,
		const std::iterator<std::random_access_iterator_tag, type>& end2,
		const _Func& fcn, const compare& cmp = std::less<type>())
	{
		auto f_begin = begin1;
		auto f_end = end1;
		auto s_begin = begin2;
		auto s_end = end2;

		if (cmp(*f_begin, *s_begin))
			f_begin = std::lower_bound(begin1, end1, *begin2, cmp);
		else if (cmp(*s_begin, *f_begin))
			s_begin = std::lower_bound(begin2, end2, *begin1, cmp);
		if (cmp(*(s_end - 1), *(f_end - 1)))
			f_end = std::upper_bound(f_begin, end1, *(s_end - 1), cmp);
		else if (cmp(*(f_end - 1), *(s_end - 1)))
			s_end = std::upper_bound(s_begin, end2, *(f_end - 1), cmp);
		// определили область возможного пересечения
		auto size1 = std::distance(f_begin, f_end);
		auto size2 = std::distance(s_begin, s_end);

		if (size1 == 0 || size2 == 0)
			return;
		auto t1 = size1*(log(size2) / log(2.0)),
			t2 = size2*(log(size1) / log(2.0));
		if (size1 + size2 <= t1 && size1 + size2 <= t2)
		{
			do {
				do {
					for (;;)
					{
						if (f_begin == f_end)
							return;
						if (cmp(*f_begin, *s_begin))
							++f_begin;
						else break;
					}
					for (;;)
					{
						if (s_begin == s_end)
							return;
						if (cmp(*s_begin, *f_begin))
							++s_begin;
						else break;
					}
				} while (cmp(*f_begin, *s_begin));
				fcn(*f_begin, *s_begin);
				++f_begin;
				++s_begin;

			} while (f_begin != f_end && s_begin != s_end);
		}
		else if (t1 < t2)
		{
			do {
				auto t = std::lower_bound(s_begin, s_end, *f_begin, cmp);
				if (!cmp(*t, *f_begin))
					fcn(*f_begin, *t);

			} while (++f_begin != f_end);
		}
		else
		{
			do {
				auto t = std::lower_bound(f_begin, f_end, *s_begin, cmp);
				if (!cmp(*t, *s_begin))
					fcn(*t, *s_begin);

			} while (++s_begin != s_end);
		}
	}

	template<class _Func>
	void foreach_equal_clusters(const std::set<default_graph::vertex_name_type>& cont1,
		const std::set<default_graph::vertex_name_type>& cont2, const _Func& fcn)
	{
		if (cont1.size() == 0 || cont2.size() == 0)
			return;
		auto f_begin = cont1.begin();
		auto f_end = cont1.end();
		auto s_begin = cont2.begin();
		auto s_end = cont2.end();

		auto cmp = cont1.key_comp();

		if (cmp(*f_begin, *s_begin))
			f_begin = cont1.lower_bound(*s_begin);
		else if (cmp(*s_begin, *f_begin))
			s_begin = cont2.lower_bound(*f_begin);

		auto f_rbegin = cont1.rbegin();
		auto s_rbegin = cont2.rbegin();

		if (cmp(*s_rbegin, *f_rbegin))
			f_end = cont1.upper_bound(*s_rbegin);
		else if (cmp(*f_rbegin, *s_rbegin))
			s_end = cont2.upper_bound(*f_rbegin);

		auto size1 = cont1.size();
		auto size2 = cont2.size();

		if (size1 == 0 || size2 == 0)
			return;
		auto t1 = size1*(log(size2) / log(2.0)),
			t2 = size2*(log(size1) / log(2.0));
		if (size1 + size2 <= t1 && size1 + size2 <= t2)
		{
			do {
				do {
					for (;;)
					{
						if (f_begin == f_end)
							return;
						if (cmp(*f_begin, *s_begin))
							++f_begin;
						else break;
					}
					for (;;)
					{
						if (s_begin == s_end)
							return;
						if (cmp(*s_begin, *f_begin))
							++s_begin;
						else break;
					}
				} while (cmp(*f_begin, *s_begin));
				fcn(*f_begin, *s_begin);
				++f_begin;
				++s_begin;

			} while (f_begin != f_end && s_begin != s_end);
		}
		else if (t1 < t2)
		{
			do {
				auto t = cont2.find(*f_begin);
				if (t != cont2.end())
					fcn(*f_begin, *t);
				++f_begin;
			} while (f_begin != cont1.end() && f_begin != f_end);
		}
		else
		{
			do {
				auto t = cont1.find(*s_begin);
				if (t != cont1.end())
					fcn(*t, *s_begin);
				++s_begin;
			} while (s_begin != cont2.end() && s_begin != s_end);
		}
	}

	template<typename edge_value_type>
	struct cmp_cl_coef
	{
		bool operator()(const std::pair<unsigned, edge_value_type>& first,
			const unsigned& second) const
		{
			return (first.first < second);
		}
		bool operator()(const unsigned& first,
			const std::pair<unsigned, edge_value_type>& second) const
		{
			return (first < second.first);
		}
		bool operator()(const unsigned& first, const unsigned& second) const
		{
			return (first < second);
		}
		bool operator()(const std::pair<unsigned, edge_value_type>& first,
			const std::pair<unsigned, edge_value_type>& second) const
		{
			return (first.first < second.first);
		}
	};

	// кластерный коэффициент для вершины графа
	double cl_coef(const default_graph& gr, const default_graph::iterator::value_type& it)
	{
		unsigned count = 0;
		std::vector<unsigned> t;
		it.second.foreach_neighbors([&t](const std::pair<unsigned, unsigned>& pr)
		{
			t.push_back(pr.first);
		});
		for (auto& i : t)
		{
			for (auto& v : gr[i].output)
			{
				if (std::binary_search(t.begin(), t.end(), v.first,
					cmp_cl_coef<default_graph::edge_type>()))
				{
					count += v.second;
				}
			}
		}
		if (it.second.deg() < 2)
			return 0;
		return 2.0*count / (it.second.deg()*(it.second.deg() - 1));
	}

	// Коэффициент сластеризации
	double cl_coef(const default_graph& gr, unsigned thread_count =
		std::thread::hardware_concurrency())
	{
		std::vector<double> v(thread_count, 0);
		parallel_for_each(gr.begin(), gr.end(), [&v, &gr](
			const default_graph::iterator::value_type& it, unsigned id)
		{
			v[id] += cl_coef(gr, it);
		}, thread_count);

		for (unsigned i = 1; i < v.size(); ++i)
			v[0] += v[i];
		return v[0] / gr.size();
	}

	std::map<default_graph::vertex_name_type, double> page_rank(
		const default_graph& gr, bool need_sort = false,
		const double& c = 0.85, const double& delta = 1e-5)
	{
		std::vector<double> rt(gr.size(), 0.5);
		std::vector<unsigned> v_s;
		if (need_sort)
		{
			top_sort srt;
			v_s = srt.sort(gr);
		}
		else
		{
			std::vector<double> clc(gr.size());
			v_s.resize(gr.size());
			for (auto& it : gr)
			{
				clc[it.first] = webgr::cl_coef(gr, it);
				v_s[it.first] = it.first;
			}
			std::function<bool(unsigned, unsigned)> comp([&clc, &gr](
				unsigned first, unsigned second)
			{
				if (abs(clc[first] - clc[second]) < 1e-5)
					return (gr[first].deg() < gr[second].deg());
				else
					return (clc[first] < clc[second]);
			});
			std::sort(v_s.begin(), v_s.end(), comp);
		}
		double dlt;
		unsigned count = 0;
		std::vector<unsigned> empt;
		for (auto& i : gr)
		{
			if (i.second.deg() - 2 * i.second.loop == 0)
				empt.push_back(i.first);
		}
		do
		{
			dlt = 0;
			for (auto& i : v_s)
			{
				double temp = 0, total = 0;
				gr[i].foreach_neighbors([&](const std::pair<unsigned,
					default_graph::edge_type>& v)
				{
					temp += rt[v.first] / (gr[v.first].deg() -
						2 * gr[v.first].loop);
				});
				/*for (auto& v : gr[i].input)
				{
					temp += rt[v.first] / (gr[v.first].out_d -
						gr[v.first].loop);
				}*/
				total += temp*c;
				temp = 0;
				for (auto& k : empt)
					temp += rt[k];
				total += c*temp / gr.size();
				total += (1 - c) / gr.size();

				temp = rt[i];
				if (abs(temp - total) > dlt)
					dlt = std::abs(temp - total);
				rt[i] = total;
			}
		} while (dlt > delta);
		std::map<default_graph::vertex_name_type, double> ret;
		dlt = rt.front();
		for (auto& i : rt)
		{
			if (i > dlt)
				dlt = i;
		}
		dlt = 1 / dlt;
		for (auto& i : gr)
			ret.insert({ i.second.name(), rt[i.first] * dlt });
		return ret;
	}

	void print(const default_graph& gr, std::ostream& os)
	{
		os << "graph gr {\n";
		for (auto& v : gr)
		{
			if (v.second.deg() == 0)
				os << v.first << '\n';
		}
		for (auto& v : gr)
		{
			for (unsigned k = 0; k < v.second.loop; ++k)
				os << v.first << " -- " << v.first << '\n';
			for (auto& i : v.second.output)
			{
				for (unsigned k = 0; k < i.second; ++k)
					os << v.first << " -- " << gr[i.first].name() << '\n';
			}
		}
		os << '}';
	}

	// Ассортативность
	double assort(const default_graph& gr)
	{
		unsigned t_count = std::thread::hardware_concurrency() * 2;
		std::vector<std::vector<double>> sum(t_count, std::vector<double>(3, 0));
		parallel_for_each(gr.begin(), gr.end(), [&sum, &gr](
			default_graph::iterator::value_type i, unsigned t_id)
		{
			for (auto& k : i.second.output)
			{
				sum[t_id][0] += i.second.deg()*gr[k.first].deg(); // 
				sum[t_id][1] += (unsigned)i.second.deg() + gr[k.first].deg();
				sum[t_id][2] += (unsigned)pow(i.second.deg(), 2) + pow(gr[k.first].deg(), 2);
			}
		});
		double sum1[3] = { 0,0,0 };
		for (auto& i : sum)
		{
			sum1[0] += i[0];
			sum1[1] += i[1];
			sum1[2] += i[2];
		}
		sum1[1] = pow(sum1[1] / 2.0, 2) / gr.edge_count();
		double num = sum1[0] - sum1[1];
		double denum = sum1[2] / 2.0 - sum1[1];
		return num / denum;
	}

	// Ассортативность для графа в виде списка рёбер
	// gr - путь к текстовому файлу списка рёбер
	double assort(const std::string& gr)
	{
		std::unordered_map<unsigned, std::pair<unsigned, unsigned>> degs;
		std::fstream fin(gr);
		unsigned in, out, l = 0;
		std::string str;
		while (getline(fin, str))
		{
			int pos = str.find(" -- ");
			if (pos == -1)
				continue;
			in = atoi(str.substr(0, pos).c_str());
			auto it = degs.find(in);
			if (it != degs.end())
				++(it->second.first);
			else
				degs.insert({ in, { 1, 0 } });
			out = atoi(str.substr(pos + 4, str.size() - pos - 4).c_str());
			it = degs.find(out);
			if (it != degs.end())
				++(it->second.second);
			else
				degs.insert({ in, { 0, 1 } });
			++l;
		}
		fin.close();

		fin.open(gr);
		double f_sum = 0, d_sum = 0, d_sqr = 0, f_sqr = 0, df_sum = 0;
		while (getline(fin, str))
		{
			int pos = str.find(" -- ");
			if (pos == -1)
				continue;
			unsigned in_it = degs[atoi(str.substr(0, pos).c_str())].first;
			unsigned out_it = degs[atoi(str.substr(pos + 4,
				str.size() - pos - 4).c_str())].second;
			d_sum += in_it;
			f_sum += out_it;
			d_sqr += in_it * in_it;
			f_sqr += out_it * out_it;
			df_sum += in_it * out_it;
		}
		double num = l * df_sum - d_sum*f_sum;
		double denum = sqrt((l*d_sqr - d_sum*d_sum)*(l*f_sqr - f_sum*f_sum));
		return num / denum;
	}

	// Диаметр графа(возвращает отображение длины пути на количество путей такой длины)
	std::vector<unsigned> diameter(const default_graph& gr, unsigned ret_size = 0)
	{
		std::mutex resize_mut;
		unsigned t_count = std::thread::hardware_concurrency() * 2;
		std::vector<std::vector<unsigned>> ret(ret_size, std::vector<unsigned>(t_count, 0));
		parallel_for_each(gr.begin(), gr.end(),
			[&](default_graph::iterator::value_type i, unsigned t_id)
		{
			std::vector<bool> used(gr.size(), false); // пройденные вершины
			used[i.first] = true;
			std::queue<unsigned> id;
			id.push(i.first);
			for (unsigned index = 0;;)
			{
				std::queue<unsigned> new_id;
				while (!id.empty())
				{
					auto& it1 = gr[id.front()];
					id.pop();
					for (auto it2 = it1.output.upper_bound(i.first);
						it2 != it1.output.end(); ++it2)
					{
						if (!used[it2->first])
						{
							new_id.push(it2->first);
							used[it2->first] = true;
						}
					}
					for (auto it2 = it1.input.upper_bound(i.first);
						it2 != it1.input.end(); ++it2)
					{
						if (!used[it2->first])
						{
							new_id.push(it2->first);
							used[it2->first] = true;
						}
					}
				}
				++index;
				if (new_id.empty())
					break;
				else
				{
					if (ret.size() <= index)
					{
						std::lock_guard<std::mutex> lc(resize_mut);
						if (ret.size() <= index)
							ret.resize(index + 1, std::vector<unsigned>(t_count, 0));
					}
					ret[index][t_id] += new_id.size();
				}
				id = std::move(new_id);
			}
			// std::cout << i.first << std::endl;
		});
		std::vector<unsigned> result(ret.size(), 0);
		for (unsigned i = 0; i < result.size(); ++i)
			for (unsigned k = 0; k < t_count; ++k)
				result[i] += ret[i][k];
		return result;
	}

	double modular(const result_graph_generator& res)
	{
		double sum1 = 0, sum2 = 0;
		double t = 1.0 / (2 * res.gr.edge_count());
		for (auto& i : res.clusters)
		{
			sum1 += i.aij_sum;
			sum2 += i.di_sum * i.di_sum;
		}
		return t*(2.0*sum1 - t*sum2);
	}

	void to_pajek(const default_graph& gr, const std::string& path)
	{
		std::ofstream fout(path);
		fout << "*Vertices " << gr.size() << std::endl;
		for (unsigned i = 0; i < gr.size(); ++i)
		{
			fout << i + 1 << ' ' << '\"' << gr[i].name() << '\"' << std::endl;
		}
		fout << "*Edges";
		for (auto& i : gr)
		{
			for (auto& v : i.second.output)
				fout << std::endl << i.first + 1 << ' ' << v.first + 1 << ' ' << v.second;
		}
		fout.close();
	}

	default_graph from_gv(const std::string& path)
	{
		std::ifstream fin(path);
		std::string str;
		default_graph gr;
		/*fin >> str;

		getline(fin, str);*/
		unsigned first, second;
		while (getline(fin, str))
		{
			if (str.find('}') != -1)
				break;
			sprclib::pars(str, " ", [&](PARS_ARGS(std::string))
			{
				if (params.index == 0)
					first = gr.new_vertex(stoi(params.str))->first;
				else if (params.index == 2)
				{
					second = gr.new_vertex(stoi(params.str))->first;
					BREAKPARS
				}
			});
			gr.connect_by_index(first, second);
		}
		return gr;
	}

	default_graph from_pajek(const std::string& path, bool need_name = true)
	{
		std::ifstream fin(path);
		std::string str; bool mod = false;
		default_graph ret;
		while (getline(fin, str))
		{
			int star_pos = str.find('*');
			if (star_pos != -1)
			{
				if (sprclib::substr_comp(str, star_pos + 1,
					star_pos + 9, "Vertices", 0, 8) == 0 ||
					sprclib::substr_comp(str, star_pos + 1,
						star_pos + 9, "vertices", 0, 8) == 0)
				{
					mod = 1;
					unsigned gr_size = stoul(str.substr(star_pos + 10,
						str.size() - star_pos - 10));
					ret._Resize(gr_size);
				}
				if (sprclib::substr_comp(str, star_pos + 1,
					star_pos + 6, "Edges", 0, 5) == 0 ||
					sprclib::substr_comp(str, star_pos + 1,
						star_pos + 6, "edges", 0, 5) == 0)
				{
					mod = 0;
				}
			}
			else if (mod)
			{
				unsigned v_index, v_name;
				sprclib::pars(str, " ", [&](PARS_ARGS(std::string))
				{
					if (params.index == 0)
						v_index = stoul(params.str) - 1;
					else if (need_name)
					{
						v_name = stol(params.str.substr(1, params.str.size() - 2));
						ret._New_Vertex(v_name, v_index);
						BREAKPARS
					}
					else
					{
						ret._New_Vertex(v_index, v_index);
						BREAKPARS
					}
				});
			}
			else
			{
				unsigned first, second;
				default_graph::edge_type weight;
				sprclib::pars(str, " ", [&first, &second, &weight](PARS_ARGS(std::string))
				{
					if (params.index == 0)
						first = stol(params.str) - 1;
					else if (params.index == 1)
						second = stol(params.str) - 1;
					else if (params.index == 2)
					{
						weight = stod(params.str);
						BREAKPARS
					}
				});
				ret.connect_by_index(first, second, weight);
			}
		}
		return ret;
	}

	// Получает разбиение на кластеры из файла membership
	std::vector<std::set<default_graph::vertex_name_type>>
		clusters_from(const std::string& membership, const std::function<
			default_graph::vertex_name_type(const std::string&)>& cnv =
			[](const std::string& str) { return std::stoll(str); })
	{
		std::ifstream fin(membership);
		unsigned count;
		std::string names;
		fin >> names;
		if (names.find('[') != -1)
		{
			fin >> count;
			std::string clusters;
			getline(fin, names); //костыль
			std::vector<std::set<default_graph::vertex_name_type>> ret(count);
			while ((bool)getline(fin, names) && (bool)getline(fin, clusters))
			{
				std::vector<default_graph::vertex_name_type> n;
				std::vector<unsigned> c;
				sprclib::pars(names, " ", [&n, &cnv](PARS_ARGS(std::string))
				{
					n.push_back(cnv(params.str));
				});
				sprclib::pars(clusters, " ", [&c, &cnv](PARS_ARGS(std::string))
				{
					c.push_back(std::stoll(params.str) - 1);
				});
				for (unsigned i = 0; i < n.size(); ++i)
					ret[c[i]].insert(n[i]);
			}
			fin.close();
			return ret;
		}
		else
		{
			std::vector<std::set<default_graph::vertex_name_type>> ret;
			getline(fin, names); // костыль
			while (getline(fin, names))
			{
				ret.push_back(std::set<default_graph::vertex_name_type>());
				sprclib::pars(names, " ", [&ret, &cnv](PARS_ARGS(std::string))
				{
					ret.back().insert(cnv(params.str));
				});
			}
			fin.close();
			return ret;
		}

	}

	// Получает разбиение на кластеры из результатов генерации графа
	std::vector<std::set<default_graph::vertex_name_type>>
		clusters_from(const result_graph_generator& res)
	{
		std::vector<std::set<default_graph::vertex_name_type>> ret(res.clusters.size());
		for (unsigned i = 0; i < ret.size(); ++i)
			ret[i].insert(res.clusters[i].g.begin(), res.clusters[i].g.end());
		return ret;
	}
}
#endif // _GRAPHUTILITY_H_
