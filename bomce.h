#pragma once
#ifndef _BOMCE_H_
#define _BOMCE_H_
#define _M_X64

#include <set>
#include <map>
#include <memory>
#include <math.h>
#include <cstring>
#include "graph.h"
#include "bomce_utils.h"
#include "graphutility.h"
//#include <sprclib\fraction.h>

#ifdef _DEBUG
#include <iostream>
#endif

#define BOMCE_TEMPLATE_DEF template<class type, class edge_value_type, class e_count_type, template<class> class compare_type, class mod_type>
#define BOMCE_TEMPLATE_ARGS type, edge_value_type, e_count_type, compare_type, mod_type
#define BOMCE_DEFAULT_GRAPH_ARGS unsigned, unsigned, unsigned, compare_ptr, float
#define BDGA_NOTMOD unsigned, unsigned, unsigned, compare_ptr

//#include "dna_alg.h"

namespace webgr
{
	template<class sum_type>
	// Кластер для алгоритма bomce
	struct bomce_cluster
	{
		sum_type e_count, loops; // число рёбер и петель
		unsigned id, v_count; // id и число вершин
		sum_type in_sum, out_sum; // суммы входящих и исходящих степеней вершин в кластере

		bomce_cluster(unsigned _id = 0) :id(_id), e_count(0),
			v_count(0), in_sum(0), out_sum(0), loops(0)
		{}
		bomce_cluster(const bomce_cluster& cl) :id(cl.id), e_count(cl.e_count),
			v_count(cl.v_count), in_sum(cl.in_sum), out_sum(cl.out_sum),
			loops(cl.loops)
		{}
		bomce_cluster(bomce_cluster&& cl) :id(cl.id), e_count(cl.e_count),
			v_count(cl.v_count), in_sum(cl.in_sum), out_sum(cl.out_sum),
			loops(cl.loops)
		{
			cl.e_count = 0;
			cl.v_count = 0;
			cl.in_sum = 0;
			cl.loops = 0;
			cl.out_sum = 0;
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
			in_sum = cl.in_sum;
			out_sum = cl.out_sum;
			loops = cl.loops;
			cl.e_count = 0;
			cl.v_count = 0;
			cl.in_sum = 0;
			cl.out_sum = 0;
			cl.loops = 0;
			return *this;
		}

		bomce_cluster& operator=(const bomce_cluster& cl)
		{
			e_count = cl.e_count;
			v_count = cl.v_count;
			id = cl.id;
			in_sum = cl.in_sum;
			out_sum = cl.out_sum;
			loops = cl.loops;
			return *this;
		}
	};

	BOMCE_TEMPLATE_DEF
		// Семейство распределений вершин по кластерам
		// В случае без ГА один энкземпляр на алгоритм
		class bomce_vdist_family;

	BOMCE_TEMPLATE_DEF
		// Абстракция для параметра q
		class q_modular
	{
	protected:
		const webgr::graph<BDGA_NOTMOD>& temp_graph;
		mod_type alpha_d_ecount;
	public:
		explicit q_modular(const webgr::graph<BDGA_NOTMOD>& temp_gr) :
			temp_graph(temp_gr), alpha_d_ecount(1) {}
		// При удалении вершины из кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>&,
			const e_count_type&, unsigned) const = 0;
		// При помещении вершины в кластер
		virtual mod_type operator()(unsigned, const bomce_cluster<e_count_type>&,
			const e_count_type&) const = 0;
		// Упрощённая(не нормированная) версия модулярности для одного кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>&) const = 0;
		void set_alpha(const mod_type& alpha)
		{
			alpha_d_ecount = alpha / temp_graph.edge_count();
		}
		virtual void set_graph(const webgr::graph<BDGA_NOTMOD>& temp_gr)
		{
			const_cast<webgr::graph<BDGA_NOTMOD>&>(temp_graph) = temp_gr;
		}
		virtual ~q_modular() {}
	};

	BOMCE_TEMPLATE_DEF
		class q_noorient :public q_modular<BOMCE_TEMPLATE_ARGS>
	{
	private:
		static const mod_type two;
		mod_type two_ecount;
	public:
		explicit q_noorient(const webgr::graph<BDGA_NOTMOD>& temp_gr) :
			q_modular<BOMCE_TEMPLATE_ARGS>(temp_gr){}
		// При удалении вершины из кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>& cl,
			const e_count_type& current_aij, unsigned gr_index) const
		{
            auto cur = this->temp_graph[gr_index]->deg();
            mod_type ret = two * current_aij - this->alpha_d_ecount*
				(cl.in_sum + cl.out_sum - cur)*cur;
			return ret;
		}
		// При помещении вершины в кластер
		virtual mod_type operator()(unsigned gr_index, const bomce_cluster<
			e_count_type>& it_cl, const e_count_type& aij_temp) const
		{
            mod_type ret = two * aij_temp - this->alpha_d_ecount*this->temp_graph[gr_index]->deg() *
				(it_cl.in_sum + it_cl.out_sum);
			return ret;
		}
		// Упрощённая(не нормированная) версия модулярности для одного кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>& it_cl) const
		{
			mod_type e_sum = it_cl.in_sum + it_cl.out_sum;
			mod_type ret = two*it_cl.e_count - e_sum*e_sum * two_ecount;
			return ret;
		}
		virtual void set_graph(const webgr::graph<BDGA_NOTMOD>& temp_gr)
		{
			q_modular<BOMCE_TEMPLATE_ARGS>::set_graph(temp_gr);
			two_ecount = mod_type(1) / (two*temp_gr.edge_count());
		}
		virtual ~q_noorient() {}
	};

	BOMCE_TEMPLATE_DEF
	const mod_type q_noorient<BOMCE_TEMPLATE_ARGS>::two(2);

	BOMCE_TEMPLATE_DEF
		class q_orient :public q_modular<BOMCE_TEMPLATE_ARGS>
	{
		mod_type i_ecount;
	public:
		explicit q_orient(const webgr::graph<BDGA_NOTMOD>& temp_gr) :
			q_modular<BOMCE_TEMPLATE_ARGS>(temp_gr) {}
		// При удалении вершины из кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>& cl,
			const e_count_type& current_aij, unsigned gr_index) const
		{
            auto* t = this->temp_graph[gr_index];
			auto in_sum = cl.in_sum - t->in_d;
			auto out_sum = cl.out_sum - t->out_d;
			auto in_d = t->in_d;
			auto out_d = t->out_d;
            mod_type ret = current_aij - (out_d*in_sum + out_sum*in_d)* this->alpha_d_ecount;
			return ret;
		}
		// При помещении вершины в кластер
		virtual mod_type operator()(unsigned gr_index, const bomce_cluster<
			e_count_type>& it_cl, const e_count_type& aij_temp) const
		{
			auto in_sum = it_cl.in_sum;
			auto out_sum = it_cl.out_sum;
            auto in_d = this->temp_graph[gr_index]->in_d;
            auto out_d = this->temp_graph[gr_index]->out_d;
            mod_type ret = aij_temp - (out_d*in_sum + out_sum*in_d)* this->alpha_d_ecount;
			return ret;
		}
		// Упрощённая(не нормированная) версия модулярности для одного кластера
		virtual mod_type operator()(const bomce_cluster<e_count_type>& it_cl) const
		{
			mod_type ret = it_cl.e_count - mod_type(it_cl.in_sum*it_cl.out_sum) *i_ecount;
			return ret;
		}
		virtual void set_graph(const webgr::graph<BDGA_NOTMOD>& temp_gr)
		{
			q_modular<BOMCE_TEMPLATE_ARGS>::set_graph(temp_gr);
            i_ecount = static_cast<mod_type>(1) / this->temp_graph.edge_count();
		}
		virtual ~q_orient() {}
	};

	template<class v_name_type>
	// описание вершины графа для содержимого кластеров
	struct cl_vertex :public std::pair<unsigned, v_name_type>
	{
		cl_vertex(unsigned _first, const v_name_type& _second) :
			std::pair<unsigned, v_name_type>(_first, _second)
		{}
        template<typename ...types>
        cl_vertex(const std::initializer_list<types...>& vals) :
            std::pair<unsigned, v_name_type>(vals)
		{}
		// Имя вершины
		v_name_type& name()
		{
			return this->second;
		}
		// Имя вершины
		const v_name_type& name() const
		{
			return this->second;
		}
		// Индекс вершины
		unsigned& index()
		{
			return this->first;
		}
		// Индекс вершины
		const unsigned& index() const
		{
			return this->first;
		}
	};

	BOMCE_TEMPLATE_DEF
		class bomce;

	BOMCE_TEMPLATE_DEF
		// Семейство распределений вершин по кластерам
		class bomce_vdist_family
	{
		friend class q_modular<BOMCE_TEMPLATE_ARGS>;
		friend q_orient<BOMCE_TEMPLATE_ARGS>;
		friend q_noorient<BOMCE_TEMPLATE_ARGS>;
		friend bomce<BOMCE_TEMPLATE_ARGS>;
		//friend dna_order_alg<BOMCE_TEMPLATE_ARGS>;
	private:
		const q_modular<BOMCE_TEMPLATE_ARGS>* q_mod; // метод вычисления оптимальности расположения вершины
		bomce_params params;
		// std::vector<e_count_type> loops; // петли исходного графа
		const bomce<BOMCE_TEMPLATE_ARGS>* _alg;
	public:
		// Распределение вершин
		class v_distribution
		{
			friend bomce<BOMCE_TEMPLATE_ARGS>;
			friend bomce_vdist_family<BOMCE_TEMPLATE_ARGS>;
			//friend dna_order_alg<BOMCE_TEMPLATE_ARGS>;

            static const std::uniform_real_distribution<double> dist_begin;
            static const std::uniform_real_distribution<double> dist_end;
		private:
			std::queue<unsigned> del_index; // индексы удалённых кластеров
			sprclib::list<bomce_cluster<e_count_type>> clusters; // кластеры
			std::vector<cluster_iterator<e_count_type>> v_to_cl; // индекс вершины на кластер
			mod_type alpha; // коэффициент масштабирования кластеров
			const bomce_vdist_family* fam; // семейство текущего распределения

			v_distribution(const bomce_vdist_family* _fam, const std::queue<unsigned>& d_index,
				const sprclib::list<bomce_cluster<e_count_type>>& _clusters,
				unsigned vtocl_size) :v_to_cl(vtocl_size),
				del_index(d_index), clusters(_clusters), fam(_fam)
			{}
		public:
			v_distribution(const bomce_vdist_family* _fam = nullptr) :fam(_fam) {}
			v_distribution(v_distribution&& v_d) :del_index(std::move(v_d.del_index)),
				v_to_cl(std::move(v_d.v_to_cl)), clusters(std::move(v_d.clusters)),
				fam(v_d.fam)
			{}
			v_distribution(const v_distribution& v_d) :fam(v_d.fam),
                del_index(v_d.del_index), v_to_cl(v_d.fam->_alg->temp_graph.size()),
				clusters(v_d.clusters)
			{
				// точная копия
				// работает за n*log(n) и служит "затычкой", используется copy

				std::map<unsigned, cluster_iterator<e_count_type>> id_reindex; // id старого -> итератор на новый
				auto it2 = v_d.clusters.begin();
				for (auto it = clusters.begin(); it != clusters.end(); ++it, ++it2)
				{
					id_reindex.insert({ it2->id, it });
				}
				for (unsigned i = 0; i < v_to_cl.size(); ++i)
					v_to_cl[i] = id_reindex[v_d.v_to_cl[i]->id];
			}
			v_distribution& operator=(const v_distribution& v_d)
			{
				del_index = v_d.del_index;
                v_to_cl.resize(family()->alg()->temp_graph.size());
				clusters = v_d.clusters;
				fam = v_d.fam;
				std::map<unsigned, cluster_iterator<e_count_type>> id_reindex; // id старого -> итератор на новый
				for (auto it = clusters.begin(), it2 = v_d.clusters.begin();
					it != clusters.end(); ++it, ++it2)
				{
					id_reindex.insert({ it2->id, it });
				}
				for (unsigned i = 0; i < v_to_cl.size(); ++i)
					v_to_cl[i] = id_reindex[v_d.v_to_cl[i]->id];
				return *this;
			}
			v_distribution& operator=(v_distribution&& v_d)
			{
				// работает линейно

				del_index = std::move(v_d.del_index);
				v_to_cl = std::move(v_d.v_to_cl);
				clusters = std::move(v_d.clusters);
				fam = v_d.fam;
				return *this;
			}
            template<class comp>
			inline void _Init(std::map <cluster_iterator<e_count_type>, unsigned, comp>&,
				const std::map<unsigned, edge_value_type>&);
			void init();
			void add_to_cl(const cluster_iterator<e_count_type>&,
				unsigned, e_count_type);
			void make_cluster(unsigned);
			void erase_from_cl(cluster_iterator<e_count_type>&,
				unsigned, e_count_type, bool need_del = true);
			mod_type move_vertex(unsigned);
			mod_type move_vertex(unsigned, const cluster_iterator<e_count_type>&);
			v_distribution copy();
			std::vector<v_distribution> copy(unsigned);
			static std::vector<v_distribution> copy(v_distribution&&, unsigned);
			const bomce_vdist_family<BOMCE_TEMPLATE_ARGS>* family()
			{
				return fam;
			}
			void next_iteration(unsigned long long seed);
			void clear()
			{
				v_to_cl.clear();
				clusters.clear();
				del_index = std::move(std::queue<unsigned>());
			}
			~v_distribution() {}

			/*
				bool operator<(const v_distribution& v_d) const
				{
					return modular < v_d.modular;
				}
				bool operator>(const v_distribution& v_d) const
				{
					return v_d < *this;
				}
				bool operator<=(const v_distribution& v_d) const
				{
					return !(v_d < *this);
				}
				bool operator>=(const v_distribution& v_d) const
				{
					return !(*this < v_d);
				}
				bool operator==(const v_distribution& v_d) const
				{
					return modular == v_d.modular;
				}
				bool operator!=(const v_distribution& v_d) const
				{
					return !(this == v_d);
				}
				*/
		};
	private:
		friend class v_distribution;
		std::vector<v_distribution> dist; // сами распределения(структура не утверждена)
	public:
		typedef typename std::vector<v_distribution>::iterator iterator;

		bomce_vdist_family(const bomce<BOMCE_TEMPLATE_ARGS>* bmc = nullptr,
			const q_modular<BOMCE_TEMPLATE_ARGS>* _q_mod = nullptr,
			unsigned dist_count = std::thread::hardware_concurrency()) :
			_alg(bmc), q_mod(_q_mod), dist(dist_count)
		{
			for (auto& i : dist)
				i.fam = this;
		}

		/*bomce_vdist_family(const graph<BDGA_NOTMOD>& gr,
						   const bomce_tmpresult_type<type>& tmp_result, typename
						   bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution& dist_init,
						   const q_modular<BOMCE_TEMPLATE_ARGS>* _q_mod = nullptr,
						   unsigned dist_count = std::thread::hardware_concurrency()) :
			temp_graph(gr), temp_result(tmp_result), q_mod(_q_mod)
		{
			dist = std::move(dist_init.copy(dist_count));
		}

		bomce_vdist_family(const graph<BDGA_NOTMOD>& gr,
						   const bomce_tmpresult_type<type>& tmp_result, typename
						   bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution&& dist_init,
						   const q_modular<BOMCE_TEMPLATE_ARGS>* _q_mod = nullptr,
						   unsigned dist_count = std::thread::hardware_concurrency()) :
			temp_graph(gr), temp_result(tmp_result), q_mod(_q_mod)
		{
			dist = std::move(copy(std::move(dist_init), dist_count));
		}*/

		bomce_vdist_family(const bomce_vdist_family&) = delete;
		bomce_vdist_family(bomce_vdist_family&&) = delete;

		void init(typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
			v_distribution& dist_init, unsigned count = std::thread::
			hardware_concurrency())
		{
			dist = std::move(dist_init.copy(count)); // делаем count копий распределения dist_init
		}

		void init(typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
			v_distribution&& dist_init, unsigned count = std::thread::
			hardware_concurrency())
		{
			// делаем count копий распределения dist_init
			dist = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
				v_distribution::copy(std::move(dist_init), count));
		}
		void init(const bomce_params&, unsigned count = std::thread::hardware_concurrency());
		iterator begin() const
		{
			return dist.begin();
		}

		iterator end() const
		{
			return dist.end();
		}

		v_distribution& operator[](unsigned index)
		{
			return dist[index];
		}

		const v_distribution& operator[](unsigned index) const
		{
			return dist[index];
		}

		unsigned size() const
		{
			return dist.size();
		}

		void resize(unsigned size)
		{
			dist.resize(size);
		}

		void resize(unsigned size, v_distribution& v_d)
		{
			if (dist.size() < size)
			{
				auto add_vd = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
					v_distribution::copy(std::move(copy(v_d)), size - dist.size()));
				unsigned old_size = dist.size();
				dist.resize(size);
				std::move(add_vd.begin(), add_vd.end(), dist.begin() + old_size);
			}
			else if (dist.size() > size)
				dist.resize(size);
		}

		const bomce<BOMCE_TEMPLATE_ARGS>* alg() const
		{
			return _alg;
		}

		~bomce_vdist_family()
		{
			delete q_mod;
		}
	};

	BOMCE_TEMPLATE_DEF
		// Диапазон возможных вариантов для начального значения индекса масштабирования кластеров
		const std::uniform_real_distribution<double>
		bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::dist_begin(0.6, 0.8);

	BOMCE_TEMPLATE_DEF
		// Диапазон возможных вариантов для конечного значения индекса масштабирования кластеров
		const std::uniform_real_distribution<double>
		bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::dist_end(1.1, 1.3);

	BOMCE_TEMPLATE_DEF
		using v_dist = typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>
		::v_distribution;

	BOMCE_TEMPLATE_DEF
		// Реиннициализация распределения вершин
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::init()
	{
		del_index = std::move(std::queue<unsigned>());
		clusters.clear();
		auto& temp_graph = family()->alg()->temp_graph;
		v_to_cl.resize(temp_graph.size());
		for (auto& i : temp_graph)
		{
			make_cluster(i.first);
		}
	}

	BOMCE_TEMPLATE_DEF
		// Добавление вершины и индексом gr_index в кластер _cl, где aij_sum_delt -
		// число рёбер(или их суммарный вес) между gr_index и другими вершинами этого кластера
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		add_to_cl(const cluster_iterator<e_count_type>& _cl, unsigned gr_index,
			e_count_type aij_sum_delt)
	{
		v_to_cl[gr_index] = _cl;
		auto& cl = *_cl;
		auto& bmc = *family()->alg();
		cl.v_count += bmc.temp_result[gr_index].size();
		auto* t = bmc.temp_graph[gr_index];
		cl.e_count += aij_sum_delt + t->loop;
		cl.in_sum += t->in_d;
		cl.out_sum += t->out_d;
		cl.loops += t->loop;
	}

	BOMCE_TEMPLATE_DEF
		// Создание кластера с единственной вершиной gr_index в нём
		// Перед вызовом удалить вершину и кластера, в котором она лежит
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		make_cluster(unsigned gr_index)
	{
		unsigned index;
		if (del_index.empty())
			index = clusters.size();
		else
		{
			index = del_index.front();
			del_index.pop();
		}
		clusters.push_front(bomce_cluster<e_count_type>(index));
		add_to_cl(clusters.begin(), gr_index, 0);
	}

	BOMCE_TEMPLATE_DEF
		// Удаление вершины current из кластера _it. aij_temp  - число рёбер
		// (или их суммарный вес) между current и другими вершинами из _it,
		// need_del - флаг, нужно ли удалять _it если он станет пустым
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		erase_from_cl(cluster_iterator<e_count_type>& _it, unsigned current,
			e_count_type aij_temp, bool need_del)
	{
		auto& it = *_it;
		it.v_count -= family()->alg()->temp_result[current].size();
		if (it.v_count == 0 && need_del)
		{
			del_index.push(it.id);
			clusters.erase(_it);
		}
		else
		{
			auto* t = family()->alg()->temp_graph[current];
			it.e_count -= aij_temp + t->loop;
			it.in_sum -= t->in_d;
			it.out_sum -= t->out_d;
			it.loops -= t->loop;
		}
		v_to_cl[current] = clusters.end();
	}

	BOMCE_TEMPLATE_DEF
		// перемещает вершину с индексом gr_index в кластер cl
		// удаление из предыдущего и пересчёт рёбер происходит 
		mod_type bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		move_vertex(unsigned gr_index, const cluster_iterator<e_count_type>& cl)
	{
		if (v_to_cl[gr_index]->id == cl->id)
			return 0.0;
		unsigned cur_e = 0, cl_e = 0;
		family()->alg()->temp_graph[gr_index].foreach_neighbors(
			[&](std::pair<unsigned, edge_value_type> pr) // По ссылке может быть UB?
		{
			if (v_to_cl[pr.first]->id == v_to_cl[gr_index]->id)
				++cur_e;
			else if (v_to_cl[pr.first]->id == cl->id)
				++cl_e;
		});
        const auto*& q_mod = family()->alg()->q_mod;
        auto ret = q_mod->operator()(gr_index, *cl, cl_e) -
			q_mod->operator()(*v_to_cl[gr_index], cur_e, gr_index);
		erase_from_cl(v_to_cl[gr_index], gr_index, cur_e);
		add_to_cl(cl, gr_index, cl_e);
		return ret;
	}

    // иннициализация cl(смежные кластеры)
    BOMCE_TEMPLATE_DEF
    template<class comp>
	void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		_Init(std::map <cluster_iterator<e_count_type>, unsigned,
            comp>& cl, const std::map<unsigned,
			edge_value_type>& graph_v)
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

	BOMCE_TEMPLATE_DEF
		// Перемещение вершины в соседний кластер на основе q_mod, где максимум - туда и перемещаем
		mod_type bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
		move_vertex(unsigned gr_index)
	{
		/*
		пройтись по смежным вершинам
		получить множество смежных кластеров
		определить предпочтительный для присоединения кластер
		выполнить либо add_ещ_сд либо make_cluster(модулярность не должна ухудшиться),
		перед этим выполнив cluster_erase если нужно
		*/
		typedef std::map<cluster_iterator<e_count_type>,
			unsigned, cl_cmp<e_count_type>> mp_type;
		mp_type cl; // смежный кластер на aij_temp(кол. рёбер между it и кластером)
		auto& temp_graph = family()->alg()->temp_graph;
		_Init(cl, temp_graph[gr_index]->input);
		_Init(cl, temp_graph[gr_index]->output);
		unsigned current_aij = 0;
		auto it = cl.find(v_to_cl[gr_index]);
		if (it != cl.end())
		{
			current_aij = it->second;
			cl.erase(v_to_cl[gr_index]);
		}
		if (cl.empty())
			return 0.0;
		typename mp_type::iterator it1 = cl.begin();
		cluster_iterator<e_count_type>& cl_it = const_cast<
			cluster_iterator<e_count_type>&>(it1->first);
		auto& qmod = *family()->q_mod;
		std::pair<mod_type, unsigned> max = { qmod(gr_index,
			*cl_it, it1->second), it1->second };
		// max - наилучшее внекластерное
		for (++it1; it1 != cl.end(); ++it1)
		{
			cluster_iterator<e_count_type>& it_tmp = const_cast<
				cluster_iterator<e_count_type>&>(it1->first);
			std::pair<mod_type, unsigned> tmp = { qmod(gr_index,
				*it_tmp, it1->second), it1->second };
			if (tmp.first > max.first)
			{
				max = tmp;
				cl_it = it_tmp;
			}
		}
		mod_type p = qmod(*v_to_cl[gr_index],
			current_aij, gr_index);
		if (p < max.first && max.first > 0) // при перемещении в друго максимезируется Q
		{
			this->erase_from_cl(v_to_cl[gr_index], gr_index, current_aij);
			this->add_to_cl(cl_it, gr_index, max.second);
			return max.first - p;
		}
		else if (p < 0) // При удалении из текущего и перемещении в новый максимезируется Q
		{
			this->erase_from_cl(v_to_cl[gr_index], gr_index, current_aij);
			make_cluster(gr_index);
			return -p;
		}
		return 0.0;
	}

	BOMCE_TEMPLATE_DEF
		// Возвращает копию текущего распределения за линейное время
		// При копировании изменяет внутренние индексы а после возвращает как было
		typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution
		bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::copy()
	{
		v_dist<BOMCE_TEMPLATE_ARGS> ret(fam, del_index, clusters, v_to_cl.size());
		std::vector<unsigned> old_index(clusters.size()); // старые индексы
		unsigned i = 0;
		for (auto& cl : clusters)
		{
			old_index[i] = cl.id;
			cl.id = i++;
		}
		// переиндексировали кластеры
		std::vector<cluster_iterator<e_count_type>> new_it(clusters.size());
		i = 0;
		for (auto it = ret.clusters.begin(); it != ret.clusters.end(); ++it)
			new_it[i++] = it;
		// new_it - индексированый контейнер итераторов на кластеры
		for (i = 0; i < v_to_cl.size(); ++i)
			ret.v_to_cl[i] = new_it[v_to_cl[i]->id];
		i = 0;
		for (auto& cl : clusters)
			cl.id = old_index[i++];
		return ret;
	}

	BOMCE_TEMPLATE_DEF
		// Возвращает вектор из copy_count копий распределения v_d
		// При копировании изменяет внутренние индексы а после возвращает как было
		std::vector<typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
		v_distribution> bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
		v_distribution::copy(typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
			v_distribution&& v_d, unsigned copy_count)
	{
		std::vector<v_dist<BOMCE_TEMPLATE_ARGS>> ret(copy_count);
		ret[0] = std::move(v_d);
		//static_assert(false, "browken begin iterator in ret[0]");
		for (unsigned i = 1; i < ret.size(); ++i)
		{
			ret[i] = std::move(v_distribution(ret[0].fam, ret[0].del_index,
				ret[0].clusters, ret[0].v_to_cl.size()));
		}
		std::vector<unsigned> old_index(ret[0].clusters.size());
		unsigned i = 0;
		for (auto& cl : ret[0].clusters)
		{
			old_index[i] = cl.id;
			cl.id = i++;
		}
		std::vector<cluster_iterator<e_count_type>> new_it(ret[0].clusters.size());
		for (unsigned k = 1; k < copy_count; ++k)
		{
			i = 0;
			for (auto it = ret[k].clusters.begin(); it != ret[k].clusters.end(); ++it)
				new_it[i++] = it;
			for (i = 0; i < ret[0].v_to_cl.size(); ++i)
				ret[k].v_to_cl[i] = new_it[ret[0].v_to_cl[i]->id];
		}
		i = 0;
		for (auto& cl : ret[0].clusters)
			cl.id = old_index[i++];
		return ret;
	}

	BOMCE_TEMPLATE_DEF
		// Возвращает вектор из copy_count копий текущего распределения
		// При копировании изменяет внутренние индексы а после возвращает как было
		std::vector<typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
		v_distribution> bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
		v_distribution::copy(unsigned copy_count)
	{
		std::vector<v_dist<BOMCE_TEMPLATE_ARGS>> ret(copy_count);
		for (auto& i : ret)
		{
			i = std::move(v_distribution(fam, del_index,
				clusters, v_to_cl.size()));
		}
		std::vector<unsigned> old_index(clusters.size());
		unsigned i = 0;
		for (auto& cl : clusters)
		{
			old_index[i] = cl.id;
			cl.id = i++;
		}
		std::vector<cluster_iterator<e_count_type>> new_it(clusters.size());
		for (unsigned k = 0; k < copy_count; ++k)
		{
			i = 0;
			for (auto it = ret[k].clusters.begin();
				it != ret[k].clusters.end(); ++it)
			{
				new_it[i++] = it;
			}
			for (i = 0; i < v_to_cl.size(); ++i)
				ret[k].v_to_cl[i] = new_it[v_to_cl[i]->id];
		}
		i = 0;
		for (auto& cl : clusters)
			cl.id = old_index[i++];
		return ret;
	}

	BOMCE_TEMPLATE_DEF
		// Серия проходов по всем вершинам графа в случайном порядке для каждого значения индекса
		// масштабирования графа
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
		v_distribution::next_iteration(unsigned long long seed)
	{
		mod_type contin;
		auto& temp_graph = family()->alg()->temp_graph;
		std::vector<unsigned> prep(temp_graph.size());
		for (unsigned i = 0; i < temp_graph.size(); ++i)
			prep[i] = i;
#ifdef _M_X64
		std::mt19937_64 gen(seed); // Генератор для индекса масштабирования графа
#else
		std::mt19937 gen(seed); // Генератор для индекса масштабирования графа
#endif
		mod_type begin = family()->params[0];// dist_begin<BOMCE_TEMPLATE_ARGS>(gen);
		mod_type end = family()->params[1];// dist_end<BOMCE_TEMPLATE_ARGS>(gen);
		mod_type step = family()->params[2];

		static const mod_type _01(0.1), _001(0.01), contin_lim(1e-4);
		if (step < _01)
		{
            throw std::string(("'step' value too small, step = " + std::to_string(step) +
				". 'step' must be >= 0.1").c_str());
		}
		for (alpha = begin; alpha <= end; )
		{
			const_cast<q_modular<BOMCE_TEMPLATE_ARGS>*>(family()->q_mod)->
				set_alpha(mod_type(alpha));
			do // обходы по графу проводятся до тех пор, пока удаётся максимезировать модулярность
			{
				contin = 0;
				std::shuffle(prep.begin(), prep.end(), gen);
				for (auto& i : prep)
				{
					contin += this->move_vertex(i);
				}
				//static_assert(false, "Method does not converge");

			} while (contin > contin_lim);

			alpha += step;
			if (alpha > end && alpha - end < step - _001)
				alpha = end;
		}
	}

	BOMCE_TEMPLATE_DEF
		// Инициализация count копиями начальных распределений
		void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::init(const bomce_params& _param,
			unsigned count)
	{
		v_distribution v_d(this);
		v_d.init();
		this->params = _param;
		dist = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
			v_distribution::copy(std::move(v_d), count));
	}

	BOMCE_TEMPLATE_DEF
		// Результат кластеризации
		struct bomce_result
	{
		// Сруктура связей кластеров
		graph<BDGA_NOTMOD> const& connections;
		// Параметры кластеров
		std::vector<cluster_iterator<e_count_type>> const& v_to_cl;
		// Состав вешин в каждом кластере
		bomce_tmpresult_type<type> const& vertexes;
		// Список кластеров
		sprclib::list<bomce_cluster<e_count_type>> const& clusters;
		// число петель в исходном графе
		e_count_type loop_count;
		mod_type modularity;

		bomce_result(graph<BDGA_NOTMOD> const& gr,
			std::vector<cluster_iterator<e_count_type>> const& _v_to_cl,
			bomce_tmpresult_type<type> const& temp_result, sprclib::list<
			bomce_cluster<e_count_type>> const& _clusters, e_count_type
			loops, const mod_type& mod) : v_to_cl(_v_to_cl), clusters(_clusters),
			connections(gr), vertexes(temp_result), loop_count(loops),
			modularity(mod)
		{}

		/*
		Запись в файл в формате:
		{модулярность}
		{число кластеров}
		{имя вершины};{номер кластера}
		*/
		void contents(const std::string& csv_file_path) const
		{
			std::ofstream fout(csv_file_path);
			fout << modularity << std::endl << clusters.size();
			for (unsigned cl = 0; cl < vertexes.size(); ++cl)
			{
				for (auto& v : vertexes[cl])
					fout << std::endl << v.name() << ';' << cl;
			}
			fout.close();
		}

		//BOMCE_TEMPLATE_DEF
			/*
			Запись в файл в формате:
			{модулярность};{число кластеров}
			{имя вершины};{номер кластера}
			*/
			/*void propertes(const std::string& csv_file_path) const
		{
			std::ofstream fout(csv_file_path);
			fout << clusters.size() << std::endl << loop_count;
			for(auto& cl : clusters)

			for (unsigned cl = 0; cl < vertexes.size(); ++cl)
			{
				for (auto& v : vertexes[cl])
					fout << std::endl << v.second << ';' << cl;
			}
			fout.close();
		}*/
	};

	enum class orientation
	{
		orient,
		no_orient
	};

	BOMCE_TEMPLATE_DEF
		// Кластеризатор
		class bomce
	{
        friend class bomce_vdist_family<BOMCE_TEMPLATE_ARGS>;
        friend class bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution;

    private:
        std::uniform_int_distribution<unsigned long long> seed_dist;

		e_count_type edge_count, loop_count; // число рёбер и число петель в графе
		graph<BDGA_NOTMOD> temp_graph; // кластеризуемый граф или его конденсация
		bomce_tmpresult_type<type> temp_result; // номер кластера на список верщин в нём
		bomce_vdist_family<BOMCE_TEMPLATE_ARGS> v_f; // Семейство распределений по кластерам
		//std::vector<e_count_type> loops; // Петли исходного графа
		std::vector<mod_type> mods; // Значения модулярности для каждого распределения
		unsigned t_count; // Число потоков

#ifdef _M_X64
        std::mt19937_64 seed_gen; // зерно генерации зёрен
#else
        std::mt19937 seed_gen; // зерно генерации зёрен
#endif

		void condense(v_dist<BOMCE_TEMPLATE_ARGS>&);
		void _Init(unsigned, const bomce_params&);
		mod_type _Cluster_Mod(unsigned index);
	public:

		explicit bomce(orientation ornt) : edge_count(0),
            v_f(this), loop_count(0), seed_dist(1ull, 1ull << 25),
            seed_gen(std::random_device()())
		{
			switch (ornt)
			{
			case(orientation::no_orient):
				v_f.q_mod = new q_noorient<BOMCE_TEMPLATE_ARGS>(temp_graph);
				break;
			case(orientation::orient):
				v_f.q_mod = new q_orient<BOMCE_TEMPLATE_ARGS>(temp_graph);
				break;
			}
		}

		void init(graph<BDGA_NOTMOD>&& gr, const bomce_params& param,
			unsigned thread_count = std::thread::hardware_concurrency())
		{
			temp_graph = std::move(gr);
			_Init(thread_count, param);
		}

		void init(const graph<BDGA_NOTMOD>& gr, const bomce_params& param,
			unsigned thread_count = std::thread::hardware_concurrency())
		{
			temp_graph = gr;
			_Init(thread_count, param);
		}

		bool next_iteration(unsigned short lim);
		bool next_iteration()
		{
			static const double lg = 1.0 / log(2.0);
			unsigned short lim = log(temp_graph.size()) * lg /*/ t_count*/ + 0.1;
			if (lim == 0)
				lim = 1;
			return next_iteration(lim);
		}

		bomce_result<BOMCE_TEMPLATE_ARGS> result() const
		{
			return bomce_result<BOMCE_TEMPLATE_ARGS>(temp_graph, v_f[0].v_to_cl,
				temp_result, v_f[0].clusters, loop_count, mods[0]);
		}

		virtual ~bomce() {}
	};

    /*BOMCE_TEMPLATE_DEF
        const std::uniform_int_distribution<unsigned long long>
        bomce<BOMCE_TEMPLATE_ARGS>::seed_dist(1ull, 1ull << 25);*/

	BOMCE_TEMPLATE_DEF
		void bomce<BOMCE_TEMPLATE_ARGS>::_Init(unsigned _tcount,
			const bomce_params& param)
	{
		edge_count = temp_graph.edge_count();
		loop_count = 0;
		const_cast<q_modular<BOMCE_TEMPLATE_ARGS>*>(v_f.q_mod)->set_graph(temp_graph);
		temp_result.resize(temp_graph.capacity());
		for (auto& i : temp_graph)
			temp_result[i.first] = { { i.first, i.second->name() } };
		t_count = _tcount;
		v_f.init(param, _tcount);
		mods.resize(t_count);
		//memset(mods.data(), 0, mods.size() * sizeof(mods[0]));
		for (auto& i : mods)
			i = 0;
	}

	BOMCE_TEMPLATE_DEF
		mod_type bomce<BOMCE_TEMPLATE_ARGS>::_Cluster_Mod(unsigned index)
	{
        static constexpr unsigned long long module = 1ull << 24;
        unsigned long long seed = seed_dist(seed_gen) * (index + 1ull) +
			time(nullptr) % module;

		v_f[index].next_iteration(seed);
		mod_type mod(0);
		for (auto& v : v_f[index].clusters)
			mod += v_f.q_mod->operator()(v);
		return mod / (2 * temp_graph.edge_count());
	}

	BOMCE_TEMPLATE_DEF
		// кластеризация, а после конденсация кластеров
		bool bomce<BOMCE_TEMPLATE_ARGS>::next_iteration(unsigned short lim)
	{
		unsigned short max_mod = v_f.size();
		mod_type max_mod_value;
		static const mod_type mod_lim(1e-4);
		std::vector<bool> bl(v_f.size(), false);

		for (unsigned short i = 0; ;)
		{
			parallel_for(0u, v_f.size(), 1u, [&](unsigned index, unsigned id)
			{
				if (index != max_mod)
				{
					mod_type mod = _Cluster_Mod(index);
					bl[index] = (mod - mods[index] > mod_lim);
					mods[index] = mod;
				}

			}, t_count);
#ifdef _DEBUG
			for (auto& i : mods)
				std::cout << i << " ";
			std::cout << std::endl;
#endif
			max_mod = 0;
			max_mod_value = mods[0];
			for (unsigned short k = 1; k < v_f.size(); ++k)
			{
				if (mods[k] > max_mod_value)
				{
					max_mod_value = mods[k];
					max_mod = k;
				}
			}
			if (++i == lim)
				break;
			auto tmp = std::move(v_f[max_mod]);
			v_f[max_mod].init();
			v_f.dist = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
				v_distribution::copy(std::move(v_f[max_mod]), v_f.size()));
			v_f[max_mod] = std::move(tmp);
			for (unsigned k = 0; k < mods.size(); ++k)
			{
				if (k != max_mod)
					mods[k] = 0;
			}
			/*parallel_for(0u, v_f.size(), 1u, [&](unsigned index, unsigned id)
			{
				if (index != max_mod)
					v_f[index].init();

			}, t_count);*/
		}
		memswap(&v_f[max_mod], &v_f[0]);
        std::swap(mods[max_mod], mods[0]);
		v_f.resize(1);
		mods.resize(1);
		condense(v_f[0]);
		return bl[max_mod];
	}

	BOMCE_TEMPLATE_DEF
		// Построение конденсации графа на основе структуры кластеров
		void bomce<BOMCE_TEMPLATE_ARGS>::condense(
			v_dist<BOMCE_TEMPLATE_ARGS>& v_d)
	{
		decltype(temp_graph) meta_gr;
		meta_gr._Resize(v_d.clusters.size());
		unsigned t = 0;
		for (auto& i : v_d.clusters)
		{
			i.id = t++;
			meta_gr._New_Vertex(i.id, i.id);
		}
		bomce_tmpresult_type<type> tmp_result(v_d.clusters.size());
		for (auto& i : temp_graph)
		{
			t = v_d.v_to_cl[i.first]->id;
			tmp_result[t].emplace_after(tmp_result[t].begin(), temp_result[i.first]);
			meta_gr.connect_by_index(t, t, i.second->loop);
			for (auto& j : i.second->output)
				meta_gr.connect_by_index(t, v_d.v_to_cl[j.first]->id, j.second);
		}
		temp_result = std::move(tmp_result);
		temp_graph = std::move(meta_gr);
		v_d.v_to_cl.resize(v_d.clusters.size());
		t = 0;
		for (auto it = v_d.clusters.begin(); it != v_d.clusters.end(); ++it)
			v_d.v_to_cl[t++] = it;
		v_d.del_index = std::move(std::queue<unsigned>());
	}

	BOMCE_TEMPLATE_DEF
		double modular(const bomce_result<BOMCE_TEMPLATE_ARGS>& res)
	{
		double ret1 = 0, ret2 = 0;
		for (auto& i : res.clusters)
		{
			ret1 += i.e_count;
			ret2 += pow(i.in_sum + i.out_sum, 2);
		}
		double t = 1.0 / (2.0 * res.connections.edge_count());
		return t*(2.0*ret1 - res.loop_count - t*ret2);
	}

	BOMCE_TEMPLATE_DEF
		double modular2(const bomce_result<BOMCE_TEMPLATE_ARGS>& res,
			const graph<BDGA_NOTMOD>& gr)
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
                    d1 += i.second->deg() * k.second->deg();
				}
			}
		}
		return (2.0*e1 - d1 / (2.0*gr.edge_count())) / (2.0*gr.edge_count());
	}

	BOMCE_TEMPLATE_DEF
		double ratio_cut(const bomce_result<BOMCE_TEMPLATE_ARGS>& res)
	{
		double ret = 0;
		for (auto& i : res.clusters)
			ret += (i.in_sum + i.out_sum - 2.0*i.e_count) / res.vertexes[i.id].size();
		return ret;
	}

	BOMCE_TEMPLATE_DEF
		double normalized_cut(const bomce_result<BOMCE_TEMPLATE_ARGS>& res)
	{
		double ret = 0;
		for (auto& i : res.clusters)
			ret += 1.0 - (2.0*i.e_count) / (i.in_sum + i.out_sum);
		return ret;
	}

	BOMCE_TEMPLATE_DEF
		std::vector<std::set<default_graph::vertex_name_type>>
		clusters_from(const bomce_result<BOMCE_TEMPLATE_ARGS>& res)
	{
		std::vector<std::set<default_graph::vertex_name_type>> ret(res.clusters.size());
		for (unsigned i = 0; i < res.vertexes.size(); ++i)
		{
			for (auto& el : res.vertexes[i])
				ret[i].insert(el.name());
		}
		return ret;
	}
}
#endif // _bomce_H_
