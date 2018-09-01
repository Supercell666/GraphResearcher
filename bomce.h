#pragma once
#ifndef _BOMCE_H_
#define _BOMCE_H_

#include <set>
#include <map>
#include <vector>
#include <memory>
#include <math.h>
#include "graphutility.h"
#include "list.h"

namespace webgr
{
template<class sum_type>
struct bomce_cluster
{
    sum_type e_count, loops;
    unsigned id, v_count;
    sum_type in_sum, out_sum;

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

template<class sum_type>
struct cl_cmp
{
    bool operator()(const typename sprclib::list<bomce_cluster<sum_type>>::
                    iterator& cl1, const typename sprclib::list<bomce_cluster<
                    sum_type>>::iterator& cl2) const
    {
        return (cl1->id < cl2->id);
    }
};

#define BOMCE_TEMPLATE_DEF template<class type, class edge_value_type, class e_count_type, template<class> class compare_type>
#define BOMCE_TEMPLATE_ARGS type, edge_value_type, e_count_type, compare_type
#define BOMCE_DEFAULT_GRAPH_ARGS unsigned, unsigned, unsigned, compare_ptr

BOMCE_TEMPLATE_DEF
class bomce_vdist_family;

BOMCE_TEMPLATE_DEF
class q_modular
{
protected:
    const bomce_vdist_family<BOMCE_TEMPLATE_ARGS>& v_d;
public:
    q_modular(const bomce_vdist_family<BOMCE_TEMPLATE_ARGS>& v_dist) :v_d(v_dist) {}

    virtual double operator()(const typename sprclib::list<bomce_cluster<e_count_type>>::
                              iterator&, e_count_type, unsigned, double) const = 0;
    virtual double operator()(unsigned, const typename sprclib::list<
                              bomce_cluster<e_count_type>>::iterator&, e_count_type, double) const = 0;
    virtual double operator()(const typename sprclib::list<
                              bomce_cluster<e_count_type>>::iterator&) const = 0;
    virtual ~q_modular() {}
};

BOMCE_TEMPLATE_DEF
class q_noorient :public q_modular<BOMCE_TEMPLATE_ARGS>
{
public:
    q_noorient(const bomce_vdist_family<BOMCE_TEMPLATE_ARGS>& v_dist) :
        q_modular<BOMCE_TEMPLATE_ARGS>(v_dist) {}

    virtual double operator()(const typename sprclib::list<bomce_cluster<
                              e_count_type>>::iterator& cl, e_count_type current_aij,
                              unsigned gr_index, double alpha) const
    {
        auto cur = this->v_d.temp_graph[gr_index].deg();
        return 2.0 * current_aij - alpha*(cl->in_sum + cl->out_sum - cur)*cur /
                this->v_d.temp_graph.edge_count();
    }

    virtual double operator()(unsigned gr_index, const typename sprclib::
                              list<bomce_cluster<e_count_type>>::iterator& it_cl,
                              e_count_type aij_temp, double alpha) const
    {
        return 2.0 * aij_temp - alpha*this->v_d.temp_graph[gr_index].deg()*
                (it_cl->in_sum + it_cl->out_sum) / this->v_d.temp_graph.edge_count();
    }

    virtual double operator()(const typename sprclib::
                              list<bomce_cluster<e_count_type>>::iterator& it_cl) const
    {
        return 2.0*it_cl->e_count - std::pow(it_cl->in_sum
                                             + it_cl->out_sum, 2) / (2.0*this->v_d.temp_graph.edge_count());
    }

    virtual ~q_noorient() {}
};

BOMCE_TEMPLATE_DEF
class q_orient :public q_modular<BOMCE_TEMPLATE_ARGS>
{

public:
    q_orient(const bomce_vdist_family<BOMCE_TEMPLATE_ARGS>& v_dist) :
        q_modular<BOMCE_TEMPLATE_ARGS>(v_dist) {}

    virtual double operator()(const typename sprclib::list<bomce_cluster<
                              e_count_type>>::iterator& cl, e_count_type current_aij,
                              unsigned gr_index, double alpha) const
    {
        auto in_sum = cl->in_sum - this->v_d.temp_graph[gr_index].in_d;
        auto out_sum = cl->out_sum -this-> v_d.temp_graph[gr_index].out_d;
        auto in_d = this->v_d.temp_graph[gr_index].in_d;
        auto out_d = this->v_d.temp_graph[gr_index].out_d;
        return current_aij - (out_d*in_sum + out_sum*in_d)* alpha /
                this->v_d.temp_graph.edge_count();
    }

    virtual double operator()(unsigned gr_index, const typename sprclib::list<
                              bomce_cluster<e_count_type>>::iterator& it_cl,
                              e_count_type aij_temp, double alpha) const
    {
        auto in_sum = it_cl->in_sum;
        auto out_sum = it_cl->out_sum;
        auto in_d = this->v_d.temp_graph[gr_index].in_d;
        auto out_d = this->v_d.temp_graph[gr_index].out_d;
        return aij_temp - (out_d*in_sum + out_sum*in_d)* alpha /
                this->v_d.temp_graph.edge_count();
    }

    virtual double operator()(const typename sprclib::
                              list<bomce_cluster<e_count_type>>::iterator& it_cl) const
    {
        return it_cl->e_count - (double)it_cl->in_sum*it_cl->out_sum /
                this->v_d.temp_graph.edge_count();
    }

    virtual ~q_orient() {}
};

template<class v_name_type>
using bomce_tmpresult_type = std::vector<sprclib::list<std::pair<unsigned,
v_name_type>>>;

BOMCE_TEMPLATE_DEF
class bomce;

BOMCE_TEMPLATE_DEF
class bomce_vdist_family
{
    friend q_modular<BOMCE_TEMPLATE_ARGS>;
    friend q_orient<BOMCE_TEMPLATE_ARGS>;
    friend q_noorient<BOMCE_TEMPLATE_ARGS>;
    friend bomce<BOMCE_TEMPLATE_ARGS>;
private:
    graph<BOMCE_TEMPLATE_ARGS> const& temp_graph;
    const q_modular<BOMCE_TEMPLATE_ARGS>* q_mod;
    bomce_tmpresult_type<type> const& temp_result;
public:
    class v_distribution
    {
        friend bomce<BOMCE_TEMPLATE_ARGS>;
        friend bomce_vdist_family<BOMCE_TEMPLATE_ARGS>;
    private:
        std::queue<unsigned> del_index;
        sprclib::list<bomce_cluster<e_count_type>> clusters;
        std::vector<typename sprclib::list<bomce_cluster<
        e_count_type>>::iterator> v_to_cl;
        std::vector<e_count_type> loops; // петли исходного графа
        double alpha;
        const bomce_vdist_family* fam;

        v_distribution(const bomce_vdist_family* _fam, const std::queue<unsigned>& d_index,
                       const sprclib::list<bomce_cluster<e_count_type>>& _clusters,
                       unsigned vtocl_size, const std::vector<
                       e_count_type>& lps) :v_to_cl(vtocl_size), loops(lps),
            del_index(d_index), clusters(_clusters), fam(_fam)
        {}
    public:
        v_distribution(const bomce_vdist_family* _fam = nullptr) :fam(_fam) {}
        v_distribution(v_distribution&& v_d) :del_index(std::move(v_d.del_index)),
            v_to_cl(std::move(v_d.v_to_cl)), clusters(std::move(v_d.clusters)),
            loops(std::move(v_d.loops)), fam(v_d.fam)
        {}
        v_distribution(const v_distribution& v_d) :fam(v_d.fam),
            del_index(v_d.del_index), v_to_cl(fam->temp_graph.size()), clusters(v_d.clusters),
            loops(v_d.loops)
        {
            // точная копия
            std::map<unsigned, typename sprclib::list<bomce_cluster<
                    e_count_type>>::iterator> id_reindex; // id старого -> итератор на новый
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
            v_to_cl.resize(temp_graph.size());
            clusters = v_d.clusters;
            loops = v_d.loops;
            fam = v_d.fam;
            std::map<unsigned, typename sprclib::list<bomce_cluster<
                    e_count_type>>::iterator> id_reindex; // id старого -> итератор на новый
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
            del_index = std::move(v_d.del_index);
            v_to_cl = std::move(v_d.v_to_cl);
            clusters = std::move(v_d.clusters);
            loops = std::move(v_d.loops);
            fam = v_d.fam;
            return *this;
        }
        template<class comp>
        void _Init(std::map<typename sprclib::list<bomce_cluster<e_count_type>>::
                   iterator, unsigned, comp>&, const std::map<unsigned, edge_value_type>&);
        void init();
        void add_to_cl(const typename sprclib::list<bomce_cluster<e_count_type>>::iterator&,
                       unsigned, e_count_type);
        void make_cluster(unsigned);
        void erase_from_cl(typename sprclib::list<bomce_cluster<e_count_type>>::iterator&,
                           unsigned, e_count_type);
        double move_vertex(unsigned);
        double move_vertex(unsigned, const typename sprclib::list<bomce_cluster<
                           e_count_type>>::iterator&);
        v_distribution copy();
        std::vector<v_distribution> copy(unsigned);
        static std::vector<v_distribution> copy(v_distribution&&, unsigned);
        bomce_vdist_family<BOMCE_TEMPLATE_ARGS>& family()
        {
            return *((bomce_vdist_family<BOMCE_TEMPLATE_ARGS>*)fam);
        }
        bool next_iteration();
        void clear()
        {
            v_to_cl.clear();
            clusters.clear();
            del_index = std::move(std::queue<unsigned>());
        }

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
    std::vector<v_distribution> dist;
public:
    typedef typename std::vector<v_distribution>::iterator iterator;

    bomce_vdist_family(const graph<BOMCE_TEMPLATE_ARGS>& gr,
                       const bomce_tmpresult_type<type>& tmp_result,
                       const q_modular<BOMCE_TEMPLATE_ARGS>* _q_mod = nullptr,
                       unsigned dist_count = std::thread::hardware_concurrency()) :
        temp_graph(gr), temp_result(tmp_result),
        q_mod(_q_mod), dist(dist_count)
    {
        for (auto& i : dist)
            i.fam = this;
    }

    /*bomce_vdist_family(const graph<BOMCE_TEMPLATE_ARGS>& gr,
                       const bomce_tmpresult_type<type>& tmp_result, typename
                       bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution& dist_init,
                       const q_modular<BOMCE_TEMPLATE_ARGS>* _q_mod = nullptr,
                       unsigned dist_count = std::thread::hardware_concurrency()) :
        temp_graph(gr), temp_result(tmp_result), q_mod(_q_mod)
    {
        dist = std::move(dist_init.copy(dist_count));
    }

    bomce_vdist_family(const graph<BOMCE_TEMPLATE_ARGS>& gr,
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
        dist = std::move(dist_init.copy(count));
    }

    void init(typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
              v_distribution&& dist_init, unsigned count = std::thread::
            hardware_concurrency())
    {
        dist = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
                         v_distribution::copy(std::move(dist_init), count));
    }
    void init(unsigned count = std::thread::hardware_concurrency());
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

    ~bomce_vdist_family()
    {
        if (q_mod != nullptr)
            delete q_mod;
    }
};

BOMCE_TEMPLATE_DEF
using v_dist = typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution;

BOMCE_TEMPLATE_DEF
void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::init()
{
    del_index = std::move(std::queue<unsigned>());
    clusters.clear();
    v_to_cl.resize(family().temp_graph.size());
    loops.resize(family().temp_graph.size());
    for (auto& i : family().temp_graph)
    {
        ((bomce_tmpresult_type<type>&)family().temp_result)[i.first] = { std::make_pair(i.first,
                                                                         i.second.name()) };
        make_cluster(i.first);
        loops[i.first] = i.second.loop;
    }
}

BOMCE_TEMPLATE_DEF
void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
add_to_cl(const typename sprclib::list<bomce_cluster<
          e_count_type>>::iterator& cl, unsigned gr_index,
          e_count_type aij_sum_delt)
{
    v_to_cl[gr_index] = cl;
    cl->v_count += family().temp_result[gr_index].size();
    cl->e_count += aij_sum_delt + family().temp_graph[gr_index].loop;
    cl->in_sum += family().temp_graph[gr_index].in_d;
    cl->out_sum += family().temp_graph[gr_index].out_d;
    cl->loops += loops[gr_index];
}

BOMCE_TEMPLATE_DEF
// Для тех случаев, когда у вершины нет смежных рёбер
// Возвращает номер созданного кластера.
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
void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
erase_from_cl(typename sprclib::list<bomce_cluster<
              e_count_type>>::iterator& it, unsigned current,
              e_count_type aij_temp)
{
    if ((it->v_count -= family().temp_result[current].size()) <= 0)
    {
        del_index.push(it->id);
        clusters.erase(it);
        //it = clusters.end();
    }
    else
    {
        it->e_count -= aij_temp + family().temp_graph[current].loop;
        it->in_sum -= family().temp_graph[current].in_d;
        it->out_sum -= family().temp_graph[current].out_d;
        it->loops -= loops[current];
    }
    v_to_cl[current] = clusters.end();
}

BOMCE_TEMPLATE_DEF
double bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
move_vertex(unsigned gr_index, const typename sprclib::list<bomce_cluster<
            e_count_type>>::iterator& cl)
{
    if (v_to_cl[gr_index]->id == cl->id)
        return 0.0;
    unsigned cur_e = 0, cl_e = 0;
    temp_graph[gr_index].foreach_neighbors([&](std::pair<unsigned,
                                           edge_value_type> pr)
    {
        if (v_to_cl[pr.first]->id == v_to_cl[gr_index]->id)
            ++cur_e;
        else if (v_to_cl[pr.first]->id == cl->id)
            ++cl_e;
    });
    double ret = q_mod->operator()(gr_index, cl, cl_e) -
            q_mod->operator()(v_to_cl[gr_index], cur_e, gr_index);
    erase_from_cl(v_to_cl[gr_index], gr_index, cur_e);
    add_to_cl(cl, gr_index, cl_e);
    return ret;
}

BOMCE_TEMPLATE_DEF
template<class comp>
        // иннициализация cl(смежные кластеры)
        void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
        _Init(std::map<typename sprclib::list<bomce_cluster<
              e_count_type>>::iterator, unsigned, comp>& cl,
              const std::map<unsigned, edge_value_type>& graph_v)
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
double bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::
move_vertex(unsigned gr_index)
{
    /*
        пройтись по смежным вершинам
        получить множество смежных кластеров и смежных вершин
        определить предпочтительный для присоединения кластер
        определить предпочтительную для присоединения вершину
        выполнить либо cluster_add либо make_cluster(модулярность не должна ухудшиться),
        перед этим выполнив cluster_erase если нужно
        */
    typedef std::map<typename sprclib::list<bomce_cluster<e_count_type>>::iterator,
            unsigned, cl_cmp<e_count_type>> mp_type;
    mp_type cl; // смежный кластер на aij_temp(кол. рёбер между it и кластером)
    _Init(cl, family().temp_graph[gr_index].input);
    _Init(cl, family().temp_graph[gr_index].output);
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
    typename sprclib::list<bomce_cluster<e_count_type>>::iterator cl_it;
    typename mp_type::iterator it1 = cl.begin();
    max = { family().q_mod->operator()(gr_index, (cl_it = it1->first),
            it1->second, alpha), it1->second };
    for (++it1; it1 != cl.end(); ++it1)
    {
        typename sprclib::list<bomce_cluster<e_count_type>>::iterator it_tmp;
        std::pair<double, unsigned> tmp = { family().q_mod->operator()(gr_index,
                                            (it_tmp = it1->first), it1->second, alpha), it1->second };
        if (tmp.first > max.first)
        {
            max = tmp;
            cl_it = it_tmp;
        }
    }
    double p = family().q_mod->operator()(v_to_cl[gr_index],
                                          current_aij, gr_index, alpha);
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

BOMCE_TEMPLATE_DEF
typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution
bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::v_distribution::copy()
{
    v_dist<BOMCE_TEMPLATE_ARGS> ret(fam, del_index, clusters, v_to_cl.size());
    std::vector<unsigned> old_index(clusters.size());
    unsigned i = 0;
    for (auto& cl : clusters)
    {
        old_index[i] = cl.id;
        cl.id = i++;
    }
    std::vector<typename sprclib::list<bomce_cluster<
            e_count_type>>::iterator> new_it(clusters.size());
    i = 0;
    for (auto it = ret.clusters.begin(); it != ret.clusters.end(); ++it)
        new_it[i++] = it;
    for (i = 0; i < v_to_cl.size(); ++i)
        ret.v_to_cl[i] = new_it[v_to_cl[i]->id];
    i = 0;
    for (auto& cl : clusters)
        cl.id = old_index[i++];
    return ret;
}

BOMCE_TEMPLATE_DEF
std::vector<typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
v_distribution> bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
v_distribution::copy(typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
                     v_distribution&& v_d, unsigned copy_count)
{
    std::vector<v_dist<BOMCE_TEMPLATE_ARGS>> ret(copy_count);
    ret[0] = std::move(v_d);
    for (unsigned i = 1; i < ret.size(); ++i)
    {
        ret[i] = std::move(v_distribution(ret[0].fam, ret[0].del_index,
            ret[0].clusters, ret[0].v_to_cl.size(), ret[0].loops));
    }
    std::vector<unsigned> old_index(ret[0].clusters.size());
    unsigned i = 0;
    for (auto& cl : ret[0].clusters)
    {
        old_index[i] = cl.id;
        cl.id = i++;
    }
    std::vector<typename sprclib::list<bomce_cluster<
            e_count_type>>::iterator> new_it(ret[0].clusters.size());
    for (unsigned k = 1; k < copy_count; ++k)
    {
        i = 0;
        for (auto it = ret[k].clusters.begin();
             it != ret[k].clusters.end(); ++it)
        {
            new_it[i++] = it;
        }
        for (i = 0; i < ret[0].v_to_cl.size(); ++i)
            ret[k].v_to_cl[i] = new_it[ret[0].v_to_cl[i]->id];
    }
    i = 0;
    for (auto& cl : ret[0].clusters)
        cl.id = old_index[i++];
    return ret;
}

BOMCE_TEMPLATE_DEF
std::vector<typename bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
v_distribution> bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
v_distribution::copy(unsigned copy_count)
{
    std::vector<v_dist<BOMCE_TEMPLATE_ARGS>> ret(copy_count);
    for(auto& i : ret)
    {
        i = std::move(v_distribution(fam, del_index,
            clusters, v_to_cl.size(), loops));
    }
    std::vector<unsigned> old_index(clusters.size());
    unsigned i = 0;
    for (auto& cl : clusters)
    {
        old_index[i] = cl.id;
        cl.id = i++;
    }
    std::vector<typename sprclib::list<bomce_cluster<
            e_count_type>>::iterator> new_it(clusters.size());
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
bool bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
v_distribution::next_iteration()
{
    double contin;
    std::vector<unsigned> prep(family().temp_graph.size());
    for (unsigned i = 0; i < family().temp_graph.size(); ++i)
        prep[i] = i;
    std::mt19937 gen(time(0) + (unsigned)std::this_thread::get_id);
    bool mooved = false;
    for (alpha = 0.7; alpha <= 1.3; alpha += 0.2) // 0.222734
    {
        do
        {
            contin = 0;
            std::shuffle(prep.begin(), prep.end(), gen);
            for (auto& i : prep)
            {
                contin += this->move_vertex(i);
            }
            mooved = mooved || (contin > 1e-5);
        } while (contin > 1e-5);
    }
    return mooved;
}

BOMCE_TEMPLATE_DEF
void bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::init(unsigned count)
{
    v_distribution v_d(this);
    v_d.init();
    dist = std::move(bomce_vdist_family<BOMCE_TEMPLATE_ARGS>::
                     v_distribution::copy(std::move(v_d), count));
}

BOMCE_TEMPLATE_DEF
struct bomce_result
{
    // Сруктура связей кластеров
    graph<BOMCE_TEMPLATE_ARGS> const& connections;
    // Параметры кластеров
    std::vector<typename sprclib::list<bomce_cluster<e_count_type>>::
    iterator> const& v_to_cl;
    // Состав вешин в каждом кластере
    bomce_tmpresult_type<type> const& vertexes;
    // Список кластеров
    sprclib::list<bomce_cluster<e_count_type>> const& clusters;
    // число петель в исходном графе
    e_count_type loop_count;

    bomce_result(graph<BOMCE_TEMPLATE_ARGS> const& gr,
                 std::vector<typename sprclib::list<bomce_cluster<e_count_type>>::
                 iterator> const& _v_to_cl, bomce_tmpresult_type<
                 type> const& temp_result, sprclib::list<bomce_cluster<
                 e_count_type>> const& _clusters, e_count_type loops) :
        v_to_cl(_v_to_cl), clusters(_clusters), connections(gr),
        vertexes(temp_result), loop_count(loops)
    {}
};

enum class orientation
{
    orient,
    no_orient
};

BOMCE_TEMPLATE_DEF
class bomce
{
private:
    e_count_type edge_count, loop_count;
    graph<BOMCE_TEMPLATE_ARGS> temp_graph;
    // номер кластера на список верщин в нём
    bomce_tmpresult_type<type> temp_result;
    bomce_vdist_family<BOMCE_TEMPLATE_ARGS> v_f;
    unsigned t_count;

    void meta_graph(v_dist<BOMCE_TEMPLATE_ARGS>&);
    void _Init(unsigned);
public:

    explicit bomce(orientation ornt) : edge_count(0), loop_count(0),
        v_f(temp_graph, temp_result)
    {
        if (ornt == orientation::no_orient)
            v_f.q_mod = new q_noorient<BOMCE_TEMPLATE_ARGS>(v_f);
        else if (ornt == orientation::no_orient)
            v_f.q_mod = new q_orient<BOMCE_TEMPLATE_ARGS>(v_f);
    }

    void init(graph<BOMCE_TEMPLATE_ARGS>&& gr,
              unsigned thread_count = std::thread::hardware_concurrency())
    {
        temp_graph = std::move(gr);
        _Init(thread_count);
    }

    void init(const graph<BOMCE_TEMPLATE_ARGS>& gr,
              unsigned thread_count = std::thread::hardware_concurrency())
    {
        temp_graph = gr;
        _Init(thread_count);
    }

    bool next_iteration();

    bomce_result<BOMCE_TEMPLATE_ARGS> result() const
    {
        return bomce_result<BOMCE_TEMPLATE_ARGS>(temp_graph, v_f[0].v_to_cl,
                temp_result, v_f[0].clusters, loop_count);
    }

    virtual ~bomce() {}
};

BOMCE_TEMPLATE_DEF
void bomce<BOMCE_TEMPLATE_ARGS>::_Init(unsigned _tcount)
{
    edge_count = temp_graph.edge_count();
    loop_count = 0;
    temp_result.resize(temp_graph.size());
    t_count = _tcount;
    v_f.init(_tcount);
}

BOMCE_TEMPLATE_DEF
bool bomce<BOMCE_TEMPLATE_ARGS>::next_iteration()
{
    static std::vector<double> mods(v_f.size(), 0);
    std::vector<bool> bl(v_f.size(), false);
    parallel_for(0u, v_f.size(), 1u, [&](unsigned index, unsigned tid)
    {
        if(!v_f[index].next_iteration())
        {
            bl[index] = false;
            return;
        }
        double mod = 0;
        for (auto it = v_f[index].clusters.begin();
             it != v_f[index].clusters.end(); ++it)
        {
            mod += v_f.q_mod->operator()(it);
        }
        bl[index] = (mod > mods[index]);
        mods[index] = mod;
        //mods[index] /= 2.0*temp_graph.edge_count();
    }, t_count);
    unsigned i_mod = 0;
    for (unsigned i = 1; i < mods.size(); ++i)
    {
        if (mods[i] > mods[i_mod])
            i_mod = i;
    }
    if (bl[i_mod])
    {
        std::swap(v_f[0], v_f[i_mod]);
        std::swap(mods[0], mods[i_mod]);
        unsigned t_count = v_f.size();
        v_f.dist.resize(1);
        meta_graph(v_f[0]);
        auto v_d = std::move(v_f[0]);
        v_f.dist.clear();
        v_f.init(std::move(v_d), t_count);
    }
    else
        v_f.init(v_f[i_mod], 1);
    return bl[i_mod];
}

BOMCE_TEMPLATE_DEF
// Создаёт метаграф
void bomce<BOMCE_TEMPLATE_ARGS>::meta_graph(
        v_dist<BOMCE_TEMPLATE_ARGS>& v_d)
{
    decltype(temp_graph) meta_gr;
    unsigned t = 0;
    for (auto& i : v_d.clusters)
    {
        i.id = t;
        meta_gr.new_vertex(t++);
    }
    std::vector<e_count_type> new_loops(v_d.clusters.size(), 0);
    for (auto& v : temp_graph)
        new_loops[v_d.v_to_cl[v.first]->id] += v_d.loops[v.first];
    v_d.loops = std::move(new_loops);
    bomce_tmpresult_type<type> tmp_result(v_d.clusters.size());
    for (auto& i : temp_graph)
    {
        t = v_d.v_to_cl[i.first]->id;
        tmp_result[t].splice(tmp_result[t].end(), temp_result[i.first]);
        meta_gr.connect_by_index(t, t, i.second.loop);
        for (auto& j : i.second.output)
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
                const graph<BOMCE_TEMPLATE_ARGS>& gr)
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
            ret[i].insert(el.second);
    }
    return ret;
}
}
#endif // _bomce_H_
