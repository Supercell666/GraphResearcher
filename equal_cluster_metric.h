#pragma once
#ifndef EQUAL_CLUSTER_METRIC
#define EQUAL_CLUSTER_METRIC
#include "graphutility.h"

namespace webgr
{
	class equal_cluster_metric
	{
	private:
		std::vector<std::vector<unsigned>> inrsct_matrix;
		std::vector<std::set<unsigned>>& clusters1, clusters2;
		std::vector<unsigned> g, g_;
		unsigned gr_size;

		double precision(unsigned first, unsigned second)
		{
			return (double)inrsct_matrix[first][second] / clusters2[second].size();
		}
		double recall(unsigned first, unsigned second)
		{
			return (double)inrsct_matrix[first][second] / clusters1[first].size();
		}
		double f1(unsigned first, unsigned second)
		{
			double p = precision(first, second), r = recall(first, second);
			return (fabs(p + r) > 1e-9) ? (2 * p*r) / (p + r) : 0;
		}
		unsigned get_g(unsigned arg)
		{
			unsigned max_j = 0;
			double max = f1(arg, 0);
			for (unsigned j = 1; j < clusters2.size(); ++j)
			{
				double t = f1(arg, j);
				if (t > max)
				{
					max = t;
					max_j = j;
				}
			}
			return max_j;
		}
		unsigned get_g_(unsigned arg)
		{
			unsigned max_i = 0;
			double max = f1(0, arg);
			for (unsigned i = 1; i < clusters1.size(); ++i)
			{
				double t = f1(i, arg);
				if (t > max)
				{
					max = t;
					max_i = i;
				}
			}
			return max_i;
		}
		double calc_h_first()
		{
			double s = 0;
			for (unsigned i = 0; i < clusters1.size(); ++i)
			{
				double t = (double)clusters1[i].size() / gr_size;
				s += t*log2(t);
			}
			return -s;
		}
		double calc_h_second()
		{
			double s = 0;
			for (unsigned i = 0; i < clusters2.size(); ++i)
			{
				double t = (double)clusters2[i].size() / gr_size;
				s += t*log2(t);
			}
			return -s;
		}
		double calc_h(unsigned);
		double calc_mi(unsigned);
		unsigned calc_matrix()
		{
			std::vector<unsigned> sizes(std::thread::hardware_concurrency(), 0);
			parallel_for(0u, (unsigned)clusters1.size(), 1u, [&](unsigned i, unsigned index)
			{
				sizes[index] += clusters1[i].size();
				for (unsigned k = 0; k < clusters2.size(); ++k)
				{
					foreach_equal_clusters(clusters1[i], clusters2[k], [&](
						unsigned first, unsigned second)
					{
						++inrsct_matrix[i][k];
					});
				}
			}, sizes.size());
			unsigned ret = 0;
			for (auto& i : sizes)
				ret += i;
			return ret;
		}
	public:
		equal_cluster_metric(const std::vector<std::set<unsigned>>& _clusters1,
			const std::vector<std::set<unsigned>>& _clusters2) :
			inrsct_matrix(_clusters1.size(), std::vector<unsigned>(_clusters2.size(), 0)),
			clusters1((std::vector<std::set<unsigned>>&)_clusters1),
			clusters2((std::vector<std::set<unsigned>>&)_clusters2), gr_size(0)
		{
			gr_size = calc_matrix();
		}
		equal_cluster_metric(const std::vector<std::set<unsigned>>& _clusters1,
			const std::vector<std::set<unsigned>>& _clusters2, unsigned _gr_size) :
			inrsct_matrix(_clusters1.size(), std::vector<unsigned>(_clusters2.size(), 0)),
			clusters1((std::vector<std::set<unsigned>>&)_clusters1),
			clusters2((std::vector<std::set<unsigned>>&)_clusters2), gr_size(_gr_size)
		{
			parallel_for(0u, (unsigned)clusters1.size(), 1u, [&](unsigned i, unsigned index)
			{
				for (unsigned k = 0; k < clusters2.size(); ++k)
				{
					foreach_equal_clusters(clusters1[i], clusters2[k], [&](
						unsigned first, unsigned second)
					{
						++inrsct_matrix[i][k];
					});
				}
			});
		}
		equal_cluster_metric(const equal_cluster_metric& val) :clusters1(val.clusters1),
			clusters2(val.clusters2), inrsct_matrix(val.inrsct_matrix), gr_size(val.gr_size) {}
		equal_cluster_metric(equal_cluster_metric&& val) :clusters1(val.clusters1),
			clusters2(val.clusters2), inrsct_matrix(std::move(val.inrsct_matrix)),
			gr_size(val.gr_size) {}
		equal_cluster_metric& operator=(const equal_cluster_metric& val)
		{
			if (this == &val)
				return *this;
			clusters1 = val.clusters1;
			clusters2=val.clusters2;
			inrsct_matrix = val.inrsct_matrix;
			gr_size = val.gr_size;
			return *this;
		}
		equal_cluster_metric& operator=(equal_cluster_metric&& val)
		{
			if (this == &val)
				return *this;
			clusters1 = std::move(val.clusters1);
			clusters2 = std::move(val.clusters2);
			inrsct_matrix = std::move(val.inrsct_matrix);
			gr_size = val.gr_size;
			return  *this;
		}

		void load(const std::vector<std::set<unsigned>> & _clusters1,
			const std::vector<std::set<unsigned>>& _clusters2)
		{
			clusters1 = _clusters1;
			clusters2 = _clusters2;
			inrsct_matrix.clear();
			inrsct_matrix = std::move(std::vector<std::vector<unsigned>>(clusters1.size(),
				std::vector<unsigned>(clusters2.size(), 0)));
			gr_size = calc_matrix();
			g.clear();
			g_.clear();
		}
		void calc_g()
		{
			g.resize(clusters1.size());
			parallel_for(0u, (unsigned)clusters1.size(), 1u, [this](unsigned i, unsigned index)
			{
				g[i] = get_g(i);
			});
		}
		void calc_g_()
		{
			g_.resize(clusters2.size());
			parallel_for(0u, (unsigned)clusters2.size(), 1u, [this](unsigned j, unsigned index)
			{
				g_[j] = get_g_(j);
			});
		}
		double average_f1();
		double average_recall();
		double average_precision();
		double nmi();
	};

	double equal_cluster_metric::average_f1()
	{
		double s1 = 0, s2 = 0;
        auto t=std::thread([&s1, this]()
		{
			for (unsigned i = 0; i < clusters1.size(); ++i)
				s1 += f1(i, g[i]);
		});
		for (unsigned j = 0; j < clusters2.size(); ++j)
			s2 += f1(g_[j], j);
		t.join();
		return (s1 / clusters1.size() + s2 / clusters2.size()) / 2;
	}

	double equal_cluster_metric::average_recall()
	{
		double s = 0;
		for (unsigned i = 0; i < clusters1.size(); ++i)
			s += recall(i, g[i]);
		return s / clusters1.size();
	}

	double equal_cluster_metric::average_precision()
	{
		double s = 0;
		for (unsigned i = 0; i < clusters1.size(); ++i)
			s += precision(i, g[i]);
		return s / clusters1.size();
	}

	double equal_cluster_metric::calc_h(unsigned th_count =
		std::thread::hardware_concurrency() * 2)
	{
		std::vector<double> sum(th_count, 0);
		parallel_for(0u, (unsigned)clusters1.size(), 1u, [this, &sum](
			unsigned i, unsigned index)
		{
			double s = 0;
			for (unsigned j = 0; j < clusters2.size(); ++j)
			{
				double t = (double)inrsct_matrix[i][j] / gr_size;
				s += (fabs(t) != 1e-12) ? t*log2(t) : 0.0;
			}
			sum[index] += s;
		}, th_count);
		double s = 0;
		for (auto& i : sum)
			s += i;
		return -s;
	}

	double equal_cluster_metric::calc_mi(unsigned th_count =
		std::thread::hardware_concurrency() * 2)
	{
		std::vector<double> sum(th_count, 0);
		parallel_for(0u, (unsigned)clusters1.size(), 1u, [this, &sum](
			unsigned i, unsigned index)
		{
			double s = 0;
			for (unsigned j = 0; j < clusters2.size(); ++j)
			{
				s += ((double)inrsct_matrix[i][j] / gr_size)*
					log2(((double)inrsct_matrix[i][j] * gr_size) / (
						clusters1[i].size()*clusters2[j].size()));
			}
			sum[index] += s;
		}, th_count);
		double s = 0;
		for (auto& i : sum)
			s += i;
		return s;
	}

	double equal_cluster_metric::nmi()
	{
		double h1, h2, h;
		std::thread t1([this, &h1]() {
			h1 = calc_h_first();
		});
		std::thread t2([this, &h2]() {
			h2 = calc_h_second();
		});
		h = calc_h(std::thread::hardware_concurrency());
		t1.join();
		t2.join();
		return std::abs((h1 + h2 - h) / ((h1*h2 > 0) ?
			sqrt(h1*h2) : sqrt(h1*h2 + 1e-16)));
	}
}
#endif
