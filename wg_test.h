#pragma once
#ifndef WG_TEST_H_
#define WG_TEST_H_

#include <fstream>
#include <iostream>
#include "equal_cluster_metric.h"
#include "bomce.h"

namespace webgr
{
	void test_by_gv()
	{
		bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl(orientation::no_orient);
		std::string path("C:\\testgr\\");
		std::ofstream fout("cl_count.csv");
		for (unsigned i = 0; i < 100; ++i)
		{
			auto gr = from_pajek(path + "testgr_" + std::to_string(i) + ".net");
			//cl.init(gr);
			while (cl.next_iteration())
				;
			std::ifstream fin(path + "default_" + std::to_string(i) + ".txt");
			unsigned count = 0; std::string str;
			while (getline(fin, str))
				++count;
			fin.close();
			fout << gr.size() << ';' << gr.edge_count() << ";" << count - 1 <<
				';' << cl.result().clusters.size() << ';' << clusters_from(
					path + "label_prop_" + std::to_string(i) + ".txt").size() << ';'
				<< clusters_from(path + "walktrap_" + std::to_string(i) + ".txt").size()
				<< ';' << clusters_from(path + "louvain_" + std::to_string(i) +
					".txt").size() << std::endl;
		}
		fout.close();
	}
	
	void nmi_test()
	{
		for (unsigned i = 0; i < 100; ++i)
		{
			auto first = clusters_from("C:\\testgr\\default_" +
				std::to_string(i) + ".txt"); // исходное разбиение
			auto second = clusters_from("C:\\testgr\\label_prop_" +
				std::to_string(i) + ".txt");
			equal_cluster_metric metric(first, second);
			std::ofstream out;
			out.open("C:\\testgr\\nmi2_lp.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			second = std::move(clusters_from("C:\\testgr\\walktrap_" +
				std::to_string(i) + ".txt"));
			metric.load(first, second);
			out.open("C:\\testgr\\nmi2_wt.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			second = std::move(clusters_from("C:\\testgr\\bomce_" +
				std::to_string(i) + ".txt"));
			metric.load(first, second);
			out.open("C:\\testgr\\nmi2_lv.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			{
				default_graph gr(from_pajek("C:\\testgr\\testgr_" +
					std::to_string(i) + ".net"));
                bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl1(orientation::no_orient);
                cl1.init(gr);
				while (cl1.next_iteration())
                    ;
                second = std::move(clusters_from(cl1.result()));
			}
			metric.load(first, second);
			out.open("C:\\testgr\\nmi2_my.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();
			std::cout << i << std::endl;
		}
	}

	void test0()
	{
		std::ofstream out("C:\\testgr\\my_alg_modular.txt");
		for (unsigned i = 0; i < 100; ++i)
		{
			default_graph gr(from_pajek("C:\\testgr\\testgr_" +
				std::to_string(i) + ".net"));
			bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl1(orientation::no_orient);
			double mod = 0;
			for (unsigned k = 0; k < 100; ++k)
			{
				// cl1 = std::move(bomce<BOMCE_DEFAULT_GRAPH_ARGS>(orientation::no_orient));
				//cl1.init(gr);
				while (cl1.next_iteration())
					;
				double mod1 = modular(cl1.result());
				if (mod1 > mod)
					mod = mod1;
			}
			out << mod << std::endl;
			std::cout << i << std::endl;
		}
		out.close();
	}

	void test1()
	{
		for (unsigned i = 0; i < 100; ++i)
		{
			auto first = clusters_from("C:\\testgr\\default_" +
				std::to_string(i) + ".txt"); // исходное разбиение
			auto second = clusters_from("C:\\testgr\\label_prop_" +
				std::to_string(i) + ".txt");
			equal_cluster_metric metric(first, second);
			metric.calc_g();
			metric.calc_g_();
			std::ofstream out("C:\\testgr\\f1_lp.txt", std::ios::app);
			out << metric.average_f1() << std::endl;
			out.close();
			out.open("C:\\testgr\\nmi_lp.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			second = std::move(clusters_from("C:\\testgr\\walktrap_" +
				std::to_string(i) + ".txt"));
			metric.load(first, second);
			metric.calc_g();
			metric.calc_g_();
			out.open("C:\\testgr\\f1_wt.txt", std::ios::app);
			out << metric.average_f1() << std::endl;
			out.close();
			out.open("C:\\testgr\\nmi_wt.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			second = std::move(clusters_from("C:\\testgr\\bomce_" +
				std::to_string(i) + ".txt"));
			metric.load(first, second);
			metric.calc_g();
			metric.calc_g_();
			out.open("C:\\testgr\\f1_lv.txt", std::ios::app);
			out << metric.average_f1() << std::endl;
			out.close();
			out.open("C:\\testgr\\nmi_lv.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();

			{
				default_graph gr(from_pajek("C:\\testgr\\testgr_" +
					std::to_string(i) + ".net"));
				bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl1(orientation::no_orient);
				double mod = 0;
				for (unsigned k = 0; k < 100; ++k)
				{
					// cl1 = std::move(bomce<BOMCE_DEFAULT_GRAPH_ARGS>(orientation::no_orient));
					//cl1.init(gr);
					while (cl1.next_iteration())
						;
					double mod1 = modular(cl1.result());
					if (mod1 > mod)
						mod = mod1;
				}
				out.open("C:\\testgr\\my_alg_modular.txt", std::ios::app);
				out << mod << std::endl;
				out.close();
				second = std::move(clusters_from(cl1.result()));
			}
			metric.load(first, second);
			metric.calc_g();
			metric.calc_g_();
			out.open("C:\\testgr\\f1_my.txt", std::ios::app);
			out << metric.average_f1() << std::endl;
			out.close();
			out.open("C:\\testgr\\nmi_my.txt", std::ios::app);
			out << metric.nmi() << std::endl;
			out.close();
			std::cout << i << std::endl;
		}
	}

	void calc_mods()
	{
		std::vector<std::vector<unsigned>> mods(5, std::vector<unsigned>(10, 0));
		std::ifstream fin; double mod;
		for (unsigned i = 0; i < 100; ++i)
		{
			fin.open("C:\\testgr\\default_" + std::to_string(i) + ".txt");
			fin >> mod;
			++mods[0][(unsigned)(mod * 10)];
			fin.close();
		}
		fin.open("C:\\testgr\\lblprop_mod.txt");
		while (fin >> mod)
		{
			++mods[1][(unsigned)(mod * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\my_alg_modular.txt");
		while (fin >> mod)
		{
			++mods[2][(unsigned)(mod * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\walktrap_mod.txt");
		while (fin >> mod)
		{
			++mods[3][(unsigned)(mod * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\bomce_mod.txt");
		while (fin >> mod)
		{
			++mods[4][(unsigned)(mod * 10)];
		}
		fin.close();
		std::ofstream out("test3.csv");
		for (unsigned i = 0; i < 10; ++i)
		{
			out << mods[0][i];
			for (unsigned k = 1; k < 5; ++k)
			{
				out << ';' << mods[k][i];
			}
			out << std::endl;
		}
		out.close();
	}

	void calc_f1()
	{
		std::vector<std::vector<unsigned>> f1s(4, std::vector<unsigned>(11, 0));
		std::ifstream fin; double val;
		fin.open("C:\\testgr\\f1_lp.txt");
		while (fin >> val)
		{
			++f1s[0][(unsigned)(val * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\f1_my.txt");
		while (fin >> val)
		{
			++f1s[1][(unsigned)(val * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\f1_wt.txt");
		while (fin >> val)
		{
			++f1s[2][(unsigned)(val * 10)];
		}
		fin.close();
		fin.open("C:\\testgr\\f1_lv.txt");
		while (fin >> val)
		{
			++f1s[3][(unsigned)(val * 10)];
		}
		fin.close();
		std::ofstream out("f1.csv");
		for (unsigned i = 0; i < 11; ++i)
		{
			out << i / 10.0;
			for (unsigned k = 0; k < 4; ++k)
			{
				out << ';' << f1s[k][i];
			}
			out << std::endl;
		}
		out.close();
	}

	void calc_nmi()
	{
		std::vector<std::vector<unsigned>> nmi(4, std::vector<unsigned>(10, 0));
		std::ifstream fin; double val;
		fin.open("C:\\testgr\\nmi_lp.txt");
		while (fin >> val)
		{
			++nmi[0][(unsigned)(val * 1)];
		}
		fin.close();
		fin.open("C:\\testgr\\nmi_my.txt");
		while (fin >> val)
		{
			++nmi[1][(unsigned)(val * 1)];
		}
		fin.close();
		fin.open("C:\\testgr\\nmi_wt.txt");
		while (fin >> val)
		{
			++nmi[2][(unsigned)(val * 1)];
		}
		fin.close();
		fin.open("C:\\testgr\\nmi_lv.txt");
		while (fin >> val)
		{
			++nmi[3][(unsigned)(val * 1)];
		}
		fin.close();
		std::ofstream out("nmi.csv");
		for (unsigned i = 0; i < 10; ++i)
		{
			out << i / 1.0;
			for (unsigned k = 0; k < 4; ++k)
			{
				out << ';' << nmi[k][i];
			}
			out << std::endl;
		}
		out.close();
	}

	void characteristics()
	{
		default_graph gr(from_pajek("Untitled.net"));
		std::cout << "Red: " << gr.size() << std::ends << gr.edge_count() << std::endl;
		std::cout << assort(gr) << std::ends << cl_coef(gr) << std::endl;
		std::ofstream fout;
		auto h = diameter(gr);
		fout.open("Untitled_diam.csv");
		for (unsigned i = 1; i < h.size(); ++i)
			fout << i << ';' << h[i] << std::endl;
		fout.close();
		
		gr = std::move(from_pajek("buddahs_only.net"));
		std::cout << "buddahs_only: " << gr.size() << std::ends << gr.edge_count() << std::endl;
		std::cout << assort(gr) << std::ends << cl_coef(gr) << std::endl;
		h = std::move(diameter(gr));
		fout.open("buddahs_only_diam.csv");
		for (unsigned i = 1; i < h.size(); ++i)
			fout << i << ';' << h[i] << std::endl;
		fout.close();

		gr = std::move(from_pajek("Itigelov_subscribers.net"));
		std::cout << "it sub: " << gr.size() << std::ends << gr.edge_count() << std::endl;
		std::cout << assort(gr) << std::ends << cl_coef(gr) << std::endl;
		h = std::move(diameter(gr));
		fout.open("Itigelov_subscribers_diam.csv");
		for (unsigned i = 1; i < h.size(); ++i)
			fout << i << ';' << h[i] << std::endl;
		fout.close();
		
		gr = std::move(from_pajek("writers.net"));
		std::cout << "writers: " << gr.size() << std::ends << gr.edge_count() << std::endl;
		std::cout << assort(gr) << std::ends << cl_coef(gr) << std::endl;
		h = std::move(diameter(gr));
		fout.open("writers_diam.csv");
		for (unsigned i = 1; i < h.size(); ++i)
			fout << i << ';' << h[i] << std::endl;
		fout.close();
		
	}

	void pr()
	{
		default_graph gr(from_pajek("Untitled.net", false));
		std::ofstream fout;
		auto h = page_rank(gr);
		fout.open("Untitled_pr.csv");
		std::vector<unsigned> avg(11, 0);
		for (auto& i : h)
		{
			fout << i.first << ';' << i.second << std::endl;
			++avg[(unsigned)(i.second * 10)];
		}
		fout.close();
		fout.open("Untitled_pr_avg.csv");
		for (unsigned i = 0; i < avg.size(); ++i)
			fout << i / 10.0 << ';' << avg[i] << std::endl;
		fout.close();

		gr = std::move(from_pajek("buddahs_only.net"));
		std::cout << gr.size() << std::ends << gr.edge_count() << std::endl;
		h = std::move(page_rank(gr));
		fout.open("buddahs_only_pr.csv");
		avg = std::vector<unsigned>(101, 0);
		for (auto& i : h)
		{
			fout << i.first << ';' << i.second << std::endl;
			++avg[(unsigned)(i.second * 100)];
		}
		fout.close();
		fout.open("buddahs_only_pr_avg.csv");
		for (unsigned i = 0; i < avg.size(); ++i)
			fout << i / 100.0 << ';' << avg[i] << std::endl;
		fout.close();

		gr = std::move(from_pajek("Itigelov_subscribers.net"));
		h = std::move(page_rank(gr));
		fout.open("Itigelov_subscribers_pr.csv");
		for (auto& i : avg)
			i = 0;
		for (auto& i : h)
		{
			fout << i.first << ';' << i.second << std::endl;
			++avg[(unsigned)(i.second * 100)];
		}
		fout.close();
		fout.open("Itigelov_subscribers_pr_avg.csv");
		for (unsigned i = 0; i < avg.size(); ++i)
			fout << i / 100.0 << ';' << avg[i] << std::endl;
		fout.close();

		gr = std::move(from_pajek("writers.net"));
		h = std::move(page_rank(gr));
		fout.open("writers_pr.csv");
		avg = std::move(std::vector<unsigned>(11, 0));
		for (auto& i : h)
		{
			fout << i.first << ';' << i.second << std::endl;
			++avg[(unsigned)(i.second * 10)];
		}
		fout.close();
		fout.open("writers_pr_avg.csv");
		for (unsigned i = 0; i < avg.size(); ++i)
			fout << i / 10.0 << ';' << avg[i] << std::endl;
		fout.close();
	}
}

void all_tests()
{

}

#endif
