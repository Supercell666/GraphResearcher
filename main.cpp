<<<<<<< HEAD
// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <iostream>
#include "bomce.h"
#include "bomce_utils.h"
=======
//  0.7315093
#include <iostream>
#include "wg_test.h"
>>>>>>> a2c79c1ed8cb811c9c12aec25dbae1cb0a5d6b14

using namespace std;
using namespace webgr;

int main(int argc, char** argv)
{
<<<<<<< HEAD
	string resultfile = (argc < 2) ? "buddahs_only.net" : argv[1];
	default_graph gr(webgr::from_file(resultfile));
	for (auto it = gr.begin(); it != gr.end(); ++it)
		cout << it->first << endl;
	if (gr.size() == 0)
	{
		std::cerr << "Graph is empty\n";
		return 1;
	}
	float mod = 0;
	unsigned pos = resultfile.rfind('.');
	resultfile.replace(pos, resultfile.size() - pos, ".csv");
	ifstream fin;
	try
	{
		fin.open(resultfile);
		fin >> mod;
		fin.close();
	}
	catch (exception& exc)
	{
		cerr << exc.what();
		fin.close();
	}
	vector<float> p = { 0.7f, 1.2f, 0.2f };
	try
	{
		fin.open("config.txt");
		fin >> p[0] >> p[1] >> p[2];
		fin.close();
	}
	catch (exception& exc)
	{
		cerr << exc.what();
		fin.close();
	}
	cout << "Begining from Q = " << mod << endl;
	bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl(orientation::no_orient);
	for (;;)
	{
		cl.init(gr, p, 2);
		while (cl.next_iteration())
			;
		auto md = cl.result().modularity;
		if (md > mod)
		{
			mod = md;
			std::cout << "\t" << md << std::ends << "Writing... ";
			cl.result().contents(resultfile);
			std::cout << "Done!\n";
		}
		else
			std::cout << md << std::endl;
	}
=======
    auto gr(from_pajek("/home/user/GraphResearcher/buddahs_only.net"));
    cout << gr.size() << endl;
    bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl(orientation::no_orient);
	cl.init(gr);
	while (cl.next_iteration())
		;
	cout << cl.result().clusters.size();
	system("pause");
>>>>>>> a2c79c1ed8cb811c9c12aec25dbae1cb0a5d6b14
	return 0;
}
