//  0.7315093
#include <iostream>
#include "wg_test.h"

using namespace std;
using namespace webgr;

int main(int argc, char** argv)
{
    auto gr(from_pajek("/home/user/GraphResearcher/buddahs_only.net"));
    cout << gr.size() << endl;
    bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl(orientation::no_orient);
	cl.init(gr);
	while (cl.next_iteration())
		;
	cout << cl.result().clusters.size();
	system("pause");
	return 0;
}
