//  0.7315093

#include "wg_test.h"

using namespace std;
using namespace webgr;

int main(int argc, char** argv)
{
    bomce<BOMCE_DEFAULT_GRAPH_ARGS> cl(orientation::no_orient);
    cl.init(from_pajek("buddahs_only.net"));
    while(cl.next_iteration())
        ;
    std::cout << cl.result().clusters.size();
    //system("pause");
    return 0;
}
