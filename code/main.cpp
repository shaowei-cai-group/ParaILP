#include "Master.hpp"

int main(int argc, char **argv)
{

    INIT_ARGS

    __global_paras.print_change();
    
    Master *master = new Master();
    master->Run();
    delete (master);
    return 0;
}