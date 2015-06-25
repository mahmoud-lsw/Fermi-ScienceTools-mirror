/**
 * @file gtburst.cxx
 * @brief Driver for running gtburst.py
 *
 * @author G. Vianello <giacomov@slac.stanford.edu>
 *
 */

#include <cstdlib>
#include <sstream>

int main(int iargc, char * argv[]) {
   std::ostringstream command;
   command << "python $FERMI_DIR/lib/python/gtburst.py";
   std::system(command.str().c_str());
}
