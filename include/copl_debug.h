
#ifndef COPL_DEBUG_H
#define COPL_DEBUG_H

#include <fstream>
#include <iostream>

using namespace std;

namespace copl_ip {
#ifndef DEBUG_TO_FILE
#define OUTPUT cout
#else
extern ofstream LOG_FILE_VARIABLE;
#define OUTPUT LOG_FILE_VARIABLE
#endif
}

#endif
