
#ifndef COPL_DEBUG_H
#define COPL_DEBUG_H
#define DEBUG_TO_FILE

#include <fstream>
#include <iostream>

using namespace std;

namespace copl_ip {
#ifdef DEBUG_TO_FILE
extern ofstream LOG_FILE_VARIABLE;
#define OUTPUT LOG_FILE_VARIABLE
#else
#define OUTPUT cout
#endif
}

#endif