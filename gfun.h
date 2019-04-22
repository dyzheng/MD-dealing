#ifndef GFUN_H
#define GFUN_H

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <cassert>
#include <ctime>

using namespace std;

//==========================================================
// GLOBAL FUNCTION :
//==========================================================
template <class T>
static void READ_VALUE(ifstream &ifs, T &v)
{
    ifs >> v;
    ifs.ignore(150, '\n');
    return;
}

bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart=1);
// Mohan warning : the last term can't be written as const bool &restart,
// I don't know why.



void TITLE(const string &class_name,const string &function_name);

void QUIT(const string &name);

template <class T>
inline void ZEROS(T a, const int &size)
{
	assert(size>0);
	for(int i=0; i<size; ++i) a[i]=0;
}

#endif
