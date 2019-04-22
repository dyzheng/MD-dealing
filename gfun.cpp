#include "gfun.h"
#include <cstdlib>

bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart)
{
    string SearchName;
    bool find = false;
    if (restart)
    {
        ifs.clear();
        ifs.seekg(0);
    }
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> SearchName;
		//cout << " " << SearchName << endl;
        if ( SearchName == TargetName)
        {
            find = true;
            break;
        }
    }
    if (!find)
    {
     //   cout <<" In SCAN_BEGIN, can't find: "<<TargetName<<" block."<<endl;
    }
    return find;
}



void TITLE(const string &class_name,const string &function_name)
{
	//return;
    //cout<<" ==> "<<class_name<<"::"<<function_name<<endl;
    return;
}

void QUIT(const string &reason)
{
	cout << " ------------------------------------------ " << endl;
	cout << " Quit because : " << reason << endl;
	exit(0);
}

