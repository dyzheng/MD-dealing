#ifndef Iprofile_H
#define Iprofile_H

#include "cell.h"

class Iprofile
{
	public: 
	
	Iprofile(){};
	~Iprofile(){};

	void Routine();
	
	private:

	// calculate the ionic density profile 
	void cal( const Cell &cel_in );

	// calculate the ionic density profile
	// using gaussian smearing
	void cal_gauss( const Cell &cel_in );

};


#endif
