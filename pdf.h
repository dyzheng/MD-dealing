#ifndef PDF_H
#define PDF_H

#include "cell.h"

class PDF
{
	public: 
	
	PDF(){};
	~PDF(){};

	void Routine();
	
	private:

	// calculate the pair distribution function.
	void cal();

	// calculate the paris under periodic boudary conditions.
	void periodic_pairs( const Cell &cell, double * gr, const int &option)const;

	// calculate the pairs for each (i+1)dr.
    void pairs( const Cell &cel, double *gr);

	// useless
	//void static_structure_factor( const Cell &cel, const double &rho_ion, const double *gr, double *sf) const;

	private:

	double dr; // delta r
	double rcut; // radius cutoff
	int nmesh; // number of mesh points.

	// useless.
	//double dg;
	//int ng;

};


#endif
