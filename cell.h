#ifndef CELL_H
#define CELL_H

#include "gfun.h"
#include "vec3.h"
#include "atoms.h"

// each file containing a set of type of atoms.
class Cell
{
	public: 
	
	Cell(){};
	~Cell(){};

	public:

	Atoms *atom; // atom class
	string coordinate; //which type of coordinate
	Vector3<double> a1,a2,a3; // lattice vectors
	double volume; // volume of cell
	string file_name;
	string system_name;
	int nat; // total atom number
	int ntype; //needed if we add new species to it.

	public:

	// transform the coordinates from direct to cartesian
	void direct2cartesian(const int &it, const int &i);

	// transform the coordinates from cartesian to direct 
	void cartesian2direct(const int &it, const int &i);	

	void add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		double &x, double &y, double &z) const ;

	void add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		float &x, float &y, float &z) const ;


	private:
};


#endif
