#ifndef ATOMS_H
#define ATOMS_H

#include "gfun.h"
#include "vec3.h"

// each element containing a set of data
class Atoms
{
	public:

	Atoms();
	~Atoms();

	public:

	string id;//atom id.
	string pot_file;//pseudopotential file name.
	int na; //number of atoms.
	Vector3<double>* pos; //atom positions in Cartesian coordinates.
	Vector3<double>* posd; //atom positions in Direct coordinates.
	Vector3<double>* vel; //atom velocities.

	public: 

	// read the 'fractional' coordinates of atoms
	// or the 'cartesian' coordinates of atoms.
	void read_pos(ifstream &ifs,bool frac);
	void read_pos_2(ifstream &ifs,bool frac);
	void read_pos_3(ifstream &ifs); // mohan add 2015-07-24 for ABACUS
	void read_vel(ifstream &ifs);
	void clear();

};

#endif
