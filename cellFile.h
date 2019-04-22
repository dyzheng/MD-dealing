#ifndef CELLFILE_H
#define CELLFILE_H

#include "gfun.h"
#include "vec3.h"
#include "cell.h"

// each file containing a set of type of atoms.
class CellFile : public Cell 
{
	public: 
	
	CellFile(){};
	~CellFile(){};

	static bool ReadGeometry( Cell &cel );
	static bool CheckGeometry( Cell &cel );
	static void WriteGeometry( Cell &cel, bool cartesian=false );

	static bool ReadVelocity( Cell &cel );

	private:

	// set this function as static because
	// not each cell need to read the geometry
//	static bool ReadGeometry_PROFESS( Cell &cel );
	static bool ReadGeometry_VASP( Cell &cel );
	static bool ReadGeometry_QE( Cell &cel );
	static bool ReadGeometry_LAMMPS( Cell &cel );
	static bool ReadGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24

//	static bool CheckGeometry_PROFESS( Cell &cel );
	static bool CheckGeometry_VASP( Cell &cel );
	static bool CheckGeometry_QE( Cell &cel );
	static bool CheckGeometry_LAMMPS( Cell &cel );
	static bool CheckGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24

//	static void WriteGeometry_PROFESS( Cell &cel, bool cartesian );
	static void WriteGeometry_VASP( Cell &cel, bool cartesian );
	static void WriteGeometry_QE( Cell &cel );
	static void WriteGeometry_LAMMPS( Cell &cel ); // mohan add 2013-09-20
	static void WriteGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24
	static void	WriteGeometry_XYZ( Cell &cel );

	static bool ReadVelocity_PROFESS( Cell &cel );
	static bool ReadVelocity_VASP( Cell &cel );

};


#endif
