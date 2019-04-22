#include "cellFile.h"
#include "input.h"



bool CellFile::CheckGeometry( Cell &cel )
{
	TITLE("CellFile","CheckGeometry");

	// (1) read in geometry file from PROFESS
	if(INPUT.geo_in_type=="PROFESS")
	{
	//	return CheckGeometry_PROFESS(cel);
	}
	// (2) read in geometry file from VASP
	else if(INPUT.geo_in_type=="VASP")
	{
		return CheckGeometry_VASP(cel);
	}
	// (3) read in geometry file from Quantum espresso.
	else if(INPUT.geo_in_type=="QE")
	{
//		CellQE::CheckGeometry(cel_in3);
	}
	// (4) read in geometry file from LAMMPS
	else if(INPUT.geo_in_type=="LAMMPS")
	{
		return CheckGeometry_LAMMPS(cel);
	}
	else if(INPUT.geo_in_type=="ABACUS")
	{
		return CheckGeometry_ABACUS(cel);
	}
	else
	{
		cout << " Error here: CheckGeometry." << endl;
		cout << " geo_in_type = " << INPUT.geo_in_type << endl;
		exit(0);
	}
	return false;
}	




bool CellFile::ReadGeometry( Cell &cel )
{
	TITLE("CellFile","ReadGeometry");

	// (1) read in geometry file from PROFESS
	if(INPUT.geo_in_type=="PROFESS")
	{
	//	return ReadGeometry_PROFESS(cel);
	}
	// (2) read in geometry file from VASP
	else if(INPUT.geo_in_type=="VASP")
	{
		return ReadGeometry_VASP(cel);
	}
	// (3) read in geometry file from Quantum espresso.
	else if(INPUT.geo_in_type=="QE")
	{
//		CellQE::ReadGeometry(cel_in3);
	}
	// (4) read in geometry file from LAMMPS.
	else if(INPUT.geo_in_type=="LAMMPS")
	{
		return ReadGeometry_LAMMPS(cel);
	}
	else if(INPUT.geo_in_type=="ABACUS")
	{
		return ReadGeometry_ABACUS(cel);
	}
	else
	{
		cout << " Error here: ReadGeometry." << endl;
		cout << " geo_in_type = " << INPUT.geo_in_type << endl;
		exit(0);
	}
	return false;
}	

void CellFile::WriteGeometry( Cell &cel, bool cartesian )
{
	ofstream ofs( INPUT.geo_out.c_str() );

	if(INPUT.geo_out_type=="PROFESS")
	{
		cout << " write geometry for PROFESS " << endl;
	//	WriteGeometry_PROFESS( cel, cartesian );
	}
	else if(INPUT.geo_out_type=="VASP")
	{
		cout << " write geometry for vasp " << endl;
        WriteGeometry_VASP( cel, cartesian );
		WriteGeometry_XYZ( cel );
	}
	else if(INPUT.geo_out_type=="LAMMPS")
	{
		WriteGeometry_LAMMPS( cel );
	}
	else
	{
		cout<< " Error here: WriteGeometry." << endl;
		exit(0);
	}
	return;
}


void CellFile::WriteGeometry_XYZ( Cell &cel )
{
	TITLE("CellFile","WriteGeometry_XYZ");
	ofstream ofs("FinalGeometry.xyz");
	
	ofs << cel.nat << endl;
	ofs << cel.system_name << endl;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		cout << " type=" << it+1 << endl;
		cout << " natom=" << cel.atom[it].na << endl;
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			ofs << cel.atom[it].id 
			<< " " << cel.atom[it].pos[ia].x
			<< " " << cel.atom[it].pos[ia].y
			<< " " << cel.atom[it].pos[ia].z << endl;
		}
	}

	ofs.close();
	return;
}



bool CellFile::ReadVelocity( Cell &cel )
{
	TITLE("CellFile","ReadVelocity");

	// (1) read in velocity file from PROFESS
	if(INPUT.velcor_in_type=="PROFESS")
	{
//		return ReadVelocity_PROFESS(cel);
	}
	// (2) read in velocity file from VASP
	else if(INPUT.velcor_in_type=="VASP")
	{
	}
	// (3) read in velocity file from Quantum espresso.
	else if(INPUT.velcor_in_type=="QE")
	{
	}
	else
	{
		cout << " Error in CellFile::ReadVelocity, please indiate velcor_in_type" << endl;
		exit(0);
	}
	return false;
}	


