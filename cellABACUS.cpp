#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_ABACUS( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_ABACUS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_ABACUS( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_ABACUS");
	const int ntype = INPUT.ntype;

	// (1) open the file.
	stringstream ss;
	ss << INPUT.geo_directory << "/";
	ss << "md_pos_";
	ss << cel.file_name;
	ss << ".cif";
	cout << " ReadGeometry : " << ss.str() << endl;


	ifstream ifs(ss.str().c_str());

	if(!ifs) return false;
	// many files does not exist, so we don't print
	// every file's name.
	cout << " File name is " << ss.str() << endl;


	cel.atom = new Atoms[ntype];
	getline(ifs, cel.system_name);
//	cout << " Name is " << cel.system_name << endl;

	string useless;
	getline(ifs, useless);
//	cout << useless << endl;
	getline(ifs, useless);
//	cout << useless << endl;
	getline(ifs, useless);
//	cout << useless << endl;

	ifs >> useless >> cel.a1.x;
	ifs >> useless >> cel.a2.y;
	ifs >> useless >> cel.a3.z;

	cel.a1.y = 0.0;
	cel.a1.z = 0.0;
	cel.a2.x = 0.0;
	cel.a2.z = 0.0;
	cel.a3.x = 0.0;
	cel.a3.y = 0.0;

	cout << " Cell: " << endl;
	cout << " " << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	cout << " " << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	cout << " " << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

	// (3) calculate the volume of the cell.
	cel.volume = cel.a1.x*cel.a2.y*cel.a3.z + cel.a1.y*cel.a2.z*cel.a3.x + cel.a1.z*cel.a2.x*cel.a3.y -
	  cel.a1.x*cel.a2.z*cel.a3.y - cel.a1.y*cel.a2.x*cel.a3.z - cel.a1.z*cel.a2.y*cel.a3.x;

	cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;

	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
	getline(ifs, useless);
//	cout << useless << endl;

	// (4) read in atom species and the pseudopotential file.
	cel.nat = INPUT.natom;
	if(ntype==1)
	{
		cel.atom[0].na=INPUT.natom;
	}
	else if(ntype==2)
	{
		cel.atom[0].na=INPUT.atom[0];
		cel.atom[1].na=INPUT.atom[1];
	}
        else if(ntype==3)
        {
                cel.atom[0].na=INPUT.atom[0];
                cel.atom[1].na=INPUT.atom[1];
                cel.atom[2].na=INPUT.atom[2];
        }
        else if(ntype==4)
        {
                cel.atom[0].na=INPUT.atom[0];
                cel.atom[1].na=INPUT.atom[1];
                cel.atom[2].na=INPUT.atom[2];
                cel.atom[3].na=INPUT.atom[3];
        }
	cout << " Total atom number is " << cel.nat << endl;

	cel.coordinate = "Direct";

	// (5) read in the atom positions.
	Vector3<double> add1,add2,add3;
	//cout << " These are direct coordinates" << endl;
	for(int it=0; it<ntype; ++it)
	{
		cel.atom[it].read_pos_3(ifs);
	    for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
		{
	  		cel.direct2cartesian(it, ia2);
		}
	}
	ifs.close();	

//	cout << " quit reading ABACUS. " << endl;
//	exit(0);


	return true;
}



void CellFile::WriteGeometry_ABACUS(Cell &cel)
{
	/*
	TITLE("CellFile","WriteGeometry_ABACUS");
	cout << " wirte out geometry for vasp." << endl;
	if( cel.ntype != INPUT.ntype )
	{
		cout << " Be careful the printted numbe of species is not equal to the input 'ntype' " << endl;
		cout << " now number of species : " << cel.ntype << endl;
		cout << " ntype from input file : " << INPUT.ntype << endl;
	}
		

    // (2.1) need the number of types of elements.
	
	ofstream ofs(INPUT.geo_out.c_str());

	ofs << "VASP_POSCAR_FORMAT" << endl;
	ofs << "1" << endl; // scaling factor
	ofs << setprecision(16);
	ofs << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	ofs << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	ofs << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

	// name of the element
	for(int it=0; it<cel.ntype; ++it)
	{
		if( cel.atom[it].id == "" ) 
		{
			cel.atom[it].id = "XX";
		}
		ofs << cel.atom[it].id << " ";
	}
	ofs << endl;

	// atom number
	for(int it=0; it<cel.ntype; ++it)
	{
		ofs << cel.atom[it].na << " ";
	}
	ofs << endl;

	ofs << "Selected Dynamics" << endl;

	if(!cartesian)
	{
		ofs << "Direct" << endl;
		// atom position
		for(int it=0; it<cel.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << cel.atom[it].posd[ia].x
					<< " " << cel.atom[it].posd[ia].y
					<< " " << cel.atom[it].posd[ia].z 
					<< " T T T" 
					<< endl;
			}
		}
	}
	else
	{
		ofs << "Cartesian" << endl;
		// atom position
		for(int it=0; it<cel.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << cel.atom[it].pos[ia].x
					<< " " << cel.atom[it].pos[ia].y
					<< " " << cel.atom[it].pos[ia].z 
					<< " T T T" 
					<< endl;
			}
		}
	}

	ofs.close();
	return;
	*/
}

