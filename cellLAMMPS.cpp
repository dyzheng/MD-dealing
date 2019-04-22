#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_LAMMPS( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_LAMMPS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_LAMMPS( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_LAMMPS");
	const int ntype = INPUT.ntype;


	cout << "BE CAREFUL !!!!!!!!! ONLY WORK FOR ONE ELEMENT NOW." << endl;

	// (1) open the file.
	stringstream ss;
	ss << INPUT.geo_directory << "/";
	ss << cel.file_name;
	cout << " ReadGeometry : " << ss.str() << endl;


	ifstream ifs(ss.str().c_str());

	if(!ifs) return false;
	// many files does not exist, so we don't print
	// every file's name.
	cout << " File name is " << ss.str() << endl;

	cel.atom = new Atoms[ntype];
	getline(ifs, cel.system_name);
	cout << " Name is " << cel.system_name << endl;

	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		READ_VALUE(ifs, cel.atom[it].na);
		cel.nat+=cel.atom[it].na;

		cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}
	cout << " Total atom number is " << cel.nat << endl;

	int tmp_ntype = 0;
	READ_VALUE(ifs, tmp_ntype);

	// (2) read in cell
	double x0, x1;
	double y0, y1;
	double z0, z1;
	double xy=0.00; // mohan add 2015-06-15
	double xz=0.00;
	double yz=0.00;

	ifs >> x0;
	READ_VALUE(ifs, x1);
	ifs >> y0;
	READ_VALUE(ifs, y1);
	ifs >> z0;
    READ_VALUE(ifs, z1);

	// mohan add 2015-06-15
	if( INPUT.triclinic == 1 ) 
	{
		ifs >> xy >> xz;
		READ_VALUE(ifs, yz);
	}

	cel.a1.x = x1-x0; cel.a1.y = 0.00;  cel.a1.z = 0.00;
	cel.a2.x = xy;    cel.a2.y = y1-y0; cel.a2.z = 0.00;
	cel.a3.x = xz;    cel.a3.y = yz;    cel.a3.z = z1-z0;

	cout << " Cell: " << endl;
	cout << " " << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	cout << " " << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	cout << " " << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;


	// (3) calculate the volume of the cell.
	cel.volume = cel.a1.x*cel.a2.y*cel.a3.z + cel.a1.y*cel.a2.z*cel.a3.x + cel.a1.z*cel.a2.x*cel.a3.y -
	  cel.a1.x*cel.a2.z*cel.a3.y - cel.a1.y*cel.a2.x*cel.a3.z - cel.a1.z*cel.a2.y*cel.a3.x;

	cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;


	// (4) read in atom species and the pseudopotential file.
	for(int it=0; it<ntype; ++it)
	{
		stringstream ss;
		ss << "ELE" << it+1;
		cel.atom[it].id = ss.str();
	}

	string tmp_name;
	READ_VALUE(ifs, tmp_name);


	for(int it=0; it<ntype; ++it)
	{
		cel.atom[it].pos = new Vector3<double>[cel.atom[it].na];
		cel.atom[it].posd = new Vector3<double>[cel.atom[it].na];
	}

	cout << " This is Cartesian coordinates" << endl;
	bool frac = false;
	int atom_index;
	int type_index;
	for(int ia=0; ia<cel.nat; ++ia)
	{
		ifs >> atom_index >> type_index;
		--atom_index; // because this starts from 1
		--type_index; // because this starts from 1
		cout << "ia=" << ia << " atom_index=" << atom_index << " type_index=" << type_index << endl;
		ifs >> cel.atom[type_index].pos[atom_index].x 
		>> cel.atom[type_index].pos[atom_index].y 
		>> cel.atom[type_index].pos[atom_index].z;
		cel.cartesian2direct(type_index, atom_index);
	}


	return true;
}


void CellFile::WriteGeometry_LAMMPS( Cell &cel )
{
	TITLE("CellFile","WriteGeometry_LAMMPS");

	ofstream ofs(INPUT.geo_out.c_str());
	ofs << setprecision(16);

	ofs << "COMMENT" << endl;
	ofs << endl;
	ofs << cel.nat << " atoms" << endl;
	ofs << endl;
	ofs << INPUT.ntype << " atom types" << endl;
	ofs << endl;
	ofs << "0 " << cel.a1.x << " xlo xhi" << endl;
	ofs << "0 " << cel.a2.y << " ylo yhi" << endl;
	ofs << "0 " << cel.a3.z << " zlo zhi" << endl;

	// mohan add 2015-06-15
	if( INPUT.triclinic == 1 ) 
	{
		ofs << cel.a2.x << " " << cel.a3.x << " " << cel.a3.y << " xy xz yz" << endl;
	} 

	assert(cel.a1.y == 0.0);
	assert(cel.a1.z == 0.0);
	assert(cel.a2.z == 0.0);

	ofs << endl;
	ofs << "Atoms" << endl;
	ofs << endl;

	bool frac = true;
	int ia2=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		cout << " printing geometry for type " << it+1 << endl;
		cout << " there are " << cel.atom[it].na << " atoms for thie type of atom" << endl;
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			cel.direct2cartesian(it, ia);
			//ofs << ia2+1 << " " << it+1 
			ofs << ia2+1 << " " << "change_me" // mohan update 2015-06-15 
			    << " " << cel.atom[it].pos[ia].x
				<< " " << cel.atom[it].pos[ia].y
				<< " " << cel.atom[it].pos[ia].z << endl;

			++ia2;
		}
	}

	ofs.close();
	return;
}
