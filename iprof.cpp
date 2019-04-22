#include "cellFile.h"
#include "input.h"
#include "iprof.h"
#include "math.h"

void Iprofile::Routine()
{
	TITLE("Iprof","Routine");
	// cel_in : input geometry file
	CellFile cel_in;

	cel_in.file_name = INPUT.geo_in;
	CellFile::ReadGeometry( cel_in );


	// the original method
//	cal( cel_in );

	// new method to calcultae the ionic density profile
	// 2014-03-09
	cal_gauss( cel_in );

	return;
}

void Iprofile::cal_gauss( const Cell &cel)
{
	TITLE("Iprofile","cal_gauss");

	// natom / v_liquid

	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();
	assert(INPUT.iprof_nr>0);

	// assume 'a3' is the longest axis
	const int nr = INPUT.iprof_nr;
	double dr = norm3/(double)nr;

	// allocate iprof
	double* r = new double[nr];
	double* iprof = new double[nr];
	double* iprof_unit = new double[nr];
	for(int ir=0; ir<nr; ++ir)
	{
		r[ir] = (ir+1)*dr;
		iprof[ir] = 0.0;
		iprof_unit[ir] = 0.0;
	}


	// search for each atom
	assert( INPUT.iprof_b > 0.0 );
	double b = INPUT.iprof_b; // coefficient
	cout << " coefficeint for gaussian broadening is " << b << endl; 

	double b2 = b * b;
	double a = 1.0/sqrt(3.1415926)/b;
	int iat=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			double posz = cel.atom[it].pos[ia].z;
			while( posz >= norm3 )
			{
				posz -= norm3;	
			}
			while( posz < 0 )
			{
				posz += norm3;
			}
			// calculate the gaussian distribution
			for(int ir=0; ir<nr; ++ir)
			{
				double dz = posz - r[ir];
				double dz2 = dz * dz;
				double Gauss = exp( -dz2/b2 ) * a; 
				iprof[ir] += Gauss;
			}
		}
	}


	// unit is made for ionic density profile
	double sum = 0.0;
	double* rab = new double[nr];
	for(int ir=0; ir<nr; ++ir)
	{
		rab[ir] = dr;
	}

	// the mesh number should be odd
	if(nr%2==0)
	{
    	Math::Simpson_Integral(nr-1, iprof, rab, sum);
	}
	else
	{
    	Math::Simpson_Integral(nr, iprof, rab, sum);
	}

	cout << " sum from ionic density profile = " << sum << endl;
	assert( sum > 0.0 );
	// unit is made
	for(int ir=0; ir<nr; ++ir)
	{
		iprof_unit[ir] = iprof[ir]/sum;
	}


	// print data
	cout << " Print data to " << INPUT.geo_out << endl;
	ofstream ofs(INPUT.geo_out.c_str());
	for(int ir=0; ir<nr; ++ir)
	{
		ofs << r[ir] << " " << iprof[ir] << " " << iprof_unit[ir] << endl; 
	}
	ofs.close();


	// finish
	delete[] iprof;
	delete[] iprof_unit;
	delete[] r;
	delete[] rab;

	return;
}


void Iprofile::cal( const Cell &cel )
{
	TITLE("Iprofile","cal");

	// natom / v_liquid
	// ionic density
	double rho0 = 0.0508791207;

	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();
	assert(INPUT.iprof_nr>0);

	// assume 'a3' is the longest axis
	double dr = norm3/(double)INPUT.iprof_nr;

	double v_each_z = norm1 * norm2 * dr;

	double *nat_z = new double[INPUT.iprof_nr+1];
	for(int i=0; i<INPUT.iprof_nr; ++i) nat_z[i]=0.0;

	cout << " Length of Lattice Vectors : " << norm1 << " " << norm2 << " " << norm3 << endl;
	cout << " dimension of nat_z along z : " << INPUT.iprof_nr << " dr = " << dr << endl;


	int iat=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			int which = (int)(cel.atom[it].pos[ia].z/dr);

			// consider periodic boundary condition
			if(which >= (INPUT.iprof_nr+1) )
			{
				which -= INPUT.iprof_nr;
				cout << " WARNING! atom out of z range (>0)." << endl;
			}
			else if(which < 0)
			{
				which += INPUT.iprof_nr;
				cout << " WARNING! atom out of z range (<0)." << endl;
			}
			nat_z[which]+=1.0;
			++iat;
		}
	}
	cout << " Total atom number is " << iat << endl;


	cout << " Print data to " << INPUT.geo_out << endl;
	ofstream ofs(INPUT.geo_out.c_str());
	for(int i=0; i<INPUT.iprof_nr; ++i)
	{
		ofs << dr*(i+1) << " " << nat_z[i] << " " << nat_z[i]/v_each_z/rho0 << endl;
	}
	ofs.close();

	delete[] nat_z;
	return;
}
