#include "cell.h"
#include "matrix3.h"

// transform the coordinates from direct to cartesian
void Cell::direct2cartesian(const int &it, const int &i)
{
	assert( atom[it].pos != NULL);
	assert( atom[it].posd != NULL);
	atom[it].pos[i].x = atom[it].posd[i].x * a1.x + atom[it].posd[i].y * a2.x + atom[it].posd[i].z * a3.x;
	atom[it].pos[i].y = atom[it].posd[i].x * a1.y + atom[it].posd[i].y * a2.y + atom[it].posd[i].z * a3.y;
	atom[it].pos[i].z = atom[it].posd[i].x * a1.z + atom[it].posd[i].y * a2.z + atom[it].posd[i].z * a3.z;
	//cout << " cartesian : " << atom[it].pos[i].x << " " << atom[it].pos[i].y << " " << atom[it].pos[i].z << endl;
	return;
}	

// transform the coordinates from cartesian to direct 
void Cell::cartesian2direct(const int &it, const int &i)
{
	assert( atom[it].pos != NULL);
	assert( atom[it].posd != NULL);
	
	Matrix3 lattice_vector, inv_lat;
	lattice_vector.e11=a1.x;
	lattice_vector.e12=a1.y;
	lattice_vector.e13=a1.z;
	lattice_vector.e21=a2.x;
	lattice_vector.e22=a2.y;
	lattice_vector.e23=a2.z;
	lattice_vector.e31=a3.x;
	lattice_vector.e32=a3.y;
	lattice_vector.e33=a3.z;

	inv_lat = lattice_vector.Inverse();
	Vector3<double> direct_vec, cartesian_vec;
	cartesian_vec.x = atom[it].pos[i].x;
	cartesian_vec.y = atom[it].pos[i].y;
	cartesian_vec.z = atom[it].pos[i].z;
	direct_vec = cartesian_vec * inv_lat;
	atom[it].posd[i].x = direct_vec.x;
	atom[it].posd[i].y = direct_vec.y;
	atom[it].posd[i].z = direct_vec.z;
	return;
}	



void Cell::add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		double &x, double &y, double &z) const
{
	x = this->atom[it].pos[ia].x + i*a1.x + j*a2.x + k*a3.x;
	y = this->atom[it].pos[ia].y + i*a1.y + j*a2.y + k*a3.y;
	z = this->atom[it].pos[ia].z + i*a1.z + j*a2.z + k*a3.z;
	return;
}

void Cell::add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		float &x, float &y, float &z) const
{
	x = this->atom[it].pos[ia].x + i*a1.x + j*a2.x + k*a3.x;
	y = this->atom[it].pos[ia].y + i*a1.y + j*a2.y + k*a3.y;
	z = this->atom[it].pos[ia].z + i*a1.z + j*a2.z + k*a3.z;
	return;
}

