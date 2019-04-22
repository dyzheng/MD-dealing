#include "atoms.h"
#include "input.h"

Atoms::Atoms()
{

}


Atoms::~Atoms()
{

}

// read the 'fractional' coordinates of atoms
// or the 'cartesian' coordinates of atoms.
void Atoms::read_pos(ifstream &ifs,bool frac)
{
	assert(na>0);
	string idtmp;
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	for(int i=0; i<na; ++i)
	{
		ifs >> idtmp;
		if(frac) // fractional coordinates
		{
			ifs >> posd[i].x >> posd[i].y >> posd[i].z;
			//cout << " pos=" << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl; 
		//	cout << ifs.good() << " " << ifs.bad() << " " << ifs.fail() << endl;
		}	
		else // cartesian coordinates
		{
			ifs >> pos[i].x >> pos[i].y >> pos[i].z;

			// right now only works for cartesian coordinates 2015-06-22
		    if(INPUT.movement_x != 0 || INPUT.movement_y !=0 || INPUT.movement_z !=0)
			{  
				pos[i].x += INPUT.movement_x;	
				pos[i].y += INPUT.movement_y;	
				pos[i].z += INPUT.movement_z;	
			}
			//cout << " pos=" << pos[i].x << " " << pos[i].y << " " << pos[i].z << endl;
		}
		if( ifs.fail() )
		{
			cout << " Reading atom " << i+1 << endl;
			QUIT("Error in reading atom positions");
		}
	}
	return;
}


// read the 'fractional' coordinates of atoms
// or the 'cartesian' coordinates of atoms.
void Atoms::read_pos_2(ifstream &ifs,bool frac)
{
	assert(na>0);
	string tmp;
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	for(int i=0; i<na; ++i)
	{
		if(frac)
		{
			ifs >> posd[i].x >> posd[i].y >> posd[i].z;
			getline(ifs, tmp);
			//cout << " posd=" << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl; 

		//	cout << ifs.good() << " " << ifs.bad() << " " << ifs.fail() << endl;
			if( ifs.fail() )
			{
				cout << " Reading atom " << i+1 << endl;
				QUIT("Error in reading atom positions 2");
			}
		}	
		else
		{
			ifs >> pos[i].x >> pos[i].y >> pos[i].z;
			getline(ifs, tmp);
		}
	}
	return;
}



void Atoms::read_pos_3(ifstream &ifs)
{
    assert(na>0);
    string tmp;
    pos = new Vector3<double>[na];
    posd = new Vector3<double>[na];
    for(int i=0; i<na; ++i)
    {
        ifs >> id >> posd[i].x >> posd[i].y >> posd[i].z;
        getline(ifs, tmp);
//		cout << id << " " << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl;
    }
    return;
}




void Atoms::read_vel(ifstream &ifs)
{
	assert(na>0);
	this->vel = new Vector3<double>[na];
//	cout << " atom number is " << na << endl;
	for(int i=0; i<na; ++i)
	{
		ifs >> vel[i].x >> vel[i].y >> vel[i].z;
		//cout << " vel=" << vel[i].x << " " << vel[i].y << " " << vel[i].z << endl;
	}
	return;
}


void Atoms::clear()
{
//	delete[] vel;
}
