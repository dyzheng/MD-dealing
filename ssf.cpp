//------------------------------------------
// STRUCTURE OF CLASS:
//   CLASS SSF
//     |_FUNCTION Routine
//     |_FUNCTION check_file_exist
//     |_FUNCTION cal
//     |_FUNCTION ssf_3D
//       |_FUNCTION sumup
//     |_FUNCTION cal_diff_norm
//     |_FUNCTION ssf_1D
//     |_FUNCTION write_ssf
//------------------------------------------
#include "cellFile.h"
#include "input.h"
#include "ssf.h"
#include "math.h"

void SSF::Routine()
{
	TITLE("SSF","Routine");
	this->cal();

	return;
}

void SSF::check_file_exist(const string &name)
{
	ifstream ifs(name.c_str());
	if(ifs)
	{
		cout << " The SSF has been calculated, quit." << endl;
		exit(0);
	}
	else
	{
		cout << " The file " << name << " is calculating now!" << endl;
	}
}


//-------------------------------------------------------------
// we will calculate the static structure factor S(g) at each
// g point, where the interval of g points is delta_g,
//-------------------------------------------------------------
void SSF::cal()
{
	TITLE("SSF","cal");

	// --> INITIALIZE <--

	assert(INPUT.struf_dgx>0);
	assert(INPUT.struf_dgy>0);
	assert(INPUT.struf_dgz>0);
	assert(INPUT.struf_ng>0);

	// number of delta_g, determined by 2pi/a 
	// where a is the lattice constant.
	this->dgx = INPUT.struf_dgx;
	this->dgy = INPUT.struf_dgy;
	this->dgz = INPUT.struf_dgz;

	// number of g points along each direction.
	this->ngx = INPUT.struf_ng;
	this->ngy = INPUT.struf_ng;
	this->ngz = INPUT.struf_ng;

	// total nuber of g points (3 dimensional).
	this->ngtot = (2*ngx+1) * (2*ngy+1) * (2*ngz+1);

	// static structure factor in 3 dimensional.
	float* sf = new float[ngtot];
	ZEROS(sf, ngtot);
	cout << " Number of total K=2pi/L*(" 
		<< ngx << "," 
		<< ngy << "," 
		<< ngz << ") = " << ngtot <<endl;

	// ---> BODY <--- 

	// circle from geo_1 to geo_2.
	// because the final static structure factor
	// is an average property, so we need to have
	// a log of samples here.
	assert(INPUT.geo_2 >= INPUT.geo_1);
	assert(INPUT.geo_interval > 0);

	int count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; igeo+=INPUT.geo_interval)
	{
		// cel_in : input geometry file.
		CellFile cel_in;

		stringstream ss; ss << igeo;
		cel_in.file_name = ss.str(); 
		// cout << " File name is " << ss.str() << endl;

		// Read in geometry.
		if( !CellFile::ReadGeometry( cel_in ) ) continue;
		++count_geometry_number;

		// (1) Calculate the static structure factor (3D).
		this->ssf_3D( cel_in, sf );
	}
	
	// do average of 3D ssf.
	assert(count_geometry_number>0);
	cout << " count_geometry_number = " << count_geometry_number << endl;
	if(count_geometry_number>0)
	{
		for(int ig=0; ig<ngtot; ++ig)
		{
			sf[ig] /= count_geometry_number;
		}
	}
	
	// (2) Prepare the G vectors. Next we will project 3D ssf
	// to 1D.
	int* nG_1D = new int[ngtot]; // Number of |G|.
	int* norm_index = new int[ngtot]; // Index between 3D G and 1D |G|.
	ZEROS( nG_1D, ngtot );
	ZEROS( norm_index, ngtot );
	this->diff_norm = cal_diff_norm( nG_1D, norm_index );

	// (3) Calculate the 1D SSF.
	float* G_1D = new float[diff_norm]; // Norm of |G|
	float* sf_1D = new float[diff_norm]; // SSF in 1D.
	ZEROS( G_1D, diff_norm );
	ZEROS( sf_1D, diff_norm );
	this->ssf_1D( G_1D, sf_1D, norm_index, nG_1D, sf );

	// (4) Output the final SSF.
	this->write_ssf( G_1D, sf_1D );


	// --> CLEAN <--

	delete[] sf;
	delete[] nG_1D;
	delete[] norm_index;
	delete[] G_1D;
	delete[] sf_1D;
	return;
}

void SSF::ssf_3D( 
	const Cell &cel, // cell information
	float* sf // final structure factor
	) const
{
	TITLE("SSF","ssf_3D");

	// --> INITIALIZE <--

	float* sum_exp = new float[ngtot];
	ZEROS( sum_exp, ngtot );

    //// Calculate the distance between atoms.
    float x2, y2, z2; // atom positions for atom 2.
    float dr[3]; // delta x,y,z between atom1, atom2.
    float dis;


	// --> BODY <--

    for(int it=0; it<INPUT.ntype; ++it)
    {
        for(int ia=0; ia<cel.atom[it].na; ++ia)
        {
			//// Double check the atom positions.
			//// cout << cel.atom[it].pos[ia].x << " "
			//// << cel.atom[it].pos[ia].y << " "
			//// << cel.atom[it].pos[ia].z << endl;
		  	int count2 = 0;

			//// it2 start from species 1, because its about 'atom pairs',
			//// so we only need to calculate once.
            for(int it2=it; it2<INPUT.ntype; ++it2)
            //for(int it2=0; it2<INPUT.ntype; ++it2)
            {
                for(int ia2=ia; ia2<cel.atom[it].na; ++ia2)
                //for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
                {
					// In ssf, for the identical atoms,
					// they will contribute 1 in the final
					// expression of ssf.
                    if(it==it2 && ia==ia2) continue;
					
                    // Search in the around cells.
					// Try to find the shortest atom distance
					// between atom 1 and atom 2,
					// then that's the distance we want!
					float shortest_distance2 = 10000.0;
					int which_i, which_j, which_k;
                    for(int i=-1; i<=1; ++i)
                    {
                        for(int j=-1; j<=1; ++j)
                        {
                            for(int k=-1; k<=1; ++k)
                            {
                                // add cell length to atom 2
                                cel.add_cell_length(it2, ia2, i, j, k, x2, y2, z2);
                                // calculate the distance between two atoms |r_1 - r_2|
                                dr[0] = cel.atom[it].pos[ia].x - x2;
                                dr[1] = cel.atom[it].pos[ia].y - y2;
                                dr[2] = cel.atom[it].pos[ia].z - z2;
								// to save the calculation, we avoid using sqrt.
								dis = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
								if(dis < shortest_distance2) 
								{
									shortest_distance2=dis;
									which_i=i;
									which_j=j;
									which_k=k;
								}
							}
						}
					}
					// Here we identify the atom in cell: (which_i, which_j, which_k)
					
					// we get the vector 'dr' again.
                    cel.add_cell_length(it2, ia2, which_i, which_j, which_k, x2, y2, z2);
                    dr[0] = cel.atom[it].pos[ia].x - x2;
                    dr[1] = cel.atom[it].pos[ia].y - y2;
                    dr[2] = cel.atom[it].pos[ia].z - z2;

					this->sumup(sum_exp, dr);
					++count2;
				}
            }
			// For check,
			//cout << " ia" << ia << " count=" << count2 << endl;
        }
    }

	// Use the final expression of structure factor.
	// s(G)=1 + 1/N * [ \sum_{ij,j!=i}exp(iGr) ] 
	for(int ig=0; ig<ngtot; ++ig)
	{
		// 2.0 conunts for the atom pairs.
		double this_sf = 0.0;
		this_sf = 2.0*sum_exp[ig];
		this_sf/= INPUT.natom;
		this_sf+= 1.0;
		sf[ig] += this_sf;

		//cout << " ig=" << ig << " sum_exp=" << sum_exp[ig] << endl;
	}

	// --> CLEAN <--

	delete[] sum_exp;

	return;
}


void SSF::sumup( float *sum_exp, const float dr[3] ) const
{
	int ik=0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				++ik;

				// because we only need to use half of the 
				// G vectors (same in gamma-only algorithm)
				if(ix<0) continue;

				// For the diagonalize part, the factor is 1,
				// and for the non-diagonalize G vectors, it's 2.
				float factor=2.0;
				if(ix==0) factor=1.0;

				// Get the G vector.
				float kx = ix * this->dgx;
				float ky = iy * this->dgy;
				float kz = iz * this->dgz;

				// exp(ik*(r_i - r_j)) appears here!
				float phase = kx * dr[0] + ky * dr[1] + kz * dr[2];
				sum_exp[ik-1] += factor*cos(phase);
			}
		}
	}
	return;
}

int SSF::cal_diff_norm(
	int *nG_1D, 
	int *norm_index
) const
{
	TITLE("SSF","cal_diff_norm");

	// --> INITIALIZE <--

	int* norm_tanker = new int[ngtot];
	ZEROS( norm_tanker, ngtot );
	int diff_norm = 0;
	for(int ig=0; ig<ngtot; ++ig)
	{
		norm_tanker[ig] = -1;
		nG_1D[ig] = 0;
		norm_index[ig] = -1;
	}

	// --> BODY <--

	int ig_global=0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				bool not_new = false;
				int norm = ix*ix+iy*iy+iz*iz;

				//// Compare to the old |G|.
				for(int i=0; i<diff_norm; ++i)
				{
					if( norm_tanker[i] == norm )
					{
						not_new = true;
						norm_index[ig_global] = i;
					}
				}

				if(not_new)
				{
					nG_1D[norm_index[ig_global]] += 1;
				}
				else
				{
					//// Find a new |G| !
					norm_tanker[diff_norm] = norm;
					nG_1D[diff_norm] += 1;
					norm_index[ig_global] = diff_norm;
					++diff_norm;
		//			cout << " The new norm is = " << sqrt((float)norm) << endl;
				}
				++ig_global;
			}
		}
	}
	cout << " Diff_norm = " << diff_norm<< endl;

	// --> CLEAN <--

	delete[] norm_tanker;
	return diff_norm;
}
	

void SSF::ssf_1D( 
	float *G_1D, 
	float *sf_1D, 
	const int *norm_index,
	const int *nG_1D,
	const float *sf	
) const
{
	TITLE("SSF","ssf_1D");

	// --> INITIALIZE <--

	for(int i=0; i<diff_norm; ++i)
	{
		G_1D[i] = -1.0;
		sf_1D[i] = 0.0;
	}

	// --> BODY <--

	int ig_global=0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				const int ig_1D = norm_index[ig_global];

                G_1D[ig_1D] = pow(ix*dgx,2)+pow(iy*dgy,2)+pow(iz*dgz,2);
                G_1D[ig_1D] = sqrt(G_1D[ig_1D]);

//				if(G_1D[ig_1D]<0.0)
//				{
//					G_1D[ig_1D] = sqrt((float)ix*ix+iy*iy+iz*iz)*dg;
//				}

				sf_1D[ig_1D] += sf[ig_global];
				++ig_global;
			}
		}	
	}


	for(int i=0; i<diff_norm; ++i)
	{
		assert(nG_1D[i]>0);
		sf_1D[i]/=nG_1D[i];
		//ofs << nor[i] << " " << nG_1D[i] << endl;
	}

	// --> CLEAN <--

	// no clean

	return;
}


void SSF::write_ssf( const float *G_1D, const float* sf_1D ) const
{
	TITLE("SSF","write_ssf");

	// --> INITIALIZE <--

	assert(this->diff_norm > 0);
	// output the static structure factor.
	ofstream ofs(INPUT.ssf_out.c_str());
	bool *visited = new bool[diff_norm];
	for(int i=0; i<diff_norm; ++i) visited[i] = false;

	
	// --> BODY <--

	for(int i=0; i<diff_norm; ++i)
	{
		int index = 0;
		float min = 10000000;
		for(int j=0; j<diff_norm; ++j)
		{
			if( visited[j] ) continue;
			else if( G_1D[j] < min )
			{
				index = j;
				min = G_1D[j];
			}
		}
		visited[index] = true;
		if(G_1D[index]>0)
		{
			ofs << G_1D[index] << " " << sf_1D[index] << endl;
		}
	}

	// --> CLEAN <--

	delete[] visited;
	return;
}

