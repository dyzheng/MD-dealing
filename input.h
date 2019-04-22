#ifndef INPUT_H
#define INPUT_H

#include "gfun.h"

class Input
{
public:
    
	Input();
    ~Input();
    
	// read in the parameters, be called in main.cpp
	void Init(const string &fn, const int &argc);

	private:

    void Read(const string &fn);
    void Default(void);
    void Check(void);

	public:

	// >>> General <<<
	string calculation; // type of calculations.
	int debug;

	// >>> For function 1 : geo <<< 
	int ext_1;
	int ext_2;
	int ext_3;
	string geo_in;
	string geo_in_type;
	string geo_out; 
	string geo_out_type;
	string geo_directory;
	int geo_1; // start geometry file for PDF, SSF, ...
	int geo_2; // last geometry file for PDF, SSF, ...
	int geo_interval; // interval between geometries
	int ntype;
	int natom;
	bool cartesian;

	// >>> For function 2 : pdf <<<
	double pdf_dr;
	double pdf_rcut;
	int pdf_ele1;
	int pdf_ele2;
	int pdf_z0; // for pdf 2d
	int pdf_z1; // for pdf 2d

	// >>> Static Structure Factor <<<
	double struf_dgx;
	double struf_dgy;
	double struf_dgz;
	int struf_ng;
	string ssf_out;

	// >>> For function 3 : vel <<< 
	string vel_in;
	string vel_out;
	int ndv; // number of delta v.


	// >>> For function 4 : vacuum <<<
	double vacuum_x1;
	double vacuum_x2;
	double vacuum_y1;
	double vacuum_y2;
	double vacuum_z1;
	double vacuum_z2;

	// >>> For function 5 : 3 dimensional density profile <<<
	int direction;
	string data_in;
	string data_out;

	// >>> For function: pseudopotential <<< 
	int pseudo_z;
	string pseudo_type;
	string pseudo_in;
	string pseudo_out;

	// >>> For function : i profile <<<
	int iprof_nr;
	double iprof_b; // parameters for gaussian broadening,
	// small b leads to sharp curve
	// larger b leads to smooth curve

	// >>> For function : electron proflie <<<
	int format3D;

	// >>> For Function : dynamics structure factor <<<
	string dsf_out;
	float dsf_dt;
	int dsf_neqi;

	// >>> For Function : velocity correlation functions.
	string velcor_in_type;
	string velcor_directory;
	int velcor_1;
	int velcor_2;
	int velcor_neqi;
	string velcor_out;
	int velcor_atom;
	int step_interval_dynamics;

	// >>> For Function : power spectra
	int ps_nv;     // number of velocity autocorrelation functions (VAF) points
	double ps_dw;  // increasement of frequency w
	int ps_nw;     // number of w points
	string ps_in;  // input file name
	string ps_out; // output file name
	double ps_dt;  // increasement of time in VAF

	// >>> For Function : intermediate scattering functions
	float isf_target_q;
	int isf_nconfig;
	int isf_ncorrelation;
	int isf_dcorrelation;
	string isf_outfile;
	int *atom;
	float isf_dg;
	int isf_ng;
	int isf_config_start;

	// >>> For Function : bond distribution function
	int bdf_nadj;
	double bdf_dtheta;
	string bdf_out;
	int bdf_movie;

	// >>> For Function : average_ion
	double movement_x;
	double movement_y;
	double movement_z;

	// >>> For Function : insert atoms
	int natom_new;
	string element_new;
	double min_dis;


	// >>> For Functio : triclinic format from LAMMPS geometry file
	int triclinic;

        // >>> For Function : msd
        int msd_replica;

	// >> CONSTANT <<
	double pi;

	public:
    
	template <class T>
    static void read_value(ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }

	void strtolower(char *sa, char *sb);
	void readbool(ifstream &ifs, bool &var);
};

extern Input INPUT;

#endif //INPUT
