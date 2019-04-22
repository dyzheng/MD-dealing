#include "input.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>

Input INPUT;


Input::Input() 
{
	triclinic = 0; // default is 0 for LAMMPS geometry file, mohan add 2015-05-15	
}

Input::~Input() {}

void Input::Init(const string &fn, const int &argc)
{
	time_t time_start = std::time(NULL);
	cout << " ------------------------------------------------------------------------------------" << endl;
	cout << "              WELCOME TO D310 " << ctime(&time_start) << endl;
	if(argc>1)
	{
		cout << "  GOO provide the following functions and                                            " << endl;
		cout << "  Here is the introduction to the required 'INPUT' file for each function            " << endl;
		cout << " " << endl;
		cout << "  1. Extend the cell" << endl << endl; 
		cout << "  PARAMETER     VALUE     EXPLNATION                                                 " << endl; 
		cout << "  calculation   ext       Extend the cell.                                           " << endl;
		cout << "  ext_1         n1        Extend cell along 1st direction for n1 times.              " << endl;
		cout << "  ext_2         n2        Extend cell along 1st direction for n1 times.              " << endl;
		cout << "  ext_3         n3        Extend cell along 1st direction for n1 times.              " << endl;
		cout << "  geo_in        Li.ion    Input name of geometry file.                               " << endl;
		cout << "  geo_in_type   PROFESS   Input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA  " << endl;
		cout << "  geo_out       Li2.ion   Output name of geometry file.                              " << endl;
		cout << "  geo_out_type  PROFESS   Output type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA " << endl;
		cout << "  ntype         1         Number of different types of atoms.                        " << endl;
		cout << "  natom         2         Total number of atoms.                                     " << endl;

		cout << " " << endl;
		cout << "  2. Plot the pair distribution function " << endl << endl;
		cout << "  PARAMETER     VALUE     EXPLNATION                                                 " << endl; 
		cout << "  calculation   pdf       Pair distribution function.                                " << endl;
		cout << "  geo_in_type   PROFESS   Input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA  " << endl;
		cout << "  geo_in        Li.ion    Input name of geometry file.                               " << endl;
		cout << "  geo_out       Li.pdf    Output name of geometry file.                              " << endl;
		cout << "  ntype         1         Number of different types of atoms.                        " << endl;
		cout << "  natom         1024      Total number of atoms.                                     " << endl;
		cout << "  pdf_dr        0.01      Delta r in PDF calculation (Angstrom).                     " << endl;

		cout << " " << endl;
		cout << "  3. Plot the distribution of velocities " << endl << endl;
		cout << "  PARAMETER     VALUE     EXPLNATION                                                 " << endl;
		cout << "  calculation   vel       Distribution of velocities                                 " << endl;
		cout << "  vel_in        Li.ion    Input name of velocity file.                               " << endl;
		cout << "  vel_out       Li.pdf    Output name of velocity file.                              " << endl;
		cout << "  natom         1024      Total number of atoms.                                     " << endl;
		cout << "  ndv           0.01      The different number of printed velocities.                " << endl;

		cout << " " << endl;
		cout << "  5. Density operation " << endl << endl;
		cout << "  PARAMETER     VALUE     EXPLNATION                                                 " << endl;
		cout << "  calculation   profile   Profiles of given 3D data.                                 " << endl;
		cout << "  direction     3         Direction of profile.                                      " << endl;
		cout << "  data_in       den.in    Input file contains the 3D data.                           " << endl;
		cout << "  data_out      den.out   Output file contains the profile data.                     " << endl;

		cout << " " << endl;
		cout << " Insert atoms into existing structures " << endl << endl;

		cout << " " << endl;
		cout << " ------------------------------------------------------------------------------------" << endl;
	}

	cout << " AUTHOR: Mohan, LAST UPDATE: 12/06/12" << endl;

	cout << setiosflags(ios::right);
	cout << setiosflags(ios::left);
	cout << setiosflags(ios::left);

	cout << " READING from INPUT file" << endl;


    this->Default();
    this->Read(fn);
    this->Check();
    time_t  time_now = time(NULL);

   

	    return;
}

void Input::Default(void)
{
	pi = 3.1415926535897932384626;
	step_interval_dynamics = -1;

	geo_interval = -1;

	format3D = 1;

	isf_dg = -1.0;
	isf_ng = -1;

	iprof_nr = 0;
	iprof_b = 0.0;
    return;
}

void Input::Read(const string &fn)
{
    ifstream ifs(fn.c_str(), ios::in);	// "in_datas/input_parameters"

    if (!ifs) 
	{
		cout << " Can't find the INPUT file." << endl;
		exit(0);
	}

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    //ifs >> setiosflags(ios::uppercase);
    ifs.rdstate();

    while (ifs.good())
    {
        ifs >> word1;
        strtolower(word1, word);
//		cout << " word = " << word << endl;

//----------------------------------------------------------
// main parameters
//----------------------------------------------------------
		// >>> General <<<
        if (strcmp("calculation", word) == 0) read_value(ifs, calculation);
        else if (strcmp("ext_1", word) == 0) read_value(ifs, ext_1);
        else if (strcmp("ext_2", word) == 0) read_value(ifs, ext_2);
        else if (strcmp("ext_3", word) == 0) read_value(ifs, ext_3);
        else if (strcmp("geo_in", word) == 0) read_value(ifs, geo_in);
        else if (strcmp("geo_in_type", word) == 0) read_value(ifs, geo_in_type);
        else if (strcmp("geo_out", word) == 0) read_value(ifs, geo_out);
        else if (strcmp("geo_out_type", word) == 0) read_value(ifs, geo_out_type);
        else if (strcmp("geo_directory", word) == 0) read_value(ifs, geo_directory);
		else if (strcmp("geo_1", word) == 0) read_value(ifs, geo_1);
		else if (strcmp("geo_2", word) == 0) read_value(ifs, geo_2);
		else if (strcmp("geo_interval", word) == 0) read_value(ifs, geo_interval);
        else if (strcmp("ntype", word) == 0) { 
                        read_value(ifs, ntype);         
                        atom=new int[10];
                        for(int ii=0;ii<10;ii++)atom[ii]=0;
        }
        else if (strcmp("natom", word) == 0) read_value(ifs, natom);
        else if (strcmp("cartesian", word) == 0) read_value(ifs, cartesian); //mohan add 2014-04-04
		// >>> Function 2 <<<
        else if (strcmp("pdf_dr", word) == 0) read_value(ifs, pdf_dr);
        else if (strcmp("pdf_rcut", word) == 0) read_value(ifs, pdf_rcut);
        else if (strcmp("pdf_ele1", word) == 0) read_value(ifs, pdf_ele1);
        else if (strcmp("pdf_ele2", word) == 0) read_value(ifs, pdf_ele2);
        else if (strcmp("pdf_z0", word) == 0) read_value(ifs, pdf_z0);
        else if (strcmp("pdf_z1", word) == 0) read_value(ifs, pdf_z1);
        else if (strcmp("struf_dgx", word) == 0) read_value(ifs, struf_dgx);
        else if (strcmp("struf_dgy", word) == 0) read_value(ifs, struf_dgy);
        else if (strcmp("struf_dgz", word) == 0) read_value(ifs, struf_dgz);
        else if (strcmp("struf_ng", word) == 0) read_value(ifs, struf_ng);
        else if (strcmp("ssf_out", word) == 0) read_value(ifs, ssf_out);
		// >>> Function 3 <<<
        else if (strcmp("vel_in", word) == 0) read_value(ifs, vel_in);
        else if (strcmp("vel_out", word) == 0) read_value(ifs, vel_out);
        else if (strcmp("ndv", word) == 0) read_value(ifs, ndv);
		// >>> Function 4 <<<
        else if (strcmp("vacuum_x1", word) == 0) read_value(ifs, vacuum_x1);
        else if (strcmp("vacuum_x2", word) == 0) read_value(ifs, vacuum_x2);
        else if (strcmp("vacuum_y1", word) == 0) read_value(ifs, vacuum_y1);
        else if (strcmp("vacuum_y2", word) == 0) read_value(ifs, vacuum_y2);
        else if (strcmp("vacuum_z1", word) == 0) read_value(ifs, vacuum_z1);
        else if (strcmp("vacuum_z2", word) == 0) read_value(ifs, vacuum_z2);
		// >>> Function 5 <<<
        else if (strcmp("direction", word) == 0) read_value(ifs, direction);
        else if (strcmp("data_in", word) == 0) read_value(ifs, data_in);
        else if (strcmp("data_out", word) == 0) read_value(ifs, data_out);
		// >>> Function : pseudopotential <<<
        else if (strcmp("pseudo_z", word) == 0) read_value(ifs, pseudo_z);
        else if (strcmp("pseudo_type", word) == 0) read_value(ifs, pseudo_type);
        else if (strcmp("pseudo_in", word) == 0) read_value(ifs, pseudo_in);
        else if (strcmp("pseudo_out", word) == 0) read_value(ifs, pseudo_out);
		// >>> Function : iprofile
        else if (strcmp("iprof_nr", word) == 0) read_value(ifs, iprof_nr);
        else if (strcmp("iprof_b", word) == 0) read_value(ifs, iprof_b);
		// >>> Function : eprofile
        else if (strcmp("format3d", word) == 0) read_value(ifs, format3D);
		// >>> Function : dynamics structure factor
        else if (strcmp("dsf_out", word) == 0) read_value(ifs, dsf_out);
        else if (strcmp("dsf_dt", word) == 0) read_value(ifs, dsf_dt);
        else if (strcmp("dsf_neqi", word) == 0) read_value(ifs, dsf_neqi);
		// >>> Function : velocity correlation function
        else if (strcmp("velcor_in_type", word) == 0) read_value(ifs, velcor_in_type);
        else if (strcmp("velcor_directory", word) == 0) read_value(ifs, velcor_directory);
		else if (strcmp("velcor_1", word) == 0) read_value(ifs, velcor_1);
		else if (strcmp("velcor_2", word) == 0) read_value(ifs, velcor_2);
		else if (strcmp("velcor_neqi", word) == 0) read_value(ifs, velcor_neqi);
		else if (strcmp("velcor_out", word) == 0) read_value(ifs, velcor_out);
		else if (strcmp("velcor_atom", word) == 0) read_value(ifs, velcor_atom);
		else if (strcmp("step_interval_dynamics", word) == 0) read_value(ifs, step_interval_dynamics);
		// >>> Function : power spectra
		else if (strcmp("ps_nv", word) == 0) read_value(ifs, ps_nv);
		else if (strcmp("ps_dw", word) == 0) read_value(ifs, ps_dw);
		else if (strcmp("ps_nw", word) == 0) read_value(ifs, ps_nw);
		else if (strcmp("ps_in", word) == 0) read_value(ifs, ps_in);
		else if (strcmp("ps_out", word) == 0) read_value(ifs, ps_out);
		else if (strcmp("ps_dt", word) == 0) read_value(ifs, ps_dt);
		// >>> Function : intermediate scattering function
		else if (strcmp("isf_target_q", word) == 0) read_value(ifs, isf_target_q);
		else if (strcmp("isf_nconfig", word) == 0) read_value(ifs, isf_nconfig);
		else if (strcmp("isf_ncorrelation", word) == 0) read_value(ifs, isf_ncorrelation);
		else if (strcmp("isf_dcorrelation", word) == 0) read_value(ifs, isf_dcorrelation);
		else if (strcmp("isf_outfile", word) == 0) read_value(ifs, isf_outfile);
		else if (strcmp("atom1", word) == 0) read_value(ifs, atom[0]);
		else if (strcmp("atom2", word) == 0) read_value(ifs, atom[1]);
                else if (strcmp("atom3", word) == 0) read_value(ifs, atom[2]);
                else if (strcmp("atom4", word) == 0) read_value(ifs, atom[3]);
                else if (strcmp("atom5", word) == 0) read_value(ifs, atom[4]);
                else if (strcmp("atom6", word) == 0) read_value(ifs, atom[5]);
                else if (strcmp("atom7", word) == 0) read_value(ifs, atom[6]);
                else if (strcmp("atom8", word) == 0) read_value(ifs, atom[7]);
                else if (strcmp("atom9", word) == 0) read_value(ifs, atom[8]);
                else if (strcmp("atom10", word) == 0) read_value(ifs, atom[9]);
                            
		else if (strcmp("isf_dg", word) == 0) read_value(ifs, isf_dg);
		else if (strcmp("isf_ng", word) == 0) read_value(ifs, isf_ng);
		else if (strcmp("isf_config_start", word) == 0) read_value(ifs, isf_config_start);
		// >>> Function : bond-angle distribution functions
		else if (strcmp("bdf_nadj", word) == 0) read_value(ifs, bdf_nadj);
		else if (strcmp("bdf_dtheta", word) == 0) read_value(ifs, bdf_dtheta);
		else if (strcmp("bdf_out", word) == 0) read_value(ifs, bdf_out);
		else if (strcmp("bdf_movie", word) == 0) read_value(ifs, bdf_movie);
		// >>> Function : movement z
		else if (strcmp("movement_x", word) == 0) read_value(ifs, movement_x);
		else if (strcmp("movement_y", word) == 0) read_value(ifs, movement_y);
		else if (strcmp("movement_z", word) == 0) read_value(ifs, movement_z);
		// >>> Function : insert atoms
		else if (strcmp("natom_new", word) == 0) read_value(ifs, natom_new);
		else if (strcmp("element_new", word) == 0) read_value(ifs, element_new);
		else if (strcmp("min_dis", word) == 0) read_value(ifs, min_dis);
        // >>> Function : triclinic
		else if (strcmp("triclinic", word) == 0) read_value(ifs, triclinic); //mohan add 2015-06-16
                //zhengdy
                else if (strcmp("msd_replica", word) == 0) read_value(ifs, msd_replica);
		// >>> Function : new
        else
        {
            cout << " THE PARAMETER NAME '" << word
               << "' IS NOT USED!" << endl;
            ifs.ignore(150, '\n');
        }

        ifs.rdstate();
        if (ifs.eof() != 0)
        {
			break;
        }
        else if (ifs.bad() != 0)
        {
			cout << " Bad input parameters. " << endl;
			exit(0);
        }
        else if (ifs.fail() != 0)
        {
			cout << " word = " << word << endl;
			cout << " Fail to read parameters. " << endl; 
            ifs.clear();
			exit(0);
        }
        else if (ifs.good() == 0)
        {
			break;
        }
    }

    return;
}//end read_parameters


void Input::Check(void)
{
	if( vacuum_x1 < 0.0 ) QUIT("vacuum_x1 < 0.0"); 
	if( vacuum_x2 < 0.0 ) QUIT("vacuum_x2 < 0.0"); 
	if( vacuum_y1 < 0.0 ) QUIT("vacuum_y1 < 0.0"); 
	if( vacuum_y2 < 0.0 ) QUIT("vacuum_y2 < 0.0"); 
	if( vacuum_z1 < 0.0 ) QUIT("vacuum_z1 < 0.0"); 
	if( vacuum_z2 < 0.0 ) QUIT("vacuum_z2 < 0.0"); 

	cout << " ------------------------------------------------------" << endl;
	if(calculation=="ssf") cout << " Category B: Static Structure Factor. " << endl;
	else if(calculation=="dsf") cout << " Category C: Dynamics Structure Factor. " << endl;
	else if(calculation=="velcor") cout << " Category C: Velocity Auto-Correlation Function. " << endl;
	cout << " ------------------------------------------------------" << endl;
        
        cout<< "check ntype and natom:" << endl;
        if(!(natom==atom[0]+atom[1]+atom[2]+atom[3]+atom[4]+atom[5]+atom[6]+atom[7]+atom[8]+atom[9]))QUIT("NOT FIT! PLEASE CHECK!");
        //cout<< "ok~" <<endl;
    return;
}


void Input::readbool(ifstream &ifs, bool &var)
{
    string str;
    ifs >> str;
    if (str == "true")
    {
        var = true;
    }
    else
    {
        var = false;
    }
    ifs.ignore(100, '\n');
    return;
}

void Input::strtolower(char *sa, char *sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

