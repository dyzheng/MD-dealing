#include "File.h"
#include "input.h"
#include "gfun.h"

bool File::CheckGeometry_ABACUS( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_ABACUS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.msd_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
    else return true;
}

bool File::ReadGeometry_ABACUS( Cell &cel )
{
       TITLE("CellFile","ReadGeometry_ABACUS");
        const int ntype = INPUT.ntype;

        // (1) open the file.
        stringstream ss;
        ss << INPUT.ABACUS_directory << "/";
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
//      cout << " Name is " << cel.system_name << endl;
        string useless;
        getline(ifs, useless);
//      cout << useless << endl;
        getline(ifs, useless);
        getline(ifs, useless);

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

        // (3) read in atom species and the pseudopotential file.
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
        cout << " Total atom number is " << cel.nat << endl;

        cel.coordinate = "Direct";

        // (4) read in the atom positions.
        Vector3<double> add1,add2,add3;
        for(int it=0; it<ntype; ++it)
        {
                cel.atom[it].read_pos_3(ifs);
            for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
                {
                        cel.direct2cartesian(it, ia2);
                }
        }
        ifs.close();

        return true;
}

void File::ReadDeltaR(Cell &cel)
{}
