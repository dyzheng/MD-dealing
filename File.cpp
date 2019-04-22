#include "File.h"
#include "input.h"


bool File::ReadGeometry( Cell &cel )
{
          TITLE("File","ReadGeometry");
          
          return ReadGeometry_ABACUS(cel);
          
          
          else
          {
                cout << " Error here: ReadGeometry." << endl;
                cout << " msd_in_type = " << INPUT.msd_in_type << endl;
                exit(0);
          }
          return false;
}

bool File::ReadDeltaR( Cell &cel )
{
          TITLE("File","ReadDeltaR");
          return ReadDeltaR_ABACUS(cel);
          return false;
}

bool File::CheckGeometry( Cell &cel)
{
          TITLE("File","CheckGeometry");
          return CheckGeometry_ABACUS(cel);
          return false;
}
