#ifndef FILE_H
#define FILE_H

#include "gfun.h"
#include "vec3.h"
#include "cell.h"

class File : public Cell
{
        public:
        
        File(){};
        ~File(){};
        
        static bool CheckGeometry(Cell &cel); 
        static bool ReadDeltaR(Cell &cel );
        static bool ReadGeometry(Cell &cel);
        
        private:

        //not each cell need to read the geometry
        static bool ReadGeometry_ABACUS( Cell &cel );
        static bool ReadDeltaR_ABACUS(Cell &cel );
        
};
        static bool CheckGeometry_ABACUS(Cell &cel);
#endif

