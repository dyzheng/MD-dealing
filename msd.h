#ifndef MSD_H
#define MSD_H

#include "cell.h"

class MSD
{
      public:

      MSD(){};
      ~MSD(){};
     
      void Routine(); 
      
      private:
      void cal();
      void save(Vector3<double>* vec, const Cell &cel);
      void delta(Vector3<double>* vec, const Cell &cel);
      double diffusion();
      void L_R(int num);
};

#endif
      
