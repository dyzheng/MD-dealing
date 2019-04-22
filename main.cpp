#include "gfun.h" 
#include "input.h"
#include "ext.h"
#include "pdf.h"
#include "pdf2d.h"
#include "msd.h"

int main(int argc, char **argv)
{
       //read in the parameters.
       //INPUT has been generated in input.cpp
       INPUT.Init("INPUT",argc);

       // The program has several subroutines.
       // Each subroutine has its unit function.
       // Until now, there are just PDF and MSD can be choosed.

       if(INPUT.calculation == "pdf"){PDF pdf; pdf.Routine();}
       else if(INPUT.calculation == "msd"){MSD msd; msd.Routine();}
       else 
       {
               cout << " calculation=" <<INPUT.calculation << endl;
               QUIT("No 'calculation' available");
       }
    
       cout << " -------------- " << endl;
       cout << "     Finish     " << endl;
       cout << " -------------- " << endl;
      
       return 0;
}
