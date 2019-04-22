#include"msd.h"
#include "cellFile.h"
#include "input.h"
#include "math.h"


void MSD::Routine()
{
         TITLE("MSD","Routine");

         cal();

         return;
}

void MSD::cal()
{
         TITLE("MSD","cal");

         // ionic density
         assert(INPUT.geo_interval>0);
         int count_geometry_number=0;
         
         double* msd_tao=new double[INPUT.ntype];

         int nstru=(INPUT.geo_2-INPUT.geo_1)/INPUT.geo_interval+1;
         Vector3<double>** dr=new Vector3<double>*[nstru];
         Vector3<double>* drtot=new Vector3<double>[INPUT.natom];
         Vector3<double>* drtot1=new Vector3<double>[INPUT.natom];
         Vector3<double>* drtot2=new Vector3<double>[INPUT.natom];
         Vector3<double> a1,a2,a3;
         int istru=0;
         for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
         {
                if(igeo%INPUT.geo_interval!=0) continue;

                CellFile cel;

                stringstream ss; ss<<igeo;

                cel.file_name = ss.str();

                // cel : input geometry file 
                if( !CellFile::ReadGeometry( cel ) ) continue;
                ++count_geometry_number;

                dr[istru]=new Vector3<double>[INPUT.natom];

                // calculate the ionic density 
                const double rho_ion =INPUT.natom /cel.volume;

                if(count_geometry_number==1)
                {
                 //      cout<< " Volume of the input cell = " << cel.volume << " A^3 " << endl;
                 //      cout << " Average ion density = " << rho_ion <<endl;
                      
                       save(dr[istru++],cel); 
                       a1=cel.a1;
                       a2=cel.a2;
                       a3=cel.a3;
                }
                else if(count_geometry_number>1)
                {
                       delta(dr[istru-1],cel);
                       save(dr[istru],cel);
                       //cout<<"ok! "<<istru<<" "<<dr[istru][100].x<<" "<<dr[istru-1][100].x<<endl;
                       istru++;
                }
         }
         cout<< "count_geometry_number: " <<count_geometry_number<<endl;

         int tao=1;
         int t,j,k,i;
         int replica=0;
         replica=INPUT.msd_replica;
         if(replica==0) replica = count_geometry_number/2;
         cout<< "replica" << replica<<endl;
         Vector3<double> dr_cart;   
         
         ofstream ofs(INPUT.geo_out.c_str());
        
         if(count_geometry_number>tao)         
                for(tao=1;tao<count_geometry_number-replica;tao++)
                {
                       for(j=0;j<INPUT.ntype;j++)
                              msd_tao[j]=0;
                       for(t=0;t < replica;t++)
                       {
                      /*
                              for(k=0;k<INPUT.natom;k++)
                              {
                                    drtot[k].x=0;
                                    drtot[k].y=0;
                                    drtot[k].z=0;
                              }*/
                              if(t==0){
                                if(tao==1)
                                    for(k=0;k<INPUT.natom;k++){
                                          drtot1[k] = dr[0][k];
                                          drtot2[k] = drtot1[k];
                                    }
                                else for(k=0;k<INPUT.natom;k++){
                                          drtot1[k] += dr[tao-1][k];
                                          drtot2[k] =drtot1[k];
                                }
                              }
                              else for(k=0;k<INPUT.natom;k++)
                                    drtot2[k]+=dr[tao+t-1][k]-dr[t-1][k];
                              /*for(j=t;j<tao+t;j++)
                                    for(k=0;k<INPUT.natom;k++)
                                          drtot[k]+=  dr[j][k];*/
                              //if(tao==2&&t==3)cout<<drtot[0].x<<" "<<drtot[0].y<<" "<<drtot[0].z<<endl;
                              k=0;
                              for(j=0;j<INPUT.ntype;j++)
                                    for(i=0;i<INPUT.atom[j];i++)
                                    {
                                          //change to cartesian
                                          dr_cart.x = drtot2[k].x * a1.x + drtot2[k].y * a2.x + drtot2[k].z * a3.x;
                                          dr_cart.y = drtot2[k].x * a1.y + drtot2[k].y * a2.y + drtot2[k].z * a3.y;
                                          dr_cart.z = drtot2[k].x * a1.z + drtot2[k].y * a2.z + drtot2[k].z * a3.z;
                                          msd_tao[j]+=pow(dr_cart.x,2)+pow(dr_cart.y,2)+pow(dr_cart.z,2);
                                          k++;
                                          //test line
                               //           if(tao==2&&t==3&&j==0)cout<<msd_tao[j]<<" "<<drtot[k-1].x<<" "<<drtot[k-1].y<<" "<<drtot[k-1].z<<endl;
                                    }
                       }
                       ofs << tao << " " ;
                       for(j=0;j<INPUT.ntype;j++)
                       {
                              msd_tao[j]=msd_tao[j]/replica/INPUT.atom[j];
                              ofs << msd_tao[j]<< " ";
                       }
                       ofs << " " << endl ;
                }
          ofs.close();
          if(count_geometry_number > 0) L_R(count_geometry_number - replica - 1);
//deallocate
          delete[] msd_tao;
          delete[] drtot;
          delete[] drtot1;
          delete[] drtot2;
          for(int i=0;i<count_geometry_number;i++){
              delete[] dr[i];
          }
          delete[] dr;

          return;
}

void MSD::save(Vector3<double>* vec, const Cell &cel)
{
          int j=0;
          for(int it=0; it<INPUT.ntype; ++it)      
                for(int ia=0; ia<cel.atom[it].na; ++ia)
                {
                        vec[j++]=cel.atom[it].posd[ia];
                }
          return;      
}

void MSD::delta(Vector3<double>* vec, const Cell &cel)
{
          int ia,it,j=0;
          for(it=0; it<INPUT.ntype; it++)
                for(ia=0; ia<cel.atom[it].na; ia++)
                {
                        vec[j].x = cel.atom[it].posd[ia].x-vec[j].x;
                        vec[j].y = cel.atom[it].posd[ia].y-vec[j].y;
                        vec[j].z = cel.atom[it].posd[ia].z-vec[j].z;
                        j++;
                }
          for(j=0;j<INPUT.natom;j++)
          {
                if(vec[j].x>0.5) vec[j].x -= 1;
                else if(vec[j].x<-0.5) vec[j].x += 1;
                if(vec[j].y>0.5) vec[j].y -= 1;
                else if(vec[j].y<-0.5) vec[j].y += 1;
                if(vec[j].z>0.5) vec[j].z -= 1;
                else if(vec[j].z<-0.5) vec[j].z += 1;
                if(vec[j].x>0.5||vec[j].y>0.5||vec[j].z>0.5||vec[j].z<-0.5||vec[j].y<-0.5||vec[j].x<-0.5) QUIT;
          }
          return;
}
void MSD::L_R(int num){
          ifstream ifs(INPUT.geo_out.c_str());
          double *xi,*y1,*y2,*y3,*y4,*y5;
          double x0 = 0,y10=0,y20=0,y30=0,y40=0,y50 = 0;
          double x1 = 0,y11 = 0,y21 = 0,y31 = 0,y41 = 0,y51 = 0;
          double y12 = 0,y22 = 0,y32 = 0,y42 = 0,y52 = 0;
          xi = new double[num];
          y1 = new double[num];
          y2 = new double[num];
          y3 = new double[num];
          y4 = new double[num];
          y5 = new double[num];
          if(INPUT.ntype==1)
           for(int i=0;i<num;i++){
             xi[i] = 0;
             y1[i] = 0;
             ifs >> xi[i] >>y1[i];
           }
          if(INPUT.ntype==2)
           for(int i=0;i<num;i++){
             xi[i] = 0;
             y1[i] = 0;
             y2[i] = 0;
             ifs >> xi[i] >> y1[i] >> y2[i];
           }
          if(INPUT.ntype==3)
           for(int i=0;i<num;i++){
             xi[i] = 0;
             y1[i] = 0;
             y2[i] = 0;
             y3[i] = 0;
             ifs >> xi[i] >> y1[i] >> y2[i] >> y3[i];
           }
          if(INPUT.ntype==4)
           for(int i=0;i<num;i++){
             xi[i] = 0;
             y1[i] = 0;
             y2[i] = 0;
             y3[i] = 0;
             y4[i] = 0;
             ifs >> xi[i] >> y1[i] >> y2[i] >>y3[i]>> y4[i];
           }
          if(INPUT.ntype==5)
           for(int i=0;i<num;i++){
             xi[i] = 0;
             y1[i] = 0;
             y2[i] = 0;
             y3[i] = 0;
             y4[i] = 0;
             y5[i] = 0;
             ifs >> xi[i] >> y1[i] >> y2[i] >>y3[i]>> y4[i] >>y5[i];
           }
          ifs.close();
         if(INPUT.ntype > 0){
           for(int i=0;i<num;i++){
              x0 += xi[i];
              y10 += y1[i];
              x1 += xi[i] * xi[i];
              y11 += xi[i] * y1[i]; 
           }
           x0 /= num;
           y10 /= num;
           y12 = (y11 - num * x0 * y10) / (x1 - num * x0 * x0);
         }
         if(INPUT.ntype > 1){
           for(int i=0;i<num;i++){
              y20 += y2[i];
              y21 += xi[i] * y2[i];
           }
           y20 /= num;
           y22 = (y21 - num * x0 * y20) / (x1 - num * x0 * x0);
         }
         if(INPUT.ntype > 2){
           for(int i=0;i<num;i++){
              y30 += y3[i];
              y31 += xi[i] * y3[i];
           }
           y30 /= num;
           y32 = (y31 - num * x0 * y30) / (x1 - num * x0 * x0);
         }
         if(INPUT.ntype > 3){
           for(int i=0;i<num;i++){
              y40 += y4[i];
              y41 += xi[i] * y4[i];
           }
           y40 /= num;
           y42 = (y41 - num * x0 * y40) / (x1 - num * x0 * x0);
         }
          if(INPUT.ntype > 4){
           for(int i=0;i<num;i++){
              y50 += y5[i];
              y51 += xi[i] * y5[i];
           }
           y50 /= num;
           y52 = (y51 - num * x0 * y50) / (x1 - num * x0 * x0);
         }
         ofstream file;
         file.open("diffusion.txt");
         if(INPUT.ntype >0){
            file<<"Diffusion_of_element1: "<<y12/6<<endl;
         }
         if(INPUT.ntype >1){
            file<<"Diffusion_of_element2: "<<y22/6<<endl;
         }
         if(INPUT.ntype >2){
            file<<"Diffusion_of_element3: "<<y32/6<<endl;
         }
         if(INPUT.ntype >3){
            file<<"Diffusion_of_element4: "<<y42/6<<endl;
         }
         if(INPUT.ntype >4){
            file<<"Diffusion_of_element5: "<<y52/6<<endl;
         }
         file.close();
         return;
}
