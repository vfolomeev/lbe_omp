#include "solver.h"
#include "QDebug"
#include "math_api.h"
#include "omp.h"
#include <time.h>
#include "QFile"
solver::solver()
{

}
void solver::iterate(double *x1,double *y1,int N){

    int const n=100,m=n; //number of lattice nodes in two dimentions
    double f[n][m][9],rho[n][m],feq,x[n],y[m],w[9],cu,rhon;
    Vector2d c[9];
    Vector2d  u[n][m],sumu;
    double dt=1.0,dy=1.0,dx=dy;
    double u0=0.1,rho0=5.;
    int i,j,k;
        int mstep=40000;//total number of steps
    x[0]=0;
    for(i=1;i<n;i++){
        x[i]=x[i-1]+dx;
    }
    y[0]=0;
    for(i=1;i<m;i++){
        y[i]=y[i-1]+dy;
    }
    double alpha=0.01;
    double Re=u0*m/alpha;
    double sum;
    //double csq=dx*dx/(dt*dt);

    double omega=1.0/(3.*alpha+0.5);
// start wall clock
   double wall_timer = omp_get_wtime();
   clock_t clock_timer = clock();

//    w[0]=4./9;
//    w[1]=w[2]=w[3]=w[4]=1./9;
//    w[5]=w[6]=w[7]=w[8]=1./36;

//    c[0]=Vector2d();
//    c[1]=Vector2d(1,0);c[2]=Vector2d(0,1);c[3]=Vector2d(-1,0);c[4]=Vector2d(0,-1);
//    c[5]=Vector2d(1,1);c[6]=Vector2d(-1,1);c[7]=Vector2d(-1,-1);c[8]=Vector2d(1,-1);
    //init c vectors
    for(i=-1;i<2;i++){
        for(j=-1;j<2;j++){
            k=ind(i,j);
            c[k]=Vector2d(i,j);
            if((c[k]*c[k])==0) w[k]=4./9;
            if(c[k]*c[k]==1) w[k]=1./9;
            if(c[k]*c[k]==2) w[k]=1./36;
            
        }
    }


omp_set_num_threads(4);

    //initiate solution
#pragma omp parallel for private(i,j,k)
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            rho[i][j]=rho0;
            u[i][j]=Vector2d();
            for(k=0;k<9;k++){
                f[i][j][k]=w[k]*rho0;// zero initial velocity

            }

        }
    }

    //main loop
    for(int t=0;t<mstep;t++){
       // qDebug()<<t;
#pragma omp parallel for private(i,j,k,sum,sumu,feq,cu)
        for(i=0;i<n;i++){
            for(j=0;j<m;j++){
                //compute rho and u
                sum=0.;
                sumu=Vector2d();
                for(k=0;k<9;k++){
                    sum=sum+f[i][j][k];
                    sumu=sumu+f[i][j][k]*c[k];
                }
                rho[i][j]=sum;
                u[i][j]=sumu*(1/rho[i][j]);
                //collision
                for(k=0;k<9;k++){
                    cu=c[k]*u[i][j];
                    feq=w[k]*rho[i][j]*(1+3*cu+4.5*cu*cu-1.5*(u[i][j]*u[i][j]));
                    f[i][j][k]=omega*feq+(1.-omega)*f[i][j][k];
               }
           }

       }

        //streaming
#pragma omp parallel for private(i,j)
        for(k=0;k<9;k++){
            if(k!=ind(0,0)){
                auto run_contraction=Vector2d(c[k].x*c[k].x,c[k].y*c[k].y);
                auto run_dir=Vector2d(c[k].x>0?-1:1,c[k].y>0?-1:1);
                auto run_start=Vector2d(c[k].x>0?n-1:0,c[k].y>0?m-1:0);
                auto run_shift=(-1)*c[k];
//                qDebug()<<"Ck   "<<c[k];
//                qDebug()<<"run_contraction  "<<run_contraction;
//                qDebug()<<"run_dir  "<<run_dir;
//                qDebug()<<"run_start  "<<run_start;
                for(i=0;i<(n-run_contraction.x);i++){
                         for(j=0;j<(m-run_contraction.y);j++){
//                             qDebug()<<(int)(i*run_dir.x+run_start.x)<<"    "<<
//                             (int)(run_start.x+run_dir.x*(i+1));
                            f[(int)(i*run_dir.x+run_start.x)][(int)(j*run_dir.y+run_start.y)][k]=
                            f[(int)(run_start.x+run_dir.x*i+run_shift.x)][(int)(run_start.y+run_dir.y*j+run_shift.y)][k];
                            }
                         }
                
                }
            
            
            
        }
        
//#pragma omp parallel for private(i,j)
//        for(i=0;i<n;i++){
//            for(j=0;j<(m-1);j++){
//                f[j][i][3]=f[j+1][i][3];
//                f[n-j-1][i][1]=f[n-j-2][i][1];
//                f[i][j][4]=f[i][j+1][4];
//                f[i][m-j-1][2]=f[i][m-j-2][2];
//            }
//        }
//        for(i=0;i<(n-1);i++){
//            for(j=0;j<(m-1);j++){
//                f[n-i-1][m-j-1][5]=f[n-i-2][m-j-2][5];
//                f[i][m-j-1][6]=f[i+1][m-j-2][6];
//                f[i][j][7]=f[i+1][j+1][7];
//                f[n-i-1][j][8]=f[n-i-2][j+1][8];
//            }
//        }
        //boundary conditions
#pragma omp parallel for private(i,j,rhon)
        for(i=0;i<m;i++){
            for(j=-1;j<2;j++){
            //x=0 bounce back
                f[0][i][ind(1,j)]=f[0][i][ind(-1,-j)];

            //x=l bounce back
                f[n-1][i][ind(-1,j)]=f[n-1][i][ind(1,-j)];

            //y=0 bounce back
                f[i][0][ind(j,1)]=f[i][0][ind(-j,-1)];
            }
            //y=l moving lid
            rhon=f[i][m-1][ind(0,0)]+f[i][m-1][ind(1,0)]+f[i][m-1][ind(-1,0)]+2.*(f[i][m-1][ind(0,1)]+f[i][m-1][ind(-1,1)]+f[i][m-1][ind(1,1)]);
            f[i][m-1][ind(-1,-1)]=f[i][m-1][ind(1,1)]-rhon*u0/6.;
            f[i][m-1][ind(0,-1)]=f[i][m-1][ind(0,1)];
            f[i][m-1][ind(1,-1)]=f[i][m-1][ind(-1,1)]+rhon*u0/6.;

        }

    }

   qDebug()<<   " time on wall: " <<  omp_get_wtime() - wall_timer << "\n";
   //writing results file
   QFile file("D:\\Development\\2d_convection_omp\\res.csv");
   file.open(QIODevice::WriteOnly);
   QTextStream results(&file);
   results<<"X  Y   Z   u   v   w \n";
   for(i=0;i<n;i++){
       for(j=0;j<m;j++){
           results<<i<<"    "<<j<<" 0   "<<u[i][j]<<"  0"<<"\n";
       }

   }
   file.flush();
   file.close();

    for (i=0; i<N; i++)//Пробегаем по всем точкам
    {


            x1[i]=x[i];
            y1[i] =u[i][50].y;// rho[i][0];//Формула нашей функции

    }

}
