#include "solver.h"
#include "QDebug"
#include "math_api.h"
#include "omp.h"
#include <time.h>
solver::solver()
{

}
void solver::iterate(double *x1,double *y1,int N){

    int const n=100,m=n; //number of lattice nodes in two dimentions
    double f[n][m][9],rho[n][m],feq,x[n],y[m],w[9],cu,rhon;
    Vector2d c[9],u[n][m],sumu;
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

    w[0]=4./9;
    w[1]=w[2]=w[3]=w[4]=1./9;
    w[5]=w[6]=w[7]=w[8]=1./36;

    c[0]=Vector2d();
    c[1]=Vector2d(1,0);c[2]=Vector2d(0,1);c[3]=Vector2d(-1,0);c[4]=Vector2d(0,-1);
    c[5]=Vector2d(1,1);c[6]=Vector2d(-1,1);c[7]=Vector2d(-1,-1);c[8]=Vector2d(1,-1);


omp_set_num_threads(4);

    //initiate solution
#pragma omp parallel for private(i,j,k)
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            rho[i][j]=rho0;
            u[i][j]=Vector2d();
            for(k=0;k<9;k++){
                f[i][j][k]=w[k]*rho0;// zero initial velocity
                  qDebug()<<f[i][j][k];
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
        for(i=0;i<n;i++){
            for(j=0;j<(m-1);j++){
                f[j][i][3]=f[j+1][i][3];
                f[n-j-1][i][1]=f[n-j-2][i][1];
                f[i][j][4]=f[i][j+1][4];
                f[i][m-j-1][2]=f[i][m-j-2][2];
            }
        }
        for(i=0;i<(n-1);i++){
            for(j=0;j<(m-1);j++){
                f[n-i-1][m-j-1][5]=f[n-i-2][m-j-2][5];
                f[i][m-j-1][6]=f[i+1][m-j-2][6];
                f[i][j][7]=f[i+1][j+1][7];
                f[n-i-1][j][8]=f[n-i-2][j+1][8];
            }
        }
        //boundary conditions
#pragma omp parallel for private(i,rhon)
        for(i=0;i<m;i++){
            //x=0 bounce back
            f[0][i][1]=f[0][i][3];
            f[0][i][5]=f[0][i][7];
            f[0][i][8]=f[0][i][6];
            //x=l bounce back
            f[n-1][i][3]=f[n-1][i][1];
            f[n-1][i][7]=f[n-1][i][5];
            f[n-1][i][6]=f[n-1][i][8];
            //y=0 bounce back
            f[i][0][5]=f[i][0][7];
            f[i][0][2]=f[i][0][4];
            f[i][0][6]=f[i][0][8];
            //y=l moving lid
            rhon=f[i][m-1][0]+f[i][m-1][1]+f[i][m-1][3]+2.*(f[i][m-1][2]+f[i][m-1][6]+f[i][m-1][5]);
            f[i][m-1][7]=f[i][m-1][5]-rhon*u0/6.;
            f[i][m-1][4]=f[i][m-1][2];
            f[i][m-1][8]=f[i][m-1][6]+rhon*u0/6.;

        }

    }

   qDebug()<<   " time on wall: " <<  omp_get_wtime() - wall_timer << "\n";
   //writing results file

    for (i=0; i<N; i++)//Пробегаем по всем точкам
    {


            x1[i]=x[i];
            y1[i] =u[i][50].y;// rho[i][0];//Формула нашей функции

    }

}
