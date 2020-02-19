
#include <stdio.h>
#include <stdlib.h>

#include "phi_mult.h"
#include "globals.h"
#include <math.h>





double ddx(double var_[N_X][N_Y],int i, int j)
{

    double res;

    if (i==0)
    {
        res=(-3.0*var_[i][j] + 4.0*var_[i+1][j] - var_[i+2][j])/(dx*2.0);
      //  res=(var_[1][j]-var_[N_X-2][j])/(2.0*dx);

    }else
    {

        if (i==N_X-1)
        {

            res=(3.0*var_[i][j] - 4.0*var_[i-1][j] + var_[i-2][j])/(dx*2.0);
           // res=(var_[1][j]-var_[N_X-2][j])/(2.0*dx);
        }else
        {

            res=(var_[i+1][j]-var_[i-1][j])/(2.0*dx);

        }


    }


    return res;



}


double ddy(double var_[N_X][N_Y],int i, int j)
{

    double res;

    if (j==0)
    {
        res=(-3.0*var_[i][j] + 4.0*var_[i][j+1] - var_[i][j+2])/(2.0*dy);

    }else
    {

        if (j==N_Y-1)
        {

            res=(3.0*var_[i][j] - 4.0*var_[i][j-1] + var_[i][j-2])/(2.0*dy);
        }else
        {

            res=(var_[i][j+1]-var_[i][j-1])/(2.0*dy);

        }


    }


    return res;



}

double field_U(double U,double out_[N_X][N_Y],double itn)
{
    int i,j,n;
    double res=0;




    for (int j=0;j<N_Y; j++)
    {
        for (int i=0; i < N_X; i++)
        {
            // out_[i][j]=1.0*(dy*(N_Y-j))*(dy*(N_Y-j))*U;
            // res= out_[i][j];
            // out_Uy[i][j]=sqrt(U*U-out_Ux[i][j]*out_Ux[i][j]);
        }

        // printf("j= %d U=%f \n",j, res);
    }

    for (int j=0;j<N_Y-1; j++)
    {

        out_[0][j]=out_[N_X-2][j];
        out_[N_X-1][j]=out_[1][j];
        for (int i=1; i < N_X-1; i++)
        {
            out_[i][j]=dy*(N_Y-j)*dy*(N_Y-j)*U;
            //res= out_[i][j];
            // out_Uy[i][j]=sqrt(U*U-out_Ux[i][j]*out_Ux[i][j]);
        }



        // printf("j= %d U=%f \n",j, res);
    }


}


double multigrid(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn)
{
    int i,j,l,k;
    double a,b_p,b_m,c_p,c_m;
    double res=0;
    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);
    static double rhs_2[N_X][N_Y];
    static double field2[N_X][N_Y];



    for (int n=0; n<itn; n++)
    {
        for (i=1; i<N_X-1; i++)
        {
            for (j=1;j<N_Y-1;j++)
            {
                if (j==1)
                {
                    if (par.s_bc_type==0)//fixed_value
                        field[i][0]=par.s_bc_val;
                    if (par.s_bc_type==1)//fixed_gradient
                        field[i][0]=field[i][1]-dy*par.s_bc_val;
                    if (par.s_bc_type==2)//cyclic
                        field[i][0]=field[i][N_Y-2];


                }

                if (j==N_Y-2)
                {
                    if (par.n_bc_type==0)//fixed_value
                        field[i][N_Y-1]=par.n_bc_val;
                    if (par.n_bc_type==1)//fixed_gradient
                        field[i][N_Y-1]=field[i][N_Y-2]+dy*par.n_bc_val;
                    if (par.n_bc_type==2)//cyclic
                        field[i][N_Y-1]=field[i][1];
                }

                if (i==1)
                {
                    if (par.w_bc_type==0)//fixed_value
                        field[0][j]=par.w_bc_val;
                    if (par.w_bc_type==1)//fixed_gradient
                        field[0][j]=field[1][j]-dx*par.w_bc_val;
                    if (par.w_bc_type==2)//cyclic
                        field[0][j]=field[N_X-2][j];
                }

                if (i==N_X-2)
                {
                    if (par.e_bc_type==0)//fixed_value
                        field[N_X-1][j]=par.e_bc_val;
                    if (par.e_bc_type==1)//fixed_gradient
                        field[N_X-1][j]=field[N_X-2][j]+dx*par.e_bc_val;
                    if (par.e_bc_type==2)//cyclic
                        field[N_X-1][j]=field[1][j];
                }

                field[i][j]=(rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))/a;

            }
        }

    }//itn1


    for (i=1; i<N_X-1; i++)
    {
        for (j=1;j<N_Y-1;j++)
        {
            rhs_2[i][j]=rhs[i][j] - (a*field[i][j]+b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]);

            field2[i][j]=0.0;
        }
    }


    for (int n=0; n<itn; n++)
    {
        for (k=1; k< (N_X-1)/2;k++)
        {
            i=k;
            for (l=1; l<(N_Y-1)/2;l++)
            {
                j=l;


                if (j==1)
                {
                    if (par.s_bc_type==0)//fixed_value
                        field2[i][j-1]=par.s_bc_val;
                    if (par.s_bc_type==1)//fixed_gradient
                        field2[i][j-1]=field2[i][1]-dy*par.s_bc_val;
                    if (par.s_bc_type==2)//cyclic
                        field2[i][j-1]=field2[i][(N_Y-1)/2-1];
                }

                if (j==(N_Y-1)/2-1)
                {
                    if (par.n_bc_type==0)//fixed_value
                        field2[i][j+1]=par.n_bc_val;
                    if (par.n_bc_type==1)//fixed_gradient
                        field2[i][j+1]=field2[i][j]+dy*par.n_bc_val;
                    if (par.n_bc_type==2)//cyclic
                        field2[i][j+1]=field2[i][1];
                }

                if (i==1)
                {
                    if (par.w_bc_type==0)//fixed_value
                        field2[i-1][j]=par.w_bc_val;
                    if (par.w_bc_type==1)//fixed_gradient
                        field2[i-1][j]=field2[1][j]-dx*par.w_bc_val;
                    if (par.w_bc_type==2)//cyclic
                        field2[i-1][j]=field2[(N_X-1)/2-1][j];
                }

                if (i==(N_X-1)/2-1)
                {
                    if (par.e_bc_type==0)//fixed_value
                        field2[i+1][j]=par.e_bc_val;
                    if (par.e_bc_type==1)//fixed_gradient
                        field2[i+1][j]=field2[i][j]+dx*par.e_bc_val;
                    if (par.e_bc_type==2)//cyclic
                        field2[i+1][j]=field2[1][j];
                }

                field2[i][j]=(rhs_2[i*2][j*2]-(b_p*field2[i+1][j]+b_m*field2[i-1][j]+c_p*field2[i][j+1]+c_m*field2[i][j-1]))/a;

            }
        }
    }

    for(k=1; k<(N_X-1)/2; k++)
    {
        i=k*2;
        for (l=1; l<(N_Y-1)/2; l++)
        {
            j=l*2;

            field[i][j]+=field2[k][l];
            field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
            field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
            field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
        }

        l=(N_Y-1)/2; j=l*2;

        field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);

    }

    k=(N_X-1)/2;
    i=k*2;

    for (l=1;l<(N_Y-1)/2;l++)
    {
        j=l*2;

        field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
    }

    l=(N_Y-1)/2; j=l*2;

    field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);


    //res+= fabs(field[i][j]-field_new_1[i][j]);

    return res;
}


double jacobi_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn,int nx,int ny)
{
    int i,j,n;
    double res=0;
    double eps =1.0;
    double a,b_p,b_m,c_p,c_m;
    a=par.a;
    b_p=par.bp;
    b_m=par.bm;
    c_p=par.cp;
    c_m=par.cm;
    for(n=0;n<itn;n++)
    {
        for (i=1; i<nx; i++)
        {
            for (j=1; j<ny; j++)
            {
                if (j==1)
                {
                    if (par.s_bc_type==0)//fixed_value
                        field[i][j-1]=par.s_bc_val;
                    if (par.s_bc_type==1)//fixed_gradient
                        field[i][j-1]=field[i][1]-dy*par.s_bc_val;
                    if (par.s_bc_type==2)//cyclic
                        field[i][j-1]=field[i][ny-1];

                    /*if (i<nx/4) field[i][j-1]=field[i][1];*/
                }

                if (j==ny-1)
                {
                    if (par.n_bc_type==0)//fixed_value
                        field[i][j+1]=par.n_bc_val;
                    if (par.n_bc_type==1)//fixed_gradient
                        field[i][j+1]=field[i][j]+dy*par.n_bc_val;
                    if (par.n_bc_type==2)//cyclic
                        field[i][j+1]=field[i][1];
                  //  if (i>nx/4) field[i][j+1]=field[i][j];
                }

                if (i==1)
                {
                    if (par.w_bc_type==0)//fixed_value
                        field[i-1][j]=par.w_bc_val;
                    if (par.w_bc_type==1)//fixed_gradient
                        field[i-1][j]=field[1][j]-dx*par.w_bc_val;
                    if (par.w_bc_type==2)//cyclic
                        field[i-1][j]=field[nx-1][j];


                     if (j<N_Y_DIEL*nx/(N_X-1)) field[i-1][j]=field[i][j];

                }

                if (i==nx-1)
                {
                    if (par.e_bc_type==0)//fixed_value
                        field[i+1][j]=par.e_bc_val;
                    if (par.e_bc_type==1)//fixed_gradient
                        field[i+1][j]=field[i][j]+dx*par.e_bc_val;
                    if (par.e_bc_type==2)//cyclic
                        field[i+1][j]=field[1][j];
                }

                double f_xm=field[i-1][j];
                double f_xp=field[i+1][j];

                //if (i==1) f_xm=field[nx-1][j];
                //if (i==nx-1) f_xp=field[1][j];

                //field[i][j]=(rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))
                  //      /a;


                field[i][j]=(rhs[i][j]-(b_p*f_xp+b_m*f_xm+c_p*field[i][j+1]+c_m*field[i][j-1]))
                        /a;


            }
        }
    }
    return 0;
}


double interp_up(double field[N_X][N_Y], double field2[N_X][N_Y], int nx,int ny)
{
    int i,j,k,l;

    for(k=1; k<(nx)/2; k++)
    {
        i=k*2;
        for (l=1; l<(ny)/2; l++)
        {
            j=l*2;
            field[i][j]+=field2[k][l];
            field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
            field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
            field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
        }

        l=(ny)/2; j=l*2;

        field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
    }

    k=(nx)/2;
    i=k*2;

    for (l=1;l<(ny)/2;l++)
    {
        j=l*2;
        field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
    }

    l=(ny)/2; j=l*2;
    field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
}

double multigrid_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn,int N)
{
    int i,j,l,k,nn;
    double a,b_p,b_m,c_p,c_m;
    double res=0;
    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);
    static double rhs_[10][N_X][N_Y];
    static double field_[10][N_X][N_Y];



    //  int N=3;

    for (i=0; i<N_X; i++)
    {
        for (j=0;j<N_Y;j++)
        {
            field_[0][i][j]=field[i][j];
            rhs_[0][i][j]=rhs[i][j]/*-rho_[i][j]*/;
        }
    }

    for (nn=0;nn<N;nn++)
    {
        jacobi_N(par,field_[nn],rhs_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));

        for (i=1; i<(N_X-1)/pow(2,nn+1); i++)
        {
            for (j=1;j<(N_Y-1)/pow(2,nn+1);j++)
            {
                k=i*2; l=j*2;
                rhs_[nn+1][i][j]=(rhs_[nn][k][l] - (a*field_[nn][k][l]+b_p*field_[nn][k+1][l]+b_m*field_[nn][k-1][l]+c_p*field_[nn][k][l+1]+c_m*field_[nn][k][l-1]))*4.0;

                field_[nn+1][i][j]=0.0;
            }
        }
    }

    jacobi_N(par,field_[N],rhs_[N],itn*3,(N_X-1)/pow(2,N),(N_Y-1)/pow(2,N));

    // interp_up(field_[0],field_[1],(N_X-1),(N_Y-1));

    for (nn=N-1;nn>=0;nn--)
    {
        interp_up(field_[nn],field_[nn+1],(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));
        jacobi_N(par,field_[nn],rhs_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));
    }


    for (i=0; i<N_X; i++)
    {
        for (j=0;j<N_Y;j++)
        {
            field[i][j]=field_[0][i][j];
        }
    }
    return res;
}

double jacobi(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    static double field_new[N_X][N_Y];

    double a,b_p,b_m,c_p,c_m;
    /*  a=((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=-1.0/(dx*dx);
    b_m=-1.0/(dx*dx);
    c_p=-1.0/(dy*dy);
    c_m=-1.0/(dy*dy);
*/

    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);


    for(n=0;n<itn;n++)
    {
        for (i=1; i<N_X-1; i++)
        {
            //field[i][0]=1;
            for (j=1; j<N_Y-1; j++)
            {
                //p[i][j]=0;


                field_new[i][j]=field[i][j];
                field[i][j]=(rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))/a;

            }

        }

        //0--fixed value, 1--fixed gradient,2 --cyclic, 3 --init
        for (int j=0; j< N_Y; j++)
        {
            if (par.w_bc_type==0)//fixed_value
                field[0][j]=par.w_bc_val;
            if (par.w_bc_type==1)//fixed_gradient
                field[0][j]=field[1][j]-dx*par.w_bc_val;
            if (par.w_bc_type==2)//cyclic
                field[0][j]=field[N_X-2][j];
            // if (par.w_bc_type==3)// init


            if (par.w_bc_type==4)//cyclic
            {
                if (j<N_Y/2)
                field[0][j]=field[1][j];

            }


            if (par.e_bc_type==0)//fixed_value
                field[N_X-1][j]=par.e_bc_val;
            if (par.e_bc_type==1)//fixed_gradient
                field[N_X-1][j]=field[N_X-2][j]+dx*par.e_bc_val;
            if (par.e_bc_type==2)//cyclic
                field[N_X-1][j]=field[1][j];

            // if (par.e_bc_type==3)// init

        }

        for (int i=0; i<N_X; i++ )
        {
            if (par.n_bc_type==0)//fixed_value
                field[i][N_Y-1]=par.n_bc_val;
            if (par.n_bc_type==1)//fixed_gradient
                field[i][N_Y-1]=field[i][N_Y-2]+dy*par.n_bc_val;
            if (par.n_bc_type==2)//cyclic
                field[i][0]=field[i][N_Y-1];
            //  if (par.n_bc_type==3)// init

            if (par.s_bc_type==0)//fixed_value
                field[i][0]=par.s_bc_val;
            if (par.s_bc_type==1)//fixed_gradient
                field[i][0]=field[i][1]-dy*par.s_bc_val;
            if (par.s_bc_type==2)//cyclic
                field[i][N_Y-1]=field[i][1];
            // if (par.s_bc_type==3)// init

        }

    }

    //for (int i=1; i<N_X-1; i++)
    for (int j=1; j<N_Y-1; j++)
    {

        res+=fabs(field[i][j]-field_new[i][j]);


    }


    // printf("res=%f \n",res);
    return res;
}


/////


double calc_poly(poly &p, double r, double x) //возвращает значение функции f(x) = x^2-2
{
    double sum=0.0;
    double xn=x;
    for (int i=0;i<p.order;++i)
    {
        sum+=xn*p.C[i];
        xn*=x;
    }

       return sum-r;
}

float calc_d_poly(poly &p, double x) //возвращает значение производной
{
    double sum=p.C[0];
    double xn=x;
    for (int i=0;i<p.order-1;++i)
    {
        sum+=xn*p.C[i+1]*(i+2);
        xn*=x;
    }
   return sum;
}

float calc_d2_poly(poly &p, double x) // значение второй производной
{
    double sum=2.0*p.C[1];
    double xn=x;
    for (int i=0;i<p.order-2;++i)
    {
        sum+=xn*p.C[i+1]*(i+2)*(i+3);
        xn*=x;
    }
   return sum;
}

double solve_poly(poly &p,double _x0, double rhs,int itn)
{
    double x0,x,xn,eps;// вычисляемые приближения для корня
    double a, b,deltax;// границы отрезка и необходимая точность

   // deltax=deltax0;
    x=_x0;
   /* a=x0-deltax;
    b=x0+deltax;

    int nn=0;
    double sucs=f(a)*f(b);
    while ((sucs>0)&&(nn<10)) // если знаки функции на краях отрезка одинаковы
    {
        deltax=deltax*2;
        a=x0-deltax;
        b=x0+deltax;
        nn++;
        sucs=f(a)*f(b);
    }

    if (sucs>0)
    {
        printf("DIVERGENCEE in POLY \n",)
        return -1e30
    }
*/
    for (int i=0; i<itn; ++i)
    {
      x = x-calc_poly(p,rhs,x)/calc_d_poly(p,x); // считаем первое приближение
    }

       // printf("x = %lf f(x)=%e \n",x*1000.0,calc_poly(p,rhs,x));

    return x;

}

int jacobi_polynomial(INPUT_PARAM par, poly pol,double field[N_X][N_Y],double rhs[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    //static poly pol_new[N_X][N_Y];

    poly poly_new;

    double a,b_p,b_m,c_p,c_m;

    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);

    static double rhs_[N_X][N_Y];

//field_x*a+....+pol*fieldx^..=rhs

    poly_new=pol;
    poly_new.C[0]+=a;



    for(n=0;n<itn;n++)
    {
        for (i=1; i<N_X-1; i++)
        {
            //field[i][0]=1;
            for (j=shift+1; j<N_Y_DIEL-1; j++)
            {
                //p[i][j]=0;


                double f_xm=field[i-1][j];
                double f_xp=field[i+1][j];

               // if (i==1) f_xm=field[N_X-2][j];
               // if (i==N_X-2) f_xp=field[1][j];

                rhs_[i][j]=rhs[i][j]-(b_p*f_xp+b_m*f_xm+c_p*field[i][j+1]+c_m*field[i][j-1]);

                field[i][j]=field[i][j]*0.7+0.3*solve_poly(poly_new,field[i][j],rhs_[i][j],3);//(rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))/a;
            }
        }

        //0--fixed value, 1--fixed gradient,2 --cyclic, 3 --init
        for (int j=shift; j< N_Y_DIEL; j++)
        {
            if (par.w_bc_type==0)//fixed_value
                field[0][j]=par.w_bc_val;
            if (par.w_bc_type==1)//fixed_gradient
                field[0][j]=field[1][j];//-dx*par.w_bc_val;
            if (par.w_bc_type==2)//cyclic
                field[0][j]=field[N_X-2][j];
            // if (par.w_bc_type==3)// init


            if (par.e_bc_type==0)//fixed_value
                field[N_X-1][j]=par.e_bc_val;
            if (par.e_bc_type==1)//fixed_gradient
                field[N_X-1][j]=field[N_X-2][j]+dx*par.e_bc_val;
            if (par.e_bc_type==2)//cyclic
                field[N_X-1][j]=field[1][j];
            // if (par.e_bc_type==3)// init

        }

        for (int i=0; i<N_X; i++ )
        {
            if (par.n_bc_type==0)//fixed_value
                field[i][N_Y_DIEL-1]=par.n_bc_val;
            if (par.n_bc_type==1)//fixed_gradient
                field[i][N_Y_DIEL-1]=field[i][N_Y_DIEL-2]+dy*par.n_bc_val;
            if (par.n_bc_type==2)//cyclic
                field[i][0]=field[i][N_Y_DIEL-1];

            if (par.n_bc_type==5)//zero div
                field[i][N_Y_DIEL-1]=field[i+1][N_Y_DIEL-1]-dx*rhs[i][N_Y_DIEL-1];
            //  if (par.n_bc_type==3)// init

            if (par.s_bc_type==0)//fixed_value
                field[i][shift]=par.s_bc_val;
            if (par.s_bc_type==1)//fixed_gradient
                field[i][shift]=field[i][shift+1]-dy*par.s_bc_val;
            if (par.s_bc_type==2)//cyclic
                field[i][N_Y_DIEL-1]=field[i][shift+1];
            // if (par.s_bc_type==3)// init


            if (par.s_bc_type==5)//zero div
                field[i][0]=field[i+1][0]-dx*rhs[i][0];

           // j=N_Y_DIEL-1;
           // Py_[i][j]=Py_[i][j-1] -(Px_[i+1][j]-Px_[i][j])/(1.0*dx)*(1.0*dy);
            //Py_[i][j]=0.5*(Py_[i][j+1]+Py_[i][j-1]);

          /* j=1;
            Py_[i][j-1]=Py_[i][j+1] +(Px_[i+1][j]-Px_[i-1][j])/(2.0*dx)*(2.0*dy);
            Py_[i][j]=0.5*(Py_[i][j+1]+Py_[i][j-1]);*/


        }

    }

/*    //for (int i=1; i<N_X-1; i++)
    for (int j=1; j<N_Y-1; j++)
    {

        res+=fabs(field[i][j]-field_new[i][j]);


    }

*/
    // printf("res=%f \n",res);
    return 0;

}


void get_div(double f_x[N_X][N_Y],double f_y[N_X][N_Y],double out[N_X][N_Y])
{
    int i,j,k;
    double divmax,lapl_phi,div_1,dudx,dvdy,dwdy;
    divmax=0.0;
    for (j=0;j<N_Y;j++)
    {
        for (i=0;i<N_X;i++)
        {
            out[i][j]=0.0;
        }
    }

    for (j=0;j<N_Y;j++)
    {
        for (i=0;i<N_X;i++)
        {
            out[i][j]= ddx(f_x,i,j)+ddy(f_y,i,j);
            if (fabs(out[i][j])>divmax)
                divmax=fabs(out[i][j]);
        }
    }

    //  printf ("div0max=%lE \n",divmax);
}


double solve_current(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    double a,b_p,b_m,c_p,c_m;
    double dt_local=0.00015;

    D=0.4;
    b=4.0;
    a=(1/dt_local+(2.0*D)/(dx*dx)+(2.0*D)/(dy*dy));
    b_p=(-1.0*D)/(dx*dx);
    b_m=(-1.0*D)/(dx*dx);
    c_p=(-1.0*D)/(dy*dy);
    c_m=(-1.0*D)/(dy*dy);
    double sigma=1.0;

    static double rho_new[N_X][N_Y];
    static double rho_in_prev[N_X][N_Y];

    for (int n=0; n< itn; n++)
    {
        /*  for (int i=0; i<N_X; i++)
        {

            for (int j=0; j<N_Y; j++)
            {

                rho_in_prev[i][j]=rho_in[i][j];

            }
        }


        get_div(Jx_in,Jy_in,RHS);*/




        for (int i=1; i<N_X-1; i++)
        {

            rho_in[i][N_Y-1]= 0.0;//rho_in[i][0]*0.9+0.1*(-phi_[i][0]*dy+0.01*D*rho_in_prev[i][1])/(b*dy*ddy(phi_,i,0)+D);

            double dPhi=(phi_[i][1]-phi_[i][0]);
            rho_in[i][0]=-(sigma*dPhi + D*rho_in[i][1])/(b*dPhi+D);

            if (rho_in[i][0]<0)
                rho_in[i][0]=0.0;
            //rho_in[i][0]=-0.001;
            for (int j=1; j<N_Y-1; j++)
            {

                double dphidx=ddx(phi_,i,j);
                double dphidy=ddy(phi_,i,j);

                double d2phidx=(phi_[i+1][j]+phi_[i-1][j]-2.0*phi_[i][j])/(dx*dx);
                double d2phidy=(phi_[i][j+1]+phi_[i][j-1]-2.0*phi_[i][j])/(dy*dy);



                a=(/*1/dt_local+*/(2.0*D)/(dx*dx)+(2.0*D)/(dy*dy)+ b*(d2phidx+d2phidy));
                b_p=(-1.0*D)/(dx*dx)- b*dphidx/(2.0*dx);
                b_m=(-1.0*D)/(dx*dx)+ b*dphidx/(2.0*dx);
                c_p=(-1.0*D)/(dy*dy)- b*dphidy/(2.0*dy);
                c_m=(-1.0*D)/(dy*dy)+ b*dphidy/(2.0*dy);


                rho_in[i][j]=rho_in[i][j]*0.9+0.1*(-(b_p*rho_in[i+1][j]+b_m*rho_in[i-1][j]+c_p*rho_in[i][j+1]+c_m*rho_in[i][j-1]))/a;


                Jx_in[i][j]=rho_in[i][j]*(b*ddx(phi_,i,j))-D*ddx(rho_in,i,j);
                Jy_in[i][j]=rho_in[i][j]*(b*ddy(phi_,i,j))-D*ddy(rho_in,i,j);



            }

        }


        /*   for (int i=1; i<N_X-1; i++)
        {

            for (int j=1; j<N_Y-1; j++)
            {
               rho_in[i][j]=rho_in_prev[i][j];

            }
        }
*/

        for (int j=0; j<N_Y; j++)
        {
            rho_in[0][j]=rho_in[N_X-2][j];
            rho_in[N_X-1][j]=rho_in[1][j];
        }


    }




    double rho_max=-1000000000;
    for (int i=1; i<N_X-1; i++)
    {

        for (int j=1; j<N_Y-1; j++)
        {
            //    printf("%d %d \n",i,j);

            //rho_in[i][j]=rho_new[i][j];


            Jx_in[i][j]=rho_in[i][j]*(b*ddx(phi_,i,j))-D*ddx(phi_,i,j);
            Jy_in[i][j]=rho_in[i][j]*(b*ddy(phi_,i,j))-D*ddy(phi_,i,j);



            double dphidx=ddx(phi_,i,j);
            double dphidy=ddy(phi_,i,j);

            double d2phidx=(phi_[i+1][j]+phi_[i-1][j]-2.0*phi_[i][j])/(dx*dx);
            double d2phidy=(phi_[i][j+1]+phi_[i][j-1]-2.0*phi_[i][j])/(dy*dy);



            a=(/*1/dt_local+*/(2.0*D)/(dx*dx)+(2.0*D)/(dy*dy)+ b*(d2phidx+d2phidy));
            b_p=(-1.0*D)/(dx*dx)- b*dphidx/(2.0*dx);
            b_m=(-1.0*D)/(dx*dx)+ b*dphidx/(2.0*dx);
            c_p=(-1.0*D)/(dy*dy)- b*dphidy/(2.0*dy);
            c_m=(-1.0*D)/(dy*dy)+ b*dphidy/(2.0*dy);

            div_J[i][j]=a*rho_in[i][j]+b_p*rho_in[i+1][j]+b_m*rho_in[i-1][j]+c_p*rho_in[i][j+1]+c_m*rho_in[i][j-1];

            if (rho_in[i][j]>rho_max)
                rho_max=rho_in[i][j];
            /*
            if (div_J[i][j]>rho_max)
                rho_max=div_J[i][j];*/

        }
    }

    // get_div(Jx_in,Jy_in,div_J);

    printf("rho_max=%f \n",rho_max);

    return res;
}

double solve_current_1(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    double a,b_p,b_m,c_p,c_m;
    double dt_local=0.001;

    D=0.4;
    b=4.0;

    double D_1=2*10e-5;//2*10e-5;
    double D_2=4.86*10e-11;
    a=(1/dt_local+(2.0*D_2)/(dx*dx)+(2.0*D_2)/(dy*dy));
    b_p=(-1.0*D_2)/(dx*dx);
    b_m=(-1.0*D_2)/(dx*dx);
    c_p=(-1.0*D_2)/(dy*dy);
    c_m=(-1.0*D_2)/(dy*dy);
    double sigma=1.0;
    double mu_1=4.0*10e-7; // mobility
    double mu_2=2.0*10e-9;

    double beta=1.0*10e-13; //ion-electron recombination
    double alpha=4.34*10e3; // ionizatin's coefficient


   static double flow_1_x[N_X][N_Y];
   static  double flow_2_y[N_X][N_Y];
    static double flow_1_y[N_X][N_Y];
    static double flow_2_x[N_X][N_Y];
    static double div_flow_1[N_X][N_Y];
    static double div_flow_2[N_X][N_Y];
     static double mod_flow_1[N_X][N_Y];
    static double n_1_new[N_X][N_Y];
    static double n_2_new[N_X][N_Y];


   // static double rho_new[N_X][N_Y];
   // static double rho_in_prev[N_X][N_Y];

    for (int ii=0; ii< N_X;ii++ )
    {
        for (int jj=0; jj< N_Y;jj++)
        {
            flow_1_y[ii][jj]=0.0;
            flow_2_x[ii][jj]=0.0;
            flow_1_y[ii][jj]=0.0;
            flow_2_y[ii][jj]=0.0;
            div_flow_1[ii][jj]=0.0;
            div_flow_2[ii][jj]=0.0;

        }
    }
    for (int n=0; n< itn; n++)
    {



        for (int i=1; i<N_X-1; i++)
        {
            //boundary conditions
            rho_in[i][N_Y-1]= 0.0;//rho_in[i][0]*0.9+0.1*(-phi_[i][0]*dy+0.01*D*rho_in_prev[i][1])/(b*dy*ddy(phi_,i,0)+D);
            double dPhi=(phi_[i][1]-phi_[i][0]);
            rho_in[i][0]=-(sigma*dPhi + D_2*rho_in[i][1])/(mu_2*dPhi+D_2);

            if (rho_in[i][0]<0)
                rho_in[i][0]=0.0;
            //rho_in[i][0]=-0.001;
            for (int j=1; j<N_Y-1; j++)
            {

                double dphidx=ddx(phi_,i,j);
                double dphidy=ddy(phi_,i,j);

                double dn1dx=ddx(n_1,i,j);
                double dn1dy=ddy(n_1,i,j);

                double dn2dx=ddx(n_2,i,j);
                double dn2dy=ddy(n_2,i,j);

                alpha=0.00482*sqrt(dphidx*dphidx + dphidy*dphidy);

                double d2n1=(n_1[i+1][j]+n_1[i-1][j]-2.0*n_1[i][j])/(dx*dx)+(n_1[i][j+1]+n_1[i][j-1]-2.0*n_1[i][j])/(dy*dy) ;
                double d2n2=(n_2[i+1][j]+n_2[i-1][j]-2.0*n_2[i][j])/(dx*dx)+(n_2[i][j+1]+n_2[i][j-1]-2.0*n_2[i][j])/(dy*dy) ;



                /*a=(1/dt_local+(2.0*D_2)/(dx*dx)+(2.0*D_2)/(dy*dy)+ mu_2*(d2phidx+d2phidy));
                b_p=(-1.0*D_2)/(dx*dx)- mu_2*dphidx/(2.0*dx);
                b_m=(-1.0*D_2)/(dx*dx)+ mu_2*dphidx/(2.0*dx);
                c_p=(-1.0*D_2)/(dy*dy)- mu_2*dphidy/(2.0*dy);
                c_m=(-1.0*D_2)/(dy*dy)+ mu_2*dphidy/(2.0*dy);*/

                flow_1_x[i][j]=mu_1*n_1[i][j]*(dphidx);//-D_1*(ddx(n_1,i,j));//thermodiffusion
                flow_1_y[i][j]=mu_1*n_1[i][j]*(dphidy);//-D_1*(ddy(n_1,i,j));
                flow_2_x[i][j]=- mu_2*n_2[i][j]*(dphidx);//-D_2*(ddx(n_2,i,j));
                flow_2_y[i][j]=- mu_2*n_2[i][j]*(dphidy);//-D_2*(ddy(n_2,i,j));
                //mod_flow_1[i][j]=(sqrt(flow_1_x[i][j]*flow_1_x[i][j]+flow_1_y[i][j]*flow_1_y[i][j]));
              mod_flow_1[i][j]=sqrt((flow_1_x[i][j]-D_1*ddx(n_1,i,j))*(flow_1_x[i][j]-D_1*ddx(n_1,i,j))+(flow_1_y[i][j]-D_1*ddy(n_1,i,j))*(flow_1_y[i][j]-D_1*ddy(n_1,i,j)));


              double lapl_phi=(phi_[i+1][j]+phi_[i-1][j]-2.0*phi_[i][j])/(dx*dx)+
                              (phi_[i][j+1]+phi_[i][j-1]-2.0*phi_[i][j])/(dy*dy) ;


                n_1_new[i][j]=n_1[i][j]+dt_local*(-mu_1*(dn1dx*dphidx+dn1dy*dphidy+n_1[i][j]*lapl_phi)+D_1*d2n1+alpha*mod_flow_1[i][j]-beta*n_1[i][j]*n_2[i][j]);
                n_2_new[i][j]=n_2[i][j]+dt_local*(mu_2*(dn2dx*dphidx+dn2dy*dphidy+n_2[i][j]*lapl_phi)+D_2*d2n2+alpha*mod_flow_1[i][j]-beta*n_1[i][j]*n_2[i][j]);

                rho_in[i][j]=q_e*(n_1_new[i][j]-n_2_new[i][j]);




            }

        }

        for (int i=1; i< N_X-1; i++)
        {
            for(int j=1; j<N_Y-1;j++)
            {
                n_1[i][j]=n_1_new[i][j];
                n_2[i][j]=n_2_new[i][j];


                if (n_1[i][j]<0) n_1[i][j]=1e-10;
     if (n_2[i][j]<0) n_2[i][j]=1e-10;
             }
        }


        for (int i=0; i< N_X; i++)
        {

           double dphidy=(phi_[i][1]-phi_[i][0])/dy;//1.0*ddy(phi_,i,0);
            // flow_1_y[i][0]=mu_1*n_1[i][0]*(dphidy)-D_1*(n_1[i][1]-n_1[i][0])/dy;

            n_1[i][0]= 0;//n_1[i][0]*0.99+0.01* (D_1*(n_1[i][1])/dy)/(mu_1*dphidy+D_1/dy);
  n_2[i][0]= 0;
             if (n_1[i][0]<0) n_1[i][0]=1e-10;
             dphidy=(phi_[i][N_Y-1]-phi_[i][N_Y-2])/dy;//1.0*ddy(phi_,i,N_Y-1);

            // flow_1_y[i][NY-1]=mu_1*n_1[i][NY_1]*(dphidy)-D_1*(n_1[i][N_Y-1]-n_1[i][N_Y-2])/dy;

            n_1[i][N_Y-1]=0;//n_1[i][N_Y-1]*0.99+0.01* (D_1*(n_1[i][N_Y-2])/dy)/(-mu_1*dphidy+D_1/dy);
 n_2[i][N_Y-1]=0;
 if (n_1[i][N_Y-1]<0) n_1[i][N_Y-1]=1e-10;
if (n_1[i][0]<0) n_1[i][0]=1e-10;
          //   n_1[i][0]=1.0e12;
            // n_1[i][N_Y-1]=1.0e12;
        }

        for (int j=0; j< N_Y; j++)
        {
   //       n_1[0][j]=n_1[1][j];
   //       n_1[N_Y-1][j]=n_1[N_Y-2][j];

           double dphidx=0.0*ddx(phi_,0,j);
           // flow_1_y[i][0]=mu_1*n_1[i][0]*(dphidy)-D_1*(n_1[i][1]-n_1[i][0])/dx;

           n_1[0][j]=  (D_1*(n_1[1][j])/dx)/(mu_1*dphidx+D_1/dx);

            dphidx=0.0*ddx(phi_,N_X-1,j);
           // flow_1_y[i][NY-1]=mu_1*n_1[i][NY_1]*(dphidy)-D_1*(n_1[i][N_Y-1]-n_1[i][N_Y-2])/dx;

           n_1[N_X-1][j]=  (D_1*(n_1[N_X-2][j])/dx)/(-mu_1*dphidx+D_1/dx);


                  n_1[0][j]=n_1[N_X-2][j];
                  n_1[N_X-1][j]=n_1[1][j];


                  n_2[0][j]=n_2[N_X-2][j];
                  n_2[N_X-1][j]=n_2[1][j];
        }

        /*   for (int i=1; i<N_X-1; i++)
        {

            for (int j=1; j<N_Y-1; j++)
            {
               rho_in[i][j]=rho_in_prev[i][j];

            }
        }
*/

        for (int j=0; j<N_Y; j++)
        {
          //  rho_in[0][j]=rho_in[N_X-2][j];
          //  rho_in[N_X-1][j]=rho_in[1][j];
        }


    }
    double n_1_max=-10e10;
    double n_2_max=-10e10;
   // double n_1_new[N_X][N_Y];
   // double n_2_new[N_X][N_Y];

    for (int i=1; i< N_X-1; i++)
    {
        for(int j=1; j<N_Y-1;j++)
        {
            //n_1[i][j]=n_1_new[i][j];
            //n_2[i][j]=n_2_new[i][j];
            if (n_1[i][j]>n_1_max)
                n_1_max=n_1[i][j];
            if (n_2[i][j]>n_2_max)
                n_2_max=n_2[i][j];

        }
    }
    printf("n_1_max= %e\n", n_1_max, "\n");
     printf("n_2_max= %e\n", n_2_max, "\n");

     /*for (int i=1; i<N_X-1;i++ )
     {
         for (int j=0; j<  N_Y-1;N_Y++)
         {
             Jx_in[i][j]=rho_in[i][j]*(mu_2*ddx(phi_,i,j))-D_2*ddx(rho_in,i,j);
             Jy_in[i][j]=rho_in[i][j]*(mu_2*ddy(phi_,i,j))-D_2*ddy(rho_in,i,j);

             div_flow_1[i][j]=-1.0(n_1_new[i][j]-n_1[i][j])/dt_local+alpha*mod_flow_1[i][j]-beta*n_1[i][j]*n_2[i][j];
             div_flow_2[i][j]=-1.0(n_2_new[i][j]-n_2[i][j])/dt_local+alpha*mod_flow_1[i][j]-beta*n_1[i][j]*n_2[i][j];
         }
     }*/

   // double rho_max=-10e23;
 /*   for (int i=1; i<N_X-1; i++)
    {

        for (int j=1; j<N_Y-1; j++)
        {
            //    printf("%d %d \n",i,j);

            //rho_in[i][j]=rho_new[i][j];


            Jx_in[i][j]=rho_in[i][j]*(mu_2*ddx(phi_,i,j))-D_2*ddx(phi_,i,j);
            Jy_in[i][j]=rho_in[i][j]*(mu_2*ddy(phi_,i,j))-D_2*ddy(phi_,i,j);



            double dphidx=ddx(phi_,i,j);
            double dphidy=ddy(phi_,i,j);

            double d2phidx=(phi_[i+1][j]+phi_[i-1][j]-2.0*phi_[i][j])/(dx*dx);
            double d2phidy=(phi_[i][j+1]+phi_[i][j-1]-2.0*phi_[i][j])/(dy*dy);



            a=(1/dt_local+(2.0*D_2)/(dx*dx)+(2.0*D_2)/(dy*dy)+ mu_2*(d2phidx+d2phidy));
            b_p=(-1.0*D_2)/(dx*dx)- mu_2*dphidx/(2.0*dx);
            b_m=(-1.0*D_2)/(dx*dx)+ mu_2*dphidx/(2.0*dx);
            c_p=(-1.0*D_2)/(dy*dy)- mu_2*dphidy/(2.0*dy);
            c_m=(-1.0*D_2)/(dy*dy)+ mu_2*dphidy/(2.0*dy);

            div_J[i][j]=a*rho_in[i][j]+b_p*rho_in[i+1][j]+b_m*rho_in[i-1][j]+c_p*rho_in[i][j+1]+c_m*rho_in[i][j-1];

            if (rho_in[i][j]>rho_max)
                rho_max=rho_in[i][j];

        }
    }
*/
    // get_div(Jx_in,Jy_in,div_J);

  //  printf("rho_max=%f \n",rho_max);

    return res;
}

double solve_current_phi(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    static double rho_new[N_X][N_Y];

    double dt_local=0.0000015;
    for (int n=0; n< itn; n++)
    {
        for (int i=0; i<N_X; i++)
        {

            for (int j=0; j<N_Y; j++)
            {
                Jx_in[i][j]=-ddx(phi_,i,j);//(out_Ux[i][j]-b*ddx(phi_,i,j))-D*ddx(rho_in,i,j);
                Jy_in[i][j]=-ddy(phi_,i,j);
                //    printf("%d %d \n",i,j);

            }
        }


        get_div(Jx_in,Jy_in,RHS);


        for (int i=1; i<N_X-1; i++)
        {
             rho_in[i][N_Y-1]=(-phi_[i][N_Y-1]*dy+D*rho_in[i][N_Y-2])/(b*dy*ddy(phi_,i,N_Y-1)+D);

               rho_in[i][0]=(-phi_[i][0]*dy+D*rho_in[i][1])/(b*dy*ddy(phi_,i,0)+D);

           // rho_in[i][N_Y-1]=1.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);

           // rho_in[i][0]=A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);


            //rho_in[i][0]=-0.001;
            // rho_in[i][N_Y-1]=rho_in[i][1];
            for (int j=1; j<N_Y-1; j++)
            {
                //    printf("%d %d \n",i,j);

               //rho_in[i][j]=-RHS[i][j]*dt_local+rho_in[i][j];
                rho_in[i][j]=-RHS[i][j]/(4*M_PI);
                // rho_in[i][j]=rho_new[i][j];
            }
        }
        for (int j=0; j<N_Y; j++)
        {
            rho_in[0][j]=rho_in[N_X-2][j];
            rho_in[N_X-1][j]=rho_in[1][j];
        }



        /*  for (int i=1; i<N_X-1; i++)
        {

            for (int j=1; j<N_Y-1; j++)
            {
                //    printf("%d %d \n",i,j);
                if (rho_in[i][j]<0)
                    rho_in[i][j]=0.0;
                //rho_in[i][j]=rho_new[i][j];
            }
        }*/


    }




    double rho_max=-1000000000;
    for (int i=1; i<N_X-1; i++)
    {

        for (int j=1; j<N_Y-1; j++)
        {
            //    printf("%d %d \n",i,j);
            if (rho_in[i][j]>rho_max)
                rho_max=rho_in[i][j];
            //rho_in[i][j]=rho_new[i][j];
        }
    }

    get_div(Jx_in,Jy_in,div_J);
    printf("rho_max=%e \n",rho_max);

    return res;
}


double solve_current_old(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;

    static double rho_new[N_X][N_Y];

    double dt_local=0.0000015;
    for (int n=0; n< itn; n++)
    {
        for (int i=0; i<N_X; i++)
        {

            for (int j=0; j<N_Y; j++)
            {
                Jx_in[i][j]=rho_in[i][j]*(/*0.001*out_Ux[i][j]*/b*ddx(phi_,i,j))-D*ddx(rho_in,i,j);//(out_Ux[i][j]-b*ddx(phi_,i,j))-D*ddx(rho_in,i,j);
                Jy_in[i][j]=rho_in[i][j]*(/*0.001*out_Uy[i][j]*/b*ddy(phi_,i,j))-D*ddy(rho_in,i,j);
                //    printf("%d %d \n",i,j);

            }
        }

        get_div(Jx_in,Jy_in,RHS);


        for (int i=1; i<N_X-1; i++)
        {
            // rho_in[i][N_Y-1]=(-phi_[i][N_Y-1]*dy+D*rho_in[i][N_Y-2])/(b*dy*ddy(phi_,i,N_Y-1)+D);

            //   rho_in[i][0]=(-phi_[i][0]*dy+D*rho_in[i][1])/(b*dy*ddy(phi_,i,0)+D);

            rho_in[i][N_Y-1]=0.0;//rho_in[i][N_Y-2];// 0.0;//rho_in[i][0]*0.9+0.1*(-phi_[i][0]*dy+0.01*D*rho_in_prev[i][1])/(b*dy*ddy(phi_,i,0)+D);
            double sigma=1.0;
            double dPhi=(phi_[i][1]-phi_[i][0]);
            rho_in[i][0]=(sigma*dPhi + D*rho_in[i][1])/(b*dPhi+D);

            if (rho_in[i][0]<0)
                rho_in[i][0]=0.0;
            //rho_in[i][0]=-0.001;
            // rho_in[i][N_Y-1]=rho_in[i][1];
            for (int j=1; j<N_Y-1; j++)
            {
                //    printf("%d %d \n",i,j);

                rho_new[i][j]=-RHS[i][j]*dt_local+rho_in[i][j];
                rho_in[i][j]=rho_new[i][j];
            }
        }
        for (int j=0; j<N_Y; j++)
        {
            rho_in[0][j]=rho_in[N_X-2][j];
            rho_in[N_X-1][j]=rho_in[1][j];
        }

        /*  for (int i=1; i<N_X-1; i++)
        {

            for (int j=1; j<N_Y-1; j++)
            {
                //    printf("%d %d \n",i,j);
                if (rho_in[i][j]<0)
                    rho_in[i][j]=0.0;
                //rho_in[i][j]=rho_new[i][j];
            }
        }*/


    }


    double rho_max=-1000000000;
    for (int i=1; i<N_X-1; i++)
    {

        for (int j=1; j<N_Y-1; j++)
        {
            //    printf("%d %d \n",i,j);
            if (0.5*rho_in[i][j]>rho_max)
                rho_max=0.5*rho_in[i][j];
            //rho_in[i][j]=rho_new[i][j];
        }
    }

    get_div(Jx_in,Jy_in,div_J);
    printf("rho_max=%f \n",rho_max);

    return res;
}

double solve_current_new(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn)
{
    int i,j,n;
    double res=0;
    double a,b_p,b_m,c_p,c_m;
    double dt_local=0.0000015;

    a=(1/dt_local+(2.0*D)/(dx*dx)+(2.0*D)/(dy*dy));
    b_p=(-1.0*D)/(dx*dx);
    b_m=(-1.0*D)/(dx*dx);
    c_p=(-1.0*D)/(dy*dy);
    c_m=(-1.0*D)/(dy*dy);
    static double rho_new[N_X][N_Y];
    static double rho_in_prev[N_X][N_Y];

    for (int n=0; n< itn; n++)
    {
        for (int i=0; i<N_X; i++)
        {

            for (int j=0; j<N_Y; j++)
            {
                Jx_in[i][j]=rho_in[i][j]*(/*0.001*out_Ux[i][j]*/-b*ddx(phi_,i,j)-D*ddx(rho_in,i,j));//(out_Ux[i][j]-b*ddx(phi_,i,j))-D*ddx(rho_in,i,j);
                Jy_in[i][j]=rho_in[i][j]*(/*0.001*out_Uy[i][j]*/-b*ddy(phi_,i,j)-D*ddy(rho_in,i,j));

                //    printf("%d %d \n",i,j);
                rho_in_prev[i][j]=rho_in[i][j];

            }
        }
        get_div(Jx_in,Jy_in,RHS);


        for (int i=1; i<N_X-1; i++)
        {
            rho_in[i][N_Y-1]=(-phi_[i][N_Y-1]*dy+D*rho_in[i][N_Y-2])/(b*dy*ddy(phi_,i,N_Y-1)+D);// Boundary conditions: Jy_[N_Y-1]=phi_[N_Y-1]
            rho_in[i][0]=(-phi_[i][0]*dy+D*rho_in[i][1])/(b*dy*ddy(phi_,i,0)+D);
            //rho_in[i][0]=-0.001;
            // rho_in[i][N_Y-1]=rho_in[i][1];
            for (int j=1; j<N_Y-1; j++)
            {
                //    printf("%d %d \n",i,j);
                rho_in[i][j]=(RHS[i][j]+rho_in_prev[i][j]/dt_local-(b_p*rho_in[i+1][j]+b_m*rho_in[i-1][j]+c_p*rho_in[i][j+1]+c_m*rho_in[i][j-1]))/a;
                //RHS[i][j]= ((rho_in[i-1][j]+rho_in[i+1][j])*b_p+(rho_in[i][j-1]+rho_in[i][j+1])*c_p)/a;
                //rho_in[i][j]=-RHS[i][j]*0.01*0.000015+rho_in_prev[i][j];

            }
        }
        for (int j=0; j<N_Y; j++)
        {
            rho_in[0][j]=rho_in[N_X-2][j];
            rho_in[N_X-1][j]=rho_in[1][j];
        }

    }

    double rho_max=-1000000000;
    for (int i=1; i<N_X-1; i++)
    {

        for (int j=1; j<N_Y-1; j++)
        {
            //    printf("%d %d \n",i,j);
            if (rho_in[i][j]>rho_max)
                rho_max=rho_in[i][j];
            //rho_in[i][j]=rho_new[i][j];
        }
    }

    printf("rho_max=%f \n",rho_max);

    return res;
}


double NavierStoks_solver(double p_in[N_X][N_Y], double out_x[N_X][N_Y],double out_y[N_X][N_Y],int itn)
{
    static double field_Ux[N_X][N_Y];
    static double field_Uy[N_X][N_Y];
    static double p_new[N_X][N_Y];

    static double arr_X[N_X][N_Y];
    static double arr_Y[N_X][N_Y];

    double a,a_,b,c;
    //double res=0;
    double u=0;
    double v=0;
    a=(1/dt+(2.0)/(dx*dx)+2.0/(dy*dy));
    a_=((2.0)/(dx*dx)+2.0/(dy*dy));
    b=1.0/(dx*dx);
    c=1.0/(dy*dy);


    for(int n=0;n<itn;n++)
    {
        for (int i=1; i<N_X-1; i++)
        {

            for (int j=1; j<N_Y-1; j++)
            {

                field_Ux[i][j]=out_x[i][j];
                field_Uy[i][j]=out_y[i][j];
                field_Ux[i][j]=(-out_x[i][j]*ddx(out_x,i,j)-out_x[i][j]*ddy(out_x,i,j)-ddx(p_in,i,j)+(out_x[i+1][j]+out_x[i-1][j])/b+(out_x[i][j+1]+out_x[i][j-1])/c)/a;
                field_Uy[i][j]=(-out_x[i][j]*ddx(out_y,i,j)-out_y[i][j]*ddy(out_y,i,j)-ddy(p_in,i,j)+(out_y[i+1][j]+out_y[i-1][j])/b+(out_y[i][j+1]+out_y[i][j-1])/c)/a;
                field_Ux[i][j]=u;
                field_Uy[i][j]=v;
            }
            // printf("Ux=%f Uy=%f\n",u ,v);
        }
        //     get_div(out_Ux, out_Uy, RHS);
    }
    /*   for (int n=0; n< itn; n++)
    {
        for (int i=0;i< N_X; i++ )
        {
            for (int j=0; j< N_Y; j++)
            {
                arr_X[i][j]=field_Ux[i][j]*ddx(field_Ux,i,j)+field_Uy[i][j]*ddy(field_Ux,i,j);
                arr_Y[i][j]=field_Ux[i][j]*ddx(field_Uy,i,j)+field_Uy[i][j]*ddy(field_Uy,i,j);
                p_in[i][j]=(ddx(arr_X,i,j)+ddy(arr_Y,i,j)+(p_in[i+1][j]+p_in[i-1][j])/(-b)+(p_in[i][j+1]+p_in[i][j-1])/(-c))/a_;
                p_new[i][j]=p_in[i][j];
               // res=p_new[i][j];
            }
            //printf("p=%f \n",res);
        }
    }
*/

    /* double coef= ;
    for (int n=0; n< itn; n++)
    {

        for(int i=0; i< N_X; i++)
        {
            for (int j=0; j<N_Y; j++)
            {
                out_Ux[i][j]=(U[i][j+1]+U[i][j-1]-RHS)*dy*dy/2.0;
            }
        }
        for (int i=1; i<N_X-1; i++)
        {
            p_in[i][j]=rho_in[i][j]*(ddx*out_Ux[i][j]+ddy*out_Ux[i][j])+ddx*(p_in[i+1][j]+p_in[i-1][j])+ddy*(p_in[i][j+1]+p_in[i][j-1]/coef) ;
        }
    }
    */
}

