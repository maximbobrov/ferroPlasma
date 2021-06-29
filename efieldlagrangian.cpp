#include "efieldlagrangian.h"
#include "globals.h"
#include <assert.h>

void eFieldLagrangian::init()
{
    m_elec_num=0;

    m_chargeNum = 0 ;

    /*vec2 p[5];

    p[0].x=w_x0+50e-6-0e-9;         p[0].y=0.5*(w_y0+w_y1);
    p[1].x=w_x0+50e-6-20e-6;        p[1].y=w_y1;
    p[2].x=w_x0+50e-6-25e-6;        p[2].y=w_y1;
    p[3].x=w_x0+50e-6-25e-6;        p[3].y=0.5*(w_y0+w_y1);
    p[4].x=p[0].x;            p[4].y=p[0].y;

    int emit_1[4] = {1,0,0,0};
    addQuad(p,2e-6,-1250,emit_1);
    printf("elecnum1 = %d\n", m_elec_num);

    p[3].x=w_x0-25e-6;        p[3].y=w_y0;
    p[2].x=w_x1;              p[2].y=w_y0;
    p[1].x=w_x1;              p[1].y=w_y0-5e-6;
    p[0].x=w_x0-25e-6;        p[0].y=w_y0-5e-6;
    p[4].x=p[0].x;            p[4].y=p[0].y;

    int emit_2[4] = {0,0,0,0};
    addQuad(p,2e-6,1250,emit_2);*/

    vec2 p[5];

    p[0].x=w_x0-0e-6;         p[0].y=0;
    p[1].x=w_x0 /*-40e-6*/;        p[1].y=w_y1;
    p[2].x=w_x0-30e-6;        p[2].y=w_y1;
    p[3].x=w_x0-30e-6;        p[3].y=0;
    p[4].x=p[0].x;            p[4].y=p[0].y;

    int emit_1[4] = {1,0,0,0};
    //double dl[5] = {0.5 * 0.125e-6, 2e-6, 4e-6, 2e-6, 0.5 * 0.25e-6};
    double dl[5] = {1. * 0.15e-6, 2e-6, 4e-6, 2e-6, 1.0 * 0.5e-6};

    //   addQuad_stabilized(p,dl, -g_phi ,emit_1, w_y0+25e-6 + 0.5 * dl_pz, 30,0);

    addQuad_noSmooth(p,dl, -g_phi ,emit_1, 0);

    printf("elecnum1 = %d\n", m_elec_num);

    p[4].x=w_x0-30e-6;        p[4].y=-dl_pz-1e-6;//w_y0;
    p[3].x=w_x1;              p[3].y=-dl_pz-1e-6;
    p[2].x=w_x1;              p[2].y=-dl_pz-10e-6;
    p[1].x=w_x0-30e-6;        p[1].y=-dl_pz-10e-6;
    p[0].x=p[4].x;            p[0].y=p[4].y;

    int emit_2[4] = {0,0,0,0};

    double coef = 1.0;
    double dl2 = coef * 2e-6;
    vec2 p2[2];

    p2[1].x=w_x1-10e-6;               p2[1].y=-dl_pz;//-1e-6;
    p2[0].x=w_x0-20e-6;         p2[0].y=-dl_pz;//-1e-6;

    addLine( p2, 3e-6,g_phi, 1e-6, 1);
    /*    double dl2[5] = {coef * 2e-6, coef * 2e-6,coef * 2e-6,coef * 2e-6,coef * 2e-6};
     * addQuad(p,dl2, g_phi,emit_2,  w_y0+25e-6 - 0.5 * dl_pz, 1,1);*/
    //addQuad_stabilized(p,dl2,20 * 12.50,emit_2,  w_y0+25e-6 - 0.5 * dl_pz, 1);
    //addQuad2Layers(p,dl,-20 * 12.50,emit_1);

    printf("elecnum2 = %d\n", m_elec_num);

    //initW_PhiE();
    initW();
    //getInv_PhiE();
    getInv();

}

double LU[3000][3000],M_[3000][3000],Inv[3000][3000];
double x_m[3000],b_m[3000],mwb[3000];
int ps[3000];

double m_W[1000][1000]; // weight of charge i at point j [i][j]

void LU_decompose(int num)
{
    int i,j,k,pivotindex;
    static double scales[3000];
    double normrow,pivot,size,biggest,mult;

    for (i=0;i<num;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<num;j++)
        {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<num-1;k++)
    {
        biggest=0;
        for (i=k; i<num;i++)
        {
            size=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size)
            {
                biggest=size;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<num;i++)
        {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<num;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
}

void m_solve(int num)
{
    int i,j;
    double dot;

    for (i=0;i<num;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=num-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<num;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}


void m_invert(int num)
{
    int i,j;
    //err_code(1);
    LU_decompose(num);

    for (j=0;j<num;j++)
    {

        for (i=0;i<num;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve(num);

        for (i=0;i<num;i++)
            Inv[i][j]=x_m[i];
    }
}

eFieldLagrangian::eFieldLagrangian()
{
    m_electrodes=new eElem[2000];
    m_rCentre=new vec2[2000];

    init();
printf("elec_num   ==  %d \n",m_elec_num);

    assert(m_elec_num<=ELEC_MAX);


}



void eFieldLagrangian::updatePhi()
{
    for (int i=0;i<m_elec_num;i++)
    {
        double x,y;
        x=m_electrodes[i].r.x;
        y=m_electrodes[i].r.y;
        m_electrodes[i].phi_ext=getPhi(x,y);
    }
}

void eFieldLagrangian::solvePhi(int itn)
{
    for (int j=0;j<itn;j++)
    {
        updatePhi();
        /* for (int i=0;i<m_elec_num;i++)
        {
            m_electrodes[i].rho1+=4e-3*(m_electrodes[i].phi_fix+m_electrodes[i].phi_fix_charges-m_electrodes[i].phi_ext);
            m_electrodes[i].rho2=m_electrodes[i].rho1;
        }*/
    }
    //  updateGridProp();
}


void  getProgrCoef(double b1, double bn, double l, int& n, double &q)
{
    bool ok = true;
    n = 4;
    double sum = 0;
    while( sum < l && n<1000)
    {
        q = pow(bn/b1, 1.0/(n-1));
        if(fabs(q - 1) < 1e-5)
            sum = b1 * n;
        else
            sum = b1 * (1-pow(q,n)) / (1-q);
        n++;
    }
    n-=1;
    q = pow(bn/b1, 1.0/(n-1));
}

void eFieldLagrangian::addQuad_stabilized(vec2 p[5], double dl[5],double phi, int emit[4], double coordYDIel, int smoothingCount,int I) //the field direction is stabilized here
{
    static double l_[4];
    static int n_[4];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    double nx, ny;
    for (int i=0;i<4;i++)
    {
        nx = p[i+1].y - p[i].y;
        ny = -p[i+1].x + p[i].x;
        l_[i]=sqrt((p[i+1].x-p[i].x)*(p[i+1].x-p[i].x)+(p[i+1].y-p[i].y)*(p[i+1].y-p[i].y));
        double q;
        double s=0;
        getProgrCoef(dl[i], dl[i+1], l_[i], n_[i], q);
        double l = dl[i];
        for (int j=0;j<n_[i];j++)
        {

            double alpha=s*1.0/(l_[i]);
            double x,y;
            x = p[i].x*(1.0-alpha)+p[i+1].x*(alpha);
            y = p[i].y*(1.0-alpha)+p[i+1].y*(alpha);

            m_electrodes[m_elec_num].INDX=I;
            m_electrodes[m_elec_num].r.x=x;
            m_electrodes[m_elec_num].r.y=y;
            m_electrodes[m_elec_num].dl=l;

            m_electrodes[m_elec_num].phi_fix=phi*2.0;
            m_electrodes[m_elec_num].phi_ext=0.0;

            if (emit[i]>0)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;
            m_electrodes[m_elec_num].nx = nx / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].ny = ny / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].eToEmit=0.0;

            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;


            l = dl[i] * pow(q,j);
            s+=l;
        }

        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    static double xCoord[1000];
    static double yCoord[1000];
    for (int f =0 ; f<smoothingCount;f++) {

        for (int i=n0;i<m_elec_num;i++)
        {
            xCoord[i] = m_electrodes[i].r.x;
            yCoord[i] = m_electrodes[i].r.y;
        }

        for (int i=n0;i<m_elec_num;i++)
        {
            int im = i-1;
            if(i == n0)
                im = m_elec_num - 1;
            int ip = i+1;
            if(i == m_elec_num - 1)
                ip = n0;

            m_electrodes[i].r.x = (xCoord[im] + xCoord[i] + xCoord[ip])/3;
            m_electrodes[i].r.y = (yCoord[im] + yCoord[i] + yCoord[ip])/3;
        }
    }

    for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;
        nx = m_electrodes[ip].r.y - m_electrodes[im].r.y;
        ny = -m_electrodes[ip].r.x + m_electrodes[im].r.x;
        m_electrodes[i].nx = nx / sqrt(nx*nx + ny*ny);
        m_electrodes[i].ny = ny / sqrt(nx*nx + ny*ny);
    }

    int el_num0=m_elec_num;
    for (int i=n0;i<el_num0/*m_elec_num*/;i++)
    {
        xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
        //additional monitoring points with 0 potentioal for stability.
        // this leads to dimminishing of the practical surface potential by 2 times but the field direction is now stable

        m_electrodes[m_elec_num].INDX=I;
        m_electrodes[m_elec_num].r.x=xCoord[i];
        m_electrodes[m_elec_num].r.y=yCoord[i];
        m_electrodes[m_elec_num].phi_fix=phi*0.0;
        m_electrodes[m_elec_num].phi_ext=0.0;
        m_electrodes[m_elec_num].canEmit=false;
        m_electrodes[m_elec_num].nx = m_electrodes[i].nx ;
        m_electrodes[m_elec_num].ny = m_electrodes[i].ny;
        m_electrodes[m_elec_num].eToEmit=0.0;
        m_elec_num++;
    }

    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;

    int n1=m_elec_num;
}





void eFieldLagrangian::addQuad(vec2 p[5], double dl[5],double phi, int emit[4], double coordYDIel, int smoothingCount, int I) //last point should coincide with the first one emit is the side number that can emit
{
    static double l_[4];
    static int n_[4];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    double nx, ny;
    for (int i=0;i<4;i++)
    {
        nx = p[i+1].y - p[i].y;
        ny = -p[i+1].x + p[i].x;
        l_[i]=sqrt((p[i+1].x-p[i].x)*(p[i+1].x-p[i].x)+(p[i+1].y-p[i].y)*(p[i+1].y-p[i].y));
        double q;
        double s=0;
        getProgrCoef(dl[i], dl[i+1], l_[i], n_[i], q);
        //printf("ni=%d\n",n_[i]);//n_[i]-=2;
        if(n_[i]<=4)
            n_[i]-=2;
        double l = dl[i];
        for (int j=0;j<n_[i];j++)
        {

            double alpha=s*1.0/(l_[i]);
            double x,y;
            x = p[i].x*(1.0-alpha)+p[i+1].x*(alpha);
            y = p[i].y*(1.0-alpha)+p[i+1].y*(alpha);

            m_electrodes[m_elec_num].INDX=I;
            m_electrodes[m_elec_num].r.x=x;
            m_electrodes[m_elec_num].r.y=y;
            m_electrodes[m_elec_num].dl=l;

            m_electrodes[m_elec_num].phi_fix=phi;
            m_electrodes[m_elec_num].phi_ext=0.0;

            if (emit[i]>0)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;
            m_electrodes[m_elec_num].nx = nx / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].ny = ny / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].eToEmit=0.0;

            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;

            l = dl[i] * pow(q,j);
            s+=l;
        }

        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    static double xCoord[1000];
    static double yCoord[1000];
    for (int f =0 ; f<smoothingCount;f++) {

        for (int i=n0;i<m_elec_num;i++)
        {
            xCoord[i] = m_electrodes[i].r.x;
            yCoord[i] = m_electrodes[i].r.y;
        }

        for (int i=n0;i<m_elec_num;i++)
        {
            int im = i-1;
            if(i == n0)
                im = m_elec_num - 1;
            int ip = i+1;
            if(i == m_elec_num - 1)
                ip = n0;

            m_electrodes[i].r.x = (xCoord[im] + xCoord[i] + xCoord[ip])/3;
            m_electrodes[i].r.y = (yCoord[im] + yCoord[i] + yCoord[ip])/3;
        }
    }

    for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;
        nx = m_electrodes[ip].r.y - m_electrodes[im].r.y;
        ny = -m_electrodes[ip].r.x + m_electrodes[im].r.x;
        m_electrodes[i].nx = nx / sqrt(nx*nx + ny*ny);
        m_electrodes[i].ny = ny / sqrt(nx*nx + ny*ny);
    }

    int el_num0=m_elec_num;
    for (int i=n0;i<el_num0/*m_elec_num*/;i++)
    {
        xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
        /* xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * m_electrodes[i].dl*1.0;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * m_electrodes[i].dl*1.0;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;*/
        //additional monitoring points with 0 potentioal for stability.
        // this leads to dimminishing of the practical surface potential by 2 times but the field direction is now stable


    }

    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;

    int n1=m_elec_num;


    /*int i0 = m_chargeNum;
    for (int i=i0;i<i0+m_elec_num;i++)
    {
        double x,y;
        x = 0.955*(m_electrodes[i-i0].r.x -c_m.x)+c_m.x;
        y = 0.955*(m_electrodes[i-i0].r.y -c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/

    /*m_charges[m_chargeNum].x=c_m.x;
    m_charges[m_chargeNum].y=c_m.y;
    m_charges[m_chargeNum].charge=0.0;
    m_chargeNum++;*/
    /*for (int i=0;i<4;i++)
    {
        double x,y;
        x = 0.985*(p[i].x-c_m.x)+c_m.x;
        y = 0.985*(p[i].y-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/
}



void eFieldLagrangian::addQuad_noSmooth(vec2 p[5], double dl[5],double phi, int emit[4],int I) //last point should coincide with the first one emit is the side number that can emit
{
    static double l_[4];
    static int n_[4];
    vec2 p_inner[5],n_inner[5];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    double nx, ny;


    for (int i=0;i<4;i++)
    {

        int ip,im;
        ip=i+1;  if (ip>=4) ip=0;
        im=i-1;  if (im<0) im=3;

        double nx_m,ny_m,nx_p, ny_p,lp,lm,l;

        nx_p = p[ip].y - p[i].y;
        ny_p = -p[ip].x + p[i].x;
        lp = sqrt(nx_p*nx_p + ny_p*ny_p);
        nx_p /= lp;
        ny_p /= lp;

        nx_m = p[i].y - p[im].y;
        ny_m = -p[i].x + p[im].x;
        lm = sqrt(nx_m*nx_m + ny_m*ny_m);
        nx_m /= lm;
        ny_m /= lm;

        nx=nx_m+nx_p;
        ny=ny_m+ny_p;
        l = sqrt(nx*nx + ny*ny);
        nx /= l;
        ny /= l;

        n_inner[i].x=nx;
        n_inner[i].y=ny;

        nx*=dl[i];
        ny*=dl[i];

        p_inner[i].x=p[i].x-nx;
        p_inner[i].y=p[i].y-ny;
    }

    p_inner[4]=p_inner[0];
    n_inner[4]=n_inner[0];

    for (int i=0;i<4;i++)
    {
        nx = p[i+1].y - p[i].y;
        ny = -p[i+1].x + p[i].x;
        l_[i]=sqrt((p[i+1].x-p[i].x)*(p[i+1].x-p[i].x)+(p[i+1].y-p[i].y)*(p[i+1].y-p[i].y));
        double q;
        double s=0;
        getProgrCoef(dl[i], dl[i+1], l_[i], n_[i], q);
        //printf("ni=%d\n",n_[i]);//n_[i]-=2;
        if(n_[i]<=4)
            n_[i]-=2;
        double l = dl[i];
        for (int j=0;j<n_[i]-1;j++)
        {

            double alpha=s*1.0/(l_[i]);
            double x,y;
            x = p[i].x*(1.0-alpha)+p[i+1].x*(alpha);
            y = p[i].y*(1.0-alpha)+p[i+1].y*(alpha);

            m_electrodes[m_elec_num].INDX=I;
            m_electrodes[m_elec_num].r.x=x;
            m_electrodes[m_elec_num].r.y=y;
            m_electrodes[m_elec_num].dl=l;

            m_electrodes[m_elec_num].phi_fix=phi*2;
            m_electrodes[m_elec_num].phi_ext=0.0;

            if (emit[i]>0)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;


            m_electrodes[m_elec_num].nx = nx / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].ny = ny / sqrt(nx*nx + ny*ny);

            if (j==0) //corner points
            {
                m_electrodes[m_elec_num].nx = n_inner[i].x;
                m_electrodes[m_elec_num].ny =  n_inner[i].y;
            }

            m_electrodes[m_elec_num].eToEmit=0.0;

            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;

            l = dl[i] * pow(q,j);
            s+=l;

            /*  //now fopr charges

            x = p_inner[i].x*(1.0-alpha)+p_inner[i+1].x*(alpha);
            y = p_inner[i].y*(1.0-alpha)+p_inner[i+1].y*(alpha);

            m_charges[m_chargeNum].x = x;
            m_charges[m_chargeNum].y = y;

            m_chargeNum++;*/
        }

        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    static double xCoord[1000];
    static double yCoord[1000];

    for (int i=n0;i<m_elec_num;i++)
    {
        xCoord[i] = m_electrodes[i].r.x;
        yCoord[i] = m_electrodes[i].r.y;
    }


    /*int el_num0=m_elec_num;
    for (int i=n0;i<el_num0;i++)
    {
        xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
    }*/


    int el_num0=m_elec_num;
    for (int i=n0;i<el_num0/*m_elec_num*/;i++)
    {
        xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
        //additional monitoring points with 0 potentioal for stability.
        // this leads to dimminishing of the practical surface potential by 2 times but the field direction is now stable

        m_electrodes[m_elec_num].INDX=I;
        m_electrodes[m_elec_num].r.x=xCoord[i];
        m_electrodes[m_elec_num].r.y=yCoord[i];
        m_electrodes[m_elec_num].phi_fix=phi*0.0;
        m_electrodes[m_elec_num].phi_ext=0.0;
        m_electrodes[m_elec_num].canEmit=false;
        m_electrodes[m_elec_num].nx = m_electrodes[i].nx ;
        m_electrodes[m_elec_num].ny = m_electrodes[i].ny;
        m_electrodes[m_elec_num].eToEmit=0.0;
        m_elec_num++;
    }


    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;
}

void eFieldLagrangian::addLine(vec2 p[2], double dl,double phi, double h, int I) //last point should coincide with the first one emit is the side number that can emit
{
    int n0=m_elec_num;
    double n = (p[1].x - p[0].x) / dl;
    for (int i=0;i< n;i++)
    {
        double alpha= i * 1.0/n;
        double x,y;
        x = p[0].x*(1.0-alpha)+p[1].x*(alpha);
        y = p[0].y*(1.0-alpha)+p[1].y*(alpha);

        m_electrodes[m_elec_num].INDX=I;
        m_electrodes[m_elec_num].r.x=x;
        m_electrodes[m_elec_num].r.y=y;
        m_electrodes[m_elec_num].dl=dl;

        m_electrodes[m_elec_num].phi_fix=phi;
        m_electrodes[m_elec_num].phi_ext=0.0;


        m_electrodes[m_elec_num].canEmit=false;
        m_electrodes[m_elec_num].nx = 0;
        m_electrodes[m_elec_num].ny = 1;
        m_electrodes[m_elec_num].eToEmit=0.0;

        m_elec_num++;
    }
    int el_num0=m_elec_num;
    for (int i=n0;i<el_num0;i++)
    {
        m_charges[m_chargeNum].x = m_electrodes[i].r.x;
        m_charges[m_chargeNum].y = m_electrodes[i].r.y - h;
        m_chargeNum++;
    }
}

void eFieldLagrangian::addLinePZ(vec2 p[2], double dl,double phi, double h, int I) //last point should coincide with the first one emit is the side number that can emit
{
    int n0=m_elec_num;
    /////////////////line similar to pz
    double _dx;
    int m_p_num=100;
    _dx= (w_x1- w_x0)/(m_p_num-1);



    for (int i=0;i< m_p_num;i++)
    {
        double x = w_x0 + _dx * i -18e-6;;
        double y = p[0].y*(1.0-alpha)+p[1].y*(alpha);

        m_electrodes[m_elec_num].INDX=I;
        m_electrodes[m_elec_num].r.x=x;
        m_electrodes[m_elec_num].r.y=y;
        m_electrodes[m_elec_num].dl=dl;

        m_electrodes[m_elec_num].phi_fix=phi;
        m_electrodes[m_elec_num].phi_ext=0.0;


        m_electrodes[m_elec_num].canEmit=false;
        m_electrodes[m_elec_num].nx = 0;
        m_electrodes[m_elec_num].ny = 1;
        m_electrodes[m_elec_num].eToEmit=0.0;

        m_elec_num++;
    }
    //////////////////////////end lnie

    int el_num0=m_elec_num;
    for (int i=n0;i<el_num0;i++)
    {
        m_charges[m_chargeNum].x = m_electrodes[i].r.x;
        m_charges[m_chargeNum].y = m_electrodes[i].r.y - h;
        m_chargeNum++;
    }
}

void eFieldLagrangian::addQuad2Layers(vec2 p[5], double dl[5],double phi, int emit[4]) //last point should coincide with the first one emit is the side number that can emit
{
    static double l_[4];
    static int n_[4];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    double nx, ny;
    for (int i=0;i<4;i++)
    {
        nx = p[i+1].y - p[i].y;
        ny = -p[i+1].x + p[i].x;
        l_[i]=sqrt((p[i+1].x-p[i].x)*(p[i+1].x-p[i].x)+(p[i+1].y-p[i].y)*(p[i+1].y-p[i].y));
        double q;
        double s=0;
        getProgrCoef(dl[i], dl[i+1], l_[i], n_[i], q);
        double l = dl[i];
        for (int j=0;j<n_[i];j++)
        {

            double alpha=s*1.0/(l_[i]);
            double x,y;
            x = p[i].x*(1.0-alpha)+p[i+1].x*(alpha);
            y = p[i].y*(1.0-alpha)+p[i+1].y*(alpha);

            m_electrodes[m_elec_num].r.x=x;
            m_electrodes[m_elec_num].r.y=y;
            m_electrodes[m_elec_num].dl=l;

            m_electrodes[m_elec_num].phi_fix=phi;
            m_electrodes[m_elec_num].phi_ext=0.0;

            if (emit[i]>0)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;
            m_electrodes[m_elec_num].nx = nx / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].ny = ny / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].eToEmit=0.0;

            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;


            l = dl[i] * pow(q,j);
            s+=l;
        }

        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    static double xCoord[1000];
    static double yCoord[1000];
    for (int i=n0;i<m_elec_num;i++)
    {
        xCoord[i] = m_electrodes[i].r.x;
        yCoord[i] = m_electrodes[i].r.y;
    }

    for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;

        m_electrodes[i].r.x = (xCoord[im] + xCoord[i] + xCoord[ip])/3;
        m_electrodes[i].r.y = (yCoord[im] + yCoord[i] + yCoord[ip])/3;
    }

    for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;
        nx = m_electrodes[ip].r.y - m_electrodes[im].r.y;
        ny = -m_electrodes[ip].r.x + m_electrodes[im].r.x;
        m_electrodes[i].nx = nx / sqrt(nx*nx + ny*ny);
        m_electrodes[i].ny = ny / sqrt(nx*nx + ny*ny);
    }

    for (int i=n0;i<m_elec_num;i+=2)
    {
        xCoord[i] = (m_electrodes[i].r.x + m_electrodes[i+1].r.x) * 0.5 - m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = (m_electrodes[i].r.y + m_electrodes[i+1].r.y) * 0.5 - m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
        /*xCoord[i] = m_electrodes[i].r.x ;//+ m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y ;//+ m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;*/
    }

    for (int i=n0;i<m_elec_num;i+=2)
    {
        xCoord[i] = (m_electrodes[i].r.x + m_electrodes[i+1].r.x) * 0.5 - m_electrodes[i].nx * m_electrodes[i].dl*1.0;
        yCoord[i] = (m_electrodes[i].r.y + m_electrodes[i+1].r.y) * 0.5 - m_electrodes[i].ny * m_electrodes[i].dl*1.0;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
        /*xCoord[i] = m_electrodes[i].r.x ;//+ m_electrodes[i].nx * m_electrodes[i].dl*0.5;
        yCoord[i] = m_electrodes[i].r.y ;//+ m_electrodes[i].ny * m_electrodes[i].dl*0.5;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;*/
    }


    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;

    int n1=m_elec_num;


    /*int i0 = m_chargeNum;
    for (int i=i0;i<i0+m_elec_num;i++)
    {
        double x,y;
        x = 0.955*(m_electrodes[i-i0].r.x -c_m.x)+c_m.x;
        y = 0.955*(m_electrodes[i-i0].r.y -c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/

    /*m_charges[m_chargeNum].x=c_m.x;
    m_charges[m_chargeNum].y=c_m.y;
    m_charges[m_chargeNum].charge=0.0;
    m_chargeNum++;*/
    /*for (int i=0;i<4;i++)
    {
        double x,y;
        x = 0.985*(p[i].x-c_m.x)+c_m.x;
        y = 0.985*(p[i].y-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/
}

void eFieldLagrangian::addQuadRegular(vec2 p[5], double dl,double phi, int emit[4]) //last point should coincide with the first one emit is the side number that can emit
{
    static double l_[4];
    static int n_[4];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    double nx, ny;
    for (int i=0;i<4;i++)
    {
        nx = p[i+1].y - p[i].y;
        ny = -p[i+1].x + p[i].x;
        l_[i]=sqrt((p[i+1].x-p[i].x)*(p[i+1].x-p[i].x)+(p[i+1].y-p[i].y)*(p[i+1].y-p[i].y));
        n_[i]=(int)(l_[i]/dl);
        double dl_l=l_[i]/n_[i];
        for (int j=0;j<n_[i];j++)
        {
            double alpha=j*1.0/(n_[i]);
            double x,y;
            x = p[i].x*(1.0-alpha)+p[i+1].x*(alpha);
            y = p[i].y*(1.0-alpha)+p[i+1].y*(alpha);

            m_electrodes[m_elec_num].r.x=x;
            m_electrodes[m_elec_num].r.y=y;
            m_electrodes[m_elec_num].dl=dl_l;

            m_electrodes[m_elec_num].phi_fix=phi;
            m_electrodes[m_elec_num].phi_ext=0.0;

            if (emit[i]>0)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;
            m_electrodes[m_elec_num].nx = nx / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].ny = ny / sqrt(nx*nx + ny*ny);
            m_electrodes[m_elec_num].eToEmit=0.0;

            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;



        }

        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    static double xCoord[1000];
    static double yCoord[1000];
    for (int i=n0;i<m_elec_num;i++)
    {
        xCoord[i] = m_electrodes[i].r.x;
        yCoord[i] = m_electrodes[i].r.y;
    }

    /* for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;

        m_electrodes[i].r.x = (xCoord[im] + xCoord[i] + xCoord[ip])/3;
        m_electrodes[i].r.y = (yCoord[im] + yCoord[i] + yCoord[ip])/3;
    }
*/
    for (int i=n0;i<m_elec_num;i++)
    {
        int im = i-1;
        if(i == n0)
            im = m_elec_num - 1;
        int ip = i+1;
        if(i == m_elec_num - 1)
            ip = n0;
        nx = m_electrodes[ip].r.y - m_electrodes[im].r.y;
        ny = -m_electrodes[ip].r.x + m_electrodes[im].r.x;
        m_electrodes[i].nx = nx / sqrt(nx*nx + ny*ny);
        m_electrodes[i].ny = ny / sqrt(nx*nx + ny*ny);
    }

    for (int i=n0;i<m_elec_num;i++)
    {
        xCoord[i] = m_electrodes[i].r.x - m_electrodes[i].nx * 1e-6;
        yCoord[i] = m_electrodes[i].r.y - m_electrodes[i].ny * 1e-6;
        m_charges[m_chargeNum].x = xCoord[i];
        m_charges[m_chargeNum].y = yCoord[i];
        m_chargeNum++;
    }

    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;

    int n1=m_elec_num;

    /*for (int i=n0;i<n1-1;i+=2)
    {
        double x,y;
        x = 0.775*((m_electrodes[i].r.x + m_electrodes[i+1].r.x)*0.5-c_m.x)+c_m.x;
        y = 0.775*((m_electrodes[i].r.y + m_electrodes[i+1].r.y)*0.5-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/

    /*for (int i=n0;i<n1-3;i+=4)
    {
        double x,y;
        x = 0.955*((m_electrodes[i].r.x + m_electrodes[i+1].r.x + m_electrodes[i+2].r.x + m_electrodes[i+3].r.x)*0.25-c_m.x)+c_m.x;
        y = 0.955*((m_electrodes[i].r.y + m_electrodes[i+1].r.y + m_electrodes[i+2].r.y + m_electrodes[i+3].r.y)*0.25-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/

    /*m_charges[m_chargeNum].x=c_m.x;
    m_charges[m_chargeNum].y=c_m.y;
    m_charges[m_chargeNum].charge=0.0;
    m_chargeNum++;*/
    /*for (int i=0;i<4;i++)
    {
        double x,y;
        x = 0.985*(p[i].x-c_m.x)+c_m.x;
        y = 0.985*(p[i].y-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }*/
}

double eFieldLagrangian::getW(double s_x, double s_y,double t_x, double t_y) //get Weight function ;//source (charge) and target (monitoring point)
{

    double sum=0.0;
    //    int i=1;
    double r;
    double q;
    double dx,dy;
    //double delta=1e-6;

    dx = s_x - t_x;
    dy = s_y - t_y;
    r=sqrt(dx*dx+dy*dy);
    //       q=qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);

    //sum-=-q*log(r+delta)/(w_z1 - w_z0);
    sum=-(qe/(eps0*pi2))*log(r+delta_phi)/(w_z1 - w_z0);
    return sum;
}



double eFieldLagrangian::getW_E(int elecNum, int chargeNum)
{
    vec2 sum(0.0, 0, 0);
    //    int i=1;
    double r;
    double q;
    double dx,dy;
    // double delta=1e-6;

    double r2;

    dx = m_charges[chargeNum].x - m_electrodes[elecNum].r.x;
    dy = m_charges[chargeNum].y - m_electrodes[elecNum].r.y;
    r2=(dx*dx+dy*dy);
    q=-qe/(eps0*pi2);

    double c=q/((r2+delta_phi*delta_phi)*(w_z1 - w_z0));

    sum.x=c*dx * m_electrodes[elecNum].ny;
    sum.y=c*dy * (-m_electrodes[elecNum].nx);
    return 0.00001 * (sum.x + sum.y);
}

double eFieldLagrangian::getPhi(double x, double y)
{
    double sum=0.0;
    for (int i=0;i<m_chargeNum;i++)
    {
        //    int i=1;
        double r;
        double q;
        double dx,dy;
        //double delta=1e-6;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (m_charges[i].charge);

        sum+=-q*log(r+delta_phi)/(w_z1 - w_z0);
    }

    return sum;
}

double eFieldLagrangian::getPhi_multiplier(double x, double y, int i)
{
    double sum=0.0;

        //    int i=1;
        double r;
        double q;
        double dx,dy;
        //double delta=1e-6;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) /* (m_charges[i].charge)*/;

        sum+=-q*log(r+delta_phi)/(w_z1 - w_z0);


    return sum;
}

void eFieldLagrangian::initW()
{
    for(int i=0;i<m_chargeNum;i++)
    {
        for(int j=0;j<m_elec_num;j++)
        {
#ifdef USE_MIRROR
            m_W[i][j]=getWMirror(i,j);//getW(m_charges[i].x,m_charges[i].y,m_electrodes[j].r.x,m_electrodes[j].r.y);
#else
            m_W[i][j]=getW(m_charges[i].x,m_charges[i].y,m_electrodes[j].r.x,m_electrodes[j].r.y);
#endif
        }
    }
}

void eFieldLagrangian::initW_PhiE()
{
    for(int i=0;i<m_chargeNum;i++)
    {
        for(int j=0;j<m_elec_num;j++)
        {
            m_W[i][j]=getW(m_charges[i].x,m_charges[i].y,m_electrodes[j].r.x,m_electrodes[j].r.y);
        }
        for(int j=m_elec_num;j<2 * m_elec_num;j++)
        {
            m_W[i][j]=getW_E(j - m_elec_num, i);
        }
    }

}



void eFieldLagrangian::solve_ls()
{
    int var_num=m_chargeNum;
    int eq_num=m_elec_num;

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            M_[i][j]=0.0;
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=m_W[i][n]*m_W[j][n];  //its mvm
            }
        }
    }

    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_electrodes[i].phi_fix-m_electrodes[i].phi_fix_charges;
    }

    for (int i=0;i<var_num;i++)
    {
        mwb[i]=0.0;
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=m_W[i][n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb
    }
    LU_decompose(var_num);
    m_solve(var_num);

    for (int i=0;i<var_num;i++)
    {

        m_charges[i].charge=x_m[i];
    }
}


void eFieldLagrangian::getInv()
{
    int var_num=m_chargeNum;
    int eq_num=m_elec_num;

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            M_[i][j]=0.0;
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=m_W[i][n]*m_W[j][n];  //its mvm
            }
        }
    }


    LU_decompose(var_num);

}

void eFieldLagrangian::getInv_PhiE()
{
    int var_num=m_chargeNum;
    int eq_num=2 * m_elec_num;

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            M_[i][j]=0.0;
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=m_W[i][n]*m_W[j][n];  //its mvm
            }
        }
    }


    LU_decompose(var_num);

}

void eFieldLagrangian::solve_ls_fast()
{
    int var_num=m_chargeNum;
    int eq_num=m_elec_num;



    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_electrodes[i].phi_fix-m_electrodes[i].phi_fix_charges;
    }

    for (int i=0;i<var_num;i++)
    {
        mwb[i]=0.0;
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=m_W[i][n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb
    }

    m_solve(var_num);

    for (int i=0;i<var_num;i++)
    {

        m_charges[i].charge=x_m[i];
    }
}

void eFieldLagrangian::solve_ls_fast_PhiE()
{
    int var_num=m_chargeNum;
    int eq_num=2 * m_elec_num;



    for (int i=0;i<m_elec_num;i++)
    {
        b_m[i]=m_electrodes[i].phi_fix-m_electrodes[i].phi_fix_charges;
    }
    for (int i=m_elec_num;i<eq_num;i++)
    {
        b_m[i]=0;
    }

    for (int i=0;i<var_num;i++)
    {
        mwb[i]=0.0;
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=m_W[i][n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb
    }

    m_solve(var_num);

    for (int i=0;i<var_num;i++)
    {

        m_charges[i].charge=x_m[i];
    }
}

vec2 eFieldLagrangian::get_E_multiplier(double x, double y, int i)
{

    vec2 sum;
    sum.x=0.0; sum.y=0.0;
    //double delta=1e-6;
    double d2=delta_phi*delta_phi;
   const double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);

          double r2;
        double q;
        double dx,dy;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r2=(dx*dx+dy*dy);
        q=-qepspi;//* (m_charges[i].charge);

        double c=q/(r2+d2);

        sum.x+=c*dx;
        sum.y+=c*dy;

    return sum;
}





vec2 eFieldLagrangian::getE(double x, double y)
{

    vec2 sum;
    sum.x=0.0; sum.y=0.0;
    //double delta=1e-6;
    double d2=delta_phi*delta_phi;
    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);

    for (int i=0;i<m_chargeNum;i++)
    {
        double r2;
        double q;
        double dx,dy;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r2=(dx*dx+dy*dy);
        q=-qepspi* (m_charges[i].charge);

        double c=q/(r2+d2);

        sum.x+=c*dx;
        sum.y+=c*dy;
    }

    return sum;
}

