#include "efieldlagrangian.h"
#include "globals.h"



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

    m_elec_num=0;
    m_electrodes=new eElem[2000];
    m_rCentre=new vec2[2000];
    m_chargeNum = 0 ;

    vec2 p[5];

    p[0].x=w_x0;        p[0].y=0.5*(w_y0+w_y1);
    p[1].x=w_x0;        p[1].y=w_y1;
    p[2].x=w_x0-50e-9; p[2].y=w_y1;
    p[3].x=w_x0-50e-9; p[3].y=0.5*(w_y0+w_y1);
    p[4].x=p[0].x;      p[4].y=p[0].y;

    addQuad(p,2e-9,-1,0);


    p[0].x=w_x0;        p[0].y=w_y0;
    p[1].x=w_x1;        p[1].y=w_y0;
    p[2].x=w_x1;        p[2].y=w_y0-10e-9;
    p[3].x=w_x0;        p[3].y=w_y0-10e-9;
    p[4].x=p[0].x;      p[4].y=p[0].y;

    addQuad(p,2e-9,1,-1);

    initW();

    getInv();

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

void eFieldLagrangian::updateGridProp()
{
    m_gridProp.NX = 30;
    m_gridProp.NY = 20;
    m_gridProp.startx = w_x0 - (w_y1-w_y0)*0.5;
    m_gridProp.starty = w_y0 ;
    m_gridProp.endx = w_x1 ;
    m_gridProp.endy = w_y1 ;
    m_gridProp.dx = (m_gridProp.endx - m_gridProp.startx) / (m_gridProp.NX);
    m_gridProp.dy = (m_gridProp.endy - m_gridProp.starty) / (m_gridProp.NY);

    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            m_gridProp.gridCenters[i][j].x = 0.0;
            m_gridProp.gridCenters[i][j].y = 0.0;
            m_gridProp.gridCenters[i][j].charge = 0.0;
            m_gridProp.gridNeighbors[i][j].clear();
        }
    for( int i=0; i<m_elec_num; i++ ){
        vec2 pos(m_electrodes[i].r.x , m_electrodes[i].r.y , 0.0/*(m_electrodes[i].rho1 + m_electrodes[i].rho2)*0.5*/);
        m_rCentre[i] = pos;
        int xIdx = fmin(int((pos.x - m_gridProp.startx) / m_gridProp.dx), m_gridProp.NX - 1);
        int yIdx = fmin(int((pos.y - m_gridProp.starty) / m_gridProp.dy), m_gridProp.NY - 1);
        m_gridProp.gridCenters[xIdx][yIdx].charge += pos.charge;
        m_gridProp.gridCenters[xIdx][yIdx].x += pos.charge * pos.x;
        m_gridProp.gridCenters[xIdx][yIdx].y += pos.charge * pos.y;
        m_gridProp.gridNeighbors[xIdx][yIdx].push_back(i);
    }
    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            if(m_gridProp.gridNeighbors[i][j].size() == 0)
                continue;
            m_gridProp.gridCenters[i][j].x /= m_gridProp.gridCenters[i][j].charge;
            m_gridProp.gridCenters[i][j].y /= m_gridProp.gridCenters[i][j].charge;
        }
}

void eFieldLagrangian::addQuad(vec2 p[5], double dl,double phi, int emit) //last point should coincide with the first one emit is the side number that can emit
{
    static double l_[4];
    static int n_[4];
    int n0=m_elec_num;
    vec2 c_m(0,0,0);
    for (int i=0;i<4;i++)
    {
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

            if (i==emit)
                m_electrodes[m_elec_num].canEmit=true;
            else
                m_electrodes[m_elec_num].canEmit=false;
            m_elec_num++;

            c_m.x+=x;
            c_m.y+=y;
            c_m.charge+=1.0;
        }
        printf("i=%d ni=%d li=%e cuur_num=%d \n",i,n_[i],l_[i],m_elec_num);
    }
    c_m.x/=c_m.charge;
    c_m.y/=c_m.charge;

    int n1=m_elec_num;

    for (int i=n0/2;i<n1/2;i++)
    {
        double x,y;
        x = 0.975*((m_electrodes[i*2].r.x + m_electrodes[i*2+1].r.x)*0.5-c_m.x)+c_m.x;
        y = 0.975*((m_electrodes[i*2].r.y + m_electrodes[i*2+1].r.y)*0.5-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }

    for (int i=n0/4;i<n1/4;i++)
    {
        double x,y;
        x = 0.75*((m_electrodes[i*4].r.x + m_electrodes[i*4+1].r.x + m_electrodes[i*4+2].r.x + m_electrodes[i*4+3].r.x)*0.25-c_m.x)+c_m.x;
        y = 0.75*((m_electrodes[i*4].r.y + m_electrodes[i*4+1].r.y + m_electrodes[i*4+2].r.y + m_electrodes[i*4+3].r.y)*0.25-c_m.y)+c_m.y;

        m_charges[m_chargeNum].x=x;
        m_charges[m_chargeNum].y=y;
        m_charges[m_chargeNum].charge=0.0;
        m_chargeNum++;
    }

    m_charges[m_chargeNum].x=c_m.x;
    m_charges[m_chargeNum].y=c_m.y;
    m_charges[m_chargeNum].charge=0.0;
    m_chargeNum++;
}

double eFieldLagrangian::getW(double s_x, double s_y,double t_x, double t_y) //get Weight function ;//source (charge) and target (monitoring point)
{

    double sum=0.0;
    //    int i=1;
    double r;
    double q;
    double dx,dy;
    double delta=1e-9;

    dx = s_x - t_x;
    dy = s_y - t_y;
    r=sqrt(dx*dx+dy*dy);
    //       q=qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);

    //sum-=-q*log(r+delta)/(w_z1 - w_z0);
    sum=-(qe/(eps0*pi2))*log(r+delta)/(w_z1 - w_z0);
    return sum;
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
        double delta=1e-9;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (m_charges[i].charge);

        sum+=-q*log(r+delta)/(w_z1 - w_z0);
    }
    return sum;
}

void eFieldLagrangian::initW()
{
    for(int i=0;i<m_chargeNum;i++)
    {
        for(int j=0;j<m_elec_num;j++)
        {
            m_W[i][j]=getW(m_charges[i].x,m_charges[i].y,m_electrodes[j].r.x,m_electrodes[j].r.y);
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



vec2 eFieldLagrangian::getEField(const vec2& iCenterPos, const vec2& iFarPos)
{
    vec2 E;
    vec2 dist;
    float invDist2;
    double delta=1e-9;

    dist.x = iCenterPos.x - iFarPos.x;
    dist.y = iCenterPos.y - iFarPos.y;

    invDist2 = iFarPos.charge / (dist.x*dist.x+dist.y*dist.y+delta*delta)/log(delta);

    E.x = dist.x*invDist2;
    E.y = dist.y*invDist2;
    return  E;
}

vec2 eFieldLagrangian::getPhiField(const vec2& iCenterPos, const vec2& iFarPos)
{
    vec2 E;
    vec2 dist;
    float invDist2;
    double delta=1e-9;

    dist.x = iCenterPos.x - iFarPos.x;
    dist.y = iCenterPos.y - iFarPos.y;

    E.x = iFarPos.charge / log(dist.x*dist.x+dist.y*dist.y + delta)/log(delta);
    return  E;
}



void eFieldLagrangian::setElectrodeAngle(double deg)
{
    double _x0,_y0;
    double _x1,_y1;
    double nx,ny;

    nx=sin(deg*M_PI/180.0);
    ny=cos(deg*M_PI/180.0);

    double l = 0.5*(w_y1-w_y0);
    _x0 = w_x0; _y0 = 0.5*(w_y0+w_y1);
    _x1 = _x0+nx*l; _y1 = _y0+ny*l;

    for (int i=0;i<m_elec_num/8;i++) //first electrode
    {
        double alpha=i*1.0/(m_elec_num/8-1);
        double x1,y1,x2,y2;
        x1 = _x0*(1.0-alpha) + _x1*alpha;
        y1 = _y0*(1.0-alpha) + _y1*alpha;;

        m_electrodes[i].r.x=x1;
        m_electrodes[i].r.y=y1;
    }

    updatePhi();
    solvePhi(100);
}

vec2 eFieldLagrangian::getE(double x, double y)
{
    /*vec2 EField;
    vec2 pos(x,y,0.0);
    //getFieldFast(pos, m_rCentre, getEField, EField);
    //return EField;
    vec2 sum;
    sum.x=0.0; sum.y=0.0;

    for (int i=0;i<m_elec_num;i++)
    {
        double r2;
        double q;
        double dxx,dyy;
        double delta=1e-9;

        dxx = m_electrodes[i].r.x  - x;
        dyy = m_electrodes[i].r.y  - y;
        r2=(dxx*dxx+dyy*dyy);
        q=0.0;//(m_electrodes[i].rho1+m_electrodes[i].rho2)*0.5;

        sum.x+=q*dxx/(r2+delta*delta)/log(delta);
        sum.y+=q*dyy/(r2+delta*delta)/log(delta);
    }
    //printf("ex = %e ex2 = %e  ey = %e ey2 = %e\n", EField.x, sum.x, EField.y, sum.y );
    return sum;*/

    vec2 EField;
    vec2 pos(x,y,0.0);
    //getFieldFast(pos, m_rCentre, getEField, EField);
    //return EField;
    vec2 sum;
    sum.x=0.0; sum.y=0.0;
    for (int i=0;i<m_chargeNum;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = m_charges[i].x - x;
        dy = m_charges[i].y - y;
        r2=(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (m_charges[i].charge);

        double c=q/((r2+delta*delta)*(w_z1 - w_z0));

        sum.x+=c*dx;
        sum.y+=c*dy;

    }
    return sum;
}

