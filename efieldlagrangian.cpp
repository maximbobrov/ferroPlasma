#include "efieldlagrangian.h"
#include "globals.h"


eFieldLagrangian::eFieldLagrangian()
{
    m_elec_num=200;
    m_electrodes=new eElem[m_elec_num];
    m_rCentre=new vec2[m_elec_num];

    m_dz=1e-7;//100 nm
    for (int i=0;i<m_elec_num/8;i++) //first electrode
    {
        double alpha=i*1.0/(m_elec_num/8-1);
        double x1,y1,x2,y2;
        x1 = w_x0;
        y1 = 0.5*(w_y0+w_y1)*(1.0-alpha)+w_y1*alpha;

        m_electrodes[i].r0.x=x1;
        m_electrodes[i].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/8-1);
        x2 = w_x0;
        y2 = 0.5*(w_y0+w_y1)*(1.0-alpha)+w_y1*alpha;

        m_electrodes[i].r1.x=x2;
        m_electrodes[i].r1.y=y2;

        m_electrodes[i].rho1=0.1;
        m_electrodes[i].rho2=0.1;

        m_electrodes[i].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

        m_electrodes[i].phi_fix=-1.0;
        m_electrodes[i].phi_fix_charges=0.0;

        if (alpha<0.5)
            m_electrodes[i].canEmit=true;
        else
            m_electrodes[i].canEmit=false;
    }

    for (int j=m_elec_num/8;j<m_elec_num;j++) //second electrode
    {
        int i=j-m_elec_num/8;

        double alpha=i*1.0/(m_elec_num*7.0/8-1);

        double x1,y1,x2,y2;
        x1 = w_x0*(1.0-alpha)+w_x1*alpha;
        y1 = w_y0;
        m_electrodes[j].r0.x=x1;
        m_electrodes[j].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/2-1);
        x2 = w_x0*(1.0-alpha)+w_x1*alpha;
        y2 = w_y0;

        m_electrodes[j].r1.x=x2;
        m_electrodes[j].r1.y=y2;

        m_electrodes[j].rho1=0.1;
        m_electrodes[j].rho2=0.1;

        m_electrodes[j].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        m_electrodes[j].phi_fix=1.0;
        m_electrodes[j].phi_fix_charges=0.0;
        m_electrodes[j].canEmit=false;
    }

    /* int elec_len=m_elec_num/2;
    for (int i=0;i<elec_len;i++) //first electrode
       {
           double alpha=i*1.0/(elec_len-1);
           double x1,y1,x2,y2;
           x1 = w_x0*(1.0-alpha)+w_x1*alpha;
           y1 = 0.5*(w_y1+w_y0)+5e-9;

           m_electrodes[i].r0.x=x1;
           m_electrodes[i].r0.z=0.0;
           m_electrodes[i].r0.y=y1-pow((rand()*1.0/RAND_MAX),10.0)*1e-9;

           alpha=(i+1)*1.0/(elec_len-1);
           x2 = w_x0*(1.0-alpha)+w_x1*alpha;
           y2 = 0.5*(w_y1+w_y0)+5e-9;

           m_electrodes[i].r1.x=x2;
           m_electrodes[i].r1.z=0.0;
           m_electrodes[i].r1.y=y2;

           m_electrodes[i].rho1=0.1;
           m_electrodes[i].rho2=0.1;

           m_electrodes[i].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

           m_electrodes[i].phi_fix=-1.0;
           m_electrodes[i].phi_fix_charges=0.0;

       }

       for (int j=m_elec_num/2;j<m_elec_num;j++) //second electrode
       {
           int i=j-m_elec_num/2;

           double alpha=i*1.0/(m_elec_num/2-1);
           double x1,y1,x2,y2;
           x1 = w_x0*(1.0-alpha)+w_x1*alpha;
           y1 = w_y0 - 5e-9;

           m_electrodes[j].r0.x=x1;
           m_electrodes[j].r0.z=0.0;
           m_electrodes[j].r0.y=y1+pow((rand()*1.0/RAND_MAX),10.0)*1e-9;//+((rand()*1.0/RAND_MAX)-0.5)*1e-8;

           alpha=(i+1)*1.0/(m_elec_num/2-1);
           x2 = w_x0*(1.0-alpha)+w_x1*alpha;
           y2 = w_y0 - 5e-9;

           m_electrodes[j].r1.x=x2;
           m_electrodes[j].r1.z=0.0;
           m_electrodes[j].r1.y=y2;

           m_electrodes[j].rho1=0.1;
           m_electrodes[j].rho2=0.1;

           m_electrodes[j].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
           m_electrodes[j].phi_fix=1.0;
           m_electrodes[i].phi_fix_charges=0.0;
       }
*/
    updateGridProp();
}

void eFieldLagrangian::updatePhi()
{
    for (int i=0;i<m_elec_num;i++)
    {
        double x,y;
        x=0.5*(m_electrodes[i].r0.x+m_electrodes[i].r1.x);
        y=0.5*(m_electrodes[i].r0.y+m_electrodes[i].r1.y);
        m_electrodes[i].phi_ext=getPhi(x,y);
    }
}

void eFieldLagrangian::solvePhi(int itn)
{
    for (int j=0;j<itn;j++)
    {
        updatePhi();
        for (int i=0;i<m_elec_num;i++)
        {
            m_electrodes[i].rho1+=1e-2*(m_electrodes[i].phi_fix+m_electrodes[i].phi_fix_charges-m_electrodes[i].phi_ext);
            m_electrodes[i].rho2=m_electrodes[i].rho1;
        }
    }
    updateGridProp();
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
        vec2 pos((m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5, (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5, (m_electrodes[i].rho1 + m_electrodes[i].rho2)*0.5);
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

double eFieldLagrangian::getPhi(double x, double y)
{
    vec2 phi;
    vec2 pos(x,y,0.0);
    //getFieldFast(pos, m_rCentre, getPhiField, phi);
    //return phi.x;
    double sum=0.0;
    for (int i=0;i<m_elec_num;i++)
    {
        //    int i=1;
        double r2;
        double q;
        double dxx,dyy;
        double delta=1e-9;

        dxx = (m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5 - x;
        dyy = (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5 - y;
        r2=sqrt(dxx*dxx+dyy*dyy);
        q=(m_electrodes[i].rho1+m_electrodes[i].rho2)*0.5;

        sum+=q*log(r2+delta)/log(delta);//q*delta/(r2+delta);
        //   sum+=q*pow(delta/(r2+delta),0.25);

    }
    return sum;
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

        m_electrodes[i].r0.x=x1;
        m_electrodes[i].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/8-1);
        x2 = _x0*(1.0-alpha) + _x1*alpha;
        y2 = _y0*(1.0-alpha) + _y1*alpha;;

        m_electrodes[i].r1.x=x2;
        m_electrodes[i].r1.y=y2;

        m_electrodes[i].rho1=0.1;
        m_electrodes[i].rho2=0.1;

        m_electrodes[i].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        m_electrodes[i].phi_fix=-1.0;
    }

    updatePhi();
    solvePhi(100);
}

vec2 eFieldLagrangian::getE(double x, double y)
{
    vec2 EField;
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

        dxx = (m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5 - x;
        dyy = (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5 - y;
        r2=(dxx*dxx+dyy*dyy);
        q=(m_electrodes[i].rho1+m_electrodes[i].rho2)*0.5;

        sum.x+=q*dxx/(r2+delta*delta)/log(delta);
        sum.y+=q*dyy/(r2+delta*delta)/log(delta);
    }
    //printf("ex = %e ex2 = %e  ey = %e ey2 = %e\n", EField.x, sum.x, EField.y, sum.y );
    return sum;
}
