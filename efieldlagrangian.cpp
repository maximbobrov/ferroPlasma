#include "efieldlagrangian.h"
#include "globals.h"

eFieldLagrangian::eFieldLagrangian()
{
    m_elec_num=160;
    m_electrodes=new eElem[m_elec_num];

    for (int i=0;i<m_elec_num/8;i++) //first electrode
    {
        double alpha=i*1.0/(m_elec_num/8-1);
        double x1,y1,x2,y2;
        x1 = w_x0;
        y1 = 0.5*(w_y0+w_y1)*(1.0-alpha)+w_y1*alpha;

        m_electrodes[i].r0.x=x1;
        m_electrodes[i].r0.z=0.0;
        m_electrodes[i].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/8-1);
        x2 = w_x0;
        y2 = 0.5*(w_y0+w_y1)*(1.0-alpha)+w_y1*alpha;

        m_electrodes[i].r1.x=x2;
        m_electrodes[i].r1.z=0.0;
        m_electrodes[i].r1.y=y2;

        m_electrodes[i].rho1=0.1;
        m_electrodes[i].rho2=0.1;

        m_electrodes[i].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

        m_electrodes[i].phi_fix=-1.0;
        m_electrodes[i].phi_fix_charges=0.0;

    }

    for (int j=m_elec_num/8;j<m_elec_num;j++) //second electrode
    {
        int i=j-m_elec_num/8;

        double alpha=i*1.0/(m_elec_num*7.0/8-1);
        double x1,y1,x2,y2;
        x1 = w_x0*(1.0-alpha)+w_x1*alpha;
        y1 = w_y0;

        m_electrodes[j].r0.x=x1;
        m_electrodes[j].r0.z=0.0;
        m_electrodes[j].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/2-1);
        x2 = w_x0*(1.0-alpha)+w_x1*alpha;
        y2 = w_y0;

        m_electrodes[j].r1.x=x2;
        m_electrodes[j].r1.z=0.0;
        m_electrodes[j].r1.y=y2;

        m_electrodes[j].rho1=0.1;
        m_electrodes[j].rho2=0.1;

        m_electrodes[j].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        m_electrodes[j].phi_fix=1.0;
        m_electrodes[i].phi_fix_charges=0.0;
    }
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
}

double eFieldLagrangian::getPhi(double x, double y)
{
    double sum=0.0;
    for (int i=0;i<m_elec_num;i++)
    {
    //    int i=1;
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = (m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5 - x;
        dy = (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5 - y;
        r2=sqrt(dx*dx+dy*dy);
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
        m_electrodes[i].r0.z=0.0;
        m_electrodes[i].r0.y=y1;

        alpha=(i+1)*1.0/(m_elec_num/8-1);
        x2 = _x0*(1.0-alpha) + _x1*alpha;
        y2 = _y0*(1.0-alpha) + _y1*alpha;;

        m_electrodes[i].r1.x=x2;
        m_electrodes[i].r1.z=0.0;
        m_electrodes[i].r1.y=y2;

        m_electrodes[i].rho1=0.1;
        m_electrodes[i].rho2=0.1;

        m_electrodes[i].dl=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        m_electrodes[i].phi_fix=-1.0;
    }
    updatePhi();
    solvePhi(100);
}

vec3<double> eFieldLagrangian::getE(double x, double y)
{
    vec3<double> sum;
    sum.x=0.0; sum.y=0.0; sum.z=0.0;

    for (int i=0;i<m_elec_num;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = (m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5 - x;
        dy = (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5 - y;
        r2=sqrt(dx*dx+dy*dy);
        q=(m_electrodes[i].rho1+m_electrodes[i].rho2)*0.5;    

        sum.x+=q*delta*dx/(r2+delta);
        sum.y+=q*delta*dy/(r2+delta);

    }
    return sum;
}
