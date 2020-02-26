#ifndef EFIELDLAGRANGIAN_H
#define EFIELDLAGRANGIAN_H
#include "sse_sum.h"

class eFieldLagrangian
{

public:
    struct eElem   //linear electrode elem
    {
        double rho1, rho2,dl,phi_ext,phi_fix;//C/m2, m, V
        vec3<double> r0,r1;
    };

    eElem * m_electrodes;
    int m_elec_num;

public:
    eFieldLagrangian();
    void updatePhi(); //get external potential on electrodes
    void solvePhi(int itn); //update potential
    vec3<double> getE(double x,double y);
    double getPhi(double x,double y);

    void setElectrodeAngle(double deg);

};



#endif // EFIELDLAGRANGIAN_H
