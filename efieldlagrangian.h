#ifndef EFIELDLAGRANGIAN_H
#define EFIELDLAGRANGIAN_H
#include "globals.h"
#include "lagrangiansolver.h"

class eFieldLagrangian : lagrangianSolver
{

public:
    struct eElem   //linear electrode elem
    {
        double phi_ext,phi_fix,phi_fix_charges;//C/m2, m, V
        vec2 r;
        bool canEmit;
        double eToEmit;
        double EdotN;
        //double charge;
        double dl;
        double eCurrent;
        double nx, ny;
        int INDX;
    };
    vec2* m_rCentre;
    eElem * m_electrodes;
    vec2 m_charges[2000];
    vec2 m_mirrorCharges[2000];
    int m_chargeNum;
    int m_elec_num;

public:

    eFieldLagrangian();
    void init();
    void updatePhi(); //get external potential on electrodes
    void solvePhi(int itn); //update potential
    vec2 getE(double x,double y);
    double getPhi(double x,double y);
    void updateGridProp();
    void addLine(vec2 p[2], double dl,double phi, double h, int I);
    void addQuad(vec2 p[4], double dl[4], double phi, int emit[4], double coordYDIel, int smoothingCount, int I);
    void addQuad2Layers(vec2 p[4], double dl[4], double phi, int emit[4]);
    void addQuadRegular(vec2 p[4], double dl, double phi, int emit[4]);
    void addQuad_stabilized(vec2 p[], double dl[], double phi, int emit[], double coordYDIel, int smoothingCount, int I);

    void initW();
    void initW_PhiE();
    double getW(double s_x, double s_y, double t_x, double t_y);
    double getWMirror(int iCharge, int iElec);
    double getW_E(int elecNum, int chargeNum);
    void solve_ls();
    void getInv();
    void getInv_PhiE();
    void solve_ls_fast();
    void solve_ls_fast_PhiE();

};
#endif // EFIELDLAGRANGIAN_H
