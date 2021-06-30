#ifndef MULTISOLVER_H
#define MULTISOLVER_H


#include "efieldlagrangian.h"
#include "pzsolver.h"

#define traj_num_max 100

class multiSolver
{   
public:
    double dt_elec;
    eFieldLagrangian* m_Esolver;
    pzSolver* m_pzSolver;


    int endPosTable[1000];
    static vec2 trajectories[traj_num_max][2002];
    static int trajectories_num[traj_num_max];
    static int traj_num;

    int points_distributed; //points are distributed over the grid

    multiSolver();

    void updateEforPz();
   // void updateEforPz_self();
    void updateEforElec();

    double calcJ(double Ein);
    void electronEmissionEndMoveToElectrode(double dt);
    void updateTrajTable(bool leapfrog, double dt0, double dl0);
    int  getEndPos(int i);
    void solve(int itn);
    void solvePzAdaptive(double dtElec);
    void solvePzDOPRI();
        void solve__(int itn); //debug purposes
    void init();

    void pzEmission(double dt);

    double getPhi_at_electrode(int n);
    double getPhi_at_pz_down(int n);
    double getPhi_at_pz_up(int n);
    void prepare_caches(bool _2D_elems);

    void fast_Fields_prepare(); //get fields on eul grid
    void fast_Fields_recalculate();
    vec2 get_fast_E(double x,double y);//simple interpolation between self

    void near_Fields_recalculate_cell(int i_, int j_);
    void slower_Fields_recalculate();
    vec2 get_slower_E(double x, double y);//smart interpolation
      vec2 get_slow_E(double x, double y);//full summation
      double get_slow_phi(double x, double y, bool _2d);

      vec2 getE_grid_ij(int i, int j);
};

#endif // MULTISOLVER_H
