#ifndef ELECTRONLAGRANGIAN_H
#define ELECTRONLAGRANGIAN_H
#include "globals.h"
#include "lagrangiansolver.h"

class electronLagrangian : lagrangianSolver
{
public:

    int m_threadIdx;
    int m_maxParticles;
    int m_numParticles;
    vec2* m_bodyAccel;
    vec2* m_bodyPos;
    vec2* m_bodyVel;
    vec2* m_bodyE;

    electronLagrangian();
    int create_electron(vec2 &pos, double Emag, double Dt, double ds); //creates electrons near the electrode form field emission
    void delete_particle(int particlesIdx);
    void wall_collision(int particlesIdx);
    double calcJ(double Ein);
    void step(double dt);
    void create_random_particles();
    vec2 getEe(double x, double y);
    void getEFromElectrons(vec2 &bodyAccel_, double x, double y, double z,  int n);
    static vec2 getEField(const vec2& iPos1, const vec2& iPos2);
    void updateGridProp();

};

#endif // ELECTRONLAGRANGIAN_H
