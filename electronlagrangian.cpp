#include "electronlagrangian.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
electronLagrangian::electronLagrangian()
{
    m_threadIdx = 0;
    m_maxParticles = 81920;
    m_numParticles = 0;
    m_bodyAccel = new vec2[m_maxParticles];
    m_bodyVel = new vec2[m_maxParticles];
    m_bodyE = new vec2[m_maxParticles];
    m_bodyPos = new vec2[m_maxParticles];
}


double electronLagrangian::calcJ(double Ein)
{

    double E=Ein/100;//from V/m to V/cm
    double t2 = 1.1;
    double B = 10.0;//145.0
    double phi = 4.0;
    double y = 3.79 * 1e-4 * sqrt(fabs(B * E)) / phi;
    double tetta = 0.95 - 1.03 * y * y;
    return 1e4*(1.54 * 1e-6 * B * B * E * E / (t2  * phi)) * exp( - 6.83 * 1e7 * pow(phi, 1.5) * tetta / fabs( B * E)); //in A/m^2
}


int electronLagrangian::create_electron(vec2 &pos, double Emag, double Dt, double ds)
{
    int num_in_pack=10.0;
    double el_to_add=1e-9*calcJ(Emag*30)*Dt*ds/(qe*num_in_pack);

    int ne=(int) el_to_add;
    ne+=((rand()*1.0)/RAND_MAX < (el_to_add - ne)); //extra electron



    printf("Emag=%e el_to_ad=%e ne=%d \n",Emag,el_to_add,ne);
    int upto=MIN(m_numParticles+ne,m_maxParticles);
    for (int n = m_numParticles; n < upto; ++n)
    {
        m_bodyPos[n].x = pos.x+(rand()*2e-9/RAND_MAX)+1e-9;
        m_bodyPos[n].y = pos.y+(rand()*2e-9/RAND_MAX-1e-9);
        m_bodyPos[n].charge = num_in_pack;
        m_bodyVel[n].x = 0.0;
        m_bodyVel[n].y = 0.0;
        m_bodyAccel[n].x = 0.0;
        m_bodyAccel[n].y = 0.0;
    }
    m_numParticles=upto;
    return num_in_pack*ne;
}


void electronLagrangian::create_random_particles()
{
    /*for (int i = (N_Y-1)/2+2; i < (N_Y-1); i+=1)
    {
        double dz = w_z1 - w_z0;
        double x_ = w_x0 + 0.01 * (w_x1-w_x0);//- N_X * dx / 2;
        double y_ =w_y0 + i * (w_y1 - w_y0) / (N_Y-1);//i * dy - N_Y * dy / 2;
        double z_ = my_rand(0) * dz - dz / 2;
        vec3<float> E;
        getEFromElectrons(E, x_, y_, z_, m_numParticles);
        //printf("Ex=%e Ey = %e\n", E.x, E.y);
        E.x = Ex[0][i];
        E.y = Ey[0][i];
        E.x = fmax(0, -E.x);
        double J =calcJ(0.02 * E.x);
        int chargeNum = 1e0;
        int curNum = m_numParticles;
        int numToAdd =  std::min(int(dt * J * dy * dz * 1e4 * 6.24151 * 1e18 / chargeNum), m_maxParticles - m_numParticles-2);
        numToAdd = numToAdd > 5 ? 5 : numToAdd;

        if(my_rand(m_threadIdx) < (dt * J * dy * dz * 1e4 * 6.24151 * 1e18 / chargeNum) - int(dt * J * dy * dz * 1e4  * 6.24151 * 1e18 / chargeNum))
            numToAdd++;
        m_numParticles += numToAdd;

        for (int n = curNum; n < curNum + numToAdd; ++n)
        {
            m_bodyPos[n].x = x_;
            m_bodyPos[n].y = y_;
            m_bodyPos[n].z = z_;
            m_bodyPos[n].w = 1.6e-19 * chargeNum;
            double angle =  0.98 * acos(1 - 2 * my_rand(m_threadIdx)) - M_PI / 2;
            //printf("Ey = %f\n", angle);
            if (y_ < w_y0 + 0.1 * (w_y1 - w_y0) && angle < 0)
                angle *= -1;
            double en = getVms_from_Ev(0.00002);
            m_bodyVel[n].x = 0.0;//en * cos(angle);
            m_bodyVel[n].y = 0.0;//en * sin(angle);
            m_bodyVel[n].z = 0;
            m_bodyAccel[n].x = 0.0;
            m_bodyAccel[n].y = 0.0;
            m_bodyAccel[n].z = 0.0;
        }
    }*/
    double x_ = w_x0 + 0.01 * (w_x1-w_x0);
    double y_ = my_rand(m_threadIdx) * (w_y1);

    int curNum = m_numParticles;
    int numToAdd =  std::min(5, m_maxParticles - m_numParticles-2);
    m_numParticles += numToAdd;

    for (int n = curNum; n < curNum + numToAdd; ++n)
    {
        m_bodyPos[n].x = x_;
        m_bodyPos[n].y = y_;
        m_bodyPos[n].charge = qe;
        double angle =  0.98 * acos(1 - 2 * my_rand(m_threadIdx)) - M_PI / 2;
        if (y_ < w_y0 + 0.1 * (w_y1 - w_y0) && angle < 0)
            angle *= -1;
        double en = getVms_from_Ev(0.00002);
        m_bodyVel[n].x = 0.0;//en * cos(angle);
        m_bodyVel[n].y = 0.0;//en * sin(angle);
        m_bodyAccel[n].x = 0.0;
        m_bodyAccel[n].y = 0.0;
    }
}

vec2 getE(float x, float y)
{
    double a0=fmax(fmin(N_X-1,(N_X-1)*(x-w_x0)/(w_x1-w_x0)),0);
    int i0=(int)a0;
    a0-=i0;

    double b0=fmax(fmin(N_Y-1,(N_Y-1)*(y-w_y0)/(w_y1-w_y0)),0);
    int j0=(int)b0;
    b0-=j0;

    vec2 ret;

    ret.x=(1.0-b0)*((1.0-a0)*Ex[i0][j0]+(a0)*Ex[i0+1][j0]) + (b0)*((1.0-a0)*Ex[i0][j0+1]+(a0)*Ex[i0+1][j0+1]);
    ret.y=(1.0-b0)*((1.0-a0)*Ey[i0][j0]+(a0)*Ey[i0+1][j0]) + (b0)*((1.0-a0)*Ey[i0][j0+1]+(a0)*Ey[i0+1][j0+1]);
    return ret;
}

void electronLagrangian::delete_particle(int particlesIdx)
{
    m_bodyPos[particlesIdx] = m_bodyPos[m_numParticles-1];
    m_bodyVel[particlesIdx] = m_bodyVel[m_numParticles-1];
    m_bodyAccel[particlesIdx] = m_bodyAccel[m_numParticles-1];
    m_numParticles -= 1;
}

void electronLagrangian::wall_collision(int particlesIdx)
{
    int xIdx = (m_bodyPos[particlesIdx].x - w_x0)/dx;
    int yIdx = (m_bodyPos[particlesIdx].y - w_y0)/dy;
    if ((xIdx>=0)&&(xIdx<N_X))
    {
        vec2 ev=getE(m_bodyPos[particlesIdx].x, dy);

        //if(ev.y>0)
        {
            q[xIdx]+=m_bodyPos[particlesIdx].charge/2.75;
            double addit=0.5;
            for (int i=xIdx-1;i>=fmax(xIdx-3,0);i--)
            {
                q[i]+=addit*m_bodyPos[particlesIdx].charge/2.75;
                addit*=0.5;
            }

            addit=0.5;
            for (int i=xIdx+1;i<=fmin(xIdx+3,N_X-1);i++)
            {
                q[i]+=addit*m_bodyPos[particlesIdx].charge/2.75;
                addit*=0.5;
            }
            delete_particle(particlesIdx);
        }
    }
}

void electronLagrangian::step(double dt)
{
    int i;
    // create_random_particles();

    for( i=0; i<m_numParticles; i++ ) {
        float magn=100.0;//qe/Me;//1e-1;
        vec2 ev=m_bodyE[i];
        m_bodyVel[i].x -= magn*(ev.x)*dt;
        m_bodyVel[i].y -= magn*(ev.y)*dt;
    }

    for( i=0; i<m_numParticles; i++ )
    {
        m_bodyPos[i].x += dt*m_bodyVel[i].x;
        m_bodyPos[i].y += dt*m_bodyVel[i].y;

        /*   if (m_bodyPos[i].y<0)
        {
            wall_collision(i);
        }
*/
        if (m_bodyPos[i].y>2.0*w_y1 || m_bodyPos[i].x < w_x0-1e-7 || m_bodyPos[i].x>w_x1)
        {
            delete_particle(i);
        }

    }
    updateGridProp();
}

void electronLagrangian::updateGridProp()
{
    m_gridProp.NX = 20;
    m_gridProp.NY = 20;
    m_gridProp.dx = (w_x1-w_x0) / (m_gridProp.NX  - 1);
    m_gridProp.dy = (w_y1-w_y0) / (m_gridProp.NY  - 1);
    m_gridProp.startx = w_x0;
    m_gridProp.starty = w_y0;

    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            m_gridProp.gridCenters[i][j].x = 0.0;
            m_gridProp.gridCenters[i][j].y = 0.0;
            m_gridProp.gridCenters[i][j].charge = 0.0;
            m_gridProp.gridNeighbors[i][j].clear();
        }
    for( int i=0; i<m_numParticles; i++ ){
        int xIdx = int((m_bodyPos[i].x - w_x0) / m_gridProp.dx);
        int yIdx = int((m_bodyPos[i].y - w_y0) / m_gridProp.dy);
        m_gridProp.gridCenters[xIdx][yIdx].charge += m_bodyPos[i].charge;
        m_gridProp.gridCenters[xIdx][yIdx].x += m_bodyPos[i].charge * m_bodyPos[i].x;
        m_gridProp.gridCenters[xIdx][yIdx].y += m_bodyPos[i].charge * m_bodyPos[i].y;
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

vec2 electronLagrangian::getEField(const vec2& iCenterPos, const vec2& iFarPos)
{
    vec2 E;
    vec2 dist;
    float invDist2;

    dist.x = iCenterPos.x - iFarPos.x;
    dist.y = iCenterPos.y - iFarPos.y;

    invDist2 = iFarPos.charge / (dist.x*dist.x+dist.y*dist.y+1e-14);

    E.x = -dist.x*invDist2;
    E.y = -dist.y*invDist2;
    return  E;
}

vec2 electronLagrangian::getEe(double x, double y)
{
    //double t0 = get_time();
    vec2 ai = {0.0, 0.0, 0.0};
    vec2 iPos(x, y, 0.0);
    vec2 EField;
    // printf("num = %d\n", m_numParticles);
    getFieldFast(iPos, m_bodyPos, getEField, EField);
    ai.x = EField.x;
    ai.y = EField.y;
    ai.x/=eps0;
    ai.y/=eps0;
    return ai;
   /* double t1 = get_time();

    vec2 dist;
    float invDist2;
    vec2 ai2 = {0.0, 0.0, 0.0};
    for( int j=0; j<m_numParticles; j++ ){
        dist.x = m_bodyPos[j].x-x;
        dist.y = m_bodyPos[j].y-y;
        invDist2 = m_bodyPos[j].charge/(dist.x*dist.x+dist.y*dist.y+1e-14);
        ai2.x -= dist.x*invDist2;
        ai2.y -= dist.y*invDist2;
    }
    ai2.x/=eps0;
    ai2.y/=eps0;
    //double t2 = get_time();
    //printf("ex = %e ex2 = %e  ey = %e ey2 = %e\n", ai.x, ai2.x, ai.y, ai2.y );
    //printf("t1 = %e t2 = %e\n",(t1-t0), (t2-t1));
    return ai;*/
}

void electronLagrangian::getEFromElectrons(vec2 &bodyAccel_, double x, double y, double z,  int n)
{
    vec2 dist;
    float invDist,invDistCube;
    vec2 ai = {0.0, 0.0, 0.0};
    for( int j=0; j<n; j++ ){
        dist.x = m_bodyPos[j].x-x;
        dist.y = m_bodyPos[j].y-y;
        invDist = 1.0/sqrtf(dist.x*dist.x+dist.y*dist.y+1e-14);
        invDistCube = m_bodyPos[j].charge*invDist*invDist*invDist;
        ai.x -= dist.x*invDistCube;
        ai.y -= dist.y*invDistCube;
    }

    bodyAccel_.x = inv4PI*ai.x/eps0;
    bodyAccel_.y = inv4PI*ai.y/eps0;
}
