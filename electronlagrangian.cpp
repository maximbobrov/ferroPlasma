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

    /*for (int i = 0; i <100;i++) {
        m_bodyPos[i].x = w_x0 + i * (w_x1 - w_x0)/100;
        m_bodyPos[i].y = 10e-9;
        m_bodyPos[i].charge = 1000;
        m_bodyVel[i].x = 0.0;
        m_bodyVel[i].y = 0.0;
        m_bodyAccel[i].x = 0.0;
        m_bodyAccel[i].y = 0.0;
    }
    */
    /*m_bodyPos[0].x = w_x0 + (w_x1 - w_x0)/2- 225e-9;
    m_bodyPos[0].y = w_y0 + (w_y1 - w_y0)/2 + 25e-9;
    m_bodyPos[0].charge = 100;
    m_bodyVel[0].x = 0.0;
    m_bodyVel[0].y = 0.0;
    m_bodyAccel[0].x = 0.0;
    m_bodyAccel[0].y = 0.0;


    m_bodyPos[1].x = w_x0 + (w_x1 - w_x0)/2+ 25e-9;
    m_bodyPos[1].y = w_y0 + (w_y1 - w_y0)/2 +25e-9;
    m_bodyPos[1].charge = 100;
    m_bodyVel[1].x = 0.0;
    m_bodyVel[1].y = 0.0;
    m_bodyAccel[1].x = 0.0;
    m_bodyAccel[1].y = 0.0;
    updateGridProp();*/
}

void electronLagrangian::init()
{
    m_threadIdx = 0;
    m_numParticles = 0;
}


double electronLagrangian::calcJ(double Ein)
{
    double E=Ein/100;//from V/m to V/cm
    double t2 = 1.1;
    double B = 145.0;//145.0
    double phi = 4.0;
    double y = 3.79 * 1e-4 * sqrt(fabs(B * E)) / phi;
    double tetta = 0.95 - 1.03 * y * y;
    return 1e4*(1.54 * 1.0e-6 * B * B * E * E / (t2  * phi)) * exp ( - 6.83 * 1.0e7 * pow(phi, 1.5) * tetta / fabs( B * E)); //in A/m^2
}


int electronLagrangian::create_electron(vec2 &pos, double Emag, double Dt, double ds)
{

    int num_in_pack=10.0;
    double el_to_add = calcJ(Emag)*Dt*ds/(fabs(qe)/**num_in_pack*/);
    num_in_pack = int(el_to_add/20+1);
    el_to_add /= num_in_pack;


    int ne=(int) el_to_add;
    ne+=((rand()*1.0)/RAND_MAX < (el_to_add - ne)); //extra electron



    // printf("Emag=%e j=%e el_to_ad=%e ne=%d \n", Emag, calcJ(Emag), el_to_add,ne);
    int upto=MIN(m_numParticles+ne,m_maxParticles-1);
    for (int n = m_numParticles; n < upto; ++n)
    {
        m_bodyPos[n].x = pos.x+(rand()*2e-9/RAND_MAX)+1e-9;
        m_bodyPos[n].y = pos.y+(rand()*2e-9/RAND_MAX-1e-9);
        m_bodyPos[n].charge = num_in_pack;
        m_bodyVel[n].x = 0.0;//1000000.0;
        m_bodyVel[n].y = 0.0;
        m_bodyAccel[n].x = 0.0;
        m_bodyAccel[n].y = 0.0;
    }
    m_numParticles=upto;
    return num_in_pack*ne;
}

void electronLagrangian::create_electrons(vec2 &pos, vec2 &vel, int num)
{
    int num_in_pack=max(num/5,500);
    int left=num;
    while ((left>num_in_pack)&&(m_numParticles<m_maxParticles-1))
    {
        m_bodyPos[m_numParticles].x = pos.x+(rand()*1.0e-6/RAND_MAX)-0.5e-6;
        m_bodyPos[m_numParticles].y = pos.y+(rand()*1.0e-6/RAND_MAX)-0.5e-6;
        m_bodyPos[m_numParticles].charge = num_in_pack;
        m_bodyVel[m_numParticles].x =vel.x+(rand()*1.0e5/RAND_MAX)-0.5e5;
        m_bodyVel[m_numParticles].y = vel.y+(rand()*1.0e5/RAND_MAX)-0.5e5;
        m_bodyAccel[m_numParticles].x = 0.0;
        m_bodyAccel[m_numParticles].y = 0.0;
        m_numParticles++;
        left-=num_in_pack;
    }
    if ((left>0)&&(m_numParticles<m_maxParticles-1))
    {
        m_bodyPos[m_numParticles].x = pos.x+(rand()*1.0e-6/RAND_MAX)-0.5e-6;
        m_bodyPos[m_numParticles].y = pos.y+(rand()*1.0e-6/RAND_MAX)-0.5e-6;
        m_bodyPos[m_numParticles].charge = left;
        m_bodyVel[m_numParticles].x =vel.x+(rand()*1.0e5/RAND_MAX)-0.5e5;
        m_bodyVel[m_numParticles].y = vel.y+(rand()*1.0e5/RAND_MAX)-0.5e5;
        m_bodyAccel[m_numParticles].x = 0.0;
        m_bodyAccel[m_numParticles].y = 0.0;
        m_numParticles++;
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
        double magn=qe/Me;//1e-1;
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
        if (m_bodyPos[i].y>2.0*w_y1 || m_bodyPos[i].x < w_x0-100e-6 || m_bodyPos[i].x>w_x1)
        {
            delete_particle(i);
        }

    }
    updateGridProp();
}

void electronLagrangian::updateGridProp()
{
    m_gridProp.NX = 20;
    m_gridProp.NY = 5;
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

void electronLagrangian::create_pz_electron(double x, double y, int q)
{
    int q_left=q;
    //while((m_numParticles<m_maxParticles)&&(q_left>0))
    {
        m_bodyPos[m_numParticles].x = x+(rand()*2.0/RAND_MAX-1.0)*1e-9;
        m_bodyPos[m_numParticles].y = y+(rand()*2.0/RAND_MAX-1.0)*1e-9;
        m_bodyPos[m_numParticles].charge = q;
        m_bodyVel[m_numParticles].x = 0.0;
        m_bodyVel[m_numParticles].y = 0.0;
        m_bodyAccel[m_numParticles].x = 0.0;
        m_bodyAccel[m_numParticles].y = 0.0;
        m_numParticles++;
        //q_left--;
    }
}

vec2 electronLagrangian::getEField(const vec2& iFarPos, const vec2& iCenterPos)
{
    vec2 E;
    vec2 dist;
    float invDist2;

    double delta=1e-9;

    dist.x = iCenterPos.x - iFarPos.x;
    dist.y = iCenterPos.y - iFarPos.y;
    double r2 = (dist.x*dist.x+dist.y*dist.y);

    double q=qe/(eps0*pi2) * (iCenterPos.charge);
    invDist2 = -q / ((r2+delta*delta)*(w_z1 - w_z0));

    E.x = dist.x*invDist2;
    E.y = dist.y*invDist2;
    return  E;
}

vec2 electronLagrangian::getEe(double x, double y)
{
    double t0 = get_time();
    /*vec2 ai = {0.0, 0.0, 0.0};
    vec2 iPos(x, y, 0.0);
    vec2 EField;
    // printf("num = %d\n", m_numParticles);
    getFieldFast(iPos, m_bodyPos, getEField, EField);
    ai.x = EField.x;
    ai.y = EField.y;
    return ai;*/
    double t1 = get_time();

    vec2 dist;
    float invDist2;
        double delta=1e-9;
    vec2 ai2 = {0.0, 0.0, 0.0};
    for( int j=0; j<m_numParticles; j++ ){
        dist.x = x-m_bodyPos[j].x;
        dist.y = y-m_bodyPos[j].y;
        double r2 = (dist.x*dist.x+dist.y*dist.y);
        double q=qe/(eps0*pi2) * (m_bodyPos[j].charge);
        invDist2 = -q / ((r2+delta*delta)*(w_z1 - w_z0));
        ai2.x -= dist.x*invDist2;
        ai2.y -= dist.y*invDist2;
    }

    double t2 = get_time();

    //printf("ex = %e ex2 = %e  ey = %e ey2 = %e\n", ai.x, ai2.x, ai.y, ai2.y );
    //printf("t1 = %e t2 = %e\n",(t1-t0), (t2-t1));
    return ai2;
}

double electronLagrangian::getEmult_dipole(double d) //calculate E in the middle two elementary charges with distance d
{
      double delta=1e-9;

      double r2 = d*d;
      double q=qe/(eps0*pi2) * 1;
      double invDist2 = q / ((r2+delta*delta)*(w_z1 - w_z0));
      return fabs(2.0*d*invDist2);
}

double electronLagrangian::getPhiSlow(double x, double y)
{
    double sum=0.0;
    for (int i=0;i<m_numParticles;i++)
    {
        //    int i=1;
        double r;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = m_bodyPos[i].x - x;
        dy = m_bodyPos[i].y - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (m_bodyPos[i].charge);

        sum-=q*log(r+delta)/(w_z1 - w_z0);


    }
    return sum;
}

