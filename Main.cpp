
#include <stdio.h>
#include <stdlib.h>

//#include  <GL/gl.h>
//#include  <GL/glu.h>
//#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/


#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <iostream>
#include <vector>
#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "multisolver.h"
#include "electronlagrangian.h"




//#include <sys/time.h>



void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

//void fmm_step(double dt);
//void sweep();
//void sweep_old();
int view=VIEW_PHI;

int redr=0;

double ck=1.0;

double cv=0.001;

double conv[100];
double conv1[100];

bool clearc=true;



std::vector <double>  avPy_;
std::vector <double>  avEy_;
std::vector <double>  area_;
std::vector <double>  T_;
FILE *file_data;


bool view_px=true;

int i_tick=0;

eFieldLagrangian* lagr_solver;
pzSolver* pz_solver;
electronLagrangian* elec_solver;
multiSolver* multi_solver;


double wall_pos=0.0;
void updateEulFields()
{
    lagr_solver->updatePhi();
    lagr_solver->solvePhi(100);

    pz_solver->get_q();
    double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);
    double phi_e0=elec_solver->getPhiSlow(w_x0,w_y0);
    //double phi_fromCharges0 = lagr_solver->getPhiFromCharges(w_x0,w_y0);

    multi_solver->updateEforPz();
    for (int i=0;i<N_X;i++)
    {
        for (int j=0;j<N_Y;j++)
        {
            double x,y;
            x=(w_x0 - (w_x1 - w_x0) * 0.05+dx*0.3*(i));
            y=(w_y0 + 0.35*(w_y1-w_y0)+dy*0.3*j);


            Py_[i][j]=pz_solver->getPhidepol(x,y) - phi_depol0 +elec_solver->getPhiSlow(x,y) - phi_e0; //-phi_fromCharges0;//-phi_depol0;

            // printf("P=%e \n",Py_[i][j]);
            phi_[i][j]=lagr_solver->getPhi(x,y)+Py_[i][j];

            vec2 Ee = elec_solver->getEe(x,y);
            vec2 Ed = lagr_solver->getE(x,y);
            vec2 Epz = pz_solver->getEdepol(x,y);

            Ex[i][j] = Ee.x + Ed.x + Epz.x;
            Ey[i][j] = Ee.y + Ed.y + Epz.y;
        }
    }
}


void display(void)
{
    double phi_max=1e-20;
    double e_max=1e-20;
    double p_max=1e-20;
    double px_max=1e-20;
    double div_max=1e-20;

    double Ec,Epc,phic,phipc,dphic,dphicp;
    for (int i=1;i<N_X-1;i++)
    {
        for (int j=1;j<N_Y-1;j++)
        {
            phi_max=fmax(phi_max,fabs(phi_[i][j]));
            p_max=fmax(p_max,fabs(Py_[i][j]));
            px_max=fmax(px_max,fabs(Px_[i][j]));
            e_max=fmax(e_max,fabs(Ey[i][j]));
            div_max=fmax(div_max,fabs(Ex[i][j]));
        }
    }

    double xc,yc,yp,ym;
    xc=(w_x0+w_x1)*0.5;
    yc=w_y0+(w_y1-w_y0)*0.25;
    yp=yc+0.5e-9;
    ym=yc-0.5e-9;

    Ec=lagr_solver->getE(xc,yc).y;
    phic=lagr_solver->getPhi(xc,yc);
    dphic=-(lagr_solver->getPhi(xc,yp)-lagr_solver->getPhi(xc,ym))/(yp-ym);


    Epc=pz_solver->getEdepol(xc,yc).y;
    phipc=pz_solver->getPhidepol(xc,yc);
    dphicp=-(pz_solver->getPhidepol(xc,yp)-pz_solver->getPhidepol(xc,ym))/(yp-ym);

    //printf("phi_max=%e Phi_p_max=%e Ey_max=%e Ep_max=%e\n",phi_max,p_max,e_max,div_max);
    //  printf("phi_c=%e Phi_p_c=%e Ey_c=%e fphic=%e Ep_c=%e ffp=%e \n",phic,phipc,Ec,dphic,Epc,dphicp);


    if (redr==1)
    {
        for (int i=0;i<10;i++)
            multi_solver->solve(1);

        // dtKoef*=1.003;

        double pzmax=0.0;
        double E_in=0.0;
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
            E_in+=pz_solver->m_p[i].E;
        }
        E_in/=pz_solver->m_p_num;

        // printf("pz_max=%e  E=%e E_in=%e \n",pzmax,E_global,E_in);
        updateEulFields();

        /*double wall_coord = pz_solver->m_p[2].r.x - pz_solver->m_p[0].r.x;
        int i;
        for( i = 2; i < pz_solver->m_p_num && pz_solver->m_p[i].p > 0; i++ )
        {
        }
        wall_coord = pz_solver->m_p[i].r.x - pz_solver->m_p[0].r.x;*/
        //   printf("e0=%e e1=%e e2=%e e3=%e wall=%e\n", pz_solver->m_p[0].p, pz_solver->m_p[1].p, pz_solver->m_p[2].p, pz_solver->m_p[3].p, wall_coord);
        //  sweep();
    }



    int i,j;//,k,l;


    double l_2;//,tx,ty,tx0,ty0,vx,vy,v0x,v0y;
    /* clear window */


    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);


    glLoadIdentity();

    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);

    glColor3f(1,1,1);

    for (i=0;i<N_X-1;i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (j=0;j<N_Y;j++)
        {
            if (view==VIEW_PHI)
                l_2= ck*(phi_[i][j]);
            if (view==VIEW_P)
                l_2=ck*(Py_[i][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f((w_x0 - (w_x1 - w_x0) * 0.05+dx*0.3*(i)),(w_y0 + 0.35*(w_y1-w_y0)+dy*0.3*j));

            if (view==VIEW_PHI)
                l_2=ck*(phi_[i+1][j]);
            if (view==VIEW_P)
                l_2=ck*(Py_[i+1][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i+1][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i+1][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i+1][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f((w_x0 - (w_x1 - w_x0) * 0.05 +dx*0.3*(i+1)),(w_y0 + 0.35*(w_y1-w_y0)+dy*0.3*j));
        }
        glEnd();
    }


    //Efields---------------------------
    /*for (i=0;i<N_X-1;i++)
    {
        glBegin(GL_LINES);
        for (j=0;j<N_Y;j++)
        {

            double x,y;
            x=1.3*(w_x0+dx*(i));
            y=1.3*(w_y0+dy*j);
            //vec2 Ee = elec_solver->getEe(x,y);
            //vec2 Ed = lagr_solver->getE(x,y);
            //vec2 Ep = pz_solver->getEdepol(x,y);
            //printf("ex=%e ey=%e \n", Ee.x, Ee.y);
            glColor3f(1,1,1);
            glVertex2f(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j));

            glColor3f(0,0,0);
            glVertex2f(1.3*(w_x0+dx*(i))+1e-16*(Ex[i][j])*ck,1.3*(w_y0+dy*j)+1e-16*(Ey[i][j])*ck);
        }
        glEnd();
    }*/
    //EoEfields

    glPointSize(3);
    glBegin(GL_LINES);

    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (lagr_solver->m_electrodes[i].canEmit)
        {
            glColor3f(1,1,1);
            double x=(lagr_solver->m_electrodes[i].r.x);
            double y=(lagr_solver->m_electrodes[i].r.y);
            glVertex3f(x,y,0.0);

            vec2 Ee = elec_solver->getEe(x,y);
            vec2 Ed = lagr_solver->getE(x,y);
            vec2 Ep = pz_solver->getEdepol(x,y);

            double E_x,E_y;
            E_x=Ee.x+Ed.x+Ep.x;
            E_y=Ee.y+Ed.y+Ep.y;

            double magn=qe/Me;//1e-1;
            double dt=multi_solver->dt_elec*100;
            E_x = magn*(E_x)*5e-25 * ck;
            E_y = magn*(E_y)*5e-25 * ck;

            glVertex3f(x + E_x,y + E_y,0.0);
        }
    }
    glEnd();

    glPointSize(3);
    glBegin(GL_POINTS);

    double rhomax=0.0;
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (fabs(lagr_solver->m_electrodes[i].phi_ext)>rhomax)
            rhomax=lagr_solver->m_electrodes[i].phi_ext;
    }


    //printf("phimax=%e \n",rhomax);
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        double c=ck*lagr_solver->m_electrodes[i].phi_ext/rhomax;
        //printf("i=%d phi=%e  phi_fix=%e\n",i,lagr_solver->m_electrodes[i].phi_ext,lagr_solver->m_electrodes[i].phi_fix);
        //glColor3f(c,c,-c);
        if (lagr_solver->m_electrodes[i].canEmit)
            glColor3f(0.25,0.25,1);
        else
            glColor3f(1,1,1);
        if(i==0)
          glColor3f(1,0,0);
        double x=(lagr_solver->m_electrodes[i].r.x);
        double y=(lagr_solver->m_electrodes[i].r.y);
        glVertex3f(x,y,0.0);
    }


    for (int i=0;i<lagr_solver->m_chargeNum;i++)
    {
        glColor3f(0,0.5,0);
        double x=(lagr_solver->m_charges[i].x);
        double y=(lagr_solver->m_charges[i].y);
        glVertex3f(x,y,0.0);
    }
    glEnd();

    /* glLineWidth(2);
    glBegin(GL_LINES);

    for (int i = 0; i < lagr_solver->m_elec_num;i++)
    {
        glColor3f(1,0,0);
        double x = (lagr_solver->m_electrodes[i].r.x);
        double y = (lagr_solver->m_electrodes[i].r.y);
        glVertex3f(x,y,0.0);
        //glVertex3f(x + 50e-6, y + 50e-6,0.0);
        glVertex3f(x + ck * 1e-6 * lagr_solver->m_electrodes[i].eToEmit * lagr_solver->m_electrodes[i].nx ,
                   y + ck * 1e-6 *lagr_solver->m_electrodes[i].eToEmit * lagr_solver->m_electrodes[i].ny,0.0);
    }
    glEnd();*/

    glPointSize(5.5);
    glBegin(GL_POINTS);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        glColor3f(1.0,0.0,1.0);
        glVertex2f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y);
    }
    glEnd();


    /*glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glBegin(GL_QUADS);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        int p_n=(int) ((elec_solver->m_bodyPos[i].x-pz_solver->m_p[0].r.x+0.5 * pz_solver->m_dx)/pz_solver->m_dx);

        glColor3f(1,1,1);
        glVertex2f(pz_solver->m_p[p_n].r.x-pz_solver->m_dx*0.5,pz_solver->m_p[p_n].r.y -pz_solver->m_p[p_n].dl*0.5);
        glVertex2f(pz_solver->m_p[p_n].r.x-pz_solver->m_dx*0.5,pz_solver->m_p[p_n].r.y +pz_solver->m_p[p_n].dl*0.5);


        glVertex2f(pz_solver->m_p[p_n].r.x+pz_solver->m_dx*0.5,pz_solver->m_p[p_n].r.y +pz_solver->m_p[p_n].dl*0.5);
        glVertex2f(pz_solver->m_p[p_n].r.x+pz_solver->m_dx*0.5,pz_solver->m_p[p_n].r.y -pz_solver->m_p[p_n].dl*0.5);


    }
    glEnd();
*/


    /*double vel_scale=1900000.0;//sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.03*(w_x1-w_x0);
    glLineWidth(4);
    glBegin(GL_LINES);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1,0,0);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex2f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y);
        glColor3f(0.0,0.0,0.0);
        glVertex2f(elec_solver->m_bodyPos[i].x+leng_sacle*elec_solver->m_bodyVel[i].x/vel_scale,
                   elec_solver->m_bodyPos[i].y+leng_sacle*elec_solver->m_bodyVel[i].y/vel_scale);
    }
    glEnd();*/


    /*   glColor3f(0.5,0.5,0.5);

    glBegin(GL_LINE_LOOP);

    glVertex3f(w_x0,w_y0,w_z0);
    glVertex3f(w_x1,w_y0,w_z0);
    glVertex3f(w_x1,w_y1,w_z0);
    glVertex3f(w_x0,w_y1,w_z0);
    glEnd();

    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);

    glVertex3f(w_x0,w_y0,w_z1);
    glVertex3f(w_x1,w_y0,w_z1);
    glVertex3f(w_x1,w_y1,w_z1);
    glVertex3f(w_x0,w_y1,w_z1);
    glEnd();
*/

    if (view_px)
    {
        //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        glBegin(GL_TRIANGLE_STRIP);
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            glColor3f(ck * pz_solver->m_p[i].p/0.26,-ck * pz_solver->m_p[i].p/0.26,0);
            glVertex2f(pz_solver->m_p[i].r.x/*-pz_solver->m_dx*0.5*/,pz_solver->m_p[i].r.y -pz_solver->m_p[i].dl*0.5);
            glVertex2f(pz_solver->m_p[i].r.x/*-pz_solver->m_dx*0.5*/,pz_solver->m_p[i].r.y +pz_solver->m_p[i].dl*0.5);
        }
        glEnd();
    }


    glBegin(GL_POINTS);
    glColor3f(1,0.5,0);
    glVertex2f(xc, yc);


    glEnd();

    glPointSize(3.0);

    glLineWidth(1.0);


    glDisable(GL_DEPTH_TEST);
    glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(0.5,0.5,1.0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-3 * scale * (pz_solver->m_p[i].q_ext+pz_solver->m_p[i].q) * (w_y1 - w_y0)-5e-6);
    }
    glEnd();
    glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(0.1,0.5,0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-3 * scale * (pz_solver->m_p[i].q_ext) * (w_y1 - w_y0)-5e-6);
    }
    glEnd();


    /* double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);
    double phi_e0=elec_solver->getPhiSlow(w_x0,w_y0);
    double phi_electr0=lagr_solver->getPhi(w_x0,w_y0);*/

    /*glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(1.0,0.0,1.0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-10 * scale * pz_solver->m_p[i].E * (w_y1 - w_y0));
    }
    glEnd();*/

    /*double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);
    double phi_e0=elec_solver->getPhiSlow(w_x0,w_y0);
    double phi_electr0=lagr_solver->getPhi(w_x0,w_y0);

    glLineWidth(3.);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {


        double phi = pz_solver->getPhidepol(pz_solver->m_p[i].r.x,pz_solver->m_p[i].r.y) - phi_depol0
                +elec_solver->getPhiSlow(pz_solver->m_p[i].r.x,pz_solver->m_p[i].r.y) - phi_e0;
                //+lagr_solver->getPhi(pz_solver->m_p[i].r.x,pz_solver->m_p[i].r.y) - phi_electr0;
        glColor3f(1.0,1.0,1.0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-1 *  scale * phi * (w_y1 - w_y0));
    }
    glEnd();*/

    /*glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,1.0,1.0);
    glVertex2f(pz_solver->m_p[0].r.x, 0);
    glVertex2f(pz_solver->m_p[pz_solver->m_p_num - 1].r.x, 0);
    glEnd();*/

    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();

}


/////////1d emission below
///
///
vec2 ch[1000];
int c_num=0;
double velFermi=1.5e6;
double E1=1.7e7;
double dt1=5e-14;
double calcJ(double Ein)
{

    double E=Ein/100;//from V/m to V/cm
    double t2 = 1.1;
    double B = 145.0;//145.0
    double phi = 4.0;
    double y = 3.79 * 1e-4 * sqrt(fabs(B * E)) / phi;
    double tetta = 0.95 - 1.03 * y * y;
    return 1e4*(1.54 * 1.0e-6 * B * B * E * E / (t2  * phi)) * exp ( - 6.83 * 1.0e7 * pow(phi, 1.5) * tetta / fabs( B * E)); //in A/m^2
}

double E_dipole(double x) //calculate E in the middle two elementary charges with distance d
{
    double sum =0.0;

    for (int i=0;i<c_num;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = ch[i].x - x;

        r2=dx*dx;
        q=-qe/(eps0*pi2) * (ch[i].charge);

        double c=q/((r2+delta*delta)*(w_z1 - w_z0));

        sum+=c*dx;


    }
    return sum;
}

void step_charges(double dt)
{
    for (int i=0;i<c_num;i++)
    {
        ch[i].x+=dt*velFermi;
    }

    for (int i=0;i<c_num;i++)
    {
        if (ch[i].x>1e-6){
            ch[i]=ch[c_num-1];
            c_num--;
        }
    }

}

double Ec=0.0;
void sim_emit(double E_1,double dt)
{
    double ds = (w_y1-w_y0)*(w_z1-w_z0);

    double phi=4; //work function
    double barrier_width=(phi-sqrt(fabs(qe*E_1/(pi4*eps0))))/E_1;

    //
    double vel=1.5e6; //fermi velocity
    double E0=2.0*fabs(E_dipole(0.0));

    Ec=0.999*Ec+ 0.001*fmax(E_1-E0,0);
    double el_to_add = calcJ(fmax(Ec,0))*dt*ds/(fabs(qe));

    vec2 q(barrier_width,0,el_to_add);

    ch[c_num]=q;
    c_num++;

    /*for (int i=0;i<c_num;i++)
     {ch[i].charge=el_to_add;

     }*/

    step_charges(dt);

    //  printf("barrier_w=%e  q=%e \n",barrier_width,el_to_add);
}

void save_file()
{

    emis_tab.ex=1.1;
    emis_tab.ey=1.1;

    emis_tab.nx=100;
    emis_tab.ny=10;

    emis_tab.x0=1e6;
    emis_tab.y0=1e-14;



    FILE* f=fopen("n.txt","w");
    fprintf (f,"0 ");
    for (int j=0;j<emis_tab.ny;j++)
    {
        fprintf(f,"%e ",1e-14*pow(1.1,j));

        emis_tab.y[j]=1e-14*pow(1.1,j);
    }
    fprintf(f,"\n");

    for (int i=0;i<emis_tab.nx;i++)
    {
        emis_tab.x[i]=1.0e6*pow(1.1,i);
    }


    double dt,E;
    for (int i=0;i<100;i++)
    {
        E=1e6*pow(1.1,i);
        fprintf(f,"%e ",E);
        for (int j=0;j<10;j++)
        {
            dt=1e-14*pow(1.1,j);

            for (int  k=0;k<4000;k++)
            {
                sim_emit(E,dt);
            }
            fprintf(f,"%e ",ch[0].charge);
            double sumCharge = 0;
            double sumX = 0;
            for (int k=0;k<c_num;k++)
            {
                sumCharge+=ch[k].charge;
                sumX+=ch[k].x;
            }
            sumX/=c_num;
            emis_tab.f[i][j]=sumCharge;//ch[0].charge;
            printf("charge = %f x = %e\n", sumCharge, sumX);
        }
        fprintf(f,"\n ");
        for (int k=0;k<4000;k++)
        {
            sim_emit(E,1e-14);
        }
        printf("i=%d \n ",i);

    }
    fclose(f);
}
int jc=0;
void display_1d(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    if (redr==1)
    {
        sim_emit(E1,dt1);
        //printf("c_num=%d \n", c_num);
    }

    glLoadIdentity();

    glPointSize(3);
    glColor3f(1,1,1);
    /*glBegin(GL_LINES);
    glVertex3f(w_x0,0.0,0.0);
    glVertex3f(w_x1,0.0,0.0);

    for (int i=0;i<100;i++)
    {
        glVertex3f(w_x0+i*1.0e-7,3.0e-8,0.0);
        glVertex3f(w_x0+i*1.0e-7,-3.0e-8,0.0);
    }
    glEnd();

    glBegin(GL_LINES);

    for (int i=0;i<c_num;i++)
    {
        glVertex3f(ch[i].x,6.0e-8+6.0e-8*0.001*ch[i].charge,0.0);
        glVertex3f(ch[i].x,6.0e-8,0.0);
    }
    glEnd();
*/
    /*glColor3f(1,0,0);
    glPointSize(3);
    glBegin(GL_POINTS);

    for (int i=0;i<c_num;i++)
    {

        glVertex3f(ch[i].x,6.0e-8,0.0);
    }
    glEnd();*/

    glColor3f(1,0,0);
    glBegin(GL_LINE_STRIP);
    for (int i=0;i<emis_tab.nx;i++)
    {
        //for (int j=0;j<emis_tab.ny;j++)
        {
            glVertex3f(w_x0+scale*emis_tab.x[i]/emis_tab.x0*1.0e-7,ck*emis_tab.f[i][jc],0.0);

        }
    }
    glEnd();

    glColor3f(0,1,1);
    glBegin(GL_POINTS);
    for (int i=0;i<10000;i++)
    {
        //for (int j=0;j<emis_tab.ny;j++)
        {
            glVertex3f(w_x0+scale*(i * 1e5 + 1e5)/emis_tab.x0*1.0e-7,ck*emis_tab.get_f(i * 1e5 + 1e5, emis_tab.y0),0.0);

        }
    }
    glEnd();

    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();
}

void m_m(int x,int y) //mouse move
{

    if (rotate==1)
    {
        rx=rx0+0.5*(x-mx0);
        ry=ry0+0.5*(y-my0);
    }
    glutPostRedisplay();
}



void m_d(int button, int state,int x, int y)  //mouse down
{

    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;

    }

    mouse_x=(1.0*x)/W_WIDTH;

    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;

    glutPostRedisplay();
}

void saveInFile()
{
    if(g_save_time2 * 1e15 >200)
    {
        g_save_time2=0;
        double q_sum = 0;
        double p_pos=0.0;
        double p_full=0.0;
        double wall_coord = pz_solver->m_p[1].r.x - pz_solver->m_p[0].r.x;
        int i;
        for( i = 1; i < pz_solver->m_p_num && pz_solver->m_p[i].p > 0; i++ )
        {
        }
        wall_coord = pz_solver->m_p[i].r.x - pz_solver->m_p[0].r.x;

        for( i = 0; i < pz_solver->m_p_num; i++ )
        {
            q_sum += pz_solver->m_p[i].q_ext;
            if (pz_solver->m_p[i].p>0) p_pos += pz_solver->m_p[i].p;
            p_full+=fabs(pz_solver->m_p[i].p);
        }
        printf("wall = %e\n", wall_coord);
        fprintf(file_data,"%lf %e %lf %lf\n",g_t * 1e15, wall_coord, q_sum,p_pos/p_full);
    }

    int time = int(g_t * 1e15);
    if(g_save_time * 1e15 > 2000)
    {
        g_save_time=0;
        char filename[64];
        sprintf(filename, "prifiles%d_%d.txt", int(g_phi), time);
        FILE *file_data_=fopen(filename,"w");
        for( int i = 0; i < pz_solver->m_p_num; i++ )
        {
            fprintf(file_data_,"%e  %lf  %lf \n",pz_solver->m_p[i].r.x - pz_solver->m_p[0].r.x, pz_solver->m_p[i].q_ext, pz_solver->m_p[i].q);
        }
        fclose(file_data_);
    }
}

double angle=0.0;
void kb(unsigned char key, int x, int y)
{
    int i,j,k,nn,n;
    double m,sum;
    double max_err=0.0;
    if (key=='.')
    {
        scale*=1.1;

        printf("scale=%e \n", scale);
    }

    if (key==',')
    {
        scale/=1.1;

        printf("scale=%e \n", scale);
    }
    if (key==']')
    {
        ck*=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='[')
    {
        ck/=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='1')
    {
        view=VIEW_PHI;
        printf("viewing PHI \n");
    }


    if (key=='2')
    {
        view=VIEW_P;
        printf("viewing phi_depol \n");
    }

    if (key=='3')
    {
        view=VIEW_E;
        printf("viewing E\n");
    }

    if (key=='4')
    {
        view=VIEW_DIV;

        printf("viewing E_depol\n");
    }


    if (key=='m')
    {
        move_particles=!move_particles;
    }

    if (key=='5')
    {
        // clearc=!clearc;
        view_px=!view_px;
    }


    if (key=='8')
    {
        //    wall_pos+=0.01;
        //    pz_solver->setWallPos(wall_pos);
        //  E_global=-E_global;
        //  printf("E_global=%e \n",E_global);

        g_emitElectrons=!g_emitElectrons;
    }

    if (key=='0')
    {
        //    wall_pos+=0.01;
        //    pz_solver->setWallPos(wall_pos);
        dtKoef*=1.1;
        printf("dtKoef=%e \n",dtKoef);
    }

    if (key=='9')
    {
        //  wall_pos-=0.01;
        // pz_solver->setWallPos(wall_pos);
        dtKoef/=1.1;
        printf("dtKoef=%e \n",dtKoef);

    }

    if (key=='=')
    {
        double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);
        printf("here i am ppp=%e \n",phi_depol0);
        for (int i=0;i<lagr_solver->m_elec_num;i++)
        {
            double x,y;
            x=lagr_solver->m_electrodes[i].r.x;
            y=lagr_solver->m_electrodes[i].r.y;

            lagr_solver->m_electrodes[i].phi_fix_charges=0;//(pz_solver->getPhidepol(x,y)-phi_depol0);
        }
        printf("start ls \n");
        lagr_solver->solve_ls_fast_PhiE();//solvePhi(20);
        printf("end ls \n");
        updateEulFields();
    }

    if (key=='7')
    {
        multi_solver->checkPotential();
    }

    if (key=='-')
    {

        updateEulFields();
    }


    if (key=='d')
    {
        //  pz_solver->solvePz(100);
        multi_solver->solve(2);
        double pzmax=0.0;
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
        }
        printf(" pzmax =%e \n",pzmax);
    }

    if (key=='g')
    {
        /* for (int i = 0;i<20;i+=2)
        {
            g_t = 0;
            g_phi = (i+1);
            char filename[64];
            sprintf(filename, "output%d.txt", i+1);
            file_data=fopen(filename,"w");
            multi_solver->init();
            while (g_t * 1e15 < (g_phi / 10 + 1) * 2e4)
            {
                multi_solver->solve(1);
                if(int(g_t * 1e15) % 100 == 0)
                {
                    double pzmax=0.0;
                    double E_in=0.0;
                    for (int i=0;i<pz_solver->m_p_num;i++)
                    {
                        if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
                        E_in+=pz_solver->m_p[i].E;
                    }
                    E_in/=pz_solver->m_p_num;
                    updateEulFields();
                    glutPostRedisplay();
                    display();
                }
                saveInFile();
            }
            fclose(file_data);
        }*/
    }

    if (key==' ')
    {
        redr=!redr;
    }

    glutPostRedisplay();
}

void kb_1d(unsigned char key, int x, int y)
{


    /*  if (key==' ')
    {
        redr=!redr;
    }


    if (key=='.')
    {
        dt1*=1.1;
        printf("dt=%e \n",dt1);
    }

    if (key==',')
    {
        dt1/=1.1;
        printf("dt=%e \n",dt1);
    }


    if (key==']')
    {
        E1*=1.1;
        printf("E=%e \n",E1);
    }

    if (key=='[')
    {
        E1/=1.1;
        printf("E=%e \n",E1);
    }*/


    if (key=='.')
    {
        jc++;
        printf("jc=%d \n", jc);
    }

    if (key==',')
    {
        jc--;
        printf("jc=%d \n", jc);
    }


    if (key==']')
    {
        ck*=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='[')
    {
        ck/=1.1;
        printf("ck=%e \n", ck);
    }


    if (key=='0')
    {
        scale*=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='9')
    {
        scale/=1.1;
        printf("ck=%e \n", ck);
    }

    glutPostRedisplay();
}





void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(w_x0-(w_x1-w_x0)*0.3, w_x1*1.3, w_y0*1.3,w_y1*1.3, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);

    lagr_solver= new eFieldLagrangian();
    lagr_solver->updatePhi();
    pz_solver= new pzSolver();

    elec_solver = new electronLagrangian();

    multi_solver = new multiSolver();
    multi_solver->m_Esolver = lagr_solver;
    multi_solver->m_pzSolver = pz_solver;
    multi_solver->m_elecSolver = elec_solver;


    static double qem[100][10]; //E 10^5 - 10^9 ; dt  10^-15 10^-8
    double ds  = lagr_solver->m_electrodes[0].dl*(w_z1-w_z0);

    //for (int i=0; i<100; i++)
    {
        int i=4;
        double E1 = 1.e6 + 1.0e7*i;
        printf("\n E=%e \n ",E1);

        for (int j=0;j<2;j++)
        {
            double d_t = pow(10.0,-8.0-j);
            double E0=elec_solver->getEmult_dipole(2.0e-6);
            double el_to_add = 0.99*E1/(E0);//elec_solver->calcJ(E1)*d_t*ds/(fabs(qe)/**num_in_pack*/);
            //el_to_add = elec_solver->calcJ(E1)*d_t*ds/(fabs(qe)/**num_in_pack*/);
            double el_to_add0=el_to_add;
            printf("__E1=%e E0=%e n=%e;  \n",E1,E0*el_to_add);
            //double J=elec_solver->calcJ(E1-E0*el_to_add)*d_t*ds/(fabs(qe));
            for (int k = 0; k<6000;k++) {
                el_to_add = el_to_add * 0.9995 + 0.0005 * elec_solver->calcJ(E1-E0*el_to_add)*d_t*ds/(fabs(qe));

                // printf("E1=%e E0=%e n=%e; \n",E1,E0*el_to_add,el_to_add);
            }
            printf(" dt=%e E1=%e E0=%e n=%e; \n",d_t,E1,E0*el_to_add,el_to_add);
            //    printf("dt=%e e=%e; ",d_t,el_to_add);
            //         printf("%e %e %e; ",el_to_add0,el_to_add,el_to_add*E0/E1);
        }
    }

    //lagr_solver->solvePhi(10);
}

void resize(int w, int h)
{
    glViewport(0, 0, w, h);
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    double h0,w0;
    w0=W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0);//w_x1*1.3-w_x0+(w_x1-w_x0)*0.3;
    h0=W_HEIGHT;//(w_y1-w_y0)*1.3;
    glOrtho(w_x0-(w_x1-w_x0)*0.11, w_x0-(w_x1-w_x0)*0.11+(w_x1*1.3-(w_x0-(w_x1-w_x0)*0.3))*w*1.0/w0*(h0*1.0/h), w_y0*1.3,w_y1*1.3, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);

}


int main(int argc, char** argv)
{
    //srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();

    save_file();
    glutMainLoop();
}
