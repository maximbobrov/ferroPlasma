
#include <stdio.h>
#include <stdlib.h>
#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#endif
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <iostream>
#include <vector>
#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "multisolver.h"



//#include <sys/time.h>

void save_fields();
void load_fields();
void display(void);
void sweep_init();
void init();
void saveStatePZ(char* fname);
int view=VIEW_PHI;


int redr=0;

double ck=1.0;

double cv=0.001;

double conv[100];
double conv1[100];

double scale_f=1.0;
bool clearc=true;


double g_x0=0;
double g_x1=1;
double g_y0=0;
double g_y1=1;

void saveInFile();

std::vector <double>  avPy_;
std::vector <double>  avEy_;
std::vector <double>  area_;
std::vector <double>  T_;
FILE *file_data;


bool view_px=true;

int i_tick=0;

eFieldLagrangian* lagr_solver;
pzSolver* pz_solver;

multiSolver* multi_solver;


double wall_pos=0.0;


static double progress = 2.0;
static int fileNum = 0;
static double timePrev = 0;


void draw_fields_2d();
void updateEulFields();
void draw_traj();
void draw_fields_pz_1d();
 void draw_fields_pz_1d_highres();
void draw_charges_pz_1d();
void draw_pz();

void draw_electrode();

void display(void)
{

    glDisable(GL_DEPTH_TEST);

    if (redr==1)
    {
        // double t0 = get_time();
        for (int j=0;j<10;j++)
        {
        multi_solver->updateTrajTable(false,1e-15,0.5e-7);
        //  double t1 = get_time();
        for (int i=0;i<5;i++)
        {
       //      multi_solver->updateTrajTable();
            multi_solver->solve(10);
        }
        //    double t2 = get_time();
        }
        //updateEulFields();
        //  double t3 =get_time();
//redr=0;
        // printf("traj_upd_time=%e solve_time=%e eul_fields_time=%e \n", t1-t0, t2-t1,t3-t2);
    }
    double wall_coord;
    if(serialRegime)
    {
        if(progress>1.0)
        {

            if(file_data)
                fclose(file_data);
            g_emitElectrons = false;
            g_t = 0;
            g_phi = -(100 + fileNum*50);
            g_phi_max = g_phi;
            fileNum++;
            char filename[64];
            sprintf(filename, "output%d.txt", fileNum);
            file_data=fopen(filename,"w");
            g_i_wall=0;
            multi_solver->init();
            progress=0;
            g_q_enable = true;
            for (int kk=0;kk<300;kk++)
                multi_solver->solve(10);

            g_phi_max =0.0;//*= -1;

            for (int kk=0;kk<300;kk++)
                multi_solver->solve(10);


            g_phi_max = (100 + fileNum*50);
            g_i_wall=g_i_wall_edge;
            g_i_wall_tmp=g_i_wall;
            for (int i=1; i<pz_solver->m_p_num;i++)
            {
                if (i<g_i_wall){
                    pz_solver->m_p[i].q_ext=0;
                    pz_solver->m_p[i].p = 0.3;
                    pz_solver->m_p[i].p_prev = 0.3;
                    pz_solver->m_p[i].q = 0;
                }
            }
            for (int kk=0;kk<100;kk++)
                multi_solver->solve(10);
            g_q_enable = false;
            g_emitElectrons = true;
            g_t=0;
            g_save_time=0;
            g_save_time2=0;
        }
        else{
            for (int kk=0;kk<10;kk++){

                multi_solver->updateTrajTable(false,1e-15,0.5e-7);
                for (int i=0;i<10;i++)
                    multi_solver->solve(10);
            }
            //for (int kk=0;kk<100;kk++)
            //    multi_solver->solve(1);
            //double pzmax=0.0;
            //double E_in=0.0;
            /*for (int i=0;i<pz_solver->m_p_num;i++)
            {
                if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
                //E_in+=pz_solver->m_p[i].E;
            }*/
            printf("progress(%d) = %d \n", fileNum, int(progress * 100));
            //E_in/=pz_solver->m_p_num;
            updateEulFields();
            saveInFile();
            wall_coord = pz_solver->m_p[1].r_top.x - pz_solver->m_p[0].r_top.x;
            int i;
            for( i = g_i_wall_edge; i < pz_solver->m_p_num && pz_solver->m_p[i].p > 0; i++ )
            {
            }
            wall_coord = pz_solver->m_p[i].r_top.x - pz_solver->m_p[0].r_top.x;

        }
        double prev=progress*10;
        progress +=fabs(g_phi_max/175.0)/2000.0 ;//(g_t / multi_solver->dt_elec)/((g_phi / 10 + 1) * 5000.0);
        int i_=(int) (progress*10.0);
        if ((progress*10.0-i_>0)&&(prev-i_<0))
        {
            char nme[1000];
            sprintf(nme, "output_%dV_%d(%e)%d.state",int(g_phi_max), fileNum,g_t,int(progress*100));
            saveStatePZ(nme);
        }


    }

    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);


    glLoadIdentity();

    glTranslatef(w_x0,0,0);
    glScalef(scale_f,scale_f,scale_f);
    glTranslatef(-w_x0,0,0);
    //glRotatef(ry,1.0,0,0);
    //glRotatef(rx,0.0,1.0,0);

    glColor3f(1,1,1);


    if (view_px)
    {
        draw_pz();
        draw_traj();
        draw_fields_pz_1d();
    }else
    {
       draw_fields_2d();
       draw_fields_pz_1d_highres();
    }

    draw_electrode();

//    draw_charges_pz_1d();





    glutSwapBuffers();
    if (redr==1 || serialRegime) glutPostRedisplay();
}



void updateEulFields()
{  
    pz_solver->get_q();
    for (int i=0;i<N_X;i++)
    {
        for (int j=0;j<N_Y;j++)
        {
            double x,y;
            x=g_x0+ i*(g_x1-g_x0)/(N_X-1);
            y=g_y0+ j*(g_y1-g_y0)/(N_Y-1);

            if ((view==VIEW_EX) || (view==VIEW_EY))
            {
                vec2 E_ = multi_solver->get_slower_E(x,y);
                Ex[i][j] = E_.x;
                Ey[i][j] = E_.y;

            }
            if (view==VIEW_PHI)
            {
                phi_[i][j]=multi_solver->get_slow_phi(x,y,false);
            }

            if (view==VIEW_P)
            {
                phi_[i][j]=multi_solver->get_slow_phi(x,y,true);
            }
        }
    }
}

void draw_fields_2d()
{
    updateEulFields();


double phi0=fabs(g_phi);
double E0=10.0*phi0/(dl_pz);

    //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    for (int i=0;i<N_X-1;i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j=0;j<N_Y;j++)
        {
            double l_2;
            double x=g_x0+ i*(g_x1-g_x0)/(N_X-1);
            double y=g_y0+ j*(g_y1-g_y0)/(N_Y-1);

            if (view==VIEW_PHI)
                l_2 = ck*phi_[i][j]/phi0;
            if (view==VIEW_P)
                l_2 = ck*phi_[i][j]/phi0;
            if (view==VIEW_EX)
                l_2=ck*(Ey[i][j])/E0;
            if (view==VIEW_EY)
                l_2=ck*(Ex[i][j])/E0;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(x,y);

            if (view==VIEW_PHI)
                l_2=ck*phi_[i+1][j]/phi0;
            if (view==VIEW_P)
                l_2 = ck*phi_[i+1][j]/phi0;
            if (view==VIEW_EX)
                l_2=ck*(Ey[i+1][j])/E0;
            if (view==VIEW_EY)
                l_2=ck*(Ex[i+1][j])/E0;

            glColor3f(l_2,l_2,-l_2);
            x=g_x0+ (i+1)*(g_x1-g_x0)/(N_X-1);
            y=g_y0+ j*(g_y1-g_y0)/(N_Y-1);

            glVertex2f(x,y);
        }
        glEnd();
    }
}



void draw_traj()
{
   // printf("traj_num=%d \n",multiSolver::traj_num);
    for (int i=0;i<multiSolver::traj_num;i++)
    {
     //   printf("i=%d j=%d \n", multiSolver::traj_num,multiSolver::trajectories_num[i]);
        vec2 r(0,0,0);
        glLineWidth(1.5);
        double cur=10.0*multiSolver::trajectories[i][0].charge;
        glColor3f(cur,0,1.0-cur);
        glBegin(GL_LINE_STRIP);
        for (int j=0;j<multiSolver::trajectories_num[i];j++)
        {
            r=multiSolver::trajectories[i][j];

            glVertex2f(r.x,r.y);
        }
        glEnd();
        glColor3f(1,0,0);
        glBegin(GL_POINTS);
        glPointSize(10);
            glVertex2f(r.x,r.y);
        glEnd();
    }
}

void draw_fields_pz_1d()
{
    static vec2 Ef[1000];
    static double  phi_f[1000],phi_f2[1000];
    //double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0,is2D);

    for (int i=0;i< pz_solver->m_p_num;i++)
    {
        vec2 r;
        r.x= pz_solver->m_p[i].r_top.x;
        r.y=0;
        phi_f[i]=pz_solver->m_p[i].E*pz_solver->m_p[i].dl-g_phi_max*2;
    }

    glLineWidth(2);
    glBegin(GL_LINE_STRIP);

    for(int i=0; i < pz_solver->m_p_num; i++ )
    {

        glColor3f(0,1,1);
        glVertex2f(pz_solver->m_p[i].r_top.x, pz_solver->m_p[i].r_top.y + scale * (1.0/g_phi_max)*1000.0*(-phi_f[i])* (w_y1 - w_y0) );
    }
    glEnd();
}


 void draw_fields_pz_1d_highres()
 {

     glLineWidth(2);
     glBegin(GL_LINE_STRIP);
     for (int i=0;i< 1000;i++)
     {
         vec2 r;
         r.x= pz_solver->m_p[g_i_wall_edge-2].r_top.x+ pz_solver->m_dx*0.03*i;
         r.y=pz_solver->m_p[g_i_wall_edge-2].r_top.y;
         double phi_f=multi_solver->get_slow_phi(r.x,r.y,is2D)+g_phi_max;

         glColor3f(0,1,1);
         glVertex2f(r.x, r.y + scale * (1.0/g_phi_max)*1000.0*(phi_f)* (w_y1 - w_y0) );
     }
     glEnd();
 }
void draw_charges_pz_1d()
{
    glLineWidth(2.5);
    glBegin(GL_LINE_STRIP);


    for(int i=0; i < pz_solver->m_p_num; i++ )
    {

        if (pz_solver->m_p[i].q_ext + pz_solver->m_p[i].q>0)
            glColor3f(0.9,0.4,0);
        else
            glColor3f(0.0,0.4,0.9);

        //glVertex2f(pz_solver->m_p[i].r_top.x, 100e-1 * scale * (pz_solver->m_p[i].q_ext + pz_solver->m_p[i].q) * (w_y1 - w_y0)-5e-6);

    }
    glEnd();

    glLineWidth(2.5);
    glBegin(GL_LINE_STRIP);

    double q_extr=0.0;
    for(int i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(0.0,0.0,1);
        glVertex2f(pz_solver->m_p[i].r_top.x, 1e-1 * scale * (pz_solver->m_p[i].q_ext) * (w_y1 - w_y0)-5e-6);
        q_extr+=pz_solver->m_p[i].q_ext;
    }
    printf(" \n qext=%e phi0=%f \n",pz_solver->m_p[0].q, pz_solver->getPhi2D(2.9e-6,0.0, 0.0,2e-6,0.0,100.0/((qe/(eps0*pi2))/(w_z1 - w_z0))));
    glEnd();
}

void draw_pz()
{
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    for (int i=0;i<g_i_wall_tmp;i++)
    {

        glColor3f(ck * pz_solver->m_p[i].p/0.26,-ck * pz_solver->m_p[i].p/0.26,0);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);

        glColor3f(ck * pz_solver->m_p[i].p/0.26,-ck * pz_solver->m_p[i].p/0.26,0);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);


    }
    glEnd();



    glBegin(GL_QUADS);
    int i=g_i_wall_tmp;
    glColor3f(1,0.5,0.5);
    double rig=(pz_solver->m_p[i].p+0.26)/0.52;
    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);
    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);


    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5 + pz_solver->m_dx*rig , pz_solver->m_p[i].r_top.y);
    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5+  pz_solver->m_dx*rig  , pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);

    glColor3f(0.5,1,0.5);
    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5+  pz_solver->m_dx*rig  , pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);
    glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5 + pz_solver->m_dx*rig , pz_solver->m_p[i].r_top.y);
    glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);
    glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);

    glEnd();

    glBegin(GL_QUADS);
    for (int i=g_i_wall_tmp+1;i<pz_solver->m_p_num;i++)
    {

        glColor3f(ck * pz_solver->m_p[i].p/0.26,-ck * pz_solver->m_p[i].p/0.26,0);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);

        glColor3f(ck * pz_solver->m_p[i].p/0.26,-ck * pz_solver->m_p[i].p/0.26,0);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);


    }
    glEnd();





    glBegin(GL_LINES);
    glLineWidth(1.0);
    for (int i=0;i<pz_solver->m_p_num;i++)
    {
        glColor3f(1,1,1);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r_top.y);
    }
    glEnd();



    glDisable(GL_LINE_SMOOTH);
     glLineWidth(7);
    glBegin(GL_LINES);

    for (int i=0;i<pz_solver->m_p_num;i++)
    {
        double q0=(0.3)*pz_solver->m_p[i].ds/qe;
        double qp=pz_solver->m_p[i].q/q0;
        double q_ext=pz_solver->m_p[i].q_ext/q0;

        if (qp>0)
        glColor3f(0.5,0,0);
        else
        glColor3f(0.0,0,0.5);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.05,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1);
        glVertex2f(pz_solver->m_p[i].r_top.x-pz_solver->m_dx*0.05,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1+
                                                                    pz_solver->m_p[i].dl*0.1*fabs(qp));

        if (q_ext>0)
        glColor3f(1.0,0.5,0);
        else
        glColor3f(0.0,0.5,1.0);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.05,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.05,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1+
                                                                    pz_solver->m_p[i].dl*0.1*fabs(q_ext));


        if (q_ext+qp>0)
        glColor3f(1.0,0.0,0);
        else
        glColor3f(0.0,0.0,1.0);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.15,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1);
        glVertex2f(pz_solver->m_p[i].r_top.x+pz_solver->m_dx*0.15,pz_solver->m_p[i].r_top.y -pz_solver->m_p[i].dl*0.1+
                                                                    pz_solver->m_p[i].dl*(q_ext+qp));

    }
    glEnd();




    glLineWidth(1.0);
}



void draw_electrode()
{
    glPointSize(3);
    glLineWidth(1.0);
    glBegin(GL_LINES);

    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (lagr_solver->m_electrodes[i].canEmit)
        {
            glColor3f(1,1,1);
            double x=(lagr_solver->m_electrodes[i].r.x);
            double y=(lagr_solver->m_electrodes[i].r.y);
            glVertex3f(x,y,0.0);

            vec2 E=multi_solver->get_slower_E(x,y);

            /*double ex=E.x;
            double ey=E.y;

            vec2 Ed = lagr_solver->getE(x,y);
            vec2 Ep = pz_solver->getEdepol(x,y);*/

            double E_x,E_y;
            E_x=E.x;//Ed.x+Ep.x;
            E_y=E.y;//Ed.y+Ep.y;

            double magn=qe/Me;//1e-1;
            E_x = magn*(E_x)*5e-25 * ck;
            E_y = magn*(E_y)*5e-25 * ck;

            double nx=lagr_solver->m_electrodes[i].nx;
            double ny=lagr_solver->m_electrodes[i].ny;
            double a=nx*E_x+ny*E_y;

            glVertex3f(x + E_x,y + E_y,0.0);
        }
    }
    /*for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (lagr_solver->m_electrodes[i].canEmit)
        {
            int ind=multi_solver->endPosTable[i];

            glColor3f(1,1,1);
            double x=(lagr_solver->m_electrodes[i].r.x);
            double y=(lagr_solver->m_electrodes[i].r.y);
            glVertex3f(x,y,0.0);

            x=pz_solver->m_p[ind].r_top.x;
            y=pz_solver->m_p[ind].r_top.y-pz_solver->m_p[ind].dl*0.05;
            glVertex3f(x,y,0.0);
        }
    }*/


    glEnd();

    glPointSize(2.5);
    glBegin(GL_POINTS);

    double rhomax=0.0;
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (fabs(lagr_solver->m_electrodes[i].phi_ext)>rhomax)
            rhomax=lagr_solver->m_electrodes[i].phi_ext;
    }

    double full_flux=0.0;
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (lagr_solver->m_electrodes[i].canEmit)
            full_flux+=lagr_solver->m_electrodes[i].eCurrent;
    }
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (lagr_solver->m_electrodes[i].canEmit)
            glColor3f(0.1*lagr_solver->m_electrodes[i].eCurrent/full_flux,0.0,1-0.1*lagr_solver->m_electrodes[i].eCurrent/full_flux);
        else
            glColor3f(1,1,1);
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
    //if(g_save_time2 * 1e15 >200)
    {
        double q_sum = 0;
        double p_pos=0.0;
        double p_full=0.0;
        double wall_coord = pz_solver->m_p[1].r_top.x - pz_solver->m_p[0].r_top.x;
        int i;
        for( i = 1; i < pz_solver->m_p_num && pz_solver->m_p[i].p > 0; i++ )
        {
        }
        wall_coord = pz_solver->m_p[i].r_top.x - pz_solver->m_p[0].r_top.x;

        for( i = 0; i < pz_solver->m_p_num; i++ )
        {
            q_sum += pz_solver->m_p[i].q_ext;
            if (pz_solver->m_p[i].p>0) p_pos += pz_solver->m_p[i].p;
            p_full+=fabs(pz_solver->m_p[i].p);
        }
        printf("wall = %e t = %e\n", wall_coord, g_t);
        fprintf(file_data,"%lf %e %lf %lf\n",g_t * 1e15, wall_coord, q_sum,p_pos/p_full);
    }
}


void saveStatePZ(char* fname)
{
    pz_solver->get_q();

    FILE *file_data_=fopen(fname,"w");

    fprintf(file_data_, "%d %e %e \n",pz_solver->m_p_num,g_phi,g_phi_max);

    for (int i=0;i<pz_solver->m_p_num;i++) {
        fprintf(file_data_,"%e %e %e %e %e %e %e\n", pz_solver->m_p[i].p,
                pz_solver->m_p[i].p_prev,
                pz_solver->m_p[i].E,
                pz_solver->m_p[i].E_prev,
                pz_solver->m_p[i].q_0,
                pz_solver->m_p[i].q,
                pz_solver->m_p[i].q_ext);
    }
    fclose(file_data_);
    printf("\n SAVED \n");

}

void loadStatePZ()
{
    pz_solver->get_q();

    FILE *file_data_=fopen("out.state","r");

    int num;
    fscanf(file_data_, "%d",&num);
    if (num==pz_solver->m_p_num)
    {
        fscanf(file_data_, " %lf %lf",&g_phi,&g_phi_max);

        for (int i=0;i<pz_solver->m_p_num;i++) {
            fscanf(file_data_,"%lf %lf %lf %lf %lf %lf %lf\n", &(pz_solver->m_p[i].p),
                   &(pz_solver->m_p[i].p_prev),
                   &(pz_solver->m_p[i].E),
                   &(pz_solver->m_p[i].E_prev),
                   &(pz_solver->m_p[i].q_0),
                   &(pz_solver->m_p[i].q),
                   &(pz_solver->m_p[i].q_ext)
                   );
        }


        multi_solver->updateTrajTable(false,1e-15,0.5e-7);
        multi_solver->solve(10);
        updateEulFields();
        printf("\n LOADED! \n");
    }
    else
    {
        printf("\n \n \n FAILED TO LOAD!!!! m_num=%d m_p_num_in=%d \n \n \n",num,pz_solver->m_p_num);
    }
    fclose(file_data_);

}


void savePotential()
{
    int nx, ny;
    double x0,x1,y0,y1;
    nx = 500;
    ny = 500;
    x0 = w_x0-80e-6;
    x1 = w_x1-250e-6;
    y0 = w_y0;
    y1 = w_y1;
    FILE *file_data_=fopen("out.dat","w");

    fprintf(file_data_, "TITLE = \" \" \n VARIABLES = \"x\" \"y\" \"z\" ");

    fprintf(file_data_, "\"phi\" ");
    fprintf(file_data_, "\"Ex\" ");
    fprintf(file_data_, "\"Ey\" ");
    fprintf(file_data_, "\n");
    fprintf(file_data_, "ZONE T=\" \" \n I=%d ,J=%d, K=%d \n", nx,ny,1);

    for (int j=0;j<ny;j++) {

        for (int i=0;i<nx;i++) {
            double x,y;
            x=x0 + (x1 - x0)/(nx-1.0)*(i);
            y=y0 + (y1 - y0)/(ny-1.0)*(j);
            vec2 E = lagr_solver->getE(x,y);
            fprintf(file_data_,"%e %e %e %e %e %e\n", x, y, 0.0, lagr_solver->getPhi(x, y), E.x, E.y);
        }
    }

    fclose(file_data_);
}

double angle=0.0;

double q_spec=2*9.5e4;


double proc=0.0;

void sim_prepare()
{
    multi_solver->dt_elec = 3e-11 / 1e6;
    double phi =550;
    g_emitElectrons = false;
    g_t = 0;
    g_phi = -phi;
    g_phi_max = g_phi;
    g_i_wall=0;
    multi_solver->init();
    progress=0;
    g_q_enable = true;
    for (int kk=0;kk<300;kk++)
        multi_solver->solve(10);
    g_phi_max =0.0;//*= -1;

    for (int kk=0;kk<300;kk++)
        multi_solver->solve(10);

    g_i_wall=g_i_wall_edge;
    g_i_wall_tmp=g_i_wall;
    g_phi_max =phi;//*= -1;
    for (int i=1; i<pz_solver->m_p_num;i++)
    {
        if (i<g_i_wall){
            pz_solver->m_p[i].q_ext=0;
            pz_solver->m_p[i].p = 0.3;
            pz_solver->m_p[i].p_prev = 0.3;
            pz_solver->m_p[i].q = 0;
        }
    }
    for (int kk=0;kk<100;kk++)
        multi_solver->solve(10);
    // g_q_enable = false;
    g_emitElectrons = true;
    g_t=0;
    g_save_time=0;
    g_save_time2=0;

    g_phi_max =phi;//*= -1;
}

void spec(int key, int x, int y)
{
    if (key==GLUT_KEY_UP)
    {
        q_spec+=0.01e4;
        for (int i=1; i<pz_solver->m_p_num;i++)
        {

            if (i<g_i_wall){
                pz_solver->m_p[i].q_ext=0;
                pz_solver->m_p[i].p = 0.26;
                pz_solver->m_p[i].p_prev = 0.26;
                pz_solver->m_p[i].q = 0;
            }
        }
        printf("q_spec= %e \n",q_spec);
    }

    if (key==GLUT_KEY_DOWN)
    {
        //q_spec-=0.01e4;
        for (int i=1; i<pz_solver->m_p_num;i++)
        {
            if (i<g_i_wall)
                pz_solver->m_p[i].q_ext+=180000.0*pow((rand()*1.0/RAND_MAX),4);
            //else
            //  pz_solver->m_p[i].q_ext=-q_spec*0.9;
        }
        printf("q_spec= %e \n",q_spec);
    }



    if (key==GLUT_KEY_PAGE_UP)
    {
        electrons_in_pack*=1.1;

        printf("electrons in pack= %f \n",electrons_in_pack);
    }

    if (key==GLUT_KEY_PAGE_DOWN)
    {
        electrons_in_pack/=1.1;

        printf("electrons in pack= %f \n",electrons_in_pack);
    }
    if (key==GLUT_KEY_INSERT)
    {
        g_i_wall=35;
        for (int i=1; i<pz_solver->m_p_num;i++)
        {
            if (i<g_i_wall){
                pz_solver->m_p[i].q_ext=0;
                pz_solver->m_p[i].p = 0.3;
                pz_solver->m_p[i].p_prev = 0.3;
                pz_solver->m_p[i].q = 0;
            }
        }
    }

    if (key==GLUT_KEY_HOME)
    {
sim_prepare();
    }
    glutPostRedisplay();
}

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
        view=VIEW_EX;
        printf("viewing Ex\n");
    }

    if (key=='3')
    {
        view=VIEW_EY;

        printf("viewing Ey\n");
    }


    if (key=='4')
    {
        view=VIEW_P;

        printf("viewing P\n");
    }


    if (key=='5')
    {
        // clearc=!clearc;
        view_px=!view_px;
    }


    if (key=='6')
    {
        // clearc=!clearc;
        saveStatePZ("out.state");
    }

    if (key=='7')
    {
        // clearc=!clearc;
        loadStatePZ();
    }



    if (key=='8')
    {
        //    wall_pos+=0.01;
        //    pz_solver->setWallPos(wall_pos);
        //  E_global=-E_global;
        //  printf("E_global=%e \n",E_global);

        g_emitElectrons=!g_emitElectrons;
    }


    if (key=='=')
    {

        multi_solver->updateTrajTable(true,1e-14,0.5e-7);
        //updateEulFields();
        // savePotential();
    }


    if (key=='-')
    {
        /*   for (int i=0;i<lagr_solver->m_elec_num;i++)
        {
            lagr_solver->m_electrodes[i].phi_fix *=-1.0;
        }*/
        g_phi_max *= -1;
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
        serialRegime =!serialRegime;
        display();
        glutPostRedisplay();
    }

    if (key==' ')
    {
        redr=!redr;



        /*   g_phi=200;
        for (int i=0;i<lagr_solver->m_elec_num;i++)
        {

            if (lagr_solver->m_electrodes[i].INDX==0) //first electrode -gphi
            {
                lagr_solver->m_electrodes[i].phi_fix=-g_phi;
            }else
            {
                lagr_solver->m_electrodes[i].phi_fix=g_phi;
            }
        }

            double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);

            for (int i=0;i<lagr_solver->m_elec_num;i++)
            {

                lagr_solver->m_electrodes[i].phi_fix_charges = pz_solver->getPhidepol(lagr_solver->m_electrodes[i].r.x,lagr_solver->m_electrodes[i].r.y)-phi_depol0;

            }
            lagr_solver->solve_ls_fast();
            updateEulFields();*/

    }



    if (key=='w')
    {
        scale_f+=0.01;
    }


    if (key=='s')
    {
        scale_f-=0.01;
    }

    if (key == 'u')
    {
       // g_use_wall=!g_use_wall;
       // printf("g_use_wall=%d \n",g_use_wall);

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


    multi_solver = new multiSolver();
    multi_solver->m_Esolver = lagr_solver;
    multi_solver->m_pzSolver = pz_solver;


    multi_solver->prepare_caches(is2D);
    multi_solver->fast_Fields_prepare();

    sim_prepare();
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
    /*glOrtho(w_x0-(w_x1-w_x0)*0.11,
            w_x0-(w_x1-w_x0)*0.11+(w_x1*1.3-(w_x0-(w_x1-w_x0)*0.3))*w*1.0/w0*(h0*1.0/h),
            w_y0*1.3,w_y1*1.3,
            -10.0, 10.0);*/

    g_x0=w_x0-(w_x1-w_x0)*0.03;
    g_x1=w_x0-(w_x1-w_x0)*0.03+0.1*(w_x1*1.3-(w_x0-(w_x1-w_x0)*0.3))*w*1.0/w0*(h0*1.0/h);
    g_y0= 0.1*w_y0*1.3;
    g_y1= 0.1*w_y1*1.3;


    glOrtho(g_x0,
            g_x1,
            g_y0,g_y1,
            -10.0, 10.0);
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
    glutSpecialUpFunc(spec);
    init();

    //save_file();
    glutMainLoop();
}
