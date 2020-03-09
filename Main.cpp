
#include <stdio.h>
#include <stdlib.h>

#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/

//#include <my_include/gl.h>
//#include <my_include/glu.h>
//#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <iostream>
#include <vector>
#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "multisolver.h"
#include "electronlagrangian.h"

#include "phi_mult.h"


//#include <sys/time.h>



void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

//void fmm_step(double dt);
void sweep();
//void sweep_old();
int view=VIEW_PHI;

int redr=0;

double ck=2.0;

double cv=0.001;

double conv[100];
double conv1[100];

bool clearc=true;
/*#define MAIN
#include "fmm.h"
#undef MAIN
*/
#include "sse_sum.h"

std::vector <double>  avPy_;
std::vector <double>  avEy_;
std::vector <double>  area_;
std::vector <double>  T_;


bool view_px=true;

int i_tick=0;

eFieldLagrangian* lagr_solver;
pzSolver* pz_solver;
electronLagrangian* elec_solver;
multiSolver* multi_solver;
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
    printf("phi_c=%e Phi_p_c=%e Ey_c=%e fphic=%e Ep_c=%e ffp=%e \n",phic,phipc,Ec,dphic,Epc,dphicp);


    /*   double bmax=0.0;
    double Emax=0.0001;
    for(int i=0; i<N_X; i++ )
    {
        if (fabs(q[i])>bmax)
        {
            bmax=fabs(q[i]);//BoundaryLayerGauss[i]);
        }
    }

    i_tick++;
    if (i_tick>10){
        double emin=1e23;
        double emax=-1e23;
        for (int i=1; i<N_X-1; i++)
        {
            for (int j=1; j<N_Y-1; j++)
            {
                if (emin>Ey[i][j]) emin=Ey[i][j];
                if (emax<Ey[i][j]) emax=Ey[i][j];
            }
        }

        printf("emin=%e emax=%e \n",emin,emax);

        emin=1e23;
        emax=-1e23;
        for (int i=1; i<N_X-1; i++)
        {
            for (int j=1; j<N_Y-1; j++)
            {

                if (emin>Py_[i][j]) emin=Py_[i][j];
                if (emax<Py_[i][j]) emax=Py_[i][j];
            }
        }
        printf("pmin=%e pmax=%e pmean=%e\n",emin,emax,avPy_[avPy_.size()-1]);

        i_tick=0;
    }*/

    if (redr==1)
    {
        // for (int i=0;i<5;i++)
        multi_solver->solve(2);
        double pzmax=0.0;
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
        }
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
                l_2=ck*(phi_[i][j])/phi_max;
            if (view==VIEW_P)
                l_2=ck*(Py_[i][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j));

            if (view==VIEW_PHI)
                l_2=ck*(phi_[i+1][j])/phi_max;
            if (view==VIEW_P)
                l_2=ck*(Py_[i+1][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i+1][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i+1][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i+1][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(1.3*(w_x0+dx*(i+1)),1.3*(w_y0+dy*j));
        }
        glEnd();
    }

    /*    glEnable(GL_BLEND);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(1.5);

    double vel_scale=0.0;
    for( i=0; i<numParticles; i++ ) {
        vel_scale+=bodyVel[i].x*bodyVel[i].x+bodyVel[i].y*bodyVel[i].y+bodyVel[i].z*bodyVel[i].z;

    }

    vel_scale=sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.01*(w_x1-w_x0);
    glBegin(GL_LINES);

    for( i=0; i<numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1.0*fabs(bodyVel[i].x)/vel_scale,4.0*fabs(bodyVel[i].y)/vel_scale,3.0*fabs(bodyVel[i].z)/vel_scale);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex3f(bodyPos[i].x,bodyPos[i].y,bodyPos[i].z);
        glColor3f(0.0,0.0,0.0);
        glVertex3f(bodyPos[i].x+leng_sacle*bodyVel[i].x/vel_scale,bodyPos[i].y+leng_sacle*bodyVel[i].y/vel_scale,bodyPos[i].z+leng_sacle*bodyVel[i].z/vel_scale);
    }
    glEnd();




    for( i=0; i<N_X; i++ )
    {
        glColor3f((q[i])/bmax,-(q[i])/bmax,0);
        j = 0;
        glBegin(GL_QUADS);
        glVertex3f(w_x0+i*dx ,w_y0+j*dy ,0);
        glColor3f((q[i+1])/bmax,-(q[i+1])/bmax,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+j*dy ,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+(j-6)*dy ,0);
        glColor3f((q[i])/bmax,-(q[i])/bmax,0);
        glVertex3f(w_x0+i*dx ,w_y0+(j-6)*dy ,0);
        glEnd();

        glBegin(GL_POINTS);
        glColor3f(1,0,0);

        double esc_max = -1e20;
        double psc_max = -1e20;

        for(int i=0; i< avEy_.size();i++)
        {
            if(fabs(avEy_[i])>esc_max) esc_max = fabs(avEy_[i]);

            if(fabs(avPy_[i])>psc_max) psc_max = fabs(avPy_[i]);
        }
        for(int i=0; i< avEy_.size();i++)
        {
            glVertex3f(0.4 * (w_x1 - w_x0)*avEy_[i]/esc_max +w_x1/2 + w_x0/2 ,  0.4 * (w_y1 - w_y0)*avPy_[i]/psc_max ,0);
        }
        glEnd();

         double tmax = -1e20;


         glBegin(GL_POINTS);
         glColor3f(0,1,0);


        for(int i=0; i< T_.size();i++)
        {
            if(fabs(log(T_[i]))>tmax) tmax = fabs(log(T_[i]));

        }
        for(int i=0; i< T_.size();i++)
        {
            glVertex3f(0.4 * (w_x1 - w_x0)*log(T_[i])/tmax +w_x1/2 + w_x0/2 ,  0.4 * (w_y1 - w_y0)*0 ,0);
            glVertex3f(0.4 * (w_x1 - w_x0)*log(T_[i])/tmax +w_x1/2 + w_x0/2 ,  0.4 * (w_y1 - w_y0)*1 ,0);
            glVertex3f(0.4 * (w_x1 - w_x0)*log(T_[i])/tmax +w_x1/2 + w_x0/2 ,  0.4 * (w_y1 - w_y0)*area_[i] ,0);
        }
        glEnd();


    }
*/

    glPointSize(5);
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
        glColor3f(1,1,1);
        double x=0.5*(lagr_solver->m_electrodes[i].r0.x+lagr_solver->m_electrodes[i].r1.x);
        double y=0.5*(lagr_solver->m_electrodes[i].r0.y+lagr_solver->m_electrodes[i].r1.y);
        glVertex3f(x,y,0.0);
    }
    glEnd();


 glPointSize(2.5);
    glBegin(GL_POINTS);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        glColor3f(0.0,0.0,1.0);
        glVertex3f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y,elec_solver->m_bodyPos[i].z);
    }
    glEnd();

    double vel_scale=10.0;//sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.01*(w_x1-w_x0);
    glBegin(GL_LINES);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1.0*fabs(elec_solver->m_bodyVel[i].x)/vel_scale,4.0*fabs(elec_solver->m_bodyVel[i].y)/vel_scale,3.0*fabs(elec_solver->m_bodyVel[i].z)/vel_scale);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex3f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y,elec_solver->m_bodyPos[i].z);
        glColor3f(0.0,0.0,0.0);
        glVertex3f(elec_solver->m_bodyPos[i].x+leng_sacle*elec_solver->m_bodyVel[i].x/vel_scale,
                   elec_solver->m_bodyPos[i].y+leng_sacle*elec_solver->m_bodyVel[i].y/vel_scale,
                   elec_solver->m_bodyPos[i].z+leng_sacle*elec_solver->m_bodyVel[i].z/vel_scale);
    }
    glEnd();






    glColor3f(0.5,0.5,0.5);

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


    if (view_px)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            glColor3f(pz_solver->m_p[i].p/0.26,-pz_solver->m_p[i].p/0.26,0);
            glVertex2f(pz_solver->m_p[i].r.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r.y -pz_solver->m_p[i].dl*0.5);
            glVertex2f(pz_solver->m_p[i].r.x-pz_solver->m_dx*0.5,pz_solver->m_p[i].r.y +pz_solver->m_p[i].dl*0.5);
        }
        glEnd();
    }


    glBegin(GL_POINTS);
    glColor3f(1,0.5,0);
    glVertex2f(xc, yc);


    glEnd();

    glPointSize(3.0);

    glLineWidth(1.0);

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

int iglob=0;



//FmmSystem tre;



vec4<float> init_pos()
{
    vec4<float> res;

    res.x = w_x0+0.001*(w_x1-w_x0);
    res.y =w_y1/2 + my_rand(0)*(w_y1)/16;// 0.0+0.0251*(w_y1)+my_rand(0)*my_rand(0)*(w_y1)*0.75;
    res.z = rand()/(float) RAND_MAX*(w_z1-w_z0)+w_z0;
    res.w = qe; //single electrons

    return res;

}

vec3<float> init_vel()
{
    vec3<float> vel;

    double angle = -M_PI/4 ;//0.98*acos(1-2*my_rand(0)) - M_PI / 2;
    double en = getVms_from_Ev(0.8);//0.1e2;
    vel.x = en *  cos(angle);
    vel.y = en * sin(angle);


    return vel;
}


void fmm_init()
{
    int i;
    rand_init();

    for (i=0; i<N_X;i++)
    {
        BoundaryLayer[i]=0.000001;
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
        ck*=1.1;
    }

    if (key==',')
    {
        ck/=1.1;
    }
    if (key==']')
    {
        // dt*=1.1;
        //  printf("dt=%e \n",dt);
        angle+=10.;
        lagr_solver->setElectrodeAngle(angle);
        //pz_solver->m_dt*=1.1;

        //printf("dt=%e \n",pz_solver->m_dt);
    }

    if (key=='[')
    {
        angle-=10.;
        lagr_solver->setElectrodeAngle(angle);
        //  dt/=1.1;
        //  printf("dt=%e \n",dt);
       // pz_solver->m_dt/=1.1;
       // printf("dt=%e \n",pz_solver->m_dt);

    }

    if (key=='1')
    {
        view=VIEW_PHI;
        printf("viewing PHI \n");
    }


    if (key=='2')
    {
        view=VIEW_P;
        printf("viewing P \n");
    }

    if (key=='3')
    {
        view=VIEW_E;
        printf("viewing E\n");
    }

    if (key=='4')
    {
        view=VIEW_DIV;

        printf("viewing E\n");
    }

    if (key=='f')
    {
        sasign*=-1.0;
        for(int i=0;i<20;i++)
            sweep();
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

    if (key=='s')
    {
        sweep();
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

    if (key==' ')
    {
        redr=!redr;
    }

    glutPostRedisplay();
}


double check_frequency()
{
    if (avPy_.size()>2)
    {
        int imax=avPy_.size()-1;
        if ((avPy_[imax]<=0)&&(avPy_[imax-1]>0))
        {
            double alpha=fabs(avPy_[imax-1])/fabs(avPy_[imax-1]-avPy_[imax]);
            double e=avEy_[imax-1]*(1.0-alpha)+avEy_[imax]*(alpha);
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!! freq= %e Ec= %e \n", frequency,e);
            FILE* f=fopen("freqs.txt","a");
            fprintf(f,"%e %e \n", frequency,e);
            fclose(f);

            tt=0.0;
            frequency/=1.3;
            // avEy_.clear();
            // avPy_.clear();
        }
    }
}

double tt_prev=0.0;
double check_area()
{
    if ((tt-tt_prev)>5e-10)
    {
        double ar_p=0.0;
        double ar_m=0.0;
        double ar_z=0.0;

        double vol=(N_X-1)*(N_Y-1);
        for (int i=1; i<N_X-1; i++)
        {
            for (int j=1; j<N_Y-1; j++)
            {
                if (fabs(Py_[i][j]) <= 0.05)
                    ar_z+=1.0/vol;

                if (Py_[i][j]>0.05)
                    ar_p+=1.0/vol;

                if (Py_[i][j]<-0.05)
                    ar_m+=1.0/vol;
            }
        }
        FILE* f=fopen("1.8area_t_p_m_z.txt","a");
        fprintf(f,"%e %e %e %e\n", tt, ar_p,ar_m,ar_z);
        fclose(f);
        area_.push_back(ar_m/(ar_p+ar_m+ar_z));
        T_.push_back(tt);
        tt_prev=tt;
    }
}


void sweep_init()
{
    int i,j,k,n;
    double App,R_2,b,m;

    double topAv = 0.0;
    double botAv = 0.0;
    for (int i=0; i<N_X; i++)
    {
        q[i]=0.0;
        Pins_top[i]=pow(rand()*1.0/RAND_MAX,20)*50.0;
        topAv += Pins_top[i]/N_X;
        Pins_bottom[i]=pow(rand()*1.0/RAND_MAX,20)*50.0;
        botAv += Pins_bottom[i]/N_X;

    }
    for (int i=0; i<N_X; i++)
    {

        Pins_top[i] -= topAv ;
        Pins_bottom[i] -= botAv;
    }

    for (int i=1; i<N_X-1; i++)
    {
        q[i]=-1e-16;

    }
    q[0]=-1e-16;
    q[1]=-1e-16;

    float gsum=0.0;
    for (int i=0; i<N_Y; i++)
    {
        gau[i]=0.0;
        if (fabs(i-(N_Y_DIEL-1))<10)
        {
            int pow2n=1<<abs(i-(N_Y_DIEL-1));

            gau[i]=1.0/pow2n;
            gsum+=gau[i];
        }
    }

    for (int i=0; i<N_Y; i++)
    {
        gau[i]/=gsum;
    }

    for (i=0;i<N_X;i++)
    {
        for (j=0;j<N_Y;j++)
        {
            double x= (i-N_X/2)*dx;
            double y= (j-N_Y/2)*dy;

            div_[i][j]=0;
            rho_[i][j]=0.0;
            n_1[i][j]=exp(-(x*x+y*y)/((N_Y*0.1*dy)*(N_Y*0.1*dy)));
            n_1_prev[i][j]=n_1[i][j];

            n_2[i][j]=10e10;//10e10;
            p_in[i][j]=0;
            //  out_Ux[i][j]=sqrt(dy*j)*U;
            phi_[i][j]=0.0;//0.02 * sin(i*6*M_PI/(N_X-1));
            Ey0[i][j]=0.0;
            if ((i==0)||(i==N_X-1)||(j==0)||(j==N_Y-1))
            {
                n_1[i][j]=0.0;
                n_2[i][j]=0.0;
            }
        }
    }

    for (i=0;i<N_X;i++)
    {
        for (j=shift;j<N_Y_DIEL;j++)
        {

            phi_[i][j]=0;
            if (i>2)
                Py_[i][j]=0.26;//0 for freq//rand()*0.52/RAND_MAX - 0.26;//0.26 * sin(6.28*i*20.0/N_X + 3.14); //* sin(6.28*i*15.0/N_X /*+ 6.28*/);//rand()*0.52/RAND_MAX - 0.26;
            else
                Py_[i][j]=-0.26;
            Py0_[i][j]=Py_[i][j];
            Px_[i][j]=0.0;
            Px0_[i][j]=0.0;
        }
    }
}


void sweep()
{
    lagr_solver->updatePhi();
    lagr_solver->solvePhi(100);

    for (int i=0;i<N_X;i++)
    {
        for (int j=0;j<N_Y;j++)
        {
            phi_[i][j]=lagr_solver->getPhi(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j));



            Ey[i][j]=lagr_solver->getE(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j)).y;

            pz_solver->get_q();

            Ex[i][j]=pz_solver->getEdepol(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j)).y;


            Py_[i][j]=pz_solver->getPhidepol(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j));
        }
    }
}

void sweep_eul()
{
    double time = get_time();
    //    printf("StartTime = %f\n", time);
    double rho_min=1e30;
    double rho_max=-1e30;
    tt+=dt;
    double vol= fabs(dx*dy*(w_z1-w_z0));//volume of a cell
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            // rho_[i][j]=-n_1[i][j];//q_e*(-n_1[i][j])/eps_0;;//q_e*(-n_1[i][j]+n_2[i][j])/eps_0;
            RHS[i][j]=(Py_[i][j+1]-Py_[i][j-1])/(2.0*dy)/eps_0
                    /*+(Px_[i+1][j]-Px_[i-1][j])/(2.0*dx)/eps_0 ;*/ - 1.0*gau[j]*q[i]/(eps_0*vol);
            div_[i][j]=RHS[i][j];
            if (div_[i][j]>rho_max) rho_max=div_[i][j];
            if (div_[i][j]<rho_min) rho_min=div_[i][j];
        }
    }
    //   printf("rho_min = %e rho_max= %e q=%e  \n", rho_min,rho_max,-0.26/dy/eps_0);
    double omega=2*3.1415926535*frequency;
    double avE;
    for (int i=0;i<N_X;i++)
    {
        //phi_[i][N_Y-1]=0.0;
        //phi_[i][N_Y-1]=-sasign * 1.8 * (6 + Pins_top[i]);//100.*sin(omega*tt);//*(1.0-i*1.0/N_X);
        phi_[i][0]=0.0;//sasign * 1.8 * (6 + Pins_bottom[i]);//-100.*sin(omega*tt);

    }
    avE = (phi_[0][N_Y-1] - phi_[0][0]) / (w_y1 -w_y0);
    avEy_.push_back(avE);

    for (int j=N_Y_DIEL;j<N_Y;j++)
    {
        phi_[0][j]=-sasign*20.5;
    }

    INPUT_PARAM par;
    par.a=((2.0)/(dx*dx)+2.0/(dy*dy));
    par.bp=-1.0/(dx*dx);
    par.bm=-1.0/(dx*dx);
    par.cp=-1.0/(dy*dy);
    par.cm=-1.0/(dy*dy);

    par.w_bc_type=3;
    par.e_bc_type=1;
    par.n_bc_type=1;
    par.s_bc_type=3;

    par.w_bc_val=0.0;
    par.e_bc_val=0.0;
    par.n_bc_val=0.0;
    par.s_bc_val=0.0;
    double time1 = get_time();
    //    printf("BeforeMultigr = %f\n", time1 - time);
    //jacobi(par,phi_,div_,20);//phi
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    time1 = get_time();
    //     printf("AfterMultigr = %f\n", time1 - time);


    double mu=1.0;
    double D=1.0;

    INPUT_PARAM par2;
    par2.a=(1.0/dt+D*2.0/(dx*dx) + D*2.0/(dy*dy));
    par2.bp=-D/(dx*dx);
    par2.bm=-D/(dx*dx);
    par2.cp=-D/(dy*dy);
    par2.cm=-D/(dy*dy);

    par2.w_bc_type=0;
    par2.e_bc_type=0;
    par2.n_bc_type=1;
    par2.s_bc_type=1;

    par2.w_bc_val=-0.26;
    par2.e_bc_val=0.26;
    par2.n_bc_val=0.0;
    par2.s_bc_val=0.0;

    for (int i=0; i<N_X; i++)
    {
        for (int j=0; j<N_Y; j++)
        {
            Ex[i][j]=-ddx(phi_,i,j);
            Ey[i][j]=-ddy(phi_,i,j);// + rand() *8e7 / RAND_MAX - 4e7;
        }
    }

    /*for (int i=0; i<N_X; i++)
    {
        for (int j=0; j<N_Y; j++)
        {

            double pow2n_top=pow(1.35,abs(j-(N_Y_DIEL-1)));//1<<abs(j-(N_Y_DIEL-1));
            double pow2n_bot=pow(1.35,abs(j-(shift)));//1<<abs(j-(shift));


            Ey[i][j]+= 5e7*(Pins_top[i]/pow2n_top+Pins_bottom[i]/pow2n_bot);
        }
    }*/



    double xix=1.17*eps_0;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y_DIEL; j++)
        {
            //Px0_[i][j]=Px_[i][j];
            Px_[i][j]=xix*Ex[i][j];
        }
    }

    double dt_poly=1000*dt;
    INPUT_PARAM par_ferr;
    double kap=1.38e-10*0.15;
    par_ferr.a=(1.0/(dt_poly))+(kap*2.0/(dx*dx) + kap*2.0/(dy*dy));
    par_ferr.bp=-kap/(dx*dx);
    par_ferr.bm=-kap/(dx*dx);
    par_ferr.cp=-kap/(dy*dy);
    par_ferr.cm=-kap/(dy*dy);

    par_ferr.w_bc_type=0;
    par_ferr.e_bc_type=0;
    par_ferr.n_bc_type=1;
    par_ferr.s_bc_type=1;

    par_ferr.w_bc_val=-0.26;
    par_ferr.e_bc_val=-0.26;
    par_ferr.n_bc_val=0.0;
    par_ferr.s_bc_val=0.0;

    poly p;

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    double x0=-0.3;
    p.order=5;
    p.C[0]=2*alp*(T-T0) *0.5;//(T-T0); //x
    p.C[1]=0.0;         //xx
    p.C[2]=-4.0*bet *0.5;   //xxx
    p.C[3]=0.0;        //x^4
    p.C[4]=6.0*gam *0.5; //o.5 from crank-nikolson

    //-(fiy*0.0033- 2.0*alp*81*P1 - 4*bet*P1^3 +6*gam*P1^5)
    double emin=1e23;
    double emax=-1e23;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            double lapl0 = kap *( - Py0_[i][j] * (2.0 / (dx * dx) + 2.0 / (dy * dy)) + 1.0 / (dx * dx) * (Py0_[i+1][j] + Py0_[i-1][j]) +  1.0 / (dy * dy) * (Py0_[i][j+1] + Py0_[i][j-1]));
            double poly0 = p.C[0] * Py0_[i][j] +
                    p.C[2] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] +
                    p.C[4] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j];
            RHS_p[i][j]= /*1.4e7*(1.0-i*1.0/N_X)*/ -0.05*(Ey[i][j]+Ey0[i][j]) + lapl0  -poly0 +Py0_[i][j]/(dt_poly);

        }
    }

    //  printf("emin=%e emax=%e \n",emin,emax);
    time1 = get_time();
    //   printf("beforePolyn = %f\n", time1 - time);
    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 4);
    time1 = get_time();
    //   printf("AfterPolyn = %f\n", time1 - time);
    emin=1e23;
    emax=-1e23;
    double avPy = 0;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            if (emin>Py_[i][j]) emin=Py_[i][j];
            if (emax<Py_[i][j]) emax=Py_[i][j];
            if (shift < j &&  j < N_Y_DIEL)
                avPy+=Py_[i][j]/((N_X-1.0)*(N_Y-2.0*shift));
        }
    }
    avPy_.push_back(avPy);
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            Py0_[i][j]=Py_[i][j];
            Ey0[i][j]=Ey[i][j];

        }
    }

    for (int j=0; j<N_Y; j++)
    {
        Py0_[0][j]=Py0_[N_X-2][j];
        Py0_[N_X-1][j]=Py0_[1][j];

        Py_[0][j]=Py_[N_X-2][j];
        Py_[N_X-1][j]=Py_[1][j];

        phi_[0][j]=phi_[N_X-2][j];
        phi_[N_X-1][j]=phi_[1][j];

    }


    if (move_particles)
    {
        //for (int i=1;i<10;i++)
        //fmm_step(dt);
    }
    //check_frequency();

    time1 = get_time();
    //    printf("EndTime = %f\n",  time1 - time);

    check_area();
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


    rand_init();
    sweep_init();

    fmm_init();

    poly p;

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;

    double x1=0;
    double x0=0;
    double d_t=1e-11;
    for (int i=0;i<2000;i++)
    {
        rh=1.4e7+x0/d_t;

        p.order=5;
        p.C[0]=2*alp*(T-T0)+1.0/d_t; //x
        p.C[1]=0.0;         //xx
        p.C[2]=-4.0*bet;   //xxx
        p.C[3]=0.0;        //x^4
        p.C[4]=6.0*gam;

        x1=solve_poly(p,x1, rh,10);
        x0=x1;
        if (i%100==0)
            printf("i=%d x=%e \n",i,x1);
    }
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
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
