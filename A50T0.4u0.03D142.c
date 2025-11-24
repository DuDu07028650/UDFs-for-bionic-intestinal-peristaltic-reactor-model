#include "udf.h"
#include "sg.h" 
#include "sg_mphase.h" 
#include "flow.h" 
#include "mem.h" 
#include "metric.h" 
#include <math.h>

#define  A   0.0015        /* amplitude of the contraction, meters*/
#define  R   0.003         /* radius of the pipe, meters*/
#define  L   0.01        /* length of the contraction, meters*/
#define  T   0.4         /*频率*/
#define pi 4.*atan(1.)
#define diam2 1.42e-4

double sum_func(double x, int n)
{
	double sum = 0.0; 
	int i;
	for(i=1; i<= n; i++)
	{
    double y = exp(-pow((x - 0.0015 * i), 2)/0.000001);
	sum += y;
	}
    return sum;
}

DEFINE_GRID_MOTION(segmentation_upperwall,domain,dt,time,dtime)
{
	Thread *tf=DT_THREAD(dt); /*动网格面指针tf*/
	face_t f;
	Node *node_p;
	real x;
	real F=0.0;
	int n;
	int m = 24;


	SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));  /* 激活动壁面的相邻区域，避免偏斜*/
	begin_f_loop(f,tf)
	{
		f_node_loop(f,tf,n)
		{
			node_p=F_NODE(f,tf,n);
			if(NODE_POS_NEED_UPDATE(node_p))    /* 判断格点是否已经更新:否为1，是为0*/
				{
					
					NODE_POS_UPDATED(node_p);
					x=NODE_X(node_p); /*  获取node格点x坐标 */
					F= sin(2*M_PI*x/L);
					NODE_Y(node_p) = 0.006-A*sin(2*M_PI/T*CURRENT_TIME)*F;    /* 修改y坐标 */
					 
				}
		}
	}
	end_f_loop(f,tf)
}

DEFINE_GRID_MOTION(segmentation_lowerwall,domain,dt,time,dtime)
{
	Thread *tf=DT_THREAD(dt); /*动网格面指针tf*/
	face_t f;
	Node *node_p;
	real x;
	real F=0.0;
	int n;
	int m = 24;


	SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));  /* 激活动壁面的相邻区域，避免偏斜*/
	begin_f_loop(f,tf)
	{
		f_node_loop(f,tf,n)
		{
			node_p=F_NODE(f,tf,n);
			if(NODE_POS_NEED_UPDATE(node_p))    /* 判断格点是否已经更新:否为1，是为0*/
				{
					
					NODE_POS_UPDATED(node_p);
					x=NODE_X(node_p);           /*  获取node格点x坐标 */
					F=sin(2*M_PI*x/L);;
					NODE_Y(node_p) = A*sin(2*M_PI/T*CURRENT_TIME)*F;    /* 修改y坐标 */
					
				}
		}
	}
	end_f_loop(f,tf)
}

DEFINE_SOURCE(enzyme_source,c,c_thread,dS,eqn)
{
	int n_face, n;
	face_t f;
	Thread *t_wall;
	real tf = RP_Get_Real("flow-time");
	real source, S, V;
	real xc[ND_ND];
	C_CENTROID(xc, c, c_thread);
	real B[ND_ND];

	n_face=C_NFACES(c,c_thread);

	
		c_face_loop(c, c_thread, n)
		{
			f=C_FACE(c,c_thread,n);
			t_wall = C_FACE_THREAD(c,c_thread,n);
			if (xc[0] <= 0.1)
	      {
			if (THREAD_TYPE(t_wall)==THREAD_F_WALL)
			{
				if (fmod(tf, 1.0) >= 0.45 && fmod(tf, 1.0) < 0.55)
			 {
				if (n_face==4)
	           {
				V=C_VOLUME(c,c_thread);
				F_AREA(B,f,t_wall);	
				S=NV_MAG(B);
				source=S*0.0168/V; /* 酶量为30g/g,流速为20mm/s时，对应的分泌质量通量*/
				dS[eqn]=0;
				break;
				}
			 else
				V=C_VOLUME(c,c_thread);
				F_AREA(B,f,t_wall);	
				S=NV_MAG(B);
				source=S*0.01992/V; /* 酶量为30g/g,流速为20mm/s时，对应的分泌质量通量*/
				dS[eqn]=0;
				break;
			  }
				else
					source=0.0;
				
			}
	else
				source=0.0;
		}
			else 
				source=0.0;
		}
	
	return source;
}


 DEFINE_EXCHANGE_PROPERTY(custom_drag,cell,mix_thread,s_col,f_col)
 {
 Thread *thread_g, *thread_s;
 real x_vel_g, x_vel_s, y_vel_g, y_vel_s, abs_v, slip_x, slip_y,
 rho_g, rho_s, mu_g, reyp, void_g, vfac, fdrgs, taup, k_g_s;
 /* find the threads for the gas (primary) */
 /* and solids (secondary phases) */
 thread_g = THREAD_SUB_THREAD(mix_thread, s_col);/* gas phase */
 thread_s = THREAD_SUB_THREAD(mix_thread, f_col);/* solid phase*/
 /* find phase velocities and properties*/
 x_vel_g = C_U(cell, thread_g);
 y_vel_g = C_V(cell, thread_g);
 x_vel_s = C_U(cell, thread_s);
 y_vel_s = C_V(cell, thread_s);
 slip_x = x_vel_g - x_vel_s;
 slip_y = y_vel_g - y_vel_s;
 rho_g = C_R(cell, thread_g); rho_s = C_R(cell, thread_s);
 mu_g = C_MU_L(cell, thread_g);
 /*compute slip*/
 abs_v = sqrt(slip_x*slip_x + slip_y*slip_y);
 /*compute Reynolds number*/
 reyp = rho_g*abs_v*diam2/mu_g;
 /* compute particle relaxation time */
 taup = rho_s*diam2*diam2/18./mu_g;
 void_g = C_VOF(cell, thread_g);/* gas vol frac*/
 /*compute drag and return drag coeff, k_g_s*/
 vfac = 0.78*(exp(-46.43*(1.-void_g))-0.0388);
 fdrgs = void_g*(pow((0.63*sqrt(reyp)/
 vfac+4.8*sqrt(vfac)/vfac),2))/24.0;
 k_g_s = (1.-void_g)*rho_s*fdrgs/taup;
 return k_g_s;
 }

 DEFINE_PROFILE(backflow_mass_fraction,t,i)
{

 int n = 0;
 face_t f;

 begin_f_loop(f,t)
 {
 	F_PROFILE(f,t,i)=F_YI(f,t,n);
 }
 end_f_loop(f,t)
}

