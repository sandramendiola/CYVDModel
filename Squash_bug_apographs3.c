/* file Squash_bug_apographs.c */
#include <R.h>
#include <math.h>
static double parms[32];

#define d_A parms[0]
#define l parms[1]
#define b parms[2]
#define K_E parms[3]
#define d_E parms[4]
#define m_E parms[5]
#define d_1 parms[6]
#define m_1 parms[7]
#define d_2 parms[8]
#define m_2 parms[9]
#define d_3 parms[10]
#define m_3 parms[11]
#define d_4 parms[12]
#define m_4 parms[13]
#define d_5 parms[14]
#define m_5 parms[15]
#define p parms[16]
#define d_2o parms[17]
#define m_2o parms[18]
#define d_3o parms[19]
#define m_3o parms[20]
#define d_4o parms[21]
#define m_4o parms[22]
#define d_5o parms[23]
#define m_5o parms[24]
#define d_Ao parms[25]
#define c parms[26]
#define B_pb parms[27]
#define c_o parms[28]
#define B_bp parms[29]
#define B_o parms[30]
#define P0 parms[31]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=32;
odeparms(&N, parms);
}

/* Derivatives and 2 output variables */
void derivs (int *neq, double *t, double *y, double *ydot)
{

double OA = y[0];
double E = y[1];
double L1 = y[2];
double L2 = y[3];
double L3 = y[4];
double L4 = y[5];
double L5 = y[6];
double A = y[7];
double L2o = y[8];
double L3o = y[9];
double L4o = y[10];
double L5o = y[11];
double Ao = y[12];

double L2I = y[13];
double L3I = y[14];
double L4I = y[15];
double L5I = y[16];
double AI = y[17];

double L2oI = y[18];
double L3oI = y[19];
double L4oI = y[20];
double L5oI = y[21];
double AoI = y[22];
double P_I = y[23];
double Apo_PI = y[24];
double Sym_PI = y[25];
double OA_PI = y[26];
double OAI = y[27];

ydot[0] = -(d_A + l)*OA + c*OAI;
ydot[1] = b*.5*(OA + OAI)*(1 - E/K_E) - (d_E + m_E)*E;
ydot[2] = m_E*E - (d_1 + m_1)*L1;

ydot[3] = p*m_1*L1 - (d_2 + m_2)*L2 + c*L2I - B_pb*L2*P_I;
ydot[4] = m_2*L2 - (d_3 + m_3)*L3 + c*L3I - B_pb*L3*P_I;
ydot[5] = m_3*L3 - (d_4 + m_4)*L4 + c*L4I - B_pb*L4*P_I;
ydot[6] = m_4*L4 - (d_5 + m_5)*L5 + c*L5I - B_pb*L5*P_I;
ydot[7] = m_5*L5 - d_A*A + c*AI - B_pb*A*P_I;


ydot[8] = (1-p)*m_1*L1 - (d_2o + m_2o)*L2o + c*c_o*L2oI - B_pb*L2o*P_I;
ydot[9] = m_2o*L2o - (d_3o + m_3o)*L3o + c*c_o*L3oI - B_pb*L3o*P_I;
ydot[10] = m_3o*L3o - (d_4o + m_4o)*L4o + c*c_o*L4oI - B_pb*L4o*P_I;
ydot[11] = m_4o*L4o - (d_5o + m_5o)*L5o + c*c_o*L5oI - B_pb*L5o*P_I;
ydot[12] = m_5o*L5o - d_Ao*Ao + c*c_o*AoI - B_pb*Ao*P_I;

ydot[13] = B_pb*L2*P_I - (d_2 + m_2 + c)*L2I;
ydot[14] = m_2*L2I + B_pb*L3*P_I - (d_3 + m_3 + c)*L3I;
ydot[15] = m_3*L3I + B_pb*L4*P_I - (d_4 + m_4 + c)*L4I;
ydot[16] = m_4*L4I + B_pb*L5*P_I - (d_5 + m_5 + c)*L5I;
ydot[17] = m_5*L5I + B_pb*A*P_I - (d_A + c)*AI;


ydot[18] = B_pb*L2o*P_I - (d_2o + m_2o + c*c_o)*L2oI;
ydot[19] = m_2o*L2oI + B_pb*L3o*P_I- (d_3o + m_3o + c*c_o)*L3oI;
ydot[20] = m_3o*L3oI + B_pb*L4o*P_I- (d_4o + m_4o + c*c_o)*L4oI;
ydot[21] = m_4o*L2oI + B_pb*L5o*P_I- (d_5o + m_5o + c*c_o)*L5oI;
ydot[22] = m_5o*L5oI + B_pb*Ao*P_I- (d_A + c*c_o)*AoI;

ydot[23] = B_bp*(L2I + L3I + L4I + L5I + AI + OAI)*(P0 - P_I) + B_o*B_bp*(L2oI + L3oI + L4oI + L5oI + AoI)*(P0 - P_I);

ydot[24] = B_bp*B_o*(L2oI + L3oI + L4oI + L5oI + AoI)*(P0 - P_I);

ydot[25] = (B_bp*(L2I + L3I + L4I + L5I + AI))*(P0 - P_I);

ydot[26] = B_bp*OAI*(P0-P_I);

ydot[27] = -(d_A + l + c)*OAI;
}

/* END file Squash_bug_apographs.c */ 
