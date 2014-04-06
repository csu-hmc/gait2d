// gait2de.h
// This file defines the data structure that holds model parameters.
// It is needed to pass model parameters to the Autolev generated C code
// in gait2de_al.c

#define NDOF 9 /* number of kinematic degrees of freedom */
#define NMOM 6 /* number of joint moments */
#define NSTICK 10 /* number of stick figure points */

typedef struct {

// Body segment parameters
	double TrunkMass, TrunkInertia, TrunkCMy;
	double ThighMass, ThighInertia, ThighCMy, ThighLen;
	double ShankMass, ShankInertia, ShankCMy, ShankLen;
	double FootMass, FootInertia, FootCMx, FootCMy;

// Parameters for the ground contact model
	double ContactY;
	double ContactHeelX, ContactToeX;
	double ContactStiff, ContactDamp, ContactV0, ContactFric;

} param_struct;

// prototype for the Autolev C function for multibody dynamics
void gait2d_al(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double mom[NDOF],
	double GRF[6],
	double Stick[NSTICK*2]);
