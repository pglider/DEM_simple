#ifndef DECLA_H
#define DECLA_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Define constants
#define PI 3.1415
#define G 9.81

//Grains physical parameters
#define R 0.001         //Radius in meters
#define RHO 2500.        //Density in kg.m^-3
#define M RHO*R*R*R*PI  //Mass in kg
#define KL 1e5          //Linear stiffness taken for a drop of 10d and a 1 percent oberlap (see Malone 2008)
#define KT 0.285*KL		 //Tangential stiffness
#define MU 0.5          //Friction coefficient
#define ELL 0.5         //Elastic restitution
#define ETA -2*log(ELL)*sqrt(M*KL)/sqrt(log(ELL)*log(ELL)+PI*PI) //DAMPING
#define MAXCONTACTS 6   //Maximum number of contacts

//Simulation parameters
#define N 120            //Total grains in simulation
#define NB 50           //Bottom grains
#define MAX_TIME 1000000    //Number of iterations
#define DT 0.05*2 * sqrt(M/KL) //Time step
#define ACTIVATION_STEP 10000 //Number of iterations before activating the next grain
#define Y0 24         //Initial height of the added grains, in number of radii

//Postscript parameters
#define EPS_STEP 5000     //Number of iterations before creating a new eps file
#define SCALE 20000      //Scale for the eps file

//Declare functions
void init();
int createps(int ieps);
void collision();
void newton_second_law();


//Declare global variables
int t,eps_count; //Time and eps file counter
int next_active_grain; //Next grain to activate
int r; //Return value for the functions
int i; //Counter for the grains
double dx; //Random horizontal speed

//Declare file pointers
FILE *pos_eps;


//Define the structure that wil handle the grains
typedef  struct grain
{
    int fixed;
	int active;

    double mass;
    double inertia;
    double r;

    double x;
	double xold;
	double dx;
	double y;
	double yold;
	double dy;
	double Oz;
	double Ozold;
	double dOz;

	double fx;
	double fy;
	double Mz;

	int n_contacts;                     //Number of contacts
	int contacts[MAXCONTACTS];          //List of contacts
	int contactsold[MAXCONTACTS];        //List of previous contacts

	//Contact integrals
	double utx[MAXCONTACTS];
	double uty[MAXCONTACTS];
	double utxold[MAXCONTACTS];
	double utyold[MAXCONTACTS];
} grain;

//Create the structure
grain disk[N];

#endif