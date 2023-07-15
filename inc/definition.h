//=================================================================
// definition.h
// Donny Chan, Billy Taj

#ifndef DEFINITION_H
#define DEFINITION_H

///////////////////////////////////////////////////////////////////
// Input File Paths
///////////////////////////////////////////////////////////////////

#define INPUT_FILE_ENV "input/env.in"
#define INPUT_FILE_DIFFUSION_RATES "input/diffusion_rates.in"
#define INPUT_FILE_PARAMS "input/params.in"
#define INPUT_FILE_ACTIVATION_NETWORK "input/activation_network.in"
#define INPUT_FILE_COMPATIBLE_ENV "input/compatible_env.in"


#define SOLID 0
#define POINT 1
#define SPHEROID 2

///////////////////////////////////////////////////////////////////
// Molecules Types
///////////////////////////////////////////////////////////////////

// Membranous Molecules
#define MOL_MEMBRANOUS_FIXED_PUMP    0
#define MOL_MEMBRANOUS_FIXED_GATE    1

// Intracellular Molecules
enum {
  MOL_ENZYME_A = 0,
  MOL_ENZYME_B,
  MOL_ENZYME_C,
  MOL_ENZYME_D,
  MOL_ENZYME_E,
  NUM_MOL_ENZYME
};

// Metabolites
enum {
  META_A = 0,
  META_B,
  META_C,
  META_D,
  META_E,
  NUM_META
};


///////////////////////////////////////////////////////////////////
// Molecule Group Name
///////////////////////////////////////////////////////////////////

#define MOL_GROUP_NAME_INTRACELLULAR "Intracellular Molecule"
#define MOL_GROUP_NAME_MEMBRANOUS "Membranous Molecule"
#define MOL_GROUP_NAME_ENZYME "Metabolic Enzyme"

#define GATE_IS_CLOSED    0
#define GATE_IS_OPEN    1


typedef int moltype;
// typedef int envtype;

//////////////////////////////////////////////////////////////////
// Variable Constants
//////////////////////////////////////////////////////////////////
#define X_DIM "X_DIM"
#define Y_DIM "Y_DIM"
#define Z_DIM "Z_DIM"
#define MAX_CYCLES "MAX_CYCLES"
#define TIMESCALE "TIMESCALE"
#define SPACESCALE "SPACESCALE"
#define INACCESSIBLE_SPACE_PERCENT "INACCESSIBLE_SPACE_PERCENT"

const static float REACTION_SPEED_FACTOR = 1.0e3;

/////////////////////////////////////////////////////////////////////////////
// LATTICE ENVIRONMENT vars
/////////////////////////////////////////////////////////////////////////////
#define LATTICE_XPLANE  0
#define LATTICE_YPLANE  1
#define LATTICE_ZPLANE  2
#define LATTICE_SOLID   3
#define LATTICE_HOLLOW  4
#define LATTICE_POINT   5
#define LATTICE_SPHERE  6

///////////////////////////////////////
// CONSTANTS
///////////////////


const static double AVO_CONST = 6.02214076e23;
const static double PI = 3.14159265;

#endif
