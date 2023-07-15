//===================================================
// LatticeSite.cpp
// Donny Chan

#include "../inc/LatticeSite.h"
#include "../inc/Compartment.h"

int LatticeSite::auto_id = 0;
extern double spacescale;
extern map<string, Compart> compartments_map;   

//////////////////////////////////////////////////////////////////////
//  Constructor for the Lattice Environment
//
LatticeSite :: LatticeSite(Point & pt){

    id = auto_id;
    auto_id ++;


    compartment_name = "default";
    
    // Position of the Lattice Site
    pos.x = pt.x;    pos.y = pt.y;    pos.z = pt.z;
}



//////////////////////////////////////////////////////////////////////
//  Destructor for the Lattice Site
//
LatticeSite :: ~LatticeSite(){}



//////////////////////////////////////////////////////////////////////
//  Get the environment of the Lattice Site
//
string LatticeSite :: get_compart_name(){

    // return the environment of the Lattice Site
    return compartment_name;
}



//////////////////////////////////////////////////////////////////////
//  Set the Lattice Site to have the specified environment
//
void LatticeSite :: set_compartment(string compart_name){
    compartment_name = compart_name;
}

// (particles / 6.02e23 particles/mol) / (m^3 * 1000 L/m) = mol/L = M
double LatticeSite::particle_to_conc(int particles) {
    double converted_conc = (particles / AVO_CONST) / (spacescale * spacescale * spacescale * 1000);
    return converted_conc;
}

// conc is mol/L. mol/L * (m^3 * 1000 L / m^3) * 6.02e23 particles/mol = particles
int LatticeSite::conc_to_particle(double conc) {
    // example: 100 uM in 1e-7 m length voxel turns into 60 particles
    int converted_particle_count = conc * (spacescale * spacescale * spacescale * 1000) * AVO_CONST; // multiplies input by about 6e5 for spacescale of 1e-7
    return converted_particle_count;
}

double LatticeSite::conc_to_mole(double conc) {
    double converted_mole = conc * (spacescale * spacescale * spacescale * 1000);
    return converted_mole;
}

double LatticeSite::mole_to_conc(double mole) {
    double converted_conc = mole / (spacescale * spacescale * spacescale * 1000);
    return converted_conc;
}

double LatticeSite::mole_to_particle(double mole) {
    double converted_particle = mole * AVO_CONST;
    return converted_particle;
}

// adds the specified amount of small molecules to the voxel while also updating the ongoing compartment count
void LatticeSite::add_mole_to_voxel(string small_mol, double mole) {
    bulk_moles_map[small_mol] += mole;
    compartments_map[compartment_name].current_bulk_count[small_mol] += mole_to_particle(mole);
}

//////////////////////////////////////////////////////////////////////
//  Display different concentrations at the Lattice Site
//
void LatticeSite :: display(){

}


int LatticeSite :: getID(){
    return id;
}

//////////////////////////////////////////////////////////////////////
//  Get the location of the Lattice Site
//
Point * LatticeSite :: getLocation(){
    return & pos;
}

//////////////////////////////////////////////////////////////////////
