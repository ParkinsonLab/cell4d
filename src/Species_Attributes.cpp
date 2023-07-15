//==========================================================
// Species_Attributes.cpp
// Donny Chan, Billy Taj

#include "../inc/Species_Attributes.h"

Species_Attributes::Species_Attributes(){

}

Species_Attributes::~Species_Attributes(){
    
}

bool Species_Attributes::is_protein(){
    // simple check to see if it's a protein
    if((species_type.compare("PROTEIN") == 0) ||(species_type.compare("ENZYME") == 0)){
        return true;
    } else return false;
}

bool Species_Attributes::is_enzyme(){
    if(species_type.compare("ENZYME") == 0){
        return true;
    } else return false;
}

bool Species_Attributes::is_simple_molecule(){
    // simple check to see if it's a simple molecule, or whatever becomes of "SIMPLE_MOLECULE"
    // This was made for convenience, and is a typical design pattern when it comes to building software
    if(species_type.compare("SIMPLE_MOLECULE") == 0){
        return true;
    } else return false;
}


double Species_Attributes::particle_to_mole(int particles) {
    double converted_mole = particles / AVO_CONST;
    return converted_mole;
}

double Species_Attributes::mole_to_particle(double mole) {
    double converted_particle_count = mole * AVO_CONST;
    return converted_particle_count;
}
