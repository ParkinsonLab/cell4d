//========================================================
// Species_Attributes.h
// Donny Chan, Billy Taj

#ifndef SPECIES_ATTRIBUTES
#define SPECIES_ATTRIBUTES

#include <vector>
#include <string>
#include <iostream>
#include "definition.h"
#include <map>

using namespace std;
// This class stores the various details contained in the XML's species annotation.  
class Species_Attributes{
    public:
        Species_Attributes();
        ~Species_Attributes();

        bool is_protein();                              // ENZYMES will return true as well, since (almost) all enzymes are proteins
        bool is_simple_molecule();  
        bool is_enzyme();                               // Subtype of protein that is capable of acting as a modifier in a reaction

        // convert between concentration and particle count
        static double particle_to_mole(int particles);
        static double mole_to_particle(double mole);

        string name;                                    // the name of the species
        string species_type;                            // what kind of species this is:  simple? complex? something in between? protein? Gene? the "speciesMoleculeType"
        string ID;                                      // the ID field of the species.  why keep 2? one's human-readable, the other is unique
        vector<string> binding_sites;                   // list of binding site names, for protein-protein binding within a container
        map<string,vector<string>> possible_states;     // site name as key, store vector of possible states. molecule can have different states, and remain the same species
        vector<string> allowable_compartments;          // molecules can exist in different compartments
        map<string, map<string, float>> compartment_emission_rates; 
        bool final_reactant = false;                    // a flag to check if a particular species is a final product (that should kill the simulation) not currently used
        int current_total_amount = 0;                   // The current number of molecules of this species within the simulation at the current timestamp. Should be updated per timestep.
        double diffusion_constant;
        float color_red = 0;
        float color_green = 0;
        float color_blue = 0;
        bool membrane_display = false;

        int initial_amount = 0;                         // a place to store the initial amount of this particular molecule (replaces numMolecule) 
        map<string, int> compartments_initial_amount;   //The initial amount also needs to be associated with a compartment.  This is for a very specific objective in genSimulation
        map<string,string> default_states;              // species will be initialized in this state
        map<string,int> default_bind;                   // species will be initialized in this state
};

#endif
