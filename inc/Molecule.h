//============================================================
// Molecule.h
// Donny Chan, Billy Taj

#ifndef MOLECULE_H
#define MOLECULE_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include "Point.h"
#include "ParameterManager.h"
#include "RandomNG.h"
#include "definition.h"
#include "functionLib.h"
#include "LatticeEnvironment.h"
#include "Vector3D.h"


extern float resolution;
extern double timescale, spacescale;    // size of a unit grid, time
extern int xDim;


// This is a struct that holds a certain combination (configuration) of states for each species. Used by Container_Config which tries to do exact matches of container species and their states
struct Species_Config {
    string molecule_name;

    map<string, string> molecule_state_map;
    map<string, int> molecule_binding_map;

    public:

    // custom comparator function for the use of comparing Species_Config object instances
    // sort by molecule name, then states, then binding.
    friend bool operator<(const Species_Config& config_1, const Species_Config& config_2) {
        if(config_1.molecule_name < config_2.molecule_name) return true;
        else if(config_1.molecule_name > config_2.molecule_name) return false;
        else if(config_1.molecule_name == config_2.molecule_name) {
            if(config_1.molecule_state_map < config_2.molecule_state_map) return true;
            else if(config_1.molecule_state_map > config_2.molecule_state_map) return false;
            else if(config_1.molecule_state_map == config_2.molecule_state_map) {
                if(config_1.molecule_binding_map < config_2.molecule_binding_map) return true;
                if(config_1.molecule_binding_map > config_2.molecule_binding_map) return false;
            }
        }
        return false;
    }
};

class Molecule{

    public:

        //////////////////////////////////////////////////////////////////////
        // Constructor and Destructor

        Molecule(string name);
        Molecule(string name, unsigned long mol_id);
        Molecule(Species_Attributes * mol_info);
        virtual ~Molecule();

        //////////////////////////////////////////////////////////////////////
        // Operations on the Molecule

        Species_Config mol_config;

        // Position/Location of the Molecule
        unsigned long molecule_id = 0;    // a unique ID to label each molecule in the simulation
        unsigned long container_id;   // the id of the container that it is currently inside of
        Molecule_Container * container_pointer; // points to the container it's inside of, allows tracing back to the container
        Point pos;
        string molecule_name; // species name

        // Oct 2, 2018 - adding "*", so all sites have the "real" state and a "*" state for reaction checking
        void set_name(string name);
        string get_name();
        
        // Get the activation status of the molecule
        bool isActivated();

        // Determine if the specified molecule is near
        bool isNear(Molecule *);

        // Determine molecule reaction probability (based on distance)
        float neighbor_reaction_probability(Molecule *);

        // Get the diffusibility of the molecule
        bool isDiffusible();

        // Set the diffusibility of the molecule
        void setDiffusible(bool);

        // Set the activation status of the molecule
        void setActivationStatus(bool);

        // Deactivate the molecule
        void deactivate();

        // Get the string representation of the postion of the molecule
        string position();

        // set molecule to the specified position
        void set_position(Point &);

        // set molecule id
        void set_molecule_id();

        // Translate the molecule by the specified (x,y,z) vector
        void translate(Vector3D &);

        double get_difc();

        static bool mol_point_comp(Molecule * a, Molecule * b);

        //////////////////////////////////////////////////////////////////////
        // Virtual Methods

        // Set activation state of the molecule to true
        virtual void activate();

        // Activate the specified molecule
        virtual void activates(Molecule *);

        // Get the activation status of the molecule
        virtual bool isActivated(Molecule *);


    protected:

        //////////////////////////////////////////////////////////////////////
        // Member variables

        // Environment Manager
        LatticeEnvironment * lEnv;

        // diffusibility of the molecule
        bool diffusible;

        // The activation status of the molecule
        bool activated;

        // The degradation timer of a molecule if applicable
        int degradationTimer;

        double diffusion_constant = 3.0e-11; // use this (should be taken from param manager) to determine the mobility, units are m^2/s. set to be 30 um^2/s at the moment

        void initialize(unsigned long mol_id = 0);
        void initialize(Species_Attributes *, unsigned long mol_id = 0);
}
;

#endif
