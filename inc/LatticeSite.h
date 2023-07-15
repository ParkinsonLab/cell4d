//==============================================================
// LatticeSite.h
// Donny Chan, Billy Taj



#ifndef LATTICESITE_H
#define LATTICESITE_H

#include <map>
#include <vector>
#include "definition.h"
#include "Point.h"
#include <string>
#include <unordered_map>

using namespace std;

class Molecule_Container;

class LatticeSite{

    public:

        // Constructor and Destructor for the Lattice Site
        LatticeSite(Point &);
        ~LatticeSite();


        // Get the location of the Lattice Site
        Point * getLocation();

        // Get the environemnt of the Lattice Site
        //moltype get_compartment();
        string get_compart_name();

        // Set the environment of the site to the specified type
        //void set_compartment(envtype);
        void set_compartment(string);

        // Display properties of the Lattice Site
        void display();

        int getID();

        // convert between concentration and particle count
        static double particle_to_conc(int particles);
        static int conc_to_particle(double conc);
        static double conc_to_mole(double conc);
        static double mole_to_conc(double mole);
        static double mole_to_particle(double mole);

        void add_mole_to_voxel(string small_mol, double mole);

        // Stores the c-voxel neighbours objects using an int key that identifies the degree of contact with original
        // July 30th 2018 - added extra layer, each bulk molecule will have its own set of neighbors
        map <string, map<int, vector<LatticeSite *>>> neighbors;
        //map <string, vector<LatticeSite*> > neighbors;
        // Aug 20 2018 - A vector of molecule containers that reside in each lattice site
        unordered_map <string, unordered_map <unsigned long, Molecule_Container *>> lattice_container_map; // container name key, stores all containers within a lattice site
        
        // string -> conc map holding the moles of bulk molecules in this lattice space (voxel)
        map <string, double> bulk_moles_map;
        // tracks how much of this voxel's small molecules change in this timestep
        map <string, double> conc_change_timestep_map;

        // property for membranes with locked 2D movement
        // first is string indicating axis, second is the 1D coordinate of the locked plane where all particles reside in
        pair<string, float> locked_axis = make_pair("n",float(0));

    private:

        static int auto_id;

        int id;

        // Positon/Location of the Lattice Site
        Point pos;

        // Environment of the Lattice Site
        //envtype environment;
        string compartment_name;
}
;

#endif
