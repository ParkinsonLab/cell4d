//============================================================
// Simulation.h
// Donny Chan, Billy Taj

#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdlib>
#include "definition.h"
#include "LatticeEnvironment.h"
#include "Molecule_Container.h"
#include "Molecule_Container_Reactions.h"
#include "Molecule.h"
#include "LatticeSite.h"
#include "ParameterManager.h"
#include "Point.h"
#include "Compartment.h"
#include "Cell_Reaction.h"
#include "Container_Rules.h"
#include "Species_Attributes.h"
#include <unordered_map>
#include <functional>


using namespace std;

class Simulation
{

    public:
        
        // Static Variables and Operations
        // The Status of the Cell Lattice
        static bool isInstantiated;
        // Pointer to the Cell Lattice
        static Simulation * sim;

        //////////////////////////////////////////////////////////////////////
        // Static Operations on the Lattice
        // Get an instance of the Simulation
        //static Simulation * getSimulation();

        //////////////////////////////////////////////////////////////////////
        // Public Operations on Simulation

        // Continue Simulation
        void continueSimulation();

        // Create Molecule within the Simulation Space (num, type, env)
        void createMoleculeContainers(int, string, string);
        // Molecule_Container * createMoleculeContainers(int, vector, env);
        void createMoleculeContainers(int, vector<Species_Attributes *>, string env_name);
        // Molecule_Container * createMoleculeContainers(int, vector, Point);
        void createMoleculeContainers(int, vector<Species_Attributes *>, Point pos);

        // Molecule_Container * createMoleculeContainer(vector, int);
        void createMoleculeContainer(vector<Species_Attributes *>, Point *);

        // void recreate_molecule_container()

        // Display/output details of Interaction events (activator, activatee, [1|0])
        // void printInteractionDetails(string activator, string activatee, string product, moltype, bool);

        // Get current number of cycles executed (Time unit)
        int getCycles();

        // Get number of molecules existing in the Simulation Space
        // int getnumContainers();

        // Get number of molecules of the specifed type in the Simulation Space (type)
        // int getnumContainers(string);

        // Access the Molecule Table that contains all molecules in the Simulation
        unordered_map <string, unordered_map <unsigned long, Molecule_Container *> > & getMoleculeTable();

        // Is Simulation completed
        virtual bool isSimulationCompleted();

        // print stats at end of simulation
        void output_completed_stats(bool pop_small_mol = true, bool pop_containers = true, bool pop_compartment = true, bool multi_compart_out = true);

        void generate_checkpoint_file(bool is_final = false);

        void get_checkpointed_counts(json & json_file);

        // poslog feature functions
        void add_molpos_to_log(json & compart_json, int log_time, string compart);

        // holds all the newly created containers of this timestep
        vector<Molecule_Container *> new_containers_vector;
        

        // container name as key, each holds a vector of particles at each timecycle.
        // vector of total container counts at every timestep
        map<string, vector<int>> Container_counts_map;
        // vector of compartment-separated container counts at every timestep
        map<string, map<string, vector<int>>> Compart_container_counts_map;

        // vector of counts of small molecules at every timestep
        map<string, vector<double>> bulk_counts_map;
        map<string, map<string, vector<double>>> Compart_bulk_counts_map;

    protected:
        //////////////////////////////////////////////////////////////////////
        // Constructor and Destructor of the Lattice
        Simulation();
        virtual ~Simulation();



        //////////////////////////////////////////////////////////////////////
        // Member variables

        // Number of Molecules in the Simulation
        // int numContainers = 0;

        // Number of cycles executed (Time unit)
        int numCycles;

        // flag for Simulation completion
        bool simulationCompleted;

        //Enzyme Type of last signalling molecule
        moltype finalSigType = -1; // not actually used.  there's no default, so we gave it one.  may 28, 2018

        // Table of manages activation events (name <-> num activated)
        map <string, int> actEventTable;

        // Environment Manager
        LatticeEnvironment * lEnv;

        //////////////////////////////////////////////////////////////////////
        // Private Operations

        // Input Data
        virtual void initialDataInput();// = 0;

        void update_molecule_counts();

        // Display data from simulation
        virtual void displayResults();// = 0;

        virtual Molecule_Container * createNewMoleculeContainer(string, Point);

        virtual Molecule_Container * createNewMoleculeContainer(vector<Species_Attributes *>, Point);

        // Periodically update Behaviour/Properties
        virtual void updateSimulation();

        // Molecules move/diffuse randomly
        void move();

        // Check molecules interaction
        void check_protein_interactions();

        // check events
        void check_events();

        // do event
        void do_event(Event_Parameters * params);

        void recursive_traverse_event_list(vector<Event_Parameters *> event_vec);

        map<unsigned long, Molecule_Container *> get_list_of_event_containers(Container_Config, vector<Point> &);


        void checkEnzymaticReactions();
        bool check_spontaneous_reactions(Molecule_Container * mol_container);
        bool check_bulk_interactions(Molecule_Container * mol_container);
        // Interact with neighboring molecules
        void interact_neighbor_containers(Molecule_Container *);

        void do_reaction(Molecule_Container * mol_0, Molecule_Container * neighbor, string reaction_name);
        void do_reaction(Molecule_Container * mol_0, string bulk_name, string reaction_name);
        void do_reaction(Molecule_Container * mol_0, string reaction_name);

        
        // This version of build neighbor has int keys for a container, gives the reactions it participates in and how many times
        map<Molecule_Container *, map<string, int>> buildNeighborList(Molecule_Container *); 
        
        // print out stats of this map
        void print_counts(map<string, vector<int>> counts_map, string file_name, bool append = false);
        void print_counts(map<string, vector<double>> counts_map, string file_name, bool append = false);

}
;

#endif
