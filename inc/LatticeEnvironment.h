//=================================================
// LatticeEnvironment.h
// Donny Chan, Billy Taj
// todo: remove statematrix

#ifndef LATTICEENVIRONMENT_H
#define LATTICEENVIRONMENT_H

#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <set>


#include "definition.h"
#include "ParameterManager.h"
#include "LatticeSite.h"
#include "Matrix.h"
#include "Point.h"
#include "functionLib.h"
#include "Compartment.h"
#include "definition.h"

using namespace std;

// struct LatticeMacroFeature {
//     int    featureType;
//     float    x1, y1, z1, x2, y2, z2;
//     int    envType;
//     string env_name;
//     LatticeMacroFeature * next;
// };


struct diffuseInto {
    float multiplier;
    vector<int> ids;
    struct diffuseInto * next;
};


class LatticeEnvironment{

    // Static members

    // Does an Environment for the Lattice already exist
    static bool isInstantiated;
    // Pointer to the Environment for the Lattice
    static LatticeEnvironment * lEnv;


    public:

        

        // Get the Pointer to the Lattice Environment
        static LatticeEnvironment * getLatticeEnv();

        // Destructor for the Lattice Environment
        ~LatticeEnvironment();

        vector <string> get_compatible_env(string molecule_name);
        // Check whether the molecule is compatible in the environment
        bool isCompatibleEnvironment(string molecule_name, string compartment_name);
        bool isCompatibleEnvironment(vector<Species_Attributes *> species_attribute_vector, string compartment_name);
        bool is_membrane(string compartment_name);

        // Get the compartment at the specified location
        string get_compart_name(Point &);

        // Display the lattice Environment
        void displayLatticeEnv();

        void diffuse();

        // Get the LatticeSite at specified location
        LatticeSite * getLatticeSite(Point &);

        // get all lattice sites of this compartment
        vector <Point *> & getCompartmentSites(string compartment_name);
        

        void reactWithMetabolite(Point & point, string reactant, string product, float Ka, float Kp, float equilibrium_constant, float maximum_rate);

        Matrix * getStateMatrix();

        Matrix * getStatisticsMatrix();

        // LatticeMacroFeature * getLatticeMacroFeatures();

        void displayState();

        Tuple * getCurStat();

        // convert between concentration and particle count
        double particle_to_conc(double particles);
        double particle_to_mole(double particles);
        double conc_to_particle(double conc);
        double mole_to_particle(double mole);
        double conc_to_mole(double conc);
        double mole_to_conc(double mole);

        void add_mole_to_voxel(string name, double input_moles, LatticeSite *);
        
        // Determine whether the provided Point (location) is contained within the lattice dimensions
        bool withinCell(Point &);

    private:

        // Member Variables

        // The 3D Lattice Space (Pointers to Lattice Sites)
        LatticeSite ****c_voxel_list;

        // Compatible Table (Molecule -> compatible envs)
        //map <moltype, vector <envtype> > compatibleTable;
        //map <string, vector<int>> compatible_compartments_map; //replaces compatibleTable -> may 30, 2018
        // outdated:  info found in species_details_map, under each molecule

        // Compatible Table (Molecule -> compatible envs)
        //map <moltype, vector <Point *> > compatiblePointsTable;
        map<string, vector<Point*>> compatible_points_map; // using species id as key, identifies all points (c-voxels) that the species is allowed to reside in
        

        ///////////////////////////////////////////////////////////
        // Diffusion

        Matrix * stateMatrix;

        Matrix * updateMatrix;

        Matrix * statisticsMatrix;

        // Tuple *** diffusionRateTable;
        map<string, double> diffusion_coefficient_table;
        map<string, double> total_timestep_out;
        // Tuple ** reactionRateTable;
        map<string, diffuseInto **> diffusionPatterns;
        Tuple ** equilibriumConstantTable;

        // reduced constants for diffusion, this value multiplied by current conc gives amount to diffuse to each face-sharing neighbor
        map<string, double> diff_multiplier;


        vector<LatticeSite *> getNeighbors(int);

        void build_voxel_neighbors();

        void buildDiffusionRateTable();

        void buildReactionRateTable();

        void buildDiffusionPatterns();

        ///////////////////////////////////////////////////////////


        // member Operations

        // Constructor for the Lattice Environment
        LatticeEnvironment();

        // Build a 3D array of lattice sites
        void buildLatticeSpace();

        // Initialize the Lattice Environment
        void buildLatticeEnvironment();

        // lock membranous particle movement
        void set_membrane_properties(Compart&);

        // deprecated
        // void buildNucleus();

        // Called from buildLatticeEnvironment, initializes the inaccessible space
        void buildInaccessibleSpace();

        // Called from buildLatticeEnvironment, initializes the microtubules
        void buildMicrotubules(bool *** accessible_nucleus_lattice, int max_nucleus_to_block);

        int buildRadialCylinder(Point cell_membrane_end,
                    double endpoint_distance_from_center,
                    int radius, int environment, int max_cells_to_block,
                    bool *** accessible_nucleus_lattice, int * max_nucleus_to_block_ptr);

        int buildCylinder(Point end1, Point end2,
                  int radius, int environment, int max_cells_to_block,
                  bool *** accessible_nucleus_lattice, int * max_nucleus_to_block_ptr);

        int buildCylinder_or_check(Point end1, Point end2,
                       int radius, int environment, int max_cells_to_block,
                       bool build, bool *** accessible_nucleus_lattice, int * max_nucleus_to_block_ptr);

        // void buildPlasmaMembrane();


        bool withinCell(int, int, int);
}
;

#endif
