//========================================================
// Molecule_Container.h
// Donny Chan, Billy Taj

// This header accompanies the Complex Class, which is is to facilitate the formation of Molecular Complexes
// This structure is meant to contain molecules, as they travel in the simulation.  However, molecules tend to bond to things, 
// forming "complexes".  

#ifndef MOL_CONTAINER
#define MOL_CONTAINER

#include "Molecule.h"
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <set>
#include "Cell_Reaction.h"
#include "Container_Rules.h"

using namespace std;

// Using the species_config struct, has a custom comparator for matching specially defined combinations of species that form containers, described in listOfSpecies. This allows certain containers
// to use custom rules for diffusion, color display, etc, that do not match the implicitly generated rules that are from standard species-defined information
struct Container_Config {
    string container_name; // the semi-colon delimited combination of molecule names
    string config_name; // the name in the xml list of species
    multiset<Species_Config> species_multiset;
    
    // custom comparator function for the use of Container_Config object instances as map keys
    friend bool operator<(const Container_Config& config_1, const Container_Config& config_2) {
        if(config_1.container_name < config_2.container_name) return true;
        else if(config_1.container_name > config_2.container_name) return false;
        else if(config_1.container_name == config_2.container_name) {
            // vector comparator
            auto const & config_map_1 = config_1.species_multiset, config_map_2 = config_2.species_multiset;
            if(config_map_1 < config_map_2) return true;
            if(config_map_1 > config_map_2) return false;
        }
        return false;
    }

    // custom comparator function for the use of Container_Config object instances as map keys
    friend bool operator==(const Container_Config& config_1, const Container_Config& config_2) {
        if(config_1.container_name == config_2.container_name) {
            // vector comparator
            auto const & config_map_1 = config_1.species_multiset, config_map_2 = config_2.species_multiset;
            if(!(config_map_1 < config_map_2) && !(config_map_1 > config_map_2)) return true;
        }
        return false;
    }
};


class Molecule_Container{
    public:
        vector<Molecule *> container_vector;
        Molecule_Container(Molecule *, Point &);
        Molecule_Container();
        Molecule_Container(vector<Molecule *>, Point &);
        virtual ~Molecule_Container();

        void initialize_container();
        void set_container_id();

        // stores the configuration info of this container: name (components), the state/binding status of the components
        // reset everytime the container changes (forces re-initialization)
        Container_Config container_config;
        bool get_config_status();

        // container position
        Point pos;
        
        // red, blue, green keys for float RGB values of this container
        map<string, float> container_color;
        string container_name = "";
        bool allowed_to_react = false;
        map<int, map<string, map<int, vector <bool>>>> bimolecular_reaction_map; // reactant container index: reaction: number of times of that reaction: vector of elements used in reaction
        map<string, map<int, vector <bool>>> unimolecular_reaction_map; // reaction: number of times of that reaction: vector of elements used in reaction
        map<string, map<int, vector <bool>>> bulk_interaction_map; // reaction: number of times of that reaction: vector of elements used in reaction        
        
        unsigned long container_id = 0; //the hash key for the individual Mol_Container_Table entry
        bool diffusible = true;
        double diffusion_constant = 3.0e-12; // use this to determine the mobility, units are m^2/s. default at 3e-12 (3 um^2/s)

        pair <bool, Point> actively_transported; // first variable indicates whether or not it is being actively transported, second is the transport destination

        bool is_small_mol = false; // flag for a container that was only created to temp house a small molecule that separated from a container

        map<int, vector<bool>> fill_reaction_rules(map<int, map<string, string>>, map<int, map<string, int>>, map<int, string>); // checks if container fulfills the input map of state rules
        static vector<vector<int>> recursive_tree_leaf_collection(int level, map<int, vector<int>> total_map); // recursive method for tree traversal for leaf nodes
        bool is_in_reaction(string reaction_name, int reactant_index);
        static bool is_in_reaction(string container_name, string reaction_name, int reactant_index);
        void set_name(string name);
        string get_name();
        void set_container_color();
        void set_position(Point&);
        Point & get_position();
        void set_mobility();
        double get_mobility();
        void set_emission_rates();
        void set_difc();
        double get_difc();
        bool is_in_membrane();
        bool show_membrane_conc_colors();
        void translate(Vector3D & v);
        void move();
        map<string, string> get_molecule_states(); // return the state map of the molecule in question
        map<string, string> get_molecule_state(int); // return the state map of the molecule in question, given the hash id of the molecule

        ///////////////////////////////////////////////////////////////////
        // //for active transport
        bool is_actively_transported();
        void initialize_active_transport(string);
        void initialize_active_transport(Point voxel_target);
        void active_move();

        ///////////////////////////////////////////////////////////////////
        
        // Determine if the specified molecule is near
        double neighbor_reaction_probability(Molecule_Container *, double react_radius); // return the probability of this container interacting with the input container
        void name_mol_container();
        void build_reaction_list();
        void get_reaction_list();
        void add_container_component_map();
        // Get the diffusibility of the molecule
        bool isDiffusible();
		bool is_unimolecular_reactant();
        bool is_bulk_interaction_reactant();

        void add_to_maps();
        void remove_from_maps(Point last_pos = Point(-1,-1,-1));

        // Set the diffusibility of the molecule
        void setDiffusible(bool);

        LatticeEnvironment * get_lat_env();

    protected:

        //////////////////////////////////////////////////////////////////////
        // Member variables

        // Environment Manager
        LatticeEnvironment * lEnv;

        // settings for whether or not this container obeys config-specific rules
        void set_container_config_status();
        bool container_rule_match = false;
		bool has_unimolecular_reactions = false;
        bool has_bulk_interactions = false;
        bool in_membrane = false;
        bool show_membrane = false;
        double mobility; // average distance travelled in 1 timestep in terms of Cell4D distance
        double binding_radius;
		
        map<string, map<string, float>> emission_rate_map; // origin compart key: target compart key: emission probability

        void membrane_emission();
        void move_recursive(float distance = -1);
        int find_lattice_intersections(double intersections [], Point pos, Vector3D translation);
        int find_lattice_intersections_one_axis(
          double intersections [],
          float pos_one_axis,
          float translation_one_axis,
          double translation_magnitude,
          int num_intersections);
        static int compare_double (const void * va, const void * vb);

        // for random movement, each container has its own mobility distribution to sample from
        // recalculated everytime mobility is changed
        // lowers number of times distribution generation needs to be done
        normal_distribution<double> diffusion_distribution;
        float sample_normal();

        static map<Container_Config, bool> config_initialized_map;
        static map<Container_Config, map<int, map<string, map<int, vector <bool>>>>> config_bimolecular_reaction_map;
        static map<Container_Config, map<string, map<int, vector <bool>>>> config_unimolecular_reaction_map; 
        static map<Container_Config, map<string, map<int, vector <bool>>>> config_bulk_interaction_map;


};

#endif
