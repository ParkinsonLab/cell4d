//============================================================
// Cell_Reaction.h
// Donny Chan, Billy Taj


#ifndef CELL_REACTION
#define CELL_REACTION
#include <vector>
#include <string>
#include <iostream>
#include "definition.h"
#include <map>
#include <algorithm>
#include <cmath>
#include "Species_Attributes.h"

using namespace std;
// This class stores the various details contained in the XML's "list of Reactions"
class Cell_Reaction{
    public:
        Cell_Reaction();
        ~Cell_Reaction();
        string reaction_name;
        map<int, vector<string>> reactants_list;  // stores the reactants associated with this reaction, int index for container number
        map<int, vector<string>> products_list;   // stores the products associated with this reaction
        vector<string> modifiers_list;  // stores the modifier associated with this reaction

        // maps each element in the reactant list to the element of product list
        // this lets you know which element of the reactants is turning into which product
        // the data type for both are int pairs that indicates index of the container, then index of the position within that container
        map<pair<int, int>, pair<int, int>> react_prod_map;

        map<int, string> product_container_destinations; // int indicating index of product that now has an active transport destination

        // these maps have reactant/product container index key, for mol index inside the container, for their sites and states
        map<int, map<int, map<string,string>>>  reactants_state_map; // requirements on the reactants' modifier state (* will ignore). 
        map<int, map<string,string>>            modifiers_state_map; // requirements on the modifiers' modifier state (* will ignore). 
        map<int, map<int, map<string,string>>>  products_state_map; // requirements on the products' modifier state (* will copy its original state)
        map<int, map<int, map<string,int>>> reactants_binding_map; // requirements on the reactants' binding site (* will ignore its state)
        map<int, map<string,int>>           modifiers_binding_map; // requirements on the modifiers' binding site (* will ignore its state)
        map<int, map<int, map<string,int>>> products_binding_map; // requirements on the products' binding site (* will copy its original state)

        // the compartments where this reaction can occur. by default, applies to all compartments
        vector<string> allowed_compartments;

		float forward_reaction_rate = 0;
		float reverse_reaction_rate = 0;
		float equilibrium_constant = 0;
        float binding_radius = 0;
        float unbind_rad_min = 0;
        float unbind_rad_max = 0;
		float maximum_rate = 0;

        double reaction_probability = 0;

        bool is_reversible = false;
        Cell_Reaction * reverse_reaction;
        int spectator_index = -1;
        // if solo bulk molecules exist as a reactant, this will be 0 or 1, referring to index that gives the bulk molecule
        int bulk_index = -1;

        void determine_interaction_type();
        void determine_bulk_interaction_type();
        string get_interaction_type();
        bool is_two_reactants();
        bool is_fission();

        void set_reaction_rate(bool AB_adjust = true);
        void set_binding_radius(bool AB_adjust = true);
        void set_unbinding_radius(double unbinding_radius);

        void set_reactant_product_pairing();
        
    private:
        string interaction_type;
        double calculate_AB_radius(double target_r_rate, double m_difc);
        double find_AB_rate_from_rstep(double m_difc, double bind_radius);

};

#endif
