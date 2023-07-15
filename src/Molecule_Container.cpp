//==================================================================
// Molecular_Container.cpp
// Donny Chan, Billy Taj

/*
    This class acts as a wrapper for molecules, to facilitate the 
    combination of a complex molecule joining together.
*/

#include "../inc/Molecule_Container.h"

extern bool graphics;

extern map<string, Cell_Reaction *> bulk_reaction_details_map;
extern map<string, Cell_Reaction *> unimolecular_reaction_details_map;
extern map<string, Cell_Reaction *> bimolecular_reaction_details_map;
extern int numContainers;
extern unordered_map <string, unordered_map <unsigned long, Molecule_Container *> > Mol_Container_Table;
extern unordered_map <string, unordered_map <unsigned long, Molecule *> > Mol_Table;
extern map<string, Compart> compartments_map;   
extern map<string, vector<string>> container_components_map; // stores the vector of molecule names in a container
extern map<string, vector<float>> container_color_map; // stores the colors of each container type
extern map<string, Species_Attributes *> species_details_map;
extern map<Container_Config, Container_Rules> container_rules_map;

map<Container_Config, bool> Molecule_Container::config_initialized_map;
map<Container_Config, map<int, map<string, map<int, vector <bool>>>>> Molecule_Container::config_bimolecular_reaction_map;
map<Container_Config, map<string, map<int, vector <bool>>>> Molecule_Container::config_unimolecular_reaction_map; 
map<Container_Config, map<string, map<int, vector <bool>>>> Molecule_Container::config_bulk_interaction_map;

Molecule_Container::Molecule_Container(Molecule *base, Point & loc){
    lEnv = LatticeEnvironment :: getLatticeEnv();
	vector<Molecule *> base_vec{base};
    Molecule_Container(base_vec, loc);
}

Molecule_Container::Molecule_Container(vector<Molecule *> base_vec, Point & loc){
    lEnv = LatticeEnvironment :: getLatticeEnv();
	container_vector = base_vec;
    set_position(loc);
	initialize_container();
    numContainers++;
}

// creating a blank container
Molecule_Container::Molecule_Container(){
    lEnv = LatticeEnvironment :: getLatticeEnv();
    set_container_id();
    pos = Point(float(0), float(0), float(0));
    numContainers++;
}

Molecule_Container::~Molecule_Container() {
    if(!container_vector.empty()) {
        for(auto & molecule : container_vector) {
            unordered_map <unsigned long, Molecule *>::iterator mol_itr_pos = Mol_Table[molecule->get_name()].find(molecule->molecule_id);
            if(mol_itr_pos != Mol_Table[molecule->get_name()].end()) {
                Mol_Table[molecule->get_name()].erase(mol_itr_pos);
            }
        }
    }

    remove_from_maps();
    numContainers--;
}

void Molecule_Container::initialize_container() {
    name_mol_container();
    set_container_config_status();
    add_container_component_map();
    get_reaction_list();
    if(graphics) {
        set_container_color();
    }
    if(container_name.empty()) return;
    if(container_id == 0) {
        set_container_id();
    }
    set_difc();
    set_mobility();
    set_emission_rates();
    if(compartments_map[lEnv->get_compart_name(pos)].is_membrane) {
        in_membrane = true;
    }
}

void Molecule_Container::set_container_config_status(){
    container_config.container_name = container_name;
    container_config.species_multiset.clear();
    container_rule_match = false;
    for(auto & mol : container_vector) {
        container_config.species_multiset.insert(mol->mol_config);
    }
    if(container_rules_map.find(container_config) != container_rules_map.end()) {
        container_rule_match = true;
    }
}

bool Molecule_Container::get_config_status() {
    return container_rule_match;
}

void Molecule_Container::name_mol_container() {
    sort(container_vector.begin(), container_vector.end(), Molecule::mol_point_comp);
	vector<string> name_vec;
   	for (int i = 0; i < container_vector.size(); i++){
   	    name_vec.push_back(container_vector[i]->molecule_name);
   	}
	sort(name_vec.begin(),name_vec.end());
	string combined_name = "";
	for(string & name_itr:name_vec) {
		combined_name += name_itr;
		combined_name += ";";
	}
    if(!combined_name.empty()) {
        combined_name.pop_back();
    }
	container_name = combined_name;
}

// when a molecule is (re)initialized, check 
void Molecule_Container::set_container_color() {
    if(container_rule_match) {
        if(container_rules_map[container_config].fields_filled["displayProperties"]) {
            container_color = container_rules_map[container_config].container_color;
            if(container_rules_map[container_config].fields_filled["membraneDisplay"]) {
                show_membrane = container_rules_map[container_config].membrane_display;
            }
            return;
        }
    }
    container_color["red"] = 0;
    container_color["green"] = 0;
    container_color["blue"] = 0;
    show_membrane = false;
    int colored_container_count = 0;
    // go through container components and average out the component colors
    for(auto & molecule : container_vector) {
        float mol_color_red = species_details_map[molecule->molecule_name]->color_red;
        float mol_color_green = species_details_map[molecule->molecule_name]->color_green;
        float mol_color_blue = species_details_map[molecule->molecule_name]->color_blue;
        if(mol_color_red > 0 || mol_color_green > 0 || mol_color_blue > 0) { // don't let invisible molecules affect color display
            container_color["red"] += mol_color_red;
            container_color["green"] += mol_color_green;
            container_color["blue"] += mol_color_blue;
            colored_container_count++;
        }
        if(species_details_map[molecule->molecule_name]->membrane_display) show_membrane = true;
    }
    // fix to undefined values if all molecules in complex are invisible
    if(colored_container_count > 0) {
        container_color["red"] /= colored_container_count;
        container_color["green"] /= colored_container_count;
        container_color["blue"] /= colored_container_count;
    }

    if(container_color_map.find(container_name) != container_color_map.end()) return;
    vector<float> rgb_vector = {container_color["red"] / 255, container_color["green"] / 255, container_color["blue"] / 255};
    container_color_map[container_name] = rgb_vector;
}

void Molecule_Container::set_container_id(){

    // using mt19937_64 (high period random generator), make a container id.
    container_id = RandomNG::randLong();

    for(auto & mol : container_vector) {
        mol->container_id = container_id;
        mol->container_pointer = this;
    }
}

void Molecule_Container::set_emission_rates() {
    emission_rate_map.clear();
    if(container_rule_match) {
        if(container_rules_map[container_config].fields_filled["compartmentEmissionRate"]) {
            emission_rate_map = container_rules_map[container_config].compartment_emission_rates;
            return;
        }
    }
    // only allow emission rates for registered molecule container configurations
    // for(auto & molecule : container_vector) {
    //     Species_Attributes * current_species = species_details_map[molecule->molecule_name];
    //     map<string, map<string, float>> species_emission_map = current_species->compartment_emission_rates;
    //     emission_rate_map.insert(species_emission_map.begin(), species_emission_map.end()); // downsteam issues in the future with conflicting emission rates, ignore for now
    // }
}

void Molecule_Container::add_to_maps() {
    // add this particle to the container table
    if(Mol_Container_Table[container_name].find(container_id) == Mol_Container_Table[container_name].end()) {
        Mol_Container_Table[container_name][container_id] = this;
    } else {

    }
    // add the new molecule to the appropriate lattice container map (if not already there)
    LatticeSite * site = get_lat_env()->getLatticeSite(pos);
    if(site->lattice_container_map[container_name].find(container_id) == site->lattice_container_map[container_name].end()) {
        site->lattice_container_map[container_name][container_id] = this;
    }
}

// remove this molecule container from the appropriate lattice map and container table
void Molecule_Container::remove_from_maps(Point last_pos) { 
    if(Mol_Container_Table[container_name].find(container_id) != Mol_Container_Table[container_name].end()) {
        Mol_Container_Table[container_name].erase(container_id);
    } else if(!is_small_mol) { // if the empty container was housing a small molecule, it won't be in the container table. just skip and leave it
        // if the container name has changed from original, search for the container in the whole system
        bool found = false;
        for(auto & container_type : Mol_Container_Table) {
            if(Mol_Container_Table[container_type.first].find(container_id) != Mol_Container_Table[container_type.first].end()) {
                Mol_Container_Table[container_type.first].erase(container_id);
                found = true;
                break;
            }
        }
        if(!found) {
            cerr << "Tried to remove a particle" << container_id << "from Mol_Container_Table" << endl;
        }
    }
    // guard for uninitialized containers being deleted 
    if(lEnv == NULL) return;

    // if old location is known, remove from that lattice map
    // again, only do this if it's not a blank container that was holding a bulk molecule, since it wouldn't be placed in the lattice map either
    if(!is_small_mol) {
        if(last_pos.x != -1) {
            LatticeSite * site = get_lat_env()->getLatticeSite(last_pos);
            if(site->lattice_container_map[container_name].find(container_id) != site->lattice_container_map[container_name].end()) {
                site->lattice_container_map[container_name].erase(container_id);
            } else {
                cerr << "Tried to remove a particle" << container_id << "from lattice list location" << site->getLocation()->x << site->getLocation()->y << site->getLocation()->z << endl;
            }
        } else {
            // if not known, try last known point
            LatticeSite * site = get_lat_env()->getLatticeSite(pos); // remove container from lattice container map
            if(site->lattice_container_map[container_name].find(container_id) != site->lattice_container_map[container_name].end()) {
                site->lattice_container_map[container_name].erase(container_id);
            } else { // if the container name has changed from original, search for the container in the whole system
                bool found = false;
                for(auto & container_type : site->lattice_container_map) {
                    if(site->lattice_container_map[container_type.first].find(container_id) != site->lattice_container_map[container_type.first].end()) {
                        site->lattice_container_map[container_type.first].erase(container_id);
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    cerr << "Tried to remove a particle" << container_id << "from lattice list location" << site->getLocation()->x << site->getLocation()->y << site->getLocation()->z << endl;
                }
            }
        }
    }


}

// building a container-specific reaction list, will take into account molecule state. this is done every time a container is generated/changed, for on-the-fly reaction network building
// bimolecular_reaction_map (and unimolecular_reaction_map if applicable) holding the types and number of reactions a container can do, at either position of the reacting containers.
// during neighbor checking, the origin will be reactant_0 and neighbors will be reactant_1
void Molecule_Container::build_reaction_list() {
    // if applicable, build unimolecular reaction table    
    unimolecular_reaction_map.clear();
    for(auto & reaction : unimolecular_reaction_details_map) {
        vector <string> reactant_container = reaction.second->reactants_list[0];
        // unimolecular reaction state/binding maps will only have 1 reactant container, so take first element in the maps
        map<int, map<string, string>> container_state_rules = reaction.second->reactants_state_map[0]; // reactant index (within the container) : site : state
        map<int, map<string, int>> container_bind_rules = reaction.second->reactants_binding_map[0]; // reactant index (within the container) : site : state
        map<int, string> container_vector_index_map;
        // create index:string map for container contents
        for(int i = 0; i < reactant_container.size(); i++) container_vector_index_map[i] = reactant_container[i];
        if(is_in_reaction(reaction.first, 0)) {
            map<int, vector<bool>> react_rule_map = fill_reaction_rules(container_state_rules, container_bind_rules, container_vector_index_map);
            if(!react_rule_map.empty()) {
                unimolecular_reaction_map[reaction.first] = react_rule_map;
            }
        }
    }


    bulk_interaction_map.clear();
    for(auto & reaction : bulk_reaction_details_map) {
        // bulk index gives position of bulk molecule in reactant map, container index gives position of container
        int container_index = abs(1 - reaction.second->bulk_index);
        vector <string> & reactant_container = reaction.second->reactants_list[container_index];
        map<int, map<string, string>> container_state_rules = reaction.second->reactants_state_map[container_index]; // reactant index (within the container) : site : state
        map<int, map<string, int>> container_bind_rules = reaction.second->reactants_binding_map[container_index]; // reactant index (within the container) : site : state
        map<int, string> container_vector_index_map;
        // create index:string map for container contents
        for(int i = 0; i < reactant_container.size(); i++) container_vector_index_map[i] = reactant_container[i];
        if(is_in_reaction(reaction.first, container_index)) {
            map<int, vector<bool>> react_rule_map = fill_reaction_rules(container_state_rules, container_bind_rules, container_vector_index_map);
            if(!react_rule_map.empty()) {
                bulk_interaction_map[reaction.first] = react_rule_map;
            }
        }
    }


    // build bimolecular reaction table
    bimolecular_reaction_map.clear();
    for(auto & reaction : bimolecular_reaction_details_map) {
        int reactant_index = 0; // 0 or 1, represents the first or second container of the reactant rules
        for(auto & reactant_container : reaction.second->reactants_list) { // will iterate twice, once for each reacting container
            map<int, map<string, string>> container_state_rules = reaction.second->reactants_state_map[reactant_index]; // reactant key, map to site:state
            map<int, map<string, int>> container_bind_rules = reaction.second->reactants_binding_map[reactant_index]; // reactant key, map to site:state
            map<int, string> container_vector_index_map;
            for(int i = 0; i < reactant_container.second.size(); i++) container_vector_index_map[i] = reactant_container.second[i];
            if(is_in_reaction(reaction.first, reactant_index)) {
                map<int, vector<bool>> react_rule_map = fill_reaction_rules(container_state_rules, container_bind_rules, container_vector_index_map);
                if(!react_rule_map.empty()) {
                    bimolecular_reaction_map[reactant_index][reaction.first] = react_rule_map;
                }
            }
            reactant_index++;
        }
    }
}

void Molecule_Container::get_reaction_list() {
    // add this configuration to the class reaction maps if it doesn't yet exist
    if(config_initialized_map.find(container_config) == config_initialized_map.end()) {
        build_reaction_list();
        config_bimolecular_reaction_map[container_config] = bimolecular_reaction_map;
        config_unimolecular_reaction_map[container_config] = unimolecular_reaction_map;
        config_bulk_interaction_map[container_config] = bulk_interaction_map;
        config_initialized_map[container_config] = true;
    } else {
        // otherwise take info from configuration maps
        bimolecular_reaction_map = config_bimolecular_reaction_map[container_config];
        unimolecular_reaction_map = config_unimolecular_reaction_map[container_config];
        bulk_interaction_map = config_bulk_interaction_map[container_config];
    }
    if(!unimolecular_reaction_map.empty()) {
        has_unimolecular_reactions = true;
    }
    if(!bulk_interaction_map.empty()) {
        has_bulk_interactions = true;
    }
}

bool Molecule_Container :: is_unimolecular_reactant() {
    return has_unimolecular_reactions;
}

bool Molecule_Container :: is_bulk_interaction_reactant() {
    return has_bulk_interactions;
}

// using formula for 1D diffusion, returns average distance travelled in 1 timestep in terms of Cell4D distance
void Molecule_Container::set_mobility() {
    mobility = sqrt(diffusion_constant * timescale * 2) / spacescale;
    diffusion_distribution = normal_distribution<double>(0, 1);
}

double Molecule_Container::get_mobility() {
    return mobility;
}

void Molecule_Container::set_difc() {
    if(container_rule_match) {
        if(container_rules_map[container_config].fields_filled["diffusionConstant"]) {
            diffusion_constant = container_rules_map[container_config].diffusion_constant;
            return;
        }
    }
    double difc_sum = 0;
    for(auto & mol : container_vector) {
        difc_sum += pow(mol->get_difc(), -3);
    }
    diffusion_constant = pow(difc_sum, double(-1)/double(3));
}

double Molecule_Container::get_difc() {
    return diffusion_constant;
}

bool Molecule_Container::is_in_reaction(string reaction_name, int reactant_index) { // whether or not the container participates in a reaction at that reactant index
    vector <string> species_names;
    map<int, vector<string>> react_components; // the list of required molecules in the reactant containers
    // try to find the reaction in the reaction maps, if not found in either, return false.
    if(bimolecular_reaction_details_map.find(reaction_name) != bimolecular_reaction_details_map.end()) {
        react_components = bimolecular_reaction_details_map[reaction_name]->reactants_list;
    } else if(unimolecular_reaction_details_map.find(reaction_name) != unimolecular_reaction_details_map.end()) {
        react_components = unimolecular_reaction_details_map[reaction_name]->reactants_list;
    } else if(bulk_reaction_details_map.find(reaction_name) != bulk_reaction_details_map.end()) {
        react_components = bulk_reaction_details_map[reaction_name]->reactants_list;
    } else return false;
    if(reactant_index > react_components.size() - 1) return false;
    for(auto & mol : container_vector) species_names.push_back(mol->molecule_name); // container names
    // if the container has the requisite reactants then return true
    // includes() require that the input ranges are sorted
    sort(species_names.begin(), species_names.end());
    sort(react_components[reactant_index].begin(), react_components[reactant_index].end());
    if(includes(species_names.begin(),species_names.end(),react_components[reactant_index].begin(),react_components[reactant_index].end())) {
        return true;
    }
    return false;
}

// static overloaded version that only uses container_name string
bool Molecule_Container::is_in_reaction(string container_name, string reaction_name, int reactant_index) { // whether or not the container participates in a reaction at that reactant index
    vector <string> species_names = container_components_map[container_name];
    map<int, vector<string>> react_components; // the list of required molecules in the reactant containers
    // try to find the reaction in the reaction maps, if not found in either, return false.
    if(bimolecular_reaction_details_map.find(reaction_name) != bimolecular_reaction_details_map.end()) {
        react_components = bimolecular_reaction_details_map[reaction_name]->reactants_list;
    } else if(unimolecular_reaction_details_map.find(reaction_name) != unimolecular_reaction_details_map.end()) {
        react_components = unimolecular_reaction_details_map[reaction_name]->reactants_list;
    } else return false;
    sort(species_names.begin(), species_names.end());
    sort(react_components[reactant_index].begin(), react_components[reactant_index].end());
    // if the container has the requisite reactants then return true
    if(includes(species_names.begin(),species_names.end(),react_components[reactant_index].begin(),react_components[reactant_index].end())) {
        return true;
    }
    return false;
}

// fills out boolean vectors that stores which container molecules fit with reactants
map<int, vector<bool>> Molecule_Container::fill_reaction_rules(map<int, map<string, string>> container_state_rules, map<int, map<string, int>> container_bind_rules, map<int, string> index_string_map) {
    vector <string> species_names, react_components, react_components_copy;
    map<int, vector<bool>> filled_rules_map;
    for(auto & reactant_name : index_string_map) react_components.push_back(reactant_name.second);
    for(auto & mol : container_vector) species_names.push_back(mol->molecule_name); // container names
    react_components_copy = react_components;
    sort(species_names.begin(), species_names.end());
    sort(react_components_copy.begin(), react_components_copy.end());
    if(includes(species_names.begin(),species_names.end(),react_components_copy.begin(),react_components_copy.end())) { // reactant container rules exist in current container
        map<int, vector<int>> reactant_index_fit_map;
        int reactant_index = 0; // counter for the reactant
        for(auto & reactant : react_components) { // do the state rules exist in the container?
            vector<int> mol_element_positions(container_vector.size(), 0);
            reactant_index_fit_map[reactant_index] = mol_element_positions;
            map<string, string> mol_states_rule = container_state_rules[reactant_index]; // modifier states of reactant i in reacting container
            map<string, int> mol_binding_rule = container_bind_rules[reactant_index]; // binding states of of reactant i in reacting container

            int container_index = 0; // this index is for filling the matrix
            for(auto & mol : container_vector) { // see if any molecule in the container match the state rules
                if(mol->get_name() != reactant) { // ensure that the name of the current mol matches the current mol in container_vec
                    container_index++;
                    continue; 
                }
                bool state_rules = true, bind_rules = true;
                for(auto & rule_site : mol_states_rule) { // ensures that every modifier state is fulfilled by the molecule
                    map<string, string>::iterator state_itr = find_if(mol->mol_config.molecule_state_map.begin(), mol->mol_config.molecule_state_map.end(), [rule_site](pair<string, string> mol_site) {
                        return (rule_site.first.compare(mol_site.first) == 0 && (rule_site.second.compare(mol_site.second) == 0 || rule_site.second.compare("*") == 0));
                    });
                    if(state_itr == mol->mol_config.molecule_state_map.end()) state_rules = false;
                }
                for(auto & rule_site : mol_binding_rule) { // ensures that every binding state is fulfilled by the molecule
                    map<string, int>::iterator bind_itr = find_if(mol->mol_config.molecule_binding_map.begin(), mol->mol_config.molecule_binding_map.end(), [rule_site](pair<string, int> mol_site) {
                        bool bind_match, site_match;
                        if(rule_site.first.compare(mol_site.first) == 0) site_match = true;
                        if(rule_site.second == -1) {
                            bind_match = true;
                        } else if(rule_site.second == mol_site.second) bind_match = true;
                        return (bind_match && site_match);
                    });
                    if(bind_itr == mol->mol_config.molecule_binding_map.end()) bind_rules = false;
                }
                if(bind_rules & state_rules) reactant_index_fit_map[reactant_index][container_index] = 1;
                container_index++;
            }
            reactant_index++;
        }
        vector<vector<int>> leaf_node_vectors = recursive_tree_leaf_collection(0, reactant_index_fit_map);
        // translate the vectors into a map of boolean vectors
        // elements are set to true at positions where the molecule is used as reactant
        // this filled_rules_map is used for neighbor list construction
        int vector_counter = 0;
        for(auto & vectors : leaf_node_vectors) {
            for(auto & element : vectors) {
                bool new_element;
                if(element == -1) new_element = true; else new_element = false;
                filled_rules_map[vector_counter].push_back(new_element);
            }
            vector_counter++;
        }
    }
    return filled_rules_map;
}

// recursive "tree" traversal, takes as input the current depth of the tree as well as the "matrix" of one-hot encoded reaction rule filling by the container
// returns a vector of vectors of values that indicate which slots are filled by the requirements
// each inner vector has "-1" at positions where the molecules are used as reactant
// the size of the outer vector indicates how many combinatorial ways there are to extract the required reactants from the given container
vector<vector<int>> Molecule_Container::recursive_tree_leaf_collection(int level, map<int, vector<int>> total_map) {
    vector<int> temp_vec = total_map[level];
    vector<vector<int>> output_vector_of_vectors;
    // base case
    if(level == total_map.size()-1) { 
        for(int i = 0; i < temp_vec.size(); i++) { // do this for all elements that fit the rule
            if(temp_vec[i] == 1) { // if the index fits the current reactant, set to -1 and add to vector of vectors
                vector<int> t = temp_vec; // make a copy, non-destructive method
                t[i] = -1;
                output_vector_of_vectors.push_back(t);
            }
        }
        return output_vector_of_vectors;
    }
    // recursive case
    for(int i = 0; i < temp_vec.size(); i++) { // do this for all elements that fit the rule
        if(temp_vec[i] == 1) {
            map<int, vector<int>> new_map = total_map;
            for(auto & map_level : new_map) {
                map_level.second[i] = -1;
            }
            vector<vector<int>> temp_vector_of_vectors = recursive_tree_leaf_collection(level+1,new_map); // not at leaf level, advance 1 layer deeper
            output_vector_of_vectors.insert(output_vector_of_vectors.end(), temp_vector_of_vectors.begin(), temp_vector_of_vectors.end()); // add results to output vector
        }
    }
    return output_vector_of_vectors;
}

void Molecule_Container::add_container_component_map() { // fill the global container_components_map that tells us what is in what container
    if(container_components_map.find(container_name) != container_components_map.end()) return;
    for(auto & mol : container_vector) {
        container_components_map[container_name].push_back(mol->molecule_name);
    }
}

void Molecule_Container::set_position(Point& base){
    pos = base;
}

string Molecule_Container::get_name(){
	return container_name;
}

Point & Molecule_Container::get_position(){
	return pos;
}

LatticeEnvironment * Molecule_Container::get_lat_env() {
    return lEnv;
}

map<string, string> Molecule_Container::get_molecule_states() { // combine and return state information of the whole container
    map<string, string> combined_container_state_map;
    for(auto & molecule : container_vector) {
        map<string, string> mol_state_map = molecule->mol_config.molecule_state_map;
        combined_container_state_map.insert(mol_state_map.begin(),mol_state_map.end());
    }
    return combined_container_state_map;
}

map<string, string> Molecule_Container::get_molecule_state(int id) { // return requested molecule's states
    vector<Molecule *>::iterator mol_pos = find_if(container_vector.begin(),container_vector.end(),[id](Molecule * mol) {
        return mol->molecule_id == id;
    });
    return (*mol_pos)->mol_config.molecule_state_map;
}

// Nov 26 new method (smoluchowski based)
double Molecule_Container::neighbor_reaction_probability(Molecule_Container * m, double react_radius) {
	Vector3D diff;
	// Get relative distances of the 2 molecules
	diff.setDiff(pos, m->pos);
	float dist = diff.magnitude();
    double threshold = react_radius / spacescale;
	// Molecules are not within the attraction distance
	if(dist > threshold) return 0;
	return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Container translates/moves by the specified vector
//
void Molecule_Container::translate(Vector3D & v){

    // Apply translation
    pos.x += v.x; pos.y += v.y; pos.z += v.z;
}

// samples a value from normal distribution of mean 0 and sdv of 1
float Molecule_Container::sample_normal() {
    return diffusion_distribution(RandomNG::get_generator());
}

// find the nearest voxel of the selected compartment space
void Molecule_Container::initialize_active_transport(string compart_name) {
    vector <Point *> compart_sites = lEnv->getCompartmentSites(compart_name);
    float distance = 100;
    Vector3D diff;
    Point * destination;
    // find the lattice site that the molecule is closest to
    for(Point* location : compart_sites) {
        diff.setDiff(*location, pos);
        if(diff.magnitude() < distance) {
            destination = location;
            distance = diff.magnitude();
        }
    }
    initialize_active_transport(*destination);
}

// overloaded version with a provided voxel target
void Molecule_Container::initialize_active_transport(Point voxel_target) {

    string compart_name = lEnv->get_compart_name(voxel_target);
    Point new_dest = voxel_target;
    // if the destination is a membrane compartment, adjust destination to the membrane plane
    if(compartments_map[compart_name].is_membrane) {
        LatticeSite * ls = lEnv->getLatticeSite(voxel_target);
        float new_locked_coord = ls->locked_axis.second;
        if(ls->locked_axis.first == "x") new_dest.x = new_locked_coord; 
        else if(ls->locked_axis.first == "y") new_dest.y = new_locked_coord; 
        else if(ls->locked_axis.first == "z") new_dest.z = new_locked_coord; 
    }
    actively_transported.first = true;
    actively_transported.second = new_dest;

    allowed_to_react = false;
}

bool Molecule_Container::is_actively_transported() {
    return actively_transported.first;
}

// This disregards all movement rules, will not check for compartment collision
void Molecule_Container::active_move() {
    // Gaussian sampling in 1D (line between current pos and destination), take absolute value since it won't move backwards
    float step_length = abs(sample_normal() * get_mobility());
    Vector3D step_vector, dest_vector;
    dest_vector.setDiff(actively_transported.second, pos);
    float distance_to_destination = dest_vector.magnitude();
    step_vector = dest_vector;
    step_vector.normalize();
    step_vector.scalarize(step_length);
    if(distance_to_destination < step_length) { // make sure it doesn't overshoot
        step_length = distance_to_destination;
        step_vector.normalize();
        step_vector.scalarize(distance_to_destination);
        actively_transported.first = false; // turn off active transport after destination is reached
        Point new_pos = Point(pos.x + step_vector.x, pos.y + step_vector.y, pos.z + step_vector.z);
        cout << "Container " << container_name << " has reached destination voxel " << new_pos.x << ", " << new_pos.y << ", " << new_pos.z << endl;
        if(compartments_map[lEnv->get_compart_name(new_pos)].is_membrane) in_membrane = true;
        else in_membrane = false;
    }
    // cout << "position x:" << pos.x << " y: " << pos.y << " z: " << pos.z << endl;
    Point old_pos = pos;
    translate(step_vector);
    if(old_pos.x != roundToInt(pos.x) || old_pos.y != roundToInt(pos.y) || old_pos.z != roundToInt(pos.z)) {
        LatticeSite * ls_0 = lEnv->getLatticeSite(old_pos);
        LatticeSite * ls_new = lEnv->getLatticeSite(pos);
        if(ls_0->lattice_container_map[container_name].find(container_id) != ls_0->lattice_container_map[container_name].end()) {
            ls_0->lattice_container_map[container_name].erase(container_id);
            ls_new->lattice_container_map[container_name][container_id] = this;
        }
    }
    // cout << "new position x:" << pos.x << " y: " << pos.y << " z: " << pos.z << endl;
}

// handle the emission of particles away from the membrane, to the voxel it is currently closest to
// if the adjacent voxel is not a valid compartment for the particle, nothing will happen.
void Molecule_Container::membrane_emission() {
    LatticeSite * cur_voxel = lEnv->getLatticeSite(pos);
    const Compart & cur_compart = compartments_map[lEnv->get_compart_name(pos)];

    // figure out which axis will be shifted in the case of emission, and in which direction
    float * changing_coord;
    Point test_point = pos;
    if(cur_compart.axis == "x") changing_coord = &test_point.x;
    else if(cur_compart.axis == "y") changing_coord = &test_point.y;
    else if(cur_compart.axis == "z") changing_coord = &test_point.z;
    if(cur_compart.face == "front") *changing_coord-= 0.00002;
    else if(cur_compart.face == "back") *changing_coord+= 0.00002;
    
    // check for validity of adjacent voxel
    LatticeSite * new_site = lEnv->getLatticeSite(test_point);
    bool compat_env = lEnv->isCompatibleEnvironment(container_name, lEnv->get_compart_name(test_point));
    if(!new_site || !lEnv->isCompatibleEnvironment(container_name, lEnv->get_compart_name(test_point))) return;

    // set new particle position and adjust lattice container maps, release from membrane
    set_position(test_point);
    LatticeSite * ls_new = lEnv->getLatticeSite(test_point);
    if(cur_voxel->lattice_container_map[container_name].find(container_id) != cur_voxel->lattice_container_map[container_name].end()) {
        cur_voxel->lattice_container_map[container_name].erase(container_id);
        ls_new->lattice_container_map[container_name][container_id] = this;
    }
    in_membrane = false;
}

void Molecule_Container::move() {
    // First determine if the particle is stuck in a not-allowed location, if so, don't allow movement
    // This is insurance against infinite recursions
    string starting_compartment = lEnv->get_compart_name(get_position());
    if(!lEnv->isCompatibleEnvironment(container_name, starting_compartment)) {
        cout << "container " << container_name << " is stuck where it shouldn't be." << endl;
        return;
    } 
    if(in_membrane) {
        // the emission rate entered in input will be an "emission constant" where, rate = constant * [conc]
        // and rate is the # of emission events per second.
        // integrating over timestep results in probability of a reaction given time step
        float emission_prob = float(1) - exp(compartments_map[starting_compartment].membrane_emission_rate * timescale * -1);
        float emission_check = RandomNG::randFloat(0,1);
        if (emission_prob >= emission_check) {
            membrane_emission();
            return;
        } 
    }
    move_recursive();
}

void Molecule_Container::move_recursive(float distance) {
    float x_comp, y_comp, z_comp;
    Vector3D translation;
    if(distance == -1) {
        x_comp = sample_normal() * mobility;
        y_comp = sample_normal() * mobility;
        z_comp = sample_normal() * mobility;
    } else {
        float x_sample = sample_normal();
        float y_sample = sample_normal();
        float z_sample = sample_normal();

        x_comp = distance * x_sample / (abs(x_sample) + abs(y_sample) + abs(z_sample));
        y_comp = distance * y_sample / (abs(x_sample) + abs(y_sample) + abs(z_sample));
        z_comp = distance * z_sample / (abs(x_sample) + abs(y_sample) + abs(z_sample));
    }

    if(in_membrane) {
        LatticeSite * cur_voxel = lEnv->getLatticeSite(pos);
        if(cur_voxel->locked_axis.first == "x") x_comp = 0;
        else if(cur_voxel->locked_axis.first == "y") y_comp = 0;
        else if(cur_voxel->locked_axis.first == "z") z_comp = 0;
    }

    float distance_magnitude = sqrt(x_comp*x_comp + y_comp*y_comp + z_comp*z_comp);
    translation.set(x_comp, y_comp, z_comp);

    // "intersections" number of elements represent the upper bound of # of boundaries 1 step can cross
    double intersections[((int)mobility + 2) * 3];
    int num_intersections;
    num_intersections = find_lattice_intersections(intersections, pos, translation);
    int cur_intersection_num;
    // distance to "peek" when crossing a boundary
    const static double INTERSECTION_EPSILON = 1.0e-3;

    // do this block if any intersections are found
    for (cur_intersection_num = 0; cur_intersection_num < num_intersections; cur_intersection_num++) {
        double intersection_distance;
        intersection_distance = intersections[cur_intersection_num];
        if (intersection_distance > distance_magnitude) {
            // An intersection should not have been detected if the particle isn't moving far enough to make the intersection.
            // This occurs because of the low precision of floats, the intersection distance calculation has a lot of operations that can cause a small error
            // past the 7th decimal point.
            break;
        }

        Vector3D intersection_translation;
        // peak at a point slightly farther than intersection point to see if there is collision
        intersection_distance += INTERSECTION_EPSILON;
        intersection_translation = Vector3D::normalize(translation); // normalizing translation vector results in unit vector of 1
        intersection_translation.scalarize(intersection_distance); // then scaled to move the distance needed to peek through intersection
        // intersection_location is a point slightly farther than the intersection
        Point intersection_location(roundToInt(pos.x + intersection_translation.x),
                                    roundToInt(pos.y + intersection_translation.y),
                                    roundToInt(pos.z + intersection_translation.z));

        string new_compartment = lEnv->get_compart_name(intersection_location);
        string current_compartment = lEnv->get_compart_name(pos);
        bool leave_compart = true, membrane_absorption = true;
        if(new_compartment != current_compartment) {
            // check to see if emission rates apply for this intersection
            // if fail, treat this intersection as impermeable
            if(emission_rate_map.find(current_compartment) != emission_rate_map.end()){
                if(emission_rate_map[current_compartment].find(new_compartment) != emission_rate_map[current_compartment].end()) {
                    float probability = emission_rate_map[current_compartment][new_compartment];
                    float rand_check = RandomNG::randFloat(0,1);
                    if(rand_check > probability) {
                        leave_compart = false;
                    } else {
                        leave_compart = true;
                    }
                }
            }
        }

        // check if the particle is allowed to attach to the membrane, if it is a membrane
        // if fail, treat this intersection as impermeable
        if(compartments_map[new_compartment].is_membrane && !in_membrane) {
            double absorb_rate = compartments_map[new_compartment].membrane_absorption_rate;
            // the absorption rate entered in input will be an "absorption constant" where, rate = constant * [conc]
            // and rate is the # of absorption events per second.
            // integrating over timestep results in probability of an event given a time step
            float absorb_prob = float(1) - exp(absorb_rate * timescale * -1);
            float absorb_check = RandomNG::randFloat(0,1);
            if (absorb_prob <= absorb_check) {
                membrane_absorption = false;
            } 
        }

        if ((new_compartment.compare("inaccessible") != 0) && lEnv->isCompatibleEnvironment(container_name, new_compartment) && leave_compart && membrane_absorption) {
            // Let molecule move within the lattice if target is a valid location, allowed for container, passes compart emission check, and passes membrane absorption check (when it applies)

        } else {
            // calculate distance to right before a collision, move there, then recalculate new distance for recursive movement
            intersection_distance = intersections[cur_intersection_num];
            unsigned long contain_id = container_id;
            if (intersection_distance > INTERSECTION_EPSILON) {
                intersection_distance -= INTERSECTION_EPSILON;
                intersection_translation = Vector3D::normalize(translation); // normalizing translation vector results in unit vector of 1
                intersection_translation.scalarize(intersection_distance); // then scaled to move the distance needed to peek through intersection
                float check = pos.x + intersection_translation.x;
                Point pre_intersect_loc(roundToInt(pos.x + intersection_translation.x),
                                        roundToInt(pos.y + intersection_translation.y),
                                        roundToInt(pos.z + intersection_translation.z));
                // this can trigger if there are multiple intersection points on the same move, where the particle can pass through the first one but not the next.
                // moving the hash of this container from old lattice site to new
                if(pre_intersect_loc.x != roundToInt(pos.x) || pre_intersect_loc.y != roundToInt(pos.y) || pre_intersect_loc.z != roundToInt(pos.z)) {
                    LatticeSite * ls_0 = lEnv->getLatticeSite(pos);
                    LatticeSite * ls_new = lEnv->getLatticeSite(pre_intersect_loc);
                    // due to rounding errors, the pre-intersect coordinate might still have crossed the boundary
                    // take incremental steps back until the bound is not crossed by this step
                    while(!ls_new) {
                        intersection_distance -= INTERSECTION_EPSILON;
                        intersection_translation = Vector3D::normalize(translation); // normalizing translation vector results in unit vector of 1
                        intersection_translation.scalarize(intersection_distance); // then scaled to move the distance needed to peek through intersection
                        Point pre_intersect_loc(roundToInt(pos.x + intersection_translation.x),
                                                roundToInt(pos.y + intersection_translation.y),
                                                roundToInt(pos.z + intersection_translation.z));
                        ls_new = lEnv->getLatticeSite(pre_intersect_loc);
                    }
                    if(ls_0->lattice_container_map[container_name].find(container_id) != ls_0->lattice_container_map[container_name].end()) {
                        ls_0->lattice_container_map[container_name].erase(container_id);
                        ls_new->lattice_container_map[container_name][container_id] = this;
                    }
                    // check to see if particle is going into a membrane, from not membrane.
                    if(ls_new->locked_axis.first != "n") {
                        if(ls_new->locked_axis.first == "x") intersection_translation.x = ls_new->locked_axis.second - pos.x;
                        else if(ls_new->locked_axis.first == "y") intersection_translation.y = ls_new->locked_axis.second - pos.y;
                        else if(ls_new->locked_axis.first == "z") intersection_translation.z = ls_new->locked_axis.second - pos.z;
                        in_membrane = true;
                    }
                }

                translate(intersection_translation);

            } else {
                intersection_distance = 0.0;
            }
            move_recursive(distance_magnitude - intersection_distance);
            return;
        }
    }

    Point new_location(roundToInt(pos.x + translation.x), roundToInt(pos.y + translation.y), roundToInt(pos.z + translation.z));
    string new_compartment = lEnv->get_compart_name(new_location);

    // if no intersections are found
    if ((new_compartment.compare("inaccessible") != 0) && lEnv->isCompatibleEnvironment(container_name, new_compartment)) {
        // Let molecule move within the lattice
        if(new_location.x != roundToInt(pos.x) || new_location.y != roundToInt(pos.y) || new_location.z != roundToInt(pos.z)) { // moving the hash of this container from old lattice site to new
            LatticeSite * ls_0 = lEnv->getLatticeSite(pos);
            LatticeSite * ls_new = lEnv->getLatticeSite(new_location);
            if(ls_0->lattice_container_map[container_name].find(container_id) != ls_0->lattice_container_map[container_name].end()) {
                ls_0->lattice_container_map[container_name].erase(container_id);
                ls_new->lattice_container_map[container_name][container_id] = this;
            }
            // check to see if particle is going into a membrane, from not membrane.
            if(ls_new->locked_axis.first != "n") {
                if(ls_new->locked_axis.first == "x") translation.x = ls_new->locked_axis.second - pos.x;
                else if(ls_new->locked_axis.first == "y") translation.y = ls_new->locked_axis.second - pos.y;
                else if(ls_new->locked_axis.first == "z") translation.z = ls_new->locked_axis.second - pos.z;
                in_membrane = true;
            }
        }
        translate(translation);
    } else {
        // Otherwise, try again.
        // This can happen even if there are no "intersections" with neighbouring lattice cells,
        // if the move would cause the particle to move right to the border.
        move_recursive(distance_magnitude);
    }

}

int Molecule_Container::find_lattice_intersections(double intersections[], Point pos, Vector3D translation) {
    int num_intersections = 0;

    double translation_magnitude;
    int cur_intersection_num;

    translation_magnitude = sqrt(translation.x * translation.x + translation.y * translation.y + translation.z * translation.z);

    // The +0.5 is because rounding is done to identify the lattice cell, not truncation.
    num_intersections += find_lattice_intersections_one_axis(intersections, pos.x + 0.5, translation.x, translation_magnitude, num_intersections);
    num_intersections += find_lattice_intersections_one_axis(intersections, pos.y + 0.5, translation.y, translation_magnitude, num_intersections);
    num_intersections += find_lattice_intersections_one_axis(intersections, pos.z + 0.5, translation.z, translation_magnitude, num_intersections);

    qsort(intersections, num_intersections, sizeof(double), compare_double);

    return num_intersections;
}

int Molecule_Container::find_lattice_intersections_one_axis(
    double intersections[],
    float pos_one_axis,
    float translation_one_axis,
    double translation_magnitude,
    int num_intersections) {

    int num_one_axis_intersections;
    float one_axis_intersection_at;
    float axis_ratio;
    int cur_intersection_num;

    if (translation_one_axis == 0.0) {
        return 0;
    }

    int num_new_intersections = 0;
    axis_ratio = translation_magnitude / fabs(translation_one_axis);

    if (translation_one_axis > 0) {
        num_one_axis_intersections = (int)(translation_one_axis + (pos_one_axis - (int)pos_one_axis));
        one_axis_intersection_at = 1 - (pos_one_axis - (int)pos_one_axis);
    } else {
        num_one_axis_intersections = (int)(1 - (translation_one_axis + (pos_one_axis - (int)pos_one_axis)));
        one_axis_intersection_at = (pos_one_axis - (int)pos_one_axis);
    }
    if (num_one_axis_intersections > 0) {
        intersections[num_intersections + num_new_intersections] = one_axis_intersection_at * axis_ratio;
        num_new_intersections++;

        for (cur_intersection_num = 1; cur_intersection_num < num_one_axis_intersections; cur_intersection_num++) {
            one_axis_intersection_at += 1.0;
            intersections[num_intersections + num_new_intersections] = one_axis_intersection_at * axis_ratio;
            num_new_intersections++;
        }
    }
    return num_new_intersections;
}

bool Molecule_Container::isDiffusible(){
    return diffusible;
}

bool Molecule_Container::is_in_membrane(){
    return in_membrane;
}

bool Molecule_Container::show_membrane_conc_colors(){
    return show_membrane;
}

int Molecule_Container::compare_double(const void *va, const void *vb) {
    double a = *((double *)va);
    double b = *((double *)vb);

    if (a > b) {
        return 1;
    } else if (a < b) {
        return -1;
    } else {
        return 0;
    }
}
