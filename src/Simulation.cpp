/*
===============================================================================

    FILE:  Simulation.cpp
    Donny Chan, Billy Taj
    
    This piece seems to be the top level.  All the global vars seems to reside here

===============================================================================
*/

#include "../inc/Simulation.h"
#include "../inc/Cell_Reaction.h"
using namespace std;

// Initialize Static Variables
bool Simulation::isInstantiated = false;
Simulation * Simulation::sim;

SBMLDocument_t *sbml_document;
Model_t *model;
//map<int, Species_t *> enzyme_map;
//map<int, Species_t *> metabolite_map;
string classInfo;
int numMetabolite=0;

Reaction_t **rxns;

int max_cycles = ParameterManager::getParameterManager()->get_int(MAX_CYCLES);
extern int xDim, yDim, zDim;
extern bool from_checkpoint;
extern json checkpoint_json;

extern string whole_bulk_out_name;
extern string compart_bulk_out_name;
extern string whole_particle_out_name;
extern string compart_particle_out_name;
extern string checkpoint_out_name;
extern int checkpoint_every_timestep;
extern bool multi_out;

// experimental poslog feature
extern bool poslog;
extern int poslog_interval;
extern string poslog_out;
extern vector<string> poslog_comparts;
// string:json map object of all poslog compartments to keep track of and then dump (experimental)
map<string, json> poslog_json;

// 2018 additions:
map<string, Cell_Reaction *> enzymatic_reaction_details_map;     // stores the names of reactions associated with a modifier, aka an enzyme
map<string, Cell_Reaction *> bulk_reaction_details_map;        // stores details of bulk/particle reactant reactions
map<string, Cell_Reaction *> unimolecular_reaction_details_map;        // stores reaction details of unimolecular reactions (product, reactants, etc etc.  will be expanded)
map<string, Cell_Reaction *> bimolecular_reaction_details_map;        // stores reaction details of bimolecular reactions (product, reactants, etc etc.  will be expanded)
map<string, Species_Attributes *> species_details_map;   // store reactant details. key: reactant ID
map<string, Compart> compartments_map;                  // store compartment details.  key: compartment ID
vector<string> small_molecule_id_list;                  // new vector of all small_molecule names
map<int, string> small_molecule_int_to_id_map;          // molecule id used as key to fetch ints used in tuple (for stats)
map<string, vector<string>> container_components_map; // stores the vector of molecule names within each container name (the key)
map<string, vector<float>> container_color_map;
map<Container_Config, Container_Rules> container_rules_map; // using configuration of a container as key, get any config-specific rules of a container 

// Table of intracellular molecule containers in the simulation grouped by type
// global variable, moved from a member var of simulation - Oct 23, 2018
unordered_map <string, unordered_map <unsigned long, Molecule_Container *> > Mol_Container_Table;
unordered_map <string, unordered_map <unsigned long, Molecule *>> Mol_Table;
int numContainers;

// container_name : container map
vector<Molecule_Container *> containers_to_remove_map;
vector<Molecule_Container *> containers_to_change_map;


//////////////////////////////////////////////////////////////////////
// Get the instance of the Cell Simulation
// ->    Only allow one instance of the cell to exist at any one time
//         If not yet exist, create a new Cell Simulation
//         otherwise return the pointer to the current existing Cell Simulation
//
//Simulation * Simulation::getSimulation(){}


//////////////////////////////////////////////////////////////////////
// Constructor for the Simulation
// ->    Create and initialize properties of the Simulation space
Simulation:: Simulation(){
    max_cycles = ParameterManager::getParameterManager()->get_int(MAX_CYCLES);
    // Cell Simulation Environment Manager
    lEnv = LatticeEnvironment::getLatticeEnv();

    // Initialize status before simulation
    simulationCompleted = 0;
    if(from_checkpoint) {
        numCycles = checkpoint_json["timesteps_completed"];
        get_checkpointed_counts(checkpoint_json);
    } else {
        numCycles = 0;
    }

    if(numCycles >= max_cycles) {
        cerr << "checkpoint file loaded is past the max timecycles of this model xml. exiting" << endl;
        exit(1);
    }
}



//////////////////////////////////////////////////////////////////////
// Destructor for the Simulation
//
Simulation:: ~Simulation(){
    // free everything
    unordered_map<string, unordered_map<unsigned long, Molecule_Container*>>::iterator outer_itr;
    unordered_map<unsigned long, Molecule_Container*>::iterator inner_itr;
    
    for(outer_itr = Mol_Container_Table.begin(); outer_itr != Mol_Container_Table.end(); outer_itr++){
        for(inner_itr = outer_itr->second.begin(); inner_itr != outer_itr->second.end(); inner_itr++){
            for(int i = 0; i < inner_itr->second->container_vector.size(); i++){
                delete inner_itr->second->container_vector[i];
            }
            delete inner_itr->second;
            outer_itr->second.erase(inner_itr);
        }
        Mol_Container_Table.erase(outer_itr);
    }

}

//////////////////////////////////////////////////////////////////////
//  Continue simulationCalvin Mo
//  ->    Events involved in each cycle
//
void Simulation:: continueSimulation(){
    const static int EQUILIBRIUM_TIME = 0;
    //  const static int METABOLITE_UPTAKE_TIME = 100;

    // Increment number of cycles (Time Units)
    numCycles ++;
    
    lEnv->diffuse();

    if (numCycles > EQUILIBRIUM_TIME) {
        // Check interaction between molecules
        checkEnzymaticReactions();
    }
    
    // Move the Molecules
    move();

    // check for reactive species reactions
    check_protein_interactions();

    // perform events
    check_events();

    // update molecule concentration/count tracking at every timestep
    update_molecule_counts();

    // do regular file checkpointing, default every 5000 timesteps
    if(numCycles % checkpoint_every_timestep == 0) {
        generate_checkpoint_file();
    }

    // experimental position logging of compartment molecules
    if(poslog) {
        // at the designated intervals, add molecule information of a compartment to the molpos log
        if(numCycles % poslog_interval == 0) {
            for(auto & compart_log : poslog_json) {
                add_molpos_to_log(compart_log.second, numCycles, compart_log.first);
            }
        }
    }

    if(numCycles==max_cycles) {
        simulationCompleted = true;
    }
    if(simulationCompleted) {
        // dump molecule counts into output stream and generate final checkpoint file
        output_completed_stats();
        return;
    }
}


//////////////////////////////////////////////////////////////////////
//  Create Molecule where following are specified
//  
// generate n molecules in given compartment
void Simulation::createMoleculeContainers(int num, string molecule_name, string compart_name){

    // Get list of available locations to place molecule
    vector <Point *> list = lEnv->getCompartmentSites(compart_name); // this returns all the sites that the molecule can sit in the compartment
    int numAvailableSites = list.size();
    cout << "Creating " << num << " of " << molecule_name << " in " << compart_name << " with " << numAvailableSites << " available sites." << endl;
    Compart & cur_comp = compartments_map[compart_name];
    if(cur_comp.is_membrane) {
        for(auto & loc : list) {
            LatticeSite * cur_voxel = lEnv->getLatticeSite(*loc);

        }
    }

    Molecule_Container * m;
    // For each molecule to be created
    for(int molNum = 0 ; molNum < num ; molNum ++) {
        // Randomly choose an allowable location
        int randomIndex = RandomNG::randInt(0,numAvailableSites-1);
        //cout << "randomIndex: " << randomIndex << endl;
        Point * location;
        location = ((list)[randomIndex]);
        Point rand_loc = *location;
        rand_loc.x = RandomNG::randFloat(rand_loc.x-0.5,rand_loc.x+0.5);
        rand_loc.y = RandomNG::randFloat(rand_loc.y-0.5,rand_loc.y+0.5);
        rand_loc.z = RandomNG::randFloat(rand_loc.z-0.5,rand_loc.z+0.5);

        // lock placement of particles in a membrane to 2D space
        if(compartments_map[lEnv->get_compart_name(*location)].is_membrane) {
            LatticeSite * cur_voxel = lEnv->getLatticeSite(*location);
            if(cur_voxel->locked_axis.first == "x") rand_loc.x = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "y") rand_loc.y = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "z") rand_loc.z = cur_voxel->locked_axis.second;
        }
        m = createNewMoleculeContainer(molecule_name, rand_loc); 
    }
}

// generating complex containers in a compartment
void Simulation::createMoleculeContainers(int num, vector<Species_Attributes *> container_info_vector, string compart_name) {
    // overloaded version
    // For each molecule to be created
    vector <Point *> list = lEnv->getCompartmentSites(compart_name); // this returns all the sites that the molecule can sit in the compartment
    // cout << "sim-no-point mol count: " << num << endl;
    int numAvailableSites = list.size();
    if(num == 0) return;
    Molecule_Container * m;
    // For each molecule to be created
    for(int molNum = 0 ; molNum < num ; molNum ++) {
        // Randomly choose an allowable location
        int randomIndex = RandomNG::randInt(0,numAvailableSites-1);
        Point * location;
        location = ((list)[randomIndex]);
        Point rand_loc = *location;
        rand_loc.x = RandomNG::randFloat(rand_loc.x-0.5,rand_loc.x+0.5);
        rand_loc.y = RandomNG::randFloat(rand_loc.y-0.5,rand_loc.y+0.5);
        rand_loc.z = RandomNG::randFloat(rand_loc.z-0.5,rand_loc.z+0.5);

        // lock placement of particles in a membrane to 2D space
        if(compartments_map[lEnv->get_compart_name(*location)].is_membrane) {
            LatticeSite * cur_voxel = lEnv->getLatticeSite(*location);
            if(cur_voxel->locked_axis.first == "x") rand_loc.x = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "y") rand_loc.y = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "z") rand_loc.z = cur_voxel->locked_axis.second;
        }
        m = createNewMoleculeContainer(container_info_vector, rand_loc); 
    }
    cout << "Created " << num << " " << m->container_name << " in " << compart_name << " with " << numAvailableSites << " available sites." << endl;
}

// overloaded version generating complex containers in a specific voxel
void Simulation::createMoleculeContainers(int num, vector<Species_Attributes *> container_info_vector, Point pos) {
    // For each molecule to be created
    Molecule_Container * m;
    // For each molecule to be created
    for(int molNum = 0 ; molNum < num ; molNum ++) {
        // Randomly choose a location within specified voxel
        Point rand_pos = pos;
        rand_pos.x = RandomNG::randFloat(pos.x-0.5,pos.x+0.5);
        rand_pos.y = RandomNG::randFloat(pos.y-0.5,pos.y+0.5);
        rand_pos.z = RandomNG::randFloat(pos.z-0.5,pos.z+0.5);

        // lock placement of particles in a membrane to 2D space
        if(compartments_map[lEnv->get_compart_name(rand_pos)].is_membrane) {
            LatticeSite * cur_voxel = lEnv->getLatticeSite(rand_pos);
            if(cur_voxel->locked_axis.first == "x") rand_pos.x = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "y") rand_pos.y = cur_voxel->locked_axis.second;
            else if(cur_voxel->locked_axis.first == "z") rand_pos.z = cur_voxel->locked_axis.second;
        }
        m = createNewMoleculeContainer(container_info_vector, rand_pos); 
    }
    cout << "Created " << num << " " << m->container_name << " in voxel " << pos.x << "," << pos.y << "," << pos.z << endl;
}


//////////////////////////////////////////////////////////////////////
//  Molecules randomly move within the Simulation
//
void Simulation:: move(){
    // Step through the Molecule table

    // Iterate through molecules of specified TYPE
    for(auto & container_type_pair : Mol_Container_Table){
        // Allow molecules of the specified TYPE to move
        for(auto & container_itr : container_type_pair.second) {
            // Get the next molecule in the list
            Molecule_Container * container = container_itr.second;
            // Move molecule
            if(container->isDiffusible()) {
                unsigned long id = container->container_id;
                if(!container->is_actively_transported()) { //checks for active transport or normal move
                    container->move();
                } else {
                    container->active_move();
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////
//  Check for Interaction between molecules
//
void Simulation :: check_protein_interactions() {
    // Iterate through the list of molecules of different types
    // this has to be iterated twice.  molecule type, each container of that type
    for (auto & container_type : Mol_Container_Table) {
		unordered_map<unsigned long, Molecule_Container *> map_of_container_type = container_type.second;
		for(auto & container_pair : map_of_container_type) {
            Molecule_Container * container = container_pair.second;
            if(container->allowed_to_react) {

                interact_neighbor_containers(container);
            }
		}
    }

    // add all of the new containers generated from reactions into the container table and the appropriate lattice map
    for(auto & container : new_containers_vector) {
        // don't work with bulk molecules, they will be deleted in the next map
        if(find(small_molecule_id_list.begin(), small_molecule_id_list.end(), container->get_name()) != small_molecule_id_list.end()) {
            continue;
        }
        container->initialize_container();
        container->add_to_maps();
    }
    new_containers_vector.clear();

    // clean out empty containers
    for(auto & container : containers_to_remove_map) {
        // if the container being removed is actually a bulk, add 1 particle equivalent back to the bulk conc
        if(find(small_molecule_id_list.begin(), small_molecule_id_list.end(), container->get_name()) != small_molecule_id_list.end()) {
            lEnv->getLatticeSite(container->get_position())->add_mole_to_voxel(container->get_name(), lEnv->particle_to_mole(1));
            delete container->container_vector[0];
            container->container_vector.pop_back();
        }
        if(!container->container_vector.empty()) {
            cout << "Container is not empty when being deleted" << endl;
        }
        delete container;
    }
    containers_to_remove_map.clear();

    for(auto & container : containers_to_change_map) { // new name : container
        container->remove_from_maps();
        container->initialize_container();
        container->add_to_maps();
    }
    containers_to_change_map.clear();

    for (auto & container_type : Mol_Container_Table) { // reset all molecules so they are allowed to react again for the next timestep
		unordered_map<unsigned long, Molecule_Container *> map_of_container_type = container_type.second;
		for(auto & container_pair : map_of_container_type) {
            if(container_pair.second->is_actively_transported()) continue;
            container_pair.second->allowed_to_react = true;
		}
    }
}

//  If this container can do unimolecular reactions, check for those prior to bimolecular reactions
bool Simulation::check_spontaneous_reactions(Molecule_Container * mol_container) {
    bool done_spontaneous_reaction = false;
    auto unimolecular_reaction_map = mol_container->unimolecular_reaction_map; // reaction name : number of combinations of reactants : reactant configuration per case
    float sum_reaction_prob = 0;
    map<int, float> reactions_probability;
    map<int, string> reaction_index_map;
    string compartment_name = lEnv->get_compart_name(mol_container->get_position());
    int reaction_index = 0;
    for(auto & element : unimolecular_reaction_map) {
        string reaction_name = element.first;
        vector <string> reaction_compartments = unimolecular_reaction_details_map[reaction_name]->allowed_compartments; // list of compartments where this reaction can take place
        if (find(reaction_compartments.begin(), reaction_compartments.end(), compartment_name) == reaction_compartments.end()) continue; // if current compartment isn't in list, do not react
        int num_of_combinations = element.second.size();
        reaction_index_map[reaction_index] = reaction_name;

        // the reaction rate entered in input will be reaction constant, rate = constant * [conc]
        // in unimolecular reactions, the constant is reactions per second per molecule
        // integrating over timestep results in probability of a reaction given time step
        float react_probability = float(1) - exp(unimolecular_reaction_details_map[reaction_name]->forward_reaction_rate * timescale * -1);
        float no_react_probability = exp(unimolecular_reaction_details_map[reaction_name]->forward_reaction_rate * timescale * -1);
        reactions_probability[reaction_index] = 1 - no_react_probability;

        reaction_index++;
    }
    // probability of reaction to this container is found by finding 1 - p(no reaction), so % of at least 1 reaction happening
    // if a reaction happens, use the vector of reaction probabilities to pick 1 reaction, weighed towards the more likely ones
    vector <float> probability_vector;
    float no_reaction_prob = RandomNG::randFloat(0,1);
    float reaction_occurring = 1;
    float no_reaction = 1;
    for(auto & reaction : reactions_probability) {
        probability_vector.push_back(reaction.second);
        no_reaction *= 1 - reaction.second;
    } 
    reaction_occurring = 1 - no_reaction;
    
    if (no_reaction_prob >= reaction_occurring) return done_spontaneous_reaction; // skip to next container if nothing will happen

    int chosen_reaction_index = RandomNG::custom_discrete_distribution(probability_vector);
    string chosen_reaction = reaction_index_map[chosen_reaction_index];
    cout << "Timestep " << getCycles() << ": reaction " << chosen_reaction << " " << mol_container->container_name << endl;
    // cout << "container " << mol_container->container_id << endl;
    do_reaction(mol_container, chosen_reaction);
    // cout << "done spontaneous reaction" << endl;
    done_spontaneous_reaction = true;
    return done_spontaneous_reaction;
}

bool Simulation::check_bulk_interactions(Molecule_Container * mol_container) {
    bool done_bulk_interaction = false;
    auto & bulk_reaction_map = mol_container->bulk_interaction_map; // reaction name : number of combinations of reactants : reactant configuration per case
    float sum_reaction_prob = 0;
    LatticeSite * current_site = lEnv->getLatticeSite(mol_container->get_position());
    map<int, float> reactions_probability;
    map<int, string> reaction_index_map;
    string compartment_name = lEnv->get_compart_name(mol_container->get_position());
    int reaction_index = 0;
    for(auto & element : bulk_reaction_map) {
        string reaction_name = element.first;
        string bulk_mol_name = bulk_reaction_details_map[reaction_name]->reactants_list[bulk_reaction_details_map[reaction_name]->bulk_index][0];
        // at least the equivalent of 1 bulk molecule must be in the voxel + 0-th degree neighbors for this to occur
        if(int(lEnv->mole_to_particle(current_site->bulk_moles_map[bulk_mol_name])) < 1) {
            double sum_particles = 0;
            for(auto & neighbor : current_site->neighbors[bulk_mol_name][0]) {
                sum_particles += neighbor->bulk_moles_map[bulk_mol_name];
            }
            sum_particles = lEnv->mole_to_particle(sum_particles);
            if(sum_particles < 1) continue;
        }

        vector <string> reaction_compartments = bulk_reaction_details_map[reaction_name]->allowed_compartments; // list of compartments where this reaction can take place
        if (find(reaction_compartments.begin(), reaction_compartments.end(), compartment_name) == reaction_compartments.end()) continue; // if current compartment isn't in list, do not react
        reaction_index_map[reaction_index] = reaction_name;

        // reaction rate will be based on the reaction radius and amount of bulk in the voxel
        // using Smoluchowski's collision equation rate = 4 * pi * mutual_difc * binding_R
        // gives rate constant of reaction, multiplied by bulk molecule amount divided by volume gives reactions/sec (rate)
        // scale rate to Cell4D timescale gives reactions per timestep, also known as reaction probability
        // this is the elementary rate const, multiply by AVO const * 100 to get the familiar M^-1s^-1 units
        double reaction_const = bulk_reaction_details_map[reaction_name]->forward_reaction_rate / (AVO_CONST * 1000);
        double reactions_per_sec = reaction_const * lEnv->mole_to_particle(current_site->bulk_moles_map[bulk_mol_name]) / (spacescale * spacescale * spacescale);
        double reactions_per_ts = reactions_per_sec * timescale;

        if(reactions_per_ts > 1) reactions_per_ts = 1;
        reactions_probability[reaction_index] = reactions_per_ts;
        reaction_index++;
    }
    // probability of reaction to this container is found by finding 1 - p(no reaction), so % of at least 1 reaction happening
    // if a reaction happens, use the vector of reaction probabilities to pick 1 reaction, weighed towards the more likely ones
    vector <float> probability_vector;
    float no_reaction_prob = RandomNG::randFloat(0,1);
    float reaction_occurring = 1;
    float no_reaction = 1;
    for(auto & reaction : reactions_probability) {
        probability_vector.push_back(reaction.second);
        no_reaction *= 1 - reaction.second;
    } 
    reaction_occurring = 1 - no_reaction;
    
    if (no_reaction_prob >= reaction_occurring) return done_bulk_interaction; // skip to next container if nothing will happen

    int chosen_reaction_index = RandomNG::custom_discrete_distribution(probability_vector);
    string chosen_reaction = reaction_index_map[chosen_reaction_index];
    string chosen_reaction_bulk_name = bulk_reaction_details_map[chosen_reaction]->reactants_list[bulk_reaction_details_map[chosen_reaction]->bulk_index][0];
    // cout << "container " << mol_container->container_id << endl;
    do_reaction(mol_container, chosen_reaction_bulk_name, chosen_reaction);
    cout << "Timestep " << getCycles() << ": reaction " << chosen_reaction << " " << mol_container->container_name << endl;
    done_bulk_interaction = true;
    return done_bulk_interaction;
}

//////////////////////////////////////////////////////////////////////
//  Check for interacting molecules of a container
//  Aug 21 2018 - This will require checking which components in the current container can conduct reactions (includes complexes within a container (AA + B -> AAB))
//  Oct 3 2018 - These will also check the states of the reactants to match requirements specified in a reaction rule
//              Note: This also allows us to implement reaction-specific reaction rates
void Simulation::interact_neighbor_containers(Molecule_Container * mol_container_0) {

    // sees if the container undergoes a unimolecular reaction. if so, skip the neighbor checking
    if(mol_container_0->is_unimolecular_reactant()) {
        bool spontaneous_done = check_spontaneous_reactions(mol_container_0);
        if(spontaneous_done) return;
    }
    if(mol_container_0->is_bulk_interaction_reactant()) {
        bool bulk_done = check_bulk_interactions(mol_container_0);
        if(bulk_done) return;
    }
    // list of surrounding containers that contain valid interaction partners (matches name and state)
    map<Molecule_Container *, map<string, int>> container_0_neighbors = buildNeighborList(mol_container_0); // reaction: how many times it does that reaction
    int num_neighboring_containers = container_0_neighbors.size();
    if(num_neighboring_containers == 0) return;
    string compartment_name = lEnv->get_compart_name(mol_container_0->get_position());
    map<int ,double> int_neighboring_container_reaction_probability_map; // Stores reaction probabilities of each neighboring container using neighbor_index int as id for the neighbor #
    map<int,Molecule_Container*> int_to_mol_container_map;
    // For each molecule that are valid reaction types of the origin molecule, calculate the probability of interaction
    float sum_prob = 0; // total probability of reacting with a (any) neighbor
    // every molecule in the origin container needs to be checked for interactions
    int neighbor_index = 0;

    for(auto & neighbor_element : container_0_neighbors) {
        Molecule_Container * neighbor = neighbor_element.first;
        int_to_mol_container_map[neighbor_index] = neighbor;
        int_neighboring_container_reaction_probability_map[neighbor_index] = 0;
        string neighbor_name = neighbor->container_name;
        
        for(auto & reaction : neighbor_element.second) {
            vector <string> reaction_compartments = bimolecular_reaction_details_map[reaction.first]->allowed_compartments; // list of compartments where this reaction can take place
            if (find(reaction_compartments.begin(), reaction_compartments.end(), compartment_name) == reaction_compartments.end()) continue; // if current compartment isn't in list, do not react
            float reaction_radius = bimolecular_reaction_details_map[reaction.first]->binding_radius;
            double base_distance_probability = mol_container_0->neighbor_reaction_probability(neighbor, reaction_radius);
            int_neighboring_container_reaction_probability_map[neighbor_index] = base_distance_probability;
            // limit it to one reaction for now
            if(base_distance_probability > 0) break;
        }
        sum_prob += int_neighboring_container_reaction_probability_map[neighbor_index];

        neighbor_index++;
    }   

    if(sum_prob <= 0) return;

    // Calculate the probability of no reaction happening for this mol_container_0 at the timestep
    float no_reaction_prob = 1;
    for(auto & container_prob_itr : int_neighboring_container_reaction_probability_map) {
        no_reaction_prob *= 1 - container_prob_itr.second;
    }
    float reaction_occurring = RandomNG::randFloat(0,1);
    // cout << "reaction prob:" << reaction_occurring << endl;
    if (no_reaction_prob > reaction_occurring) return; // skip to next container if nothing will happen

    // using the Mersenne twister random number generator, pick a random neighboring container from the list with the previous reaction probabilities as weights.
    vector <float> react_prob_vector;
    for(auto & itr : int_neighboring_container_reaction_probability_map) react_prob_vector.push_back(itr.second);
    int chosen_neighbor_index = RandomNG::custom_discrete_distribution(react_prob_vector);
    Molecule_Container * chosen_neighbor = int_to_mol_container_map[chosen_neighbor_index];

    vector <string> chosen_reactions;

    // add shared reactions between this container and neighbor, and randomly pick one. 
    // do not allow reactions that cannot occur in current compartment
    int reaction_index = 0;
    for(auto & reaction_pair : container_0_neighbors[chosen_neighbor]) {
        vector <string> reaction_compartments = bimolecular_reaction_details_map[reaction_pair.first]->allowed_compartments; // list of compartments where this reaction can take place
        if(find(reaction_compartments.begin(), reaction_compartments.end(), compartment_name) == reaction_compartments.end()) continue; // if current compartment isn't in list, do not react
        chosen_reactions.push_back(reaction_pair.first);
    }

    // Oct 3 2018 - these should take the relative pReact of each reaction and pick, since right now everything is equally likely, just randomly choose one
    int reaction_int = RandomNG::randInt(0,chosen_reactions.size()-1);
    string chosen_reaction = chosen_reactions[reaction_int];

    cout << "Timestep " << getCycles() << ": reaction " << chosen_reaction << " " << mol_container_0->container_name << " + " << chosen_neighbor->container_name << endl;
    do_reaction(mol_container_0, chosen_neighbor, chosen_reaction);    
}

void Simulation::do_reaction(Molecule_Container * mol_0, Molecule_Container * neighbor, string reaction_name) { // for standard bi-molecular interactions
    string reaction_type = bimolecular_reaction_details_map[reaction_name]->get_interaction_type();
    string placehold_string = "@@";
    map<int, vector <string>> reactant_names = bimolecular_reaction_details_map[reaction_name]->reactants_list; // reactant order follows reaction rule
    map<int, vector <string>> product_names = bimolecular_reaction_details_map[reaction_name]->products_list;

    // pick 1 set of reactants from origin and neighbor containers
    map<int, vector<bool>> mol_0_reactant_vectors = mol_0->bimolecular_reaction_map[0][reaction_name]; 
    vector<bool> mol_0_choice = mol_0_reactant_vectors[RandomNG::randInt(0, mol_0_reactant_vectors.size()-1)];
    map<int, vector<bool>> mol_1_reactant_vectors = neighbor->bimolecular_reaction_map[1][reaction_name];
    vector<bool> mol_1_choice = mol_1_reactant_vectors[RandomNG::randInt(0, mol_1_reactant_vectors.size()-1)];

    vector <Molecule *> react_0_container, react_1_container; // the species in the containers that are undergoing a reaction
    int mol_container_index = 0;
    for(auto element : mol_0_choice) {
        if(element) react_0_container.push_back(mol_0->container_vector[mol_container_index]);
        mol_container_index++;
    }
    mol_container_index = 0;
    for(auto element : mol_1_choice) {
        if(element) react_1_container.push_back(neighbor->container_vector[mol_container_index]);
        mol_container_index++;
    }

    if(reaction_type.compare("state_change") == 0) { // this includes modifier state as well as binding state
        bool in_order = false;
        if(reactant_names[0] == product_names[0] && reactant_names[1] == product_names[1]) {
            in_order = true;
        } else if(reactant_names[0] == product_names[1] && reactant_names[1] == product_names[0]) {
            in_order = false;
        } else {
            cout << "error in reaction " << reaction_name << " labelling." << endl;
            exit(1);
        }
        map<int, map<string, string>> prod_container_0_state_map;
        map<int, map<string, string>> prod_container_1_state_map;
        // depending on if the input reactant/product order was shuffled, adjust the maps accordingly
        if(in_order) {
            prod_container_0_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map[0];
            prod_container_1_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map[1];
        } else {
            prod_container_0_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map[1];
            prod_container_1_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map[0];
        }
        for(auto & mol : react_0_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator mol_pos = find_if(reactant_names[0].begin(), reactant_names[0].end(), [mol](string reactant) {
                return mol->molecule_name.compare(reactant) == 0;
            });
            auto pos = mol_pos - reactant_names[0].begin(); // this pos is the index of the matching reactant name in reactant_names
            Molecule_Container_Reactions::molecule_state_change(mol, prod_container_0_state_map[pos]);
            reactant_names[0].erase(mol_pos);
            prod_container_0_state_map.erase(pos); // erasing the same entry on reactant_names and prod_container should prevent it from being found next round, and also not change indexes
        }
        for(auto & mol : react_1_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator mol_pos = find_if(reactant_names[1].begin(), reactant_names[1].end(), [mol](string reactant) {
                return mol->molecule_name.compare(reactant) == 0;
            });
            auto pos = mol_pos - reactant_names[1].begin(); // this pos is the index of the matching reactant name in reactant_names
            Molecule_Container_Reactions::molecule_state_change(mol, prod_container_1_state_map[pos]);
            reactant_names[1].erase(mol_pos);
            prod_container_1_state_map.erase(pos); // erasing the same entry on reactant_names and prod_container should prevent it from being found next round, and also not change indexes
        }

        // use unbind parameter to shuffle product molecule location


        mol_0->initialize_container();
        neighbor->initialize_container();

        // If any products of reaction have an active transport destination, change the container status to allow for that movement
        if(!bimolecular_reaction_details_map[reaction_name]->product_container_destinations.empty()) {
            if(in_order) {
                for(auto const & prod_destination : bimolecular_reaction_details_map[reaction_name]->product_container_destinations) {
                    if(prod_destination.first == 0) {
                        mol_0->initialize_active_transport(prod_destination.second);
                    } else if(prod_destination.first == 1) {
                        neighbor->initialize_active_transport(prod_destination.second);
                    }
                }
            } else {
                for(auto const & prod_destination : bimolecular_reaction_details_map[reaction_name]->product_container_destinations) {
                    if(prod_destination.first == 0) {
                        neighbor->initialize_active_transport(prod_destination.second);
                    } else if(prod_destination.first == 1) {
                        mol_0->initialize_active_transport(prod_destination.second);
                    }
                }
            }
        }

    } else if(reaction_type.compare("fusion") == 0) { // assumes no species change
        map<Molecule *, pair<int, int>> mol_react_pos_map; // match each molecule to a int pair of original container:reactant index
        auto & reactant_bind_map = bimolecular_reaction_details_map[reaction_name]->reactants_binding_map;
        auto & reactant_state_map = bimolecular_reaction_details_map[reaction_name]->reactants_state_map;
        map<int, map<string, string>> & prod_container_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map[0];
        map<int, map<string, int>> & prod_container_bind_map = bimolecular_reaction_details_map[reaction_name]->products_binding_map[0];
        map<Molecule *, int> mol_to_prod_index_map; // mapping each mol to the product index position, not paired values since only 1 product container
        auto reactant_prod_map = bimolecular_reaction_details_map[reaction_name]->react_prod_map; // pair maps that matches reactant index to prod index
        vector<string> reactants_0 = reactant_names[0], reactants_1 = reactant_names[1], products = product_names[0];

        for(auto & mol : react_0_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(reactants_0.begin(), reactants_0.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == reactants_0.end()) continue;
            auto pos_reactant = reactant_mol_pos - reactants_0.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[0][pos_reactant], reactant_bind_map[0][pos_reactant])) continue;
            reactants_0[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            // combined_mol.push_back(mol);
            mol_react_pos_map[mol] = make_pair(0, pos_reactant);
        }
        for(auto & mol : react_1_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(reactants_1.begin(), reactants_1.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == reactants_1.end()) continue;
            auto pos_reactant = reactant_mol_pos - reactants_1.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[1][pos_reactant], reactant_bind_map[1][pos_reactant])) continue;

            reactants_1[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            // combined_mol.push_back(mol);
            mol_react_pos_map[mol] = make_pair(1, pos_reactant);
        }

        for(auto & mol_pos_pair : mol_react_pos_map) {  // change the state info of the reactants
            Molecule * reacting_mol = mol_pos_pair.first;
            pair<int, int> react_pos_pair = mol_pos_pair.second;
            mol_to_prod_index_map[reacting_mol] = reactant_prod_map[react_pos_pair].second; // since only 1 product container
            int prod_container_index = reactant_prod_map[react_pos_pair].first;
            int prod_pos_index = reactant_prod_map[react_pos_pair].second;
            
            map<string, string> current_new_states = prod_container_state_map[prod_pos_index];
            map<string, int> current_new_binding = prod_container_bind_map[prod_pos_index];
            Molecule_Container_Reactions::molecule_state_change(reacting_mol, current_new_states);
            Molecule_Container_Reactions::molecule_bind_change(reacting_mol, current_new_binding);
        }

        Molecule_Container * new_container = Molecule_Container_Reactions::fusion(mol_0, neighbor);
        containers_to_change_map.push_back(new_container);
        // cout << "added " << new_container.first << " " << new_container.second->container_id << " to map" << endl;
        // if the fusion product is actively transported, set that flag and find destination
        if(!bimolecular_reaction_details_map[reaction_name]->product_container_destinations.empty()) {
            new_container->initialize_active_transport(bimolecular_reaction_details_map[reaction_name]->product_container_destinations[0]);
        }

    } else if(reaction_type.compare("multi_product") == 0) { // no species change
        map<Molecule *, pair<int, int>> mol_react_pos_map; // match each molecule to a int pair of original container:reactant index
        auto & reactant_bind_map = bimolecular_reaction_details_map[reaction_name]->reactants_binding_map;
        auto & reactant_state_map = bimolecular_reaction_details_map[reaction_name]->reactants_state_map;
        auto & prod_bind_map = bimolecular_reaction_details_map[reaction_name]->products_binding_map;
        auto & prod_state_map = bimolecular_reaction_details_map[reaction_name]->products_state_map;
        vector<string> reactants_0 = reactant_names[0], reactants_1 = reactant_names[1];
        map<Molecule *, pair<int, int>> mol_to_prod_index_map; // mol pointer key, gives pair that indicates the container and the position inside that container
        auto & reactant_prod_map = bimolecular_reaction_details_map[reaction_name]->react_prod_map; // pair maps that matches reactant index to prod index
        // shuffle all reacting molecules into mol_react_pos_map, and also match each reacting molecule with a product
        for(auto & mol : react_0_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(reactants_0.begin(), reactants_0.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == reactants_0.end()) continue;
            auto pos_reactant = reactant_mol_pos - reactants_0.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[0][pos_reactant], reactant_bind_map[0][pos_reactant])) continue;
            reactants_0[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            mol_react_pos_map[mol] = make_pair(0, pos_reactant);
        }
        for(auto & mol : react_1_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(reactants_1.begin(), reactants_1.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == reactants_1.end()) continue;
            auto pos_reactant = reactant_mol_pos - reactants_1.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[1][pos_reactant], reactant_bind_map[1][pos_reactant])) continue;

            reactants_1[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            mol_react_pos_map[mol] = make_pair(1, pos_reactant);
        }

        // create mapping of reacting molecules to the product container placement
        for(auto & mol_pos_pair : mol_react_pos_map) { 
            Molecule * reacting_mol = mol_pos_pair.first;
            pair<int, int> react_pos_pair = mol_pos_pair.second;
            mol_to_prod_index_map[reacting_mol] = reactant_prod_map[react_pos_pair];
        }

        // when all species in the reacting containers are involved in the multi-product reaction
        if(reactants_0.size() == mol_0->container_vector.size() && reactants_1.size() == neighbor->container_vector.size()) {
            // handle the separation of molecules into new containers, and deletes the origin empty container
            // returns a vector of all the new containers to be inserted
            vector<Molecule_Container *> new_containers = Molecule_Container_Reactions::fission(mol_to_prod_index_map, bimolecular_reaction_details_map[reaction_name]);
            // change the state info of the reactants after containers are shifted
            for(auto & mol_pos_pair : mol_react_pos_map) {  
                Molecule * reacting_mol = mol_pos_pair.first;
                pair<int, int> react_pos_pair = mol_pos_pair.second;
                int prod_container_index = reactant_prod_map[react_pos_pair].first;
                int prod_pos_index = reactant_prod_map[react_pos_pair].second;
                
                map<string, string> current_new_states = prod_state_map[prod_container_index][prod_pos_index];
                map<string, int> current_new_binding = prod_bind_map[prod_container_index][prod_pos_index];
                Molecule_Container_Reactions::molecule_state_change(reacting_mol, current_new_states);
                Molecule_Container_Reactions::molecule_bind_change(reacting_mol, current_new_binding);
            }
            // if applicable, initiate active transport for products
            for(auto & container_dest_pair : bimolecular_reaction_details_map[reaction_name]->product_container_destinations) {
                new_containers[container_dest_pair.first]->initialize_active_transport(container_dest_pair.second);
            }
            new_containers_vector.insert(new_containers_vector.end(), new_containers.begin(), new_containers.end());

        } else if(reactants_0.size() == mol_0->container_vector.size()) {
            cerr << "not all molecules in reactant participating in fission reaction" << endl;
            exit(1);
        } else if(reactants_1.size() == neighbor->container_vector.size()) {
            cerr << "not all molecules in reactant participating in fission reaction" << endl;
            exit(1);
        } else {
            cerr << "not all molecules in reactant participating in fission reaction" << endl;
            exit(1);
        }

    } else if(reaction_type.compare("modify") == 0) {
        cout << numCycles << ": final reaction";
        cout << endl;
    }
}

void Simulation::do_reaction(Molecule_Container * mol_0, string reaction_name) { // overloaded for unimolecular reactions
    Cell_Reaction * reaction_details = unimolecular_reaction_details_map[reaction_name];
    string reaction_type = reaction_details->get_interaction_type();
    string placehold_string = "@@";
    vector <string> reactant_names = reaction_details->reactants_list[0]; // reactant order follows reaction rule
    map<int, vector <string>> product_names = reaction_details->products_list;

    // pick 1 set of reactants from origin and neighbor containers
    map<int, vector<bool>> mol_0_reactant_vectors = mol_0->unimolecular_reaction_map[reaction_name]; 
    vector<bool> mol_0_choice = mol_0_reactant_vectors[RandomNG::randInt(0, mol_0_reactant_vectors.size()-1)];

    vector <Molecule *> react_container;
    int mol_container_index = 0;
    for(auto element : mol_0_choice) {
        if(element) react_container.push_back(mol_0->container_vector[mol_container_index]);
        mol_container_index++;
    }

    if(reaction_type.compare("state_change") == 0) { // this includes modifier state as well as binding state
        map<int, map<string, string>> prod_container_0_state_map = reaction_details->products_state_map[0];
        auto & reactant_state_map = reaction_details->reactants_state_map[0];
        auto & reactant_bind_map = reaction_details->reactants_binding_map[0];

        for(auto & mol : react_container) {
            //  find the position of the current molecule in the reactant list's container
            //  first try to match the names, then check to see if the states match as well
            vector<string>::iterator mol_pos = find_if(reactant_names.begin(), reactant_names.end(), [mol](string reactant) {
                return mol->molecule_name.compare(reactant) == 0;
            });
            auto pos = mol_pos - reactant_names.begin(); // this pos is the index of the matching reactant name in reactant_names
            auto mol_binding_map = mol->mol_config.molecule_binding_map;
            auto mol_state_map = mol->mol_config.molecule_state_map;
            bool paired = false;
            for(int i = 0; i < reaction_details->reactants_list.size(); i++) {
                if(includes(mol_binding_map.begin(),mol_binding_map.end(),reactant_bind_map[i].begin(),reactant_bind_map[i].end())) {
                if(includes(mol_state_map.begin(),mol_state_map.end(),reactant_state_map[i].begin(),reactant_state_map[i].end())) {
                    paired = true;
                    break;
                }
                }
            }
            if(!paired) continue;

            Molecule_Container_Reactions::molecule_state_change(mol, prod_container_0_state_map[pos]);
            reactant_names.erase(mol_pos);
            prod_container_0_state_map.erase(pos); // erasing the same entry on reactant_names and prod_container should prevent it from being found next round, and also not change indexes
        }

        // optional: if unbind parameter exists, place molecule somewhere else
        // use unbind parameter to shuffle product molecule location
        // this is code from fission in reactions to separate reversible products
        if(reaction_details->unbind_rad_min != 0) {
            // flag for whether or not this is a membrane reaction
            pair<bool, Point> membrane_reaction = make_pair(false, Point(0, 0, 0));
            if(mol_0->is_in_membrane()) {
                membrane_reaction.first = true;
                membrane_reaction.second = mol_0->get_position();
            }
            double rand_radius_min = reaction_details->unbind_rad_min;
            double rand_radius_range = reaction_details->unbind_rad_max - rand_radius_min;

            Point old_loc = mol_0->pos;
            Point rand_loc;
            // instead of using RandomNG float generator multiple times, create distribution once and sample from it
            uniform_real_distribution<float> radius_range_dist(-rand_radius_range, rand_radius_range);
            rand_loc.x = radius_range_dist(RandomNG::get_generator());
            rand_loc.x = rand_loc.x > 0 ? old_loc.x + rand_radius_range + rand_radius_min : old_loc.x - rand_radius_range - rand_radius_min;
            rand_loc.y = radius_range_dist(RandomNG::get_generator());
            rand_loc.y = rand_loc.y > 0 ? old_loc.y + rand_radius_range + rand_radius_min : old_loc.y - rand_radius_range - rand_radius_min;
            rand_loc.z = radius_range_dist(RandomNG::get_generator());
            rand_loc.z = rand_loc.z > 0 ? old_loc.z + rand_radius_range + rand_radius_min : old_loc.z - rand_radius_range - rand_radius_min;

            string new_env = lEnv->get_compart_name(rand_loc);
            int counter = 0;
            while(!lEnv->isCompatibleEnvironment(mol_0->container_name, new_env) && counter < 10) {
                rand_loc.x = radius_range_dist(RandomNG::get_generator());
                rand_loc.x = rand_loc.x > 0 ? old_loc.x + rand_radius_range + rand_radius_min : old_loc.x - rand_radius_range - rand_radius_min;
                rand_loc.y = radius_range_dist(RandomNG::get_generator());
                rand_loc.y = rand_loc.y > 0 ? old_loc.y + rand_radius_range + rand_radius_min : old_loc.y - rand_radius_range - rand_radius_min;
                rand_loc.z = radius_range_dist(RandomNG::get_generator());
                rand_loc.z = rand_loc.z > 0 ? old_loc.z + rand_radius_range + rand_radius_min : old_loc.z - rand_radius_range - rand_radius_min;
                new_env = lEnv->get_compart_name(rand_loc);
                counter++;
            }
            // if we can't find a good destination around the radius, find closest voxel that is a valid compartment and drop it in there
            if(counter >= 10) {
                vector <string> allowed_comparts = lEnv->get_compatible_env(mol_0->container_name);
                vector <Point *> list_of_sites;
                for(auto & compartment : allowed_comparts) {
                    vector <Point *> compart_sites = lEnv->getCompartmentSites(compartment);
                    list_of_sites.insert(list_of_sites.end(), compart_sites.begin(), compart_sites.end());
                }
                float distance = 100;
                Vector3D diff;
                Point dest_voxel = *(list_of_sites.front());
                for(Point* location : list_of_sites) {
                    diff.setDiff(*location, dest_voxel);
                    if(diff.magnitude() < distance) {
                        dest_voxel = *location;
                        distance = diff.magnitude();
                    }
                }
                if(!lEnv->isCompatibleEnvironment(mol_0->container_name, lEnv->get_compart_name(dest_voxel))) {
                    cerr << "error when attempting to place molecule in destination voxel." << endl;
                    exit(1);
                }
                if(membrane_reaction.first && lEnv->isCompatibleEnvironment(mol_0->container_name, lEnv->get_compart_name(membrane_reaction.second))) {
                    LatticeSite * membrane_voxel = lEnv->getLatticeSite(membrane_reaction.second);
                    if(membrane_voxel->locked_axis.first == "x") dest_voxel.x = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "y") dest_voxel.y = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "z") dest_voxel.z = membrane_voxel->locked_axis.second;
                }
                mol_0->remove_from_maps(old_loc);
                mol_0->set_position(dest_voxel);
                mol_0->add_to_maps();
            } else {
                if(membrane_reaction.first && lEnv->isCompatibleEnvironment(mol_0->container_name, lEnv->get_compart_name(membrane_reaction.second))) {
                    LatticeSite * membrane_voxel = lEnv->getLatticeSite(membrane_reaction.second);
                    if(membrane_voxel->locked_axis.first == "x") rand_loc.x = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "y") rand_loc.y = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "z") rand_loc.z = membrane_voxel->locked_axis.second;
                }
                mol_0->remove_from_maps(old_loc);
                mol_0->set_position(rand_loc);
                mol_0->add_to_maps();
            }
        }
        
        mol_0->initialize_container();

        // If any products of reaction have an active transport destination, change the container status to allow for that movement
        if(!reaction_details->product_container_destinations.empty()) {
            for(auto const & prod_destination : reaction_details->product_container_destinations) {
                if(prod_destination.first == 0) {
                    mol_0->initialize_active_transport(prod_destination.second);
                }
            }
        }
        
    } else if(reaction_type.compare("fission") == 0) {
        map<Molecule *, int> mol_react_pos_map; // match each molecule to a int pair of original container:reactant index
        vector<Molecule *> combined_mol;
        auto & reactant_bind_map = reaction_details->reactants_binding_map;
        auto & reactant_state_map = reaction_details->reactants_state_map;
        auto & prod_bind_map = reaction_details->products_binding_map;
        auto & prod_state_map = reaction_details->products_state_map;
        vector<string> reactants_0 = reactant_names;
        map<Molecule *, pair<int, int>> mol_to_prod_index_map; // mol pointer key, gives pair that indicates the container and the position inside that container of the products
        auto & reactant_prod_map = reaction_details->react_prod_map; // pair maps that matches reactant index to prod index

        // shuffle all reacting molecules into mol_react_pos_map, and also match each reacting molecule with a product
        for(auto & mol : react_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(reactants_0.begin(), reactants_0.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == reactants_0.end()) continue;
            auto pos_reactant = reactant_mol_pos - reactants_0.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[0][pos_reactant], reactant_bind_map[0][pos_reactant])) continue;

            reactants_0[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            mol_react_pos_map[mol] = pos_reactant;
        }

        // create mapping of reacting molecules to the product container placement
        for(auto & mol_pos_pair : mol_react_pos_map) {
            Molecule * reacting_mol = mol_pos_pair.first;
            pair<int, int> react_pos = make_pair(0, mol_pos_pair.second);
            mol_to_prod_index_map[reacting_mol] = reactant_prod_map[react_pos];
        }

        // handle the separation of molecules into new containers, and deletes the origin empty container
        // returns a vector of all the new containers to be inserted
        vector<Molecule_Container *> new_containers = Molecule_Container_Reactions::fission(mol_to_prod_index_map, reaction_details); 
        // change the state info of the reactants after containers are shifted
        for(auto & mol_pos_pair : mol_react_pos_map) {
            Molecule * reacting_mol = mol_pos_pair.first;
            pair<int, int> react_pos = make_pair(0, mol_pos_pair.second);
            int prod_container_index = reactant_prod_map[react_pos].first;
            int prod_pos_index = reactant_prod_map[react_pos].second;
            
            map<string, string> current_new_states = prod_state_map[prod_container_index][prod_pos_index];
            map<string, int> current_new_binding = prod_bind_map[prod_container_index][prod_pos_index];
            Molecule_Container_Reactions::molecule_state_change(reacting_mol, current_new_states);
            Molecule_Container_Reactions::molecule_bind_change(reacting_mol, current_new_binding);
        }
        for(auto & container_dest_pair : reaction_details->product_container_destinations) {
            new_containers[container_dest_pair.first]->initialize_active_transport(container_dest_pair.second);
        }
        new_containers_vector.insert(new_containers_vector.end(), new_containers.begin(), new_containers.end());
    }
}

void Simulation::do_reaction(Molecule_Container * mol_0, string bulk_name, string reaction_name) { // overloaded for bulk_particle interactions
    string reaction_type = bulk_reaction_details_map[reaction_name]->get_interaction_type();
    int bulk_index = bulk_reaction_details_map[reaction_name]->bulk_index; // position of bulk molecule in reactant (0 or 1)
    int container_index = abs(1 - bulk_index); // position of the container that reacts with bulk mol (0 or 1)
    LatticeSite * current_voxel = lEnv->getLatticeSite(mol_0->get_position());
    string placehold_string = "@@";
    vector <string> particle_reactant_names = bulk_reaction_details_map[reaction_name]->reactants_list[container_index]; // reactant order follows reaction rule
    // pick 1 set of reactants from origin and neighbor containers
    map<int, vector<bool>> mol_0_reactant_vectors = mol_0->bulk_interaction_map[reaction_name]; 
    vector<bool> mol_0_choice = mol_0_reactant_vectors[RandomNG::randInt(0, mol_0_reactant_vectors.size()-1)];

    vector <Molecule *> react_container;
    int mol_container_index = 0;
    for(auto element : mol_0_choice) {
        if(element) react_container.push_back(mol_0->container_vector[mol_container_index]);
        mol_container_index++;
    }

    if(reaction_type.compare("fusion") == 0) { // assumes no species change
        map<Molecule *, int> mol_react_pos_map; // match each molecule to an int reactant index
        auto & reactant_bind_map = bulk_reaction_details_map[reaction_name]->reactants_binding_map[container_index];
        auto & reactant_state_map = bulk_reaction_details_map[reaction_name]->reactants_state_map[container_index];
        map<int, map<string, string>> & prod_container_state_map = bulk_reaction_details_map[reaction_name]->products_state_map[0];
        map<int, map<string, int>> & prod_container_bind_map = bulk_reaction_details_map[reaction_name]->products_binding_map[0];
        auto & reactant_prod_map = bulk_reaction_details_map[reaction_name]->react_prod_map; // pair maps that matches reactant index to prod index

        // fill the mol_react_pos_map, which maps position of reacting molecules to the respective products
        for(auto & mol : react_container) {
            // find the position of the current molecule in the reactant list's container
            vector<string>::iterator reactant_mol_pos = find_if(particle_reactant_names.begin(), particle_reactant_names.end(), [mol](string reactant) {
                return (mol->molecule_name.compare(reactant) == 0);
            });
            if(reactant_mol_pos == particle_reactant_names.end()) continue;
            auto pos_reactant = reactant_mol_pos - particle_reactant_names.begin(); // this pos is the index of the matching reactant name in reactant_names
            if(!Molecule_Container_Reactions::molecule_fills_rule(mol, reactant_state_map[pos_reactant], reactant_bind_map[pos_reactant])) continue;
            particle_reactant_names[pos_reactant] = placehold_string; // placeholder symbol to overwrite already-found elements
            // combined_mol.push_back(mol);
            mol_react_pos_map[mol] = pos_reactant;
        }

        // iterate through each reactant molecule inside container and change the state info
        for(auto & mol_pos_pair : mol_react_pos_map) { 
            Molecule * reacting_mol = mol_pos_pair.first;
            pair<int,int> react_pos = make_pair(container_index,mol_pos_pair.second);
            int prod_pos_index = reactant_prod_map[react_pos].second;
            
            map<string, string> current_new_states = prod_container_state_map[prod_pos_index];
            map<string, int> current_new_binding = prod_container_bind_map[prod_pos_index];
            Molecule_Container_Reactions::molecule_state_change(reacting_mol, current_new_states);
            Molecule_Container_Reactions::molecule_bind_change(reacting_mol, current_new_binding);
        }

        // once states of reactants are done, create new molecule object from bulk reactant
        // remove corresponding amount of bulk from the bulk map of this voxel
        // if current voxel doesnt have enough (< 1) then take from neighbors too
        if(current_voxel->bulk_moles_map[bulk_name] < lEnv->particle_to_mole(1)) {
            double remainder_moles = lEnv->particle_to_mole(1) - current_voxel->bulk_moles_map[bulk_name];
            current_voxel->add_mole_to_voxel(bulk_name, -current_voxel->bulk_moles_map[bulk_name]);
            for(auto & neighbor : current_voxel->neighbors[bulk_name][0]) {
                if(remainder_moles > neighbor->bulk_moles_map[bulk_name]) {
                    remainder_moles -= neighbor->bulk_moles_map[bulk_name];
                    neighbor->add_mole_to_voxel(bulk_name, -neighbor->bulk_moles_map[bulk_name]);
                } else {
                    neighbor->add_mole_to_voxel(bulk_name, -remainder_moles);
                    remainder_moles = 0;
                }
                if(remainder_moles == 0) break;
            }
            if(remainder_moles != 0) {
                cerr << "expanded bulk_mol removal in bulk reaction failed." << endl;
                exit(1);
            }
        } else {
            current_voxel->add_mole_to_voxel(bulk_name, -lEnv->particle_to_mole(1));
        }
        Molecule * new_bulk_molecule = new Molecule(bulk_name);

        Molecule_Container * new_container = Molecule_Container_Reactions::bulk_fusion(mol_0, new_bulk_molecule);
        containers_to_change_map.push_back(new_container);
        // if the fusion product is actively transported, set that flag and find destination
        if(!bulk_reaction_details_map[reaction_name]->product_container_destinations.empty()) {
            new_container->initialize_active_transport(bulk_reaction_details_map[reaction_name]->product_container_destinations[0]);
        }

    }
}


//////////////////////////////////////////////////////////////////////
//  Build a list of possible interacting molecules using the container's internal reaction list
//  returns a "map" (really a list of pairs) with containers as keys, and information about how many times (possibilities, combinations) this container can do that reaction
map<Molecule_Container *, map<string, int>> Simulation::buildNeighborList(Molecule_Container * mol_container_0){

    // Initialize a list that holds containers neighboring mol_container_0
    // vector <Molecule_Container *> neighborList;
    map<Molecule_Container *, map<string, int>> neighborList;

    // Get the list of possible interacting molecules that are near
    LatticeSite * ls = lEnv->getLatticeSite(mol_container_0->get_position());

    unordered_map<string, unordered_map<unsigned long, Molecule_Container *>> latt_map = ls->lattice_container_map; // all containers in current c-voxel
    
    // list of reactions the origin molecule can do, and the number of times they happen
    // this map is filled at container creation, and re-created every time it goes through a reaction
    map<string, map<int, vector<bool>>> origin_reaction_node_map = mol_container_0->bimolecular_reaction_map[0];
    vector <string> origin_reaction_names;
    for(auto & reaction : origin_reaction_node_map) origin_reaction_names.push_back(reaction.first);

    // construct interacting_neighbors, a list of container names that share 1 or more reactions with the origin container
    set<string> interacting_neighbors; // list of container names that hold valid interaction partners for mol_container_0
    for(auto & latt_container_itr : latt_map) { //iterate neighbor container names
        string container_name = latt_container_itr.first;
        for(auto & reaction : origin_reaction_node_map) {
            string origin_reaction = reaction.first;
            if(Molecule_Container::is_in_reaction(container_name,origin_reaction,1)) {
                interacting_neighbors.insert(container_name);
                break;
            }
        }
    }

    // if this is a membrane molecule, add adjacent voxel molecules to neighbors
    if(mol_container_0->is_in_membrane()) {
        Point vox_pos = mol_container_0->get_position();
        Compart & mem_compart = compartments_map[lEnv->get_compart_name(vox_pos)];
        // adjust target adjacent voxel depending on which way the membrane is facing and which axis it is aligned with
        if(mem_compart.face == "front") {
            if(mem_compart.axis == "x") vox_pos.x -= 0.5;
            else if(mem_compart.axis == "y") vox_pos.y -= 0.5;
            else if(mem_compart.axis == "z") vox_pos.z -= 0.5;
        } else if(mem_compart.face == "back") {
            if(mem_compart.axis == "x") vox_pos.x += 0.5;
            else if(mem_compart.axis == "y") vox_pos.y += 0.5;
            else if(mem_compart.axis == "z") vox_pos.z += 0.5;
        } else {
            cerr << "membrane compartment definition error." << endl;
            exit(1);
        }
        // add the molecules of the adjacent voxel into the current latt_map
        unordered_map<string, unordered_map<unsigned long, Molecule_Container *>> adj_voxel_latt_map = lEnv->getLatticeSite(vox_pos)->lattice_container_map;
        for(auto & container_name_map_pair : adj_voxel_latt_map) {
            string container_name = container_name_map_pair.first;
            // if the new container type already exists in lattice map, insert new collection of the same containers to map
            // otherwise, add it as a new type of container in the lattice map and also register them as valid neighbors if reactions are shared
            if(latt_map.find(container_name) == latt_map.end()) {
                latt_map[container_name] = container_name_map_pair.second;
                for(auto & reaction : origin_reaction_node_map) {
                    string origin_reaction = reaction.first;
                    if(Molecule_Container::is_in_reaction(container_name,origin_reaction,0) || Molecule_Container::is_in_reaction(container_name,origin_reaction,1)) {
                        interacting_neighbors.insert(container_name);
                        break;
                    }
                }
            } else {
                latt_map[container_name].insert(container_name_map_pair.second.begin(), container_name_map_pair.second.end());
            }
        }
    }

    // iterate over only the neighboring containers that share reactions
    for(auto & neighbor_name : interacting_neighbors) {
        // neighbors organized by container name, this is a copy to manipulate (delete itself from the neighbor selection)
        unordered_map <unsigned long, Molecule_Container *> container_list = latt_map[neighbor_name];
        if(neighbor_name == mol_container_0->container_name) {
            //delete itself from latt_container_map so it isn't its own neighbor
            if(container_list.find(mol_container_0->container_id) != container_list.end()) {
                container_list.erase(mol_container_0->container_id);
            }
        }

        for(auto & container_pair : container_list) { // iterate through all containers of that name
            // unsigned long container_key = container_pair.first;
            Molecule_Container * current_container = Mol_Container_Table[neighbor_name][container_pair.second->container_id];
            // ignore molecules that are not permitted to react, do not place them into neighbor list
            if(!current_container->allowed_to_react) continue;
            // list of reactions the neighbor molecule can do, and the number of times they happen
            map<string, map<int, vector<bool>>> neighbor_reaction_node_map = current_container->bimolecular_reaction_map[1];
            vector <string> neighbor_reaction_names;
            for(auto & reaction : neighbor_reaction_node_map) neighbor_reaction_names.push_back(reaction.first);

            // intersecting the two lists will give the reactions that both are qualified for
            vector<string> common_reactions;
            set_intersection(origin_reaction_names.begin(), origin_reaction_names.end(), neighbor_reaction_names.begin(), neighbor_reaction_names.end(), back_inserter(common_reactions));

            for(auto & reaction : common_reactions) {
                neighborList[current_container][reaction] = min(origin_reaction_node_map[reaction].size(), neighbor_reaction_node_map[reaction].size());
            }
        }
    }

    // a list of all containers that share reactions with origin container within the same c-voxel that have the correct modifier states
    // container pointer: neighbor's reaction_map
    return neighborList;
}


//////////////////////////////////////////////////////////////////////
//  Return the current number of cycles executed in the simulation
//  ->    A cycle is the basic time unit in the simulation
//
int Simulation :: getCycles(){
    // Return the current number of cycles
    return numCycles;
}

// //////////////////////////////////////////////////////////////////////
// //  Get number of molecules with the specified type in the Simulation
// //
// int Simulation :: getnumContainers(string molecule_name){
//     // Return the number of molecules of the specified type in the Simulation Space
//     return Mol_Container_Table[molecule_name].size();
// }



//////////////////////////////////////////////////////////////////////
//  Is simulation completed
//  ->    Is the termination condition reached?
//
bool Simulation :: isSimulationCompleted(){
    // Return whether the simulation is completed
    return simulationCompleted;
}



//////////////////////////////////////////////////////////////////////
//  Get access to the Molecule table
//
unordered_map <string, unordered_map <unsigned long, Molecule_Container *> > & Simulation :: getMoleculeTable(){
    // Return a pointer to the Molecule table that manages all
    // existing molecules in the Simulation Space
    return Mol_Container_Table;
}


//////////////////////////////////////////////////////////////////////
//  Perform periodic modifications to the Simulation Space
//  ->    Temporal events updating behaviour/Properties of the Simulation
//
void Simulation :: updateSimulation(){

}



// recursively adds all events in list (if fulfills conditions) into the event queue
void Simulation::recursive_traverse_event_list(vector<Event_Parameters *> event_vec) {
    map<int, vector<Event_Parameters *>> & coming_events_queue = ParameterManager::getParameterManager()->event_queue;

    while(!event_vec.empty()) {
        for(auto triggered_event_itr = event_vec.begin(), next_itr = triggered_event_itr; triggered_event_itr != event_vec.end(); triggered_event_itr = next_itr) {
            ++next_itr;
            auto triggered_event = *triggered_event_itr;
            bool skip = false;
            if(triggered_event->event_done) skip = true; // check if event is already done
            if(RandomNG::randFloat(0, 1) >= triggered_event->event_probability) skip = true; // check if event will occur, prob=1 will always happen

            // if make event go, add to queue, then check that one's event list
            if(!skip) {
                coming_events_queue[numCycles + triggered_event->event_interval].push_back(triggered_event);
                recursive_traverse_event_list(triggered_event->events_to_trigger);
            } 
            event_vec.erase(triggered_event_itr);
            if(event_vec.empty()) break;
        }
    }
}

map<unsigned long, Molecule_Container *> Simulation::get_list_of_event_containers(Container_Config config, vector<Point> & voxel_locations) {

}

void Simulation::check_events() {
    map<string, Event_Parameters *> & PM_event_map = ParameterManager::getParameterManager()->event_map;
    map<int, vector<Event_Parameters *>> & coming_events_queue = ParameterManager::getParameterManager()->event_queue;

    // check event dependencies of this cycle, checking only independent events
    for(auto & event : PM_event_map) {
        Event_Parameters * & event_params = event.second;
        if(event_params->event_done) continue; // check if event is already done
        if(RandomNG::randFloat(0, 1) >= event_params->event_probability) continue; // check if event will occur, prob=1 will always happen
        // first, see if time triggered
        if(event_params->event_trigger == "time") {
            if(event_params->event_initial_time > numCycles) continue;
            // if it's after the final time, don't trigger event again
            if(event_params->event_final_time < numCycles && event_params->event_final_time != -1) {
                event_params->event_done = true;
                continue;
            }

            // if trigger_calc is 0 (from hitting interval of repeatable events, or if it is the initial time), do event
            int trigger_calc;
            if(event_params->event_interval == 0) trigger_calc = numCycles - event_params->event_initial_time;
            else if(event_params->event_interval > 0) trigger_calc = (numCycles - event_params->event_initial_time) % event_params->event_interval;
            if(trigger_calc == 0) {
                // do event
                do_event(event_params);
                recursive_traverse_event_list(event_params->events_to_trigger);
                if(!event_params->event_repeat) {
                    event_params->event_done = true;
                }
            }
            
        } else if(event_params->event_trigger == "event") {
            continue;

        } else if(event_params->event_trigger == "state") {
            auto & species_name_attr_map = ParameterManager::getParameterManager()->list_of_species_map;
            Event_Parameters * & event_params = event.second;
            if(event_params->event_done) continue; // check if event is already done
            // for intervals
            int trigger_calc;
            if(event_params->event_interval == 0) trigger_calc = numCycles - event_params->event_initial_time;
            else if(event_params->event_interval > 0) trigger_calc = (numCycles - event_params->event_initial_time) % event_params->event_interval;
            if(trigger_calc != 0) continue;
            
            // if it's after the final time, don't trigger event again
            if(event_params->event_final_time < numCycles && event_params->event_final_time != -1) {
                event_params->event_done = true;
                continue;
            }
            if(RandomNG::randFloat(0, 1) >= event_params->event_probability) continue; // check if event will occur, prob=1 will always happen

            // determine container name, and whether the triggering mol is a bulk or container
            auto species_attr_vec = species_name_attr_map[event_params->state_trigger_molecule];
            string trigger_mol_name;
            Container_Config event_container_config;
            bool trigger_mol_is_bulk = false;
            if(species_attr_vec.size() == 1 && species_attr_vec[0]->is_simple_molecule()) {
                trigger_mol_name = species_attr_vec[0]->ID;
                trigger_mol_is_bulk = true;
            } else {
                // figure out what container is checked as the trigger
                for(auto & config_rule_pair : container_rules_map) {
                    if(event_params->state_trigger_molecule == config_rule_pair.first.config_name) {
                        trigger_mol_name = config_rule_pair.first.container_name;
                        event_container_config = config_rule_pair.first;
                        break;
                    }
                }
                if(trigger_mol_name.empty()) {
                    cerr << "can't find state trigger containers. exiting" << endl;
                    exit(1);
                }
                trigger_mol_is_bulk = false;
            }

            // check locations and sum up matching molecules to see if event threshold is met
            int threshold_amount = event_params->state_trigger_mol_amount;
            double mol_counter = 0;
            vector<Point> & list_of_location_voxels = event_params->state_trigger_voxels;

            if(trigger_mol_is_bulk) {
                for(auto & compart_voxel_coord : list_of_location_voxels) {
                    LatticeSite * ls = lEnv->getLatticeSite(compart_voxel_coord);
                    mol_counter += lEnv->mole_to_particle(ls->bulk_moles_map[trigger_mol_name]);
                }

            } else {
                map<unsigned long, Molecule_Container *> event_trigger_mol_map;
                for(auto & compart_voxel_coord : list_of_location_voxels) {
                    LatticeSite * ls = lEnv->getLatticeSite(compart_voxel_coord);
                    for(auto & container : ls->lattice_container_map[trigger_mol_name]) {
                        // if the container matches the config specified in the event info, add to grand_totals
                        if(container.second->container_config == event_container_config) {
                            mol_counter++;
                        }
                    }
                }
            }
            // depending on if the state trigger is greater or less than threshold
            if(event_params->event_triggering_condition == 1) {
                if(mol_counter > threshold_amount) {
                    // do event
                    do_event(event_params);
                    recursive_traverse_event_list(event_params->events_to_trigger);
                    if(!event_params->event_repeat) {
                        event_params->event_done = true;
                    }
                }
            } else if(event_params->event_triggering_condition == 0) {
                if(mol_counter < threshold_amount) {
                    // do event
                    do_event(event_params);
                    recursive_traverse_event_list(event_params->events_to_trigger);
                    if(!event_params->event_repeat) {
                        event_params->event_done = true;
                    }
                }
            } else {
                cerr << "event condition somehow undefined. exiting" << endl;
                exit(1);
            }

        }
    }

    // do any queue'd up events
    if(coming_events_queue.find(numCycles) != coming_events_queue.end()) {
        vector<Event_Parameters *> event_list = coming_events_queue[numCycles];
        for(auto & event : event_list) {
            // do event
            do_event(event);
            if(!event->event_repeat) {
                event->event_done = true;
            }
        }
        // clean up after it's done
        coming_events_queue.erase(numCycles);
    }
}

void Simulation::do_event(Event_Parameters * params) {
    string event_type = params->event_type;
    auto & species_name_attr_map = ParameterManager::getParameterManager()->list_of_species_map;
    // make sure the event species actually exists in the simulation and defined within listOfSpecies in XML
    if(species_name_attr_map.find(params->event_container_species_name) == species_name_attr_map.end()) exit(1);

    if(event_type == "add_mols") {
        vector<Species_Attributes *> & attr_vec = species_name_attr_map[params->event_container_species_name];
        vector<pair<Point, int>> & container_placements = params->event_molecule_placement_pairs;
        for(auto & location_pair : container_placements) {
            auto & amount = location_pair.second;
            auto & location = location_pair.first;
            // check if generating bulk or containers
            if(attr_vec.size() == 1 && attr_vec[0]->is_simple_molecule()) {
                lEnv->getLatticeSite(location)->add_mole_to_voxel(attr_vec[0]->ID, lEnv->particle_to_mole(amount));
                cout << "added " << amount << " of " << attr_vec[0]->ID << " to " << location.x << "," << location.y << "," << location.z << endl;
            } else {
                // number of particles, attributes, position of voxel
                createMoleculeContainers(amount, attr_vec, location);
            }
        }

    } else if(event_type == "remove_mols") {
        vector<Species_Attributes *> & attr_vec = species_name_attr_map[params->event_container_species_name];
        // for bulk molecule removal
        // tries to remove molecules equally across all designated voxels
        if(attr_vec.size() == 1 && attr_vec[0]->is_simple_molecule()) {
            string bulk_name = attr_vec[0]->ID;
            double total_bulk_particle_amount = 0;
            if(!params->event_compartment.empty() && compartments_map[params->event_compartment].list_of_voxels.size() == params->event_molecule_placement_pairs.size()) {
                total_bulk_particle_amount = compartments_map[params->event_compartment].current_bulk_count[bulk_name];
            } else {
                for(auto & voxel_location : params->event_placement_voxels) {
                    total_bulk_particle_amount += lEnv->mole_to_particle(lEnv->getLatticeSite(voxel_location)->bulk_moles_map[bulk_name]);
                }
            }
            double amount_to_remove = params->event_container_amount_total;
            if(amount_to_remove > total_bulk_particle_amount) amount_to_remove = total_bulk_particle_amount;
            double proportion_to_remove = double(amount_to_remove)/double(total_bulk_particle_amount);
            for(auto & voxel_location : params->event_placement_voxels) {
                double mole_to_remove = lEnv->getLatticeSite(voxel_location)->bulk_moles_map[bulk_name] * proportion_to_remove;
                double particles_to_remove = lEnv->mole_to_particle(mole_to_remove);
                lEnv->getLatticeSite(voxel_location)->add_mole_to_voxel(bulk_name, -mole_to_remove);
            }
            cout << "removed " << amount_to_remove << " of " << bulk_name << endl;
            
        } else { // for molecule containers
            // picks random eligible containers within the area
            map<unsigned long, Molecule_Container *> grand_total_container_map;
            string container_name;
            Container_Config event_container_config;
            // figure out what container is being removed
            for(auto & config_rule_pair : container_rules_map) {
                if(params->event_container_species_name == config_rule_pair.first.config_name) {
                    container_name = config_rule_pair.first.container_name;
                    event_container_config = config_rule_pair.first;
                    break;
                }
            }
            if(container_name.empty()) {
                cerr << "can't find container to remove. exiting" << endl;
                exit(1);
            }
            vector<Point> & list_of_location_voxels = params->event_placement_voxels;
            int total_amount = params->event_container_amount_total;

            // fill up Grand_Total container map by going through the containers in each relevant voxel
            for(auto & compart_voxel_coord : list_of_location_voxels) {
                LatticeSite * ls = lEnv->getLatticeSite(compart_voxel_coord);
                for(auto & container : ls->lattice_container_map[container_name]) {
                    // if the container matches the config specified in the event info, add to grand_totals
                    if(container.second->container_config == event_container_config) {
                        grand_total_container_map[container.second->container_id] = container.second;
                    }
                }
            }
            if(total_amount > grand_total_container_map.size()) total_amount = grand_total_container_map.size();
            
            // randomly select the molecules to delete by accessing container map using a shuffled index
            vector<int> shuffled_indices = RandomNG::create_shuffled_index_vector(grand_total_container_map.size());
            int counter = 0;
            for(auto & index : shuffled_indices) {
                if(counter == total_amount) break;
                auto itAtOffset = std::next(grand_total_container_map.begin(), index);
                delete itAtOffset->second;
                counter++;
            }
        }

    } else if(event_type == "transport_mols") {
        vector<Species_Attributes *> & attr_vec = species_name_attr_map[params->event_container_species_name];
        // for bulk molecule removal
        if(attr_vec.size() == 1 && attr_vec[0]->is_simple_molecule()) {
            cerr << "transport events not allowed for bulk molecules. exiting" << endl;
            exit(1);
        } else { // for molecule containers
            map<unsigned long, Molecule_Container *> grand_total_container_map; // map of containers that currently fulfills requirements in applicable space which can be moved
            string container_name;
            Container_Config event_container_config;
            // figure out what container is being transported
            for(auto & config_rule_pair : container_rules_map) {
                if(params->event_container_species_name == config_rule_pair.first.config_name) {
                    container_name = config_rule_pair.first.container_name;
                    event_container_config = config_rule_pair.first;
                    break;
                }
            }
            if(container_name.empty()) {
                cerr << "can't find container to remove. exiting" << endl;
                exit(1);
            }
            vector<Point> & list_of_location_voxels = params->event_placement_voxels; // from the "location" node in XML
            int total_amount = params->event_container_amount_total;

            // fill up Grand_Total container map by going through the containers in each relevant voxel
            for(auto & compart_voxel_coord : list_of_location_voxels) {
                LatticeSite * ls = lEnv->getLatticeSite(compart_voxel_coord);
                for(auto & container : ls->lattice_container_map[container_name]) {
                    // if the container matches the config specified in the event info, add to grand_totals
                    if(container.second->container_config == event_container_config) {
                        grand_total_container_map[container.second->container_id] = container.second;
                    }
                }
            }
            if(total_amount > grand_total_container_map.size()) total_amount = grand_total_container_map.size();
            
            // randomly select the molecules to delete by accessing container map using a shuffled index of the container map
            vector<int> shuffled_indices = RandomNG::create_shuffled_index_vector(grand_total_container_map.size());
            // create vector of sampled ints that can access index of potential target voxel destinations, size of vec is number of molecules being transported
            vector<Point> & list_of_destination_voxels = params->transport_destination_voxels;
            vector<int> rand_voxel_index = RandomNG::randIntVec(0, list_of_destination_voxels.size() - 1, shuffled_indices.size());
            int counter = 0;
            // iterating through each of the randomly selected containers, set their transport destination to a random possible voxel
            for(auto & index : shuffled_indices) {
                if(counter == total_amount) break;
                auto itAtOffset = std::next(grand_total_container_map.begin(), index);
                itAtOffset->second->initialize_active_transport(list_of_destination_voxels[rand_voxel_index[counter]]);
                counter++;
            }
        }

    }
}

void Simulation::output_completed_stats(bool pop_small_mol, bool pop_containers, bool pop_compartment, bool multi_compart_out) {
    vector<string> small_mol_list = small_molecule_id_list;
    // These stats are for the diffusible small molecules.
    if(!small_molecule_id_list.empty()) {
        print_counts(bulk_counts_map, whole_bulk_out_name);
        if(pop_compartment) {
            if(multi_out) { 
                // if creating separate tsv for every compartment
                for(auto & compartment_pair : Compart_bulk_counts_map) {
                    if(compartment_pair.first == "default" || compartment_pair.first == "inaccessible") continue; // don't print these compartments that are meant to be empty
                    string file_out_name;
                    // insert compartment name to output file name before extension, probably tsv
                    if(compart_bulk_out_name.find(".") != string::npos) {
                        size_t ext_pos = string::npos;
                        ext_pos = compart_bulk_out_name.find(".");
                        // original up until extension + _compartment + extension
                        // if no extension then just append to the end
                        file_out_name = compart_bulk_out_name.substr(0, ext_pos) + "_" + compartment_pair.first + compart_bulk_out_name.substr(ext_pos);
                    } else {
                        file_out_name = compart_bulk_out_name + "_" + compartment_pair.first + ".tsv";
                    }

                    auto compart_bulk_count_map = compartment_pair.second;
                    print_counts(compart_bulk_count_map, file_out_name);
                }
            } else {
                // create single row concatenated TSV compartment output
                string spacer = "";
                // clear existing file
                ofstream compart_bulk_count(compart_bulk_out_name);
                compart_bulk_count.close();
                for(auto & compartment_pair : Compart_bulk_counts_map) {
                    if(compartment_pair.first == "default" || compartment_pair.first == "inaccessible") continue; // don't print these compartments that are meant to be empty
                    auto compart_bulk_count_map = compartment_pair.second;
                    ofstream compart_bulk_count(compart_bulk_out_name, ios::app);
                    compart_bulk_count << spacer << "cell4d_compartment " << compartment_pair.first << endl;
                    compart_bulk_count.close();
                    print_counts(compart_bulk_count_map, compart_bulk_out_name, true);
                    spacer = "\n";
                }
            }
        }
    }

    // print all container names
    if(!Container_counts_map.empty() && pop_containers) {
        // make copy of counts map, deleting bulk containers
        auto counts_map_copy = Container_counts_map;
        vector<string> delete_keys;
        for(auto element : counts_map_copy) {
            if(find(small_mol_list.begin(), small_mol_list.end(), element.first) != small_mol_list.end()) {
                delete_keys.push_back(element.first);
            }
        }
        for(auto element : delete_keys) counts_map_copy.erase(element);

        // output stream to designated file name, default is "total_container_stats.tsv"
        print_counts(counts_map_copy, whole_particle_out_name);
    }

    // sort container outputs into compartments
    if(!compartments_map.empty() && pop_compartment) {
        if(multi_out) {
            for(auto & compartment_pair : Compart_container_counts_map) {
                if(compartment_pair.first == "default" || compartment_pair.first == "inaccessible") continue; // don't print these compartments that are meant to be empty
                string file_out_name;
                // insert compartment name to output file name before extension, probably tsv
                if(compart_particle_out_name.find(".") != string::npos) {
                    size_t ext_pos = string::npos;
                    ext_pos = compart_particle_out_name.find(".");
                    // original up until extension + _compartment + extension, if no extension then just append to the end
                    file_out_name = compart_particle_out_name.substr(0, ext_pos) + "_" + compartment_pair.first + compart_particle_out_name.substr(ext_pos);
                } else {
                    file_out_name = compart_particle_out_name + "_" + compartment_pair.first + ".tsv";
                }


                // container name : vector of counts at each timestep
                auto compartment_container_map = compartment_pair.second;
                // remove entries that are actually bulk molecules
                vector<string> delete_keys;
                for(auto element : compartment_container_map) {
                    if(find(small_mol_list.begin(), small_mol_list.end(), element.first) != small_mol_list.end()) {
                        delete_keys.push_back(element.first);
                    }
                }
                for(auto element : delete_keys) compartment_container_map.erase(element);

                print_counts(compartment_container_map, file_out_name);
            }
        } else {
            string spacer = "";
            // clear existing file
            ofstream compartment_container_counts(compart_particle_out_name);
            compartment_container_counts.close();
            // separate stats into compartment-specific counts
            for(auto & compartment_pair : Compart_container_counts_map) {
                if(compartment_pair.first == "default" || compartment_pair.first == "inaccessible") continue; // don't print these compartments that are meant to be empty
                auto compartment_container_map = compartment_pair.second;
                // remove entries that are actually bulk molecules
                vector<string> delete_keys;
                for(auto element : compartment_container_map) {
                    if(find(small_mol_list.begin(), small_mol_list.end(), element.first) != small_mol_list.end()) {
                        delete_keys.push_back(element.first);
                    }
                }
                for(auto element : delete_keys) compartment_container_map.erase(element);

                ofstream compartment_container_counts(compart_particle_out_name, ios::app);
                compartment_container_counts << spacer << "cell4d_compartment " << compartment_pair.first << endl;
                compartment_container_counts.close();
                print_counts(compartment_container_map, compart_particle_out_name, true);
                spacer = "\n";
            }
        }
    }
    generate_checkpoint_file(true);
}

void Simulation::print_counts(map<string, vector<int>> counts_map, string file_name, bool append) {
    string delimiter = "";
    ofstream counts_out_stream;
    if(append) {
        counts_out_stream.open(file_name, ios::app);
    } else {
        counts_out_stream.open(file_name);
    }

    for(auto const & mol_type : counts_map) {
        if(mol_type.first.compare("") == 0) continue;
        counts_out_stream << delimiter << mol_type.first;
        delimiter = "\t";
    }
    counts_out_stream << endl;
    delimiter = "";
    for(int timestep = 0; timestep < max_cycles; timestep++) {
        for(auto const & mol_type : counts_map) {
            if(mol_type.first.compare("") == 0) continue;
            counts_out_stream << delimiter << mol_type.second[timestep];
            delimiter = "\t";
        }
        counts_out_stream << "\n";
        delimiter = "";
    }
    counts_out_stream.close();
}

void Simulation::print_counts(map<string, vector<double>> counts_map, string file_name, bool append) {
    string delimiter = "";
    ofstream counts_out_stream;
    if(append) {
        counts_out_stream.open(file_name, ios::app);
    } else {
        counts_out_stream.open(file_name);
    }

    for(auto const & mol_type : counts_map) {
        if(mol_type.first.compare("") == 0) continue;
        counts_out_stream << delimiter << mol_type.first;
        delimiter = "\t";
    }
    counts_out_stream << endl;
    delimiter = "";
    for(int timestep = 0; timestep < max_cycles; timestep++) {
        for(auto const & mol_type : counts_map) {
            if(mol_type.first.compare("") == 0) continue;
            counts_out_stream << delimiter << mol_type.second[timestep];
            delimiter = "\t";
        }
        counts_out_stream << "\n";
        delimiter = "";
    }
    counts_out_stream.close();
}

//////////////////////////////////////////////////////////////////////
//  Create and return Molecule of the type TYPE
// this isn't actually used, and it's overloaded, but still needs a body
// Overloaded in GenSimulation
Molecule_Container * Simulation :: createNewMoleculeContainer(string type, Point loc){
    return NULL;
}

Molecule_Container * Simulation :: createNewMoleculeContainer(vector<Species_Attributes *>, Point loc){
    return NULL;
}
//////////////////////////////////////////////////////////////////////
//  Input of molecules to the Simulation Space
//
void Simulation::initialDataInput(){}


void Simulation::update_molecule_counts() {
    // update the counts of small molecules at each timestep
    for(auto & small_mol : small_molecule_id_list) {
        double current_mol_total = 0;
        // if this metabolite wasn't seen before, fill with 0 up to current timestep
        if(bulk_counts_map.count(small_mol) == 0) {
            bulk_counts_map[small_mol].reserve(max_cycles);
            vector <double> new_vec(numCycles-1, 0);
            bulk_counts_map[small_mol].swap(new_vec);
        }
        // add up total counts from each compartment
        for(auto & compartment : compartments_map) {
            current_mol_total += compartment.second.current_bulk_count[small_mol];
        }
        bulk_counts_map[small_mol].push_back(current_mol_total);
    }

    for(auto & compart_pair : compartments_map) {
        string compart_name = compart_pair.first;
        Compart & compartment = compart_pair.second;

        for(auto & small_mol : small_molecule_id_list) {
            // if this metabolite wasn't seen before, fill with 0 up to current timestep
            if(Compart_bulk_counts_map[compart_name].count(small_mol) == 0) {
                Compart_bulk_counts_map[compart_name][small_mol].reserve(max_cycles);
                vector <double> new_vec(numCycles-1, 0);
                Compart_bulk_counts_map[compart_name][small_mol].swap(new_vec);
            }

            Compart_bulk_counts_map[compart_name][small_mol].push_back(compartment.current_bulk_count[small_mol]);
        }
    }

    // Periodically updating the simulation with total counts of each container at each timestep
    for(auto const & container_pair : Mol_Container_Table) {
        string container_name = container_pair.first;
        int container_count = container_pair.second.size();
        if(Container_counts_map.count(container_name) == 0) {
            Container_counts_map[container_name].reserve(max_cycles);
            vector <int> new_vec(numCycles-1, 0);
            Container_counts_map[container_name].swap(new_vec);
        }
        Container_counts_map[container_name].push_back(container_count);
    }

    
    // do the same for compartment-specific counts
    for(auto const & compart_pair : compartments_map) {
        string compart_name = compart_pair.first;
        Compart const & compartment = compart_pair.second;
        map<string, map<string, int>>  compart_container_current_counts;
        // iterate through all voxels of this compartment and sum up the container counts
        for(auto const & compart_voxel : compartment.list_of_voxels) {
            unordered_map <string, unordered_map <unsigned long, Molecule_Container *>> voxel_map = lEnv->getLatticeSite(*compart_voxel)->lattice_container_map;
            for(auto const & container_pair : voxel_map) { // for each container type in this voxel, tally containers
                string voxel_container_name = container_pair.first;
                int voxel_container_count = container_pair.second.size();
                compart_container_current_counts[compart_name][voxel_container_name] += voxel_container_count;
            }
        }
        // fills the containers that exist, but not in the current compartment, with a count of 0
        for(auto const & container_pair : Mol_Container_Table) {
            if(compart_container_current_counts[compart_name].find(container_pair.first) == compart_container_current_counts[compart_name].end()) {
                compart_container_current_counts[compart_name][container_pair.first] = 0;
            }
        }
        for(auto & container_pair : compart_container_current_counts[compart_name]) {
            string container_name = container_pair.first;
            int compart_container_count = container_pair.second;
            if(Compart_container_counts_map[compart_name].count(container_name) == 0) {
                Compart_container_counts_map[compart_name][container_name].reserve(max_cycles);
                vector <int> new_vec(numCycles-1, 0);
                Compart_container_counts_map[compart_name][container_name].swap(new_vec);
            }
            Compart_container_counts_map[compart_name][container_name].push_back(compart_container_count);
        }
    }

}

//////////////////////////////////////////////////////////
// load counts data from checkpoint
// loads whatever it finds, so compartment/whole particle/bulk_Mols
void Simulation::get_checkpointed_counts(json & json_file) {
    if(json_file.find("whole_container_counts") != json_file.end()) {
        for(auto & container_counts : json_file["whole_container_counts"]) {
            string container_name = container_counts["container_name"];
            vector<int> counts = container_counts["timestep_counts"];
            Container_counts_map[container_name] = counts;
        }
    }

    if(json_file.find("compart_container_counts") != json_file.end()) {
        // use iterator to grab the variable compartment name key
        for(auto & compartment_counts : json_file["compart_container_counts"].items()) {
            string compart_name = compartment_counts.key();
            for(auto & container_counts : compartment_counts.value()) {
                string container_name = container_counts["container_name"];
                vector<int> counts = container_counts["timestep_counts"];
                Compart_container_counts_map[compart_name][container_name] = counts;
            }
        }
    }

    if(json_file.find("whole_bulk_counts") != json_file.end()) {
        for(auto & bulk_counts : json_file["whole_bulk_counts"]) {
            string bulk_name = bulk_counts["bulk_name"];
            vector<double> counts = bulk_counts["timestep_counts"];
            bulk_counts_map[bulk_name] = counts;
        }
    }

    if(json_file.find("compart_bulk_counts") != json_file.end()) {
        // use iterator to grab the variable compartment name key
        for(auto & compartment_counts : json_file["compart_bulk_counts"].items()) {
            string compart_name = compartment_counts.key();
            for(auto & bulk_counts : compartment_counts.value()) {
                string bulk_name = bulk_counts["bulk_name"];
                vector<double> counts = bulk_counts["timestep_counts"];
                Compart_bulk_counts_map[compart_name][bulk_name] = counts;
            }
        }
    }

}



//////////////////////////////////////////////////////////////////////
// Generate checkpoint file
// right now it doesn't actually care if this is the final timestep or not, but leaving the flag in.
//
void Simulation :: generate_checkpoint_file(bool is_final) {
    json checkpoint_output;
    checkpoint_output["timesteps_completed"] = numCycles;

    // set the "molecule_list" entry, with all molecule positions and state
    json molecule_list;
    for(auto & container_type_pair : Mol_Container_Table){
        for(auto & container_itr : container_type_pair.second) {
            // Get the next molecule in the list
            Molecule_Container * container = container_itr.second;
            json container_object;
            map<string, double> container_pos;
            for(auto & molecule : container->container_vector) {
                json molecule_json;
                molecule_json["molecule_name"] = molecule->molecule_name;
                molecule_json["molecule_id"] = molecule->molecule_id;
                molecule_json["binding_sites"] = molecule->mol_config.molecule_binding_map;
                molecule_json["mod_sites"] = molecule->mol_config.molecule_state_map;
                container_object["molecules"].push_back(molecule_json);
            }
            container_object["container_name"] = container->container_name;
            container_object["container_id"] = container->container_id;

            container_pos["x"] = container->get_position().x;
            container_pos["y"] = container->get_position().y;
            container_pos["z"] = container->get_position().z;
            container_object["pos"] = container_pos;
            container_object["actively_transported"] = container->is_actively_transported();
            json transport_destination;
            if(container->is_actively_transported()) {
                transport_destination["x"] = container->actively_transported.second.x;
                transport_destination["y"] = container->actively_transported.second.y;
                transport_destination["z"] = container->actively_transported.second.z;
                container_object["transport_destination"] = transport_destination;
            }
            molecule_list.push_back(container_object);
        }
    }
    checkpoint_output["molecule_list"] = molecule_list;

    // beyond just saving the state of checkpoint, also save history of counts throughout compartments and whole
    // set "whole_container_counts" entry, with an array of container_name and counts array saved as each entry
    json whole_container_counts;
    for(auto & container_counts_entry : Container_counts_map) {
        json container_counts;
        container_counts["container_name"] = container_counts_entry.first;
        container_counts["timestep_counts"] = container_counts_entry.second;
        whole_container_counts.push_back(container_counts);
    }
    checkpoint_output["whole_container_counts"] = whole_container_counts;

    // set "compart_container_counts" entry, same structure as above but each compartment info is saved as the value with key being compartment name
    // compart_name:compart_info (which includes container name, array of counts)
    json compart_container_counts;
    for(auto & compart_container_counts_entry : Compart_container_counts_map) {
        json compartment_counts;
        string compartment_name = compart_container_counts_entry.first;
        if(compartment_name == "default" || compartment_name == "inaccessible") continue; // don't print these compartments that are meant to be empty
        for(auto & compartment : compart_container_counts_entry.second) {
            json container_counts;
            container_counts["container_name"] = compartment.first;
            container_counts["timestep_counts"] = compartment.second;
            compartment_counts.push_back(container_counts);
        }
        compart_container_counts[compartment_name] = compartment_counts;
    }
    checkpoint_output["compart_container_counts"] = compart_container_counts;

    // saving bulk molecules into the checkpoint file (if they exist)
    // this includes bulk mol placement upon resuming as well as bulk mol counts history
    if(!small_molecule_id_list.empty()) {
        // 
        json bulk_mole_info;
        Point p = Point(0, 0, 0);
        for(p.x = 0; lEnv->withinCell(p); p.x++) {
        for(p.y = 0; lEnv->withinCell(p); p.y++) {
        for(p.z = 0; lEnv->withinCell(p); p.z++) {
            vector<float> point_coord{p.x, p.y, p.z};
            LatticeSite * ls = lEnv->getLatticeSite(p);
            json voxel_info;
            voxel_info["voxel_coordinate"] = point_coord;
            json voxel_bulk_moles;
            for(auto & bulk_mol : ls->bulk_moles_map) {
                voxel_bulk_moles[bulk_mol.first] = bulk_mol.second;
            }
            voxel_info["voxel_bulk_moles"] = voxel_bulk_moles;
            bulk_mole_info.push_back(voxel_info);
        }
        p.z = 0;
        }
        p.y = 0;
        }
        checkpoint_output["bulk_mol_positions"] = bulk_mole_info;

        // set "whole_bulk_counts" entry, with an array of bulk_name and counts array saved as each entry
        json whole_bulk_counts;
        for(auto & bulk_counts_entry : bulk_counts_map) {
            json bulk_counts;
            bulk_counts["bulk_name"] = bulk_counts_entry.first;
            bulk_counts["timestep_counts"] = bulk_counts_entry.second;
            whole_bulk_counts.push_back(bulk_counts);
        }
        checkpoint_output["whole_bulk_counts"] = whole_bulk_counts;

        // set "compart_bulk_counts" entry, same structure as above but each compartment info is saved as the value with key being compartment name
        // compart_name:compart_info (which includes bulk name, array of counts)
        json compart_bulk_counts;
        for(auto & compart_bulk_counts_entry : Compart_bulk_counts_map) {
            json compartment_counts;
            string compartment_name = compart_bulk_counts_entry.first;
            if(compartment_name == "default" || compartment_name == "inaccessible") continue; // don't print these compartments that are meant to be empty
            for(auto & compartment : compart_bulk_counts_entry.second) {
                json bulk_counts;
                bulk_counts["bulk_name"] = compartment.first;
                bulk_counts["timestep_counts"] = compartment.second;
                compartment_counts.push_back(bulk_counts);
            }
            compart_bulk_counts[compartment_name] = compartment_counts;
        }
        checkpoint_output["compart_bulk_counts"] = compart_bulk_counts;
    }


    // take the trimmed checkpoint out name, add number of cycles, add .json extension
    string current_checkpoint_name = checkpoint_out_name + "_" + to_string(numCycles) + ".json";
    ofstream checkpoint_out(current_checkpoint_name);
    checkpoint_out << checkpoint_output.dump(4) << endl;

    if(is_final) {
        // if poslog is enabled, output the json file
        if(poslog) {
            for(auto poslog_compart : poslog_json) {
                // add back the json extension to the end, so name format is <poslog>_<compart>.json
                string poslog_out_full = poslog_out + "_" + poslog_compart.first + ".json";

                // open out stream with name poslog_out_full
                ofstream poslog_out_json(poslog_out_full);

                // dump poslog_json json object with recorded mol positions across timesteps out the stream, with tab spacing of 4
                poslog_out_json << poslog_compart.second.dump(4) << endl;
            }
        }
    }
}

// triggered at poslog_interval timesteps, adding state of compartment molecules to the json log
void Simulation::add_molpos_to_log(json & compart_json, int log_time, string compart) {

    // set the "molecule_list" entry, with all molecule positions and state
    json current_mol_list;

    Compart & logged_compart = compartments_map[compart];
    // iterate through all voxels of this compartment and add molecules in them to the json log
    // this is copied from the checkpointing code
    for(auto const & compart_voxel : logged_compart.list_of_voxels) {
        unordered_map <string, unordered_map <unsigned long, Molecule_Container *>> voxel_map = lEnv->getLatticeSite(*compart_voxel)->lattice_container_map;
        for(auto & container_type_pair : voxel_map){
            for(auto & container_itr : container_type_pair.second) {
                // Get the next molecule in the list
                Molecule_Container * container = container_itr.second;
                json container_object;
                map<string, double> container_pos;
                for(auto & molecule : container->container_vector) {
                    json molecule_json;
                    molecule_json["molecule_name"] = molecule->molecule_name;
                    molecule_json["molecule_id"] = molecule->molecule_id;
                    molecule_json["binding_sites"] = molecule->mol_config.molecule_binding_map;
                    molecule_json["mod_sites"] = molecule->mol_config.molecule_state_map;
                    container_object["molecules"].push_back(molecule_json);
                }
                container_object["container_name"] = container->container_name;
                container_object["container_id"] = container->container_id;

                container_pos["x"] = container->get_position().x;
                container_pos["y"] = container->get_position().y;
                container_pos["z"] = container->get_position().z;
                container_object["pos"] = container_pos;
                container_object["actively_transported"] = container->is_actively_transported();
                json transport_destination;
                if(container->is_actively_transported()) {
                    transport_destination["x"] = container->actively_transported.second.x;
                    transport_destination["y"] = container->actively_transported.second.y;
                    transport_destination["z"] = container->actively_transported.second.z;
                    container_object["transport_destination"] = transport_destination;
                }
                current_mol_list.push_back(container_object);
            }
        }
    }

    // add molecule details within the selected compartment to json, key is timestep of logged molecules
    compart_json[to_string(log_time)] = current_mol_list;
}

//////////////////////////////////////////////////////////////////////
//  Display data collected from the Simulation
//
void Simulation :: displayResults(){}



void Simulation :: checkEnzymaticReactions(){
    // Nov 1 2018 - We'll assume that enzymatic reactions only have one modifier
    // Iterate through the list of enzymatic reacions
    for(auto & enz_reaction_itr : enzymatic_reaction_details_map) {
        string reaction_name = enz_reaction_itr.first;
        string modifier_name = enz_reaction_itr.second->modifiers_list[0];
        auto mol_name = Mol_Table.find(modifier_name); // look for the enzyme in the molecule table
        if(mol_name == Mol_Table.end()) continue;
        unordered_map<unsigned long, Molecule *> reacting_molecules = Mol_Table[modifier_name];
        float Ka                    = enz_reaction_itr.second->forward_reaction_rate; // This is fine, the Ka for enzyme kinetics is saved in this slot from PM
        float Kp                    = enz_reaction_itr.second->reverse_reaction_rate; // same as above
        float equilibriumConstant   = enz_reaction_itr.second->equilibrium_constant;
        float maximumRate           = enz_reaction_itr.second->maximum_rate;
        string reaction_reactant = enz_reaction_itr.second->reactants_list[0][0];
        string reaction_product = enz_reaction_itr.second->products_list[0][0];
        // iterate through all enzymes of this name
        for(auto & enzyme_itr : reacting_molecules) {
            lEnv->reactWithMetabolite(enzyme_itr.second->container_pointer->get_position(), reaction_reactant, reaction_product, Ka, Kp, equilibriumConstant, maximumRate);
        }
    }
}
