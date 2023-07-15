//============================================================
// GenSimulation.cpp
// Donny Chan, Billy Taj

#include "../inc/GenSimulation.h"

using namespace std;

static GenSimulation * simulation = NULL;

extern map<int, Species_t *> final_signal_enzyme_map;
extern map<string, Cell_Reaction *> bimolecular_reaction_details_map;
extern map<string, Cell_Reaction *> unimolecular_reaction_details_map;
extern map<string, Cell_Reaction *> bulk_reaction_details_map;
extern unordered_map <string, unordered_map <unsigned long, Molecule *>> Mol_Table;
extern map<string, Compart> compartments_map;
extern int numContainers;
extern bool from_checkpoint;
extern json checkpoint_json;
//////////////////////////////////////////////////////////////////////
// Get the instance of the Cell Simulation
// ->    Only allow one instance of the cell to exist at any one time
//         If not yet exist, create a new Cell Simulation
//         otherwise return the pointer to the current existing Cell Simulation
//
GenSimulation * GenSimulation::getSimulation(){
    // If the simulation cell Simulation not yet exist
    if(!isInstantiated){
        // Create new Cell Simulation
        // also calls constructor of simulation
        simulation = new GenSimulation();
        // New Cell Simulation Created
        isInstantiated = true;
    }
    // Return the pointer to the current Cell Simulation
    return simulation;
}



//////////////////////////////////////////////////////////////////////
// Constructor for the GenSimulation
// ->    Create and initialize properties of the GenSimulation space
GenSimulation::GenSimulation():Simulation(){

    // Initialize Simulation
    initialDataInput();
}



//////////////////////////////////////////////////////////////////////
// Destructor for the GenSimulation
//
GenSimulation:: ~GenSimulation(){}


//////////////////////////////////////////////////////////////////////
//  Input of molecules to the GenSimulation Space
//  (number of molecules, molecule type, environment)
//
void GenSimulation::initialDataInput(){
    
    // insert particles into simulation space
    if(from_checkpoint) {
        // if resuming simulation with json-notated molecule placements
        create_checkpoint_molecules();
    } else {
        // create all molecular containers specified in listOfSpecies using initial_species_map
        auto initial_containers_info = ParameterManager::getParameterManager()->initial_species_map;
        for(auto & initial_details_itr : initial_containers_info) {
            map<string, int> compartment_info = initial_details_itr.first;
            vector<Species_Attributes *> container_components = initial_details_itr.second;
            // if there is more than 1 species in vector, then it is a particle
            if(container_components.size() > 1 || container_components[0]->is_protein()) {
                for(auto & compartments_init_itr : compartment_info) {
                    int compartment_init_amount = compartments_init_itr.second;
                    string compartment_name = compartments_init_itr.first;
                    // createMoleculeContainers needs
                    // 1) the initial amount to create
                    // 2) the name of the molecule
                    // 3) the environment name
                    createMoleculeContainers(compartment_init_amount, container_components, compartment_name);
                }
            }
        }
    }

    // insert bulk molecules into simulation space
    if(from_checkpoint) {
        load_bulk_molecules();
    } else {
        create_bulk_molecules();
    }
    // create_bulk_molecules();

    // determine the reaction type of all reactions
    // for reversible multi-product reactions, meaning that the reverse reaction is bimolecular, an unbinding radius needs to be set to prevent geminate recombinations
    for(auto & reaction : bimolecular_reaction_details_map) {
        reaction.second->determine_interaction_type();
        cout << "reaction " << reaction.first << " type is " << reaction.second->get_interaction_type() << endl;
        // call functions with AB_adjust set to true
        reaction.second->set_binding_radius(true);
        reaction.second->set_reaction_rate(true);
    }
    for(auto & reaction : unimolecular_reaction_details_map) {
        reaction.second->determine_interaction_type();
        cout << "reaction " << reaction.first << " type is " << reaction.second->get_interaction_type() << endl;
        cout << "Reaction " << reaction.first << " rate constant is " << reaction.second->forward_reaction_rate << endl;
    }
    for(auto & reaction : bulk_reaction_details_map) {
        reaction.second->determine_bulk_interaction_type();
        cout << "reaction " << reaction.first << " type is " << reaction.second->get_interaction_type() << endl;
        // do not use AB adjustment for bulk, smoluchowski rates are fine
        reaction.second->set_binding_radius(false);
        reaction.second->set_reaction_rate(false);
    }

    for(auto & reaction : bimolecular_reaction_details_map) {
        if(reaction.second->is_fission() && reaction.second->is_reversible) {
            if(reaction.second->unbind_rad_min == 0) {
                reaction.second->set_unbinding_radius(reaction.second->reverse_reaction->binding_radius);
            }
        }
    }
    for(auto & reaction : unimolecular_reaction_details_map) {
        if(reaction.second->is_fission() && reaction.second->is_reversible) {
            if(reaction.second->unbind_rad_min == 0) {
                reaction.second->set_unbinding_radius(reaction.second->reverse_reaction->binding_radius);
            }
        }
    }

}

// using the parsed json object from param manager, iterate through each entry to create all the saved molecules with correct properties
void GenSimulation::create_checkpoint_molecules() {
    // loop through each json array molecule container entry
    for(auto & json_mol_container : checkpoint_json["molecule_list"]) {
        // retrieve molecule position, name, components, id, binding/mod sites
        Point loc;
        loc.x = json_mol_container["pos"]["x"];
        loc.y = json_mol_container["pos"]["y"];
        loc.z = json_mol_container["pos"]["z"];
        vector<Molecule *> molecule_vector;
        for(auto & molecule : json_mol_container["molecules"]) {
            unsigned long old_mol_id = molecule["molecule_id"];
            Molecule * new_mol = new Molecule(string(molecule["molecule_name"]), old_mol_id);
            map<string, int> bind_sites = molecule["binding_sites"];
            new_mol->mol_config.molecule_binding_map = bind_sites;
            map<string, string> mod_sites = molecule["mod_sites"];
            new_mol->mol_config.molecule_state_map = mod_sites;
            molecule_vector.push_back(new_mol);
        }
        unsigned long old_container_id = json_mol_container["container_id"];
        Molecule_Container * new_container = createNewMoleculeContainer(molecule_vector, loc, old_container_id);

        // ensure that molecules that were actively moving when the checkpoint was created maintains its previous destination
        if(json_mol_container["actively_transported"].get<bool>()) {
            Point dest_pos;
            dest_pos.x = json_mol_container["transport_destination"]["x"];
            dest_pos.y = json_mol_container["transport_destination"]["y"];
            dest_pos.z = json_mol_container["transport_destination"]["z"];
            new_container->actively_transported = make_pair(true, dest_pos);
        }
        // cout << "created container " << new_container->container_name << " from json." << endl;
    }
}

Molecule_Container * GenSimulation::createNewMoleculeContainer(string molecule_name, Point loc) {
    Molecule * new_molecule = new Molecule(molecule_name);

    Molecule_Container* new_molecule_container = new Molecule_Container(new_molecule, loc);
    // Add the created molecule to the Mol_Container_Table
    new_molecule_container->add_to_maps();
    new_molecule->container_id = new_molecule_container->container_id;
    new_molecule->container_pointer = new_molecule_container;

    return new_molecule_container;
}

// create molecule container using already created molecules, meant for json molecule list initialization
Molecule_Container * GenSimulation::createNewMoleculeContainer(vector<Molecule *> molecule_vector, Point loc, unsigned long old_container_id) {
    Molecule_Container* new_molecule_container = new Molecule_Container(molecule_vector, loc);
    // assign the previous container id
    new_molecule_container->container_id = old_container_id;
    // Add the created molecule to the Mol_Container_Table
    new_molecule_container->add_to_maps();
    for(auto & mol : new_molecule_container->container_vector) {
        mol->container_id = new_molecule_container->container_id;
        mol->container_pointer = new_molecule_container;
    }

    return new_molecule_container;
}

Molecule_Container * GenSimulation::createNewMoleculeContainer(vector<Species_Attributes *> molecules_stats, Point loc) {
    vector<Molecule *> mol_vec;
    for(auto & mol_info : molecules_stats) {
        string molecule_name = mol_info->ID;
        Molecule * new_molecule = new Molecule(mol_info);
        mol_vec.push_back(new_molecule);
    }
    
    Molecule_Container* new_molecule_container = new Molecule_Container(mol_vec, loc);
    // Add the created molecule to the Mol_Container_Table
    new_molecule_container->add_to_maps();
    for(auto & mol : new_molecule_container->container_vector) {
        mol->container_id = new_molecule_container->container_id;
        mol->container_pointer = new_molecule_container;
    }

    return new_molecule_container;
}


//////////////////////////////////////////////////////////////////////
//  Perform periodic modifications to the GenSimulation Space
//  ->    Temporal events updating behaviour/Properties of the GenSimulation
//
void GenSimulation :: updateSimulation(){


}

// if this is loaded from a checkpoint, grab the bulk molecules from json instead of using the xml
void GenSimulation::load_bulk_molecules() {
    for(auto & current_c_voxel : checkpoint_json["bulk_mol_positions"]) {
        vector<double> coord = current_c_voxel["voxel_coordinate"];
        Point point_coord = Point(coord);
        LatticeSite * current_voxel = lEnv->getLatticeSite(point_coord);
        for(auto & bulk_mol : current_c_voxel["voxel_bulk_moles"].items()) {
            string bulk_mol_name = bulk_mol.key();
            double bulk_mol_value = bulk_mol.value();
            current_voxel->add_mole_to_voxel(bulk_mol_name, bulk_mol_value); 
        }
    }
}


// create all small molecules at simulation start
void GenSimulation::create_bulk_molecules() {

    // "pump in" all the small molecules
    auto initial_containers_info = ParameterManager::getParameterManager()->initial_species_map;
    // pair < map<compartment:initial amount>, comtainer components info>
    for(auto & initial_details_itr : initial_containers_info) {
        map<string, int> compartment_info = initial_details_itr.first;
        vector<Species_Attributes *> species_vector = initial_details_itr.second;
        // can only be bulk molecule if 1 species exists
        if (species_vector[0]->is_simple_molecule() && species_vector.size() == 1) {
            string small_mol_name = species_vector[0]->ID;
            cout << "reactant:" << small_mol_name << endl;
            for(auto & compartments_init_itr : compartment_info) {
                string compartment_name = compartments_init_itr.first;
                int compart_init_particles = compartments_init_itr.second;
                vector <Point *> points_list = lEnv->getCompartmentSites(compartment_name);
                if (points_list.size() == 0 || compartments_init_itr.second == 0) { // skip invalid compartments
                    continue;
                }

                double compart_init_moles = Species_Attributes::particle_to_mole(compart_init_particles);

                // evenly spread out the amount of bulk molecules to pump out
                double moles_added = compart_init_moles / points_list.size();
                cout << "initial:" << compart_init_particles << " converted moles:" << compart_init_moles << " spread amount:" << moles_added << " into " << points_list.size() << " voxels:" << endl;
                // for each c-voxel, insert the determined amount of small molecules
                for (Point * compart_point : points_list) {
                    lEnv->add_mole_to_voxel(small_mol_name, moles_added, lEnv->getLatticeSite(*compart_point));
                }

            }
        }
    }
}


bool GenSimulation :: isSimulationCompleted(){

    return simulationCompleted;
}




//////////////////////////////////////////////////////////////////////
//  Display data collected from the Simulation
//
void GenSimulation :: displayResults(){
    map< string , int >::const_iterator itr = actEventTable.begin();

    for (itr = actEventTable.begin(); itr != actEventTable.end(); itr ++)
    {
        printf("~\t%s\t%d\n", (itr->first).c_str(), itr->second);
    }

    // Display total number of cycles to complete simulation
    printf("#\tTotal num of cycles\t%d\n",numCycles);
}
