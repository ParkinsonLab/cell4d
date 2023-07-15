//============================================================
// Molecular_Container_Reactions.cpp
// Donny Chan, Billy Taj


#include "../inc/Molecule_Container_Reactions.h"
extern unordered_map <string, unordered_map <unsigned long, Molecule_Container *> > Mol_Container_Table;
extern map<string, Species_Attributes *> species_details_map; 
extern vector<string> small_molecule_id_list;                               // new vector of all small_molecule names
extern map<string, Cell_Reaction *> bulk_reaction_details_map;                // stores details of bulk/particle reactant reactions
extern map<string, Cell_Reaction *> unimolecular_reaction_details_map;        // stores reaction details of unimolecular reactions (product, reactants, etc etc.  will be expanded)
extern map<string, Cell_Reaction *> bimolecular_reaction_details_map;         // stores reaction details of bimolecular reactions (product, reactants, etc etc.  will be expanded)

extern vector<Molecule_Container *> containers_to_remove_map;
extern int numContainers;

/*
    This class dictates what happens when 2 or more Molecule Container objects interact.  
    It won't have any internal pieces itself, but just serves to manipulate the Molecule Containers where appropriate
*/

Molecule_Container_Reactions::Molecule_Container_Reactions(){

}

Molecule_Container_Reactions::~Molecule_Container_Reactions(){

}

vector<Molecule_Container *> Molecule_Container_Reactions::create_new_empty_containers(int n) {
    vector<Molecule_Container *> output_vector;
    for(int count = 0; count < n; count++) {
        Molecule_Container* new_molecule_container = new Molecule_Container();
        output_vector.push_back(new_molecule_container);
    }
    return output_vector;
}

// remove molecule in second argument from container in first argument
void Molecule_Container_Reactions::remove_molecule_from_container(Molecule_Container * container, Molecule * molecule) {
    if(find_if(container->container_vector.begin(), container->container_vector.end(), [molecule](Molecule* mol) { return mol == molecule; }) != container->container_vector.end()) {
        container->container_vector.erase(remove(container->container_vector.begin(), container->container_vector.end(), molecule), container->container_vector.end());
        // container->initialize_container();
    } else {
        cerr << "did not find required molecule in container, exiting." << endl;
        exit(1);
    }
    molecule->container_id = 0;
    molecule->container_pointer = NULL;
    container->allowed_to_react = false;
}

void Molecule_Container_Reactions::add_molecule_to_container(Molecule_Container * container, Molecule * molecule) {
    container->container_vector.push_back(molecule);
    // container->initialize_container();
    molecule->container_id = container->container_id;
    molecule->container_pointer = container;
    container->allowed_to_react = false;
}

bool Molecule_Container_Reactions::molecule_fills_rule(Molecule * mol, const map<string, string> & state_rules, const map<string, int> & bind_rules) {
    if(includes(mol->mol_config.molecule_binding_map.begin(),mol->mol_config.molecule_binding_map.end(),bind_rules.begin(),bind_rules.end())) {
        if(includes(mol->mol_config.molecule_state_map.begin(),mol->mol_config.molecule_state_map.end(),state_rules.begin(),state_rules.end())) {
            return true;
        }
    }
    return false;
}

Molecule_Container * Molecule_Container_Reactions::bulk_fusion(Molecule_Container * head, Molecule * bulk){
    add_molecule_to_container(head, bulk);
    // head->initialize_container();
    head->allowed_to_react = false;
    return head;
}

Molecule_Container * Molecule_Container_Reactions::fusion(Molecule_Container * head, Molecule_Container * tail){
    // delete both the old molecules from the lattice map and container table
    // remove_container_from_maps(head);
    // remove_container_from_maps(tail);
    // put all in tail.  erase head
    for(Molecule * mol : head->container_vector) {
        add_molecule_to_container(tail, mol);
    }

    head->container_vector.clear();
    containers_to_remove_map.push_back(head);

    head->allowed_to_react = false;

    // string old_name = tail->container_name;
    // tail->name_mol_container();
    // string new_name = tail->container_name;
    // tail->container_name = old_name;

    // tail->initialize_container();
    tail->allowed_to_react = false;
    return tail;
}

// split the input "reacting" molecules into however many containers needed for the reaction
// can be a uni or bimolecular reaction
vector<Molecule_Container *> Molecule_Container_Reactions::fission(map<Molecule *, pair<int, int>> mol_to_prod_map, Cell_Reaction * reaction_details) {
    map<int, vector <string>> products = reaction_details->products_list;
    int spectator_placement = reaction_details->spectator_index;
    Molecule_Container * reactant_1 = mol_to_prod_map.begin()->first->container_pointer;
    LatticeEnvironment * lEnv = reactant_1->get_lat_env();
    map<unsigned long, vector<Molecule *>> reactant_container_count_map; // container id key for a vector of reactant molecules that share the same container
    for(auto & mol : mol_to_prod_map) {
        reactant_container_count_map[(mol.first)->container_id].push_back(mol.first);
    }

    // flag for whether or not this is a membrane reaction
    pair<bool, Point> membrane_reaction = make_pair(false, Point(0, 0, 0));
    // calculate the midpoint of the bimolecular reaction, useful for placing new containers into the simulation
    vector <Point> original_points;
    for(auto & origin_container : reactant_container_count_map) {
        original_points.push_back(origin_container.second[0]->container_pointer->get_position());
        // see if this is a membrane reaction
        if(origin_container.second[0]->container_pointer->is_in_membrane()) {
            membrane_reaction.first = true;
            membrane_reaction.second = origin_container.second[0]->container_pointer->get_position();
        }
    }
    Point mid_point;
    if(original_points.size() > 1) {
        mid_point = Point::mid_point(original_points[0], original_points[1]);
    } else if(original_points.size() == 1) {
        mid_point = original_points[0];
    }

    
    // set up new blank containers for the fission products
    vector<Molecule_Container *> new_mol_containers = create_new_empty_containers(products.size());

    // find spectating molecules of reaction (molecules that are in reacting containers, but don't react)
    set<Molecule_Container *> old_container_set;
    vector<Molecule *> spectator_molecules, reacting_molecules;
    for(auto & old_container : mol_to_prod_map) {
        old_container_set.insert(old_container.first->container_pointer);
    }
    for(auto & reacting_mol : mol_to_prod_map) reacting_molecules.push_back(reacting_mol.first);
    vector<Molecule_Container *> add_to_removed_containers;
    for(Molecule_Container * container : old_container_set) {
        vector <Molecule *> to_remove;
        for(Molecule * container_mol : container->container_vector) {
            // if molecule not found in reacting molecules, then it is a spectator. remove and place into spectator vector
            if(find_if(reacting_molecules.begin(), reacting_molecules.end(), [container_mol](Molecule* mol) { return mol == container_mol; }) == reacting_molecules.end()) {
                spectator_molecules.push_back(container_mol);
                to_remove.push_back(container_mol);
            }
        }
        for(auto & spectator : to_remove) {
            remove_molecule_from_container(container, spectator);
        }
        for(auto & reacting_mol : reacting_molecules) {
            if(reacting_mol->container_pointer == container) {
                remove_molecule_from_container(container, reacting_mol);
            }
        }
        // now all reacting and non-reacting molecules are accounted for, set up old containers for disposal
        add_to_removed_containers.push_back(container);
    }

    for (auto & old_container : add_to_removed_containers) { // don't allow the old empty containers to react
        old_container->allowed_to_react = false;
    }
    containers_to_remove_map.insert(containers_to_remove_map.end(),add_to_removed_containers.begin(), add_to_removed_containers.end());


    // store vector of molecules grouped by the new container organization
    map<int, vector<Molecule *>> product_vectors; 
    for(auto & mol_pair : mol_to_prod_map) {
        pair<int, int> mol_position = mol_pair.second;
        product_vectors[mol_position.first].push_back(mol_pair.first);
    }

    int new_prod_counter_index = 0;
    bool spectators = false, spects_placed = false;
    if(!spectator_molecules.empty()) spectators = true;
    // with all the extracted molecules, put them into the correct bins of empty containers, then initialize them
    for(auto & temp_container_vector : product_vectors) {
        auto & molecule_vector = temp_container_vector.second;
        Molecule_Container * blank_container = new_mol_containers[new_prod_counter_index];
        for(auto & moving_molecule : molecule_vector) { // add molecules into the new blank container
            add_molecule_to_container(blank_container, moving_molecule);
            // place spectator molecules in the correct designated containers
            // intention: spectators should remain bound to the molecules it was with before, so stick them with reactants that have filled sites
            // so that it doesn't get glued to a solo product container
            if(spectators && !spects_placed) { // if there are spectator molecules 
                if(spectator_placement != -1) { // if spectator placement defined, put them in the correct product container
                    if(spectator_placement == new_prod_counter_index) {
                        for(auto & molecule : spectator_molecules) {
                            add_molecule_to_container(blank_container, molecule);
                        }
                        spects_placed = true;
                    }
                } else { // if placement index undefined, put them in the first product container that has occupied sites
                    for(auto & site : moving_molecule->mol_config.molecule_binding_map) {
                        if(site.second == 1) {
                            for(auto & molecule : spectator_molecules) {
                                add_molecule_to_container(blank_container, molecule);
                            }
                            spects_placed = true;
                            break;
                        }
                    }
                }
            }
        }
        blank_container->name_mol_container();

        // appropriately deal with solo containers that are actually bulk
        if(molecule_vector.size() == 1 && find(small_molecule_id_list.begin(), small_molecule_id_list.end(), molecule_vector.front()->molecule_name) != small_molecule_id_list.end()) {
            blank_container->set_position(mid_point);
            containers_to_remove_map.push_back(blank_container);
            blank_container->is_small_mol = true;
        }
        new_prod_counter_index++;
    }

    // clear out the unused blank containers
    for(int i = new_prod_counter_index; i < new_mol_containers.size(); i++) {
        if(!new_mol_containers[i]->container_vector.empty()) {
            cerr << "Container is not empty when being deleted" << endl;
        }
        delete new_mol_containers[i];
    }

    // set the new containers' positions (randomized around mid-point of reaction, technically cube)
    double rand_radius_min = reaction_details->unbind_rad_min;
    double rand_radius_range = reaction_details->unbind_rad_max - rand_radius_min;
    for(auto & new_container : new_mol_containers) {
        // set the first container product as the "midpoint"
        if(new_container == new_mol_containers.front()) {
            if(membrane_reaction.first) {
                if(lEnv->isCompatibleEnvironment(new_container->container_name, lEnv->get_compart_name(membrane_reaction.second))) {
                    LatticeSite * membrane_voxel = lEnv->getLatticeSite(membrane_reaction.second);
                    if(membrane_voxel->locked_axis.first == "x") mid_point.x = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "y") mid_point.y = membrane_voxel->locked_axis.second;
                    else if(membrane_voxel->locked_axis.first == "z") mid_point.z = membrane_voxel->locked_axis.second;
                }
            }
            new_container->set_position(mid_point);
            // new_container->initialize_container();
            continue;
        }

        Point rand_loc;
        // instead of using RandomNG float generator multiple times, create distribution once and sample from it
        uniform_real_distribution<float> radius_range_dist(-rand_radius_range, rand_radius_range);
        rand_loc.x = radius_range_dist(RandomNG::get_generator());
        rand_loc.x = rand_loc.x > 0 ? mid_point.x + rand_radius_range + rand_radius_min : mid_point.x - rand_radius_range - rand_radius_min;
        rand_loc.y = radius_range_dist(RandomNG::get_generator());
        rand_loc.y = rand_loc.y > 0 ? mid_point.y + rand_radius_range + rand_radius_min : mid_point.y - rand_radius_range - rand_radius_min;
        rand_loc.z = radius_range_dist(RandomNG::get_generator());
        rand_loc.z = rand_loc.z > 0 ? mid_point.z + rand_radius_range + rand_radius_min : mid_point.z - rand_radius_range - rand_radius_min;

        string new_env = lEnv->get_compart_name(rand_loc);
        int counter = 0;
        while(!lEnv->isCompatibleEnvironment(new_container->container_name, new_env) && counter < 10) {
            rand_loc.x = radius_range_dist(RandomNG::get_generator());
            rand_loc.x = rand_loc.x > 0 ? mid_point.x + rand_radius_range + rand_radius_min : mid_point.x - rand_radius_range - rand_radius_min;
            rand_loc.y = radius_range_dist(RandomNG::get_generator());
            rand_loc.y = rand_loc.y > 0 ? mid_point.y + rand_radius_range + rand_radius_min : mid_point.y - rand_radius_range - rand_radius_min;
            rand_loc.z = radius_range_dist(RandomNG::get_generator());
            rand_loc.z = rand_loc.z > 0 ? mid_point.z + rand_radius_range + rand_radius_min : mid_point.z - rand_radius_range - rand_radius_min;
            new_env = lEnv->get_compart_name(rand_loc);
            counter++;
        }
        // if we can't find a good destination around the radius, find closest voxel that is a valid compartment and drop it in there
        if(counter >= 10) {
            vector <string> allowed_comparts = lEnv->get_compatible_env(new_container->container_name);
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
            if(!lEnv->isCompatibleEnvironment(new_container->container_name, lEnv->get_compart_name(dest_voxel))) {
                cerr << "can't place molecule after fission reaction." << endl;
                exit(1);
            }
            if(membrane_reaction.first && lEnv->isCompatibleEnvironment(new_container->container_name, lEnv->get_compart_name(membrane_reaction.second))) {
                LatticeSite * membrane_voxel = lEnv->getLatticeSite(membrane_reaction.second);
                if(membrane_voxel->locked_axis.first == "x") dest_voxel.x = membrane_voxel->locked_axis.second;
                else if(membrane_voxel->locked_axis.first == "y") dest_voxel.y = membrane_voxel->locked_axis.second;
                else if(membrane_voxel->locked_axis.first == "z") dest_voxel.z = membrane_voxel->locked_axis.second;
            }
            new_container->set_position(dest_voxel);
        } else {
            if(membrane_reaction.first && lEnv->isCompatibleEnvironment(new_container->container_name, lEnv->get_compart_name(membrane_reaction.second))) {
                LatticeSite * membrane_voxel = lEnv->getLatticeSite(membrane_reaction.second);
                if(membrane_voxel->locked_axis.first == "x") rand_loc.x = membrane_voxel->locked_axis.second;
                else if(membrane_voxel->locked_axis.first == "y") rand_loc.y = membrane_voxel->locked_axis.second;
                else if(membrane_voxel->locked_axis.first == "z") rand_loc.z = membrane_voxel->locked_axis.second;
            }
            new_container->set_position(rand_loc);
        }
        // new_container->initialize_container();
    }
    return new_mol_containers;
}

void Molecule_Container_Reactions::molecule_state_change(Molecule * mol_component, map <string, string> new_states) {
    if(new_states.empty()) return;
    for(auto & changed_site : new_states) { 
        mol_component->mol_config.molecule_state_map[changed_site.first] = changed_site.second;
    }
}

void Molecule_Container_Reactions::molecule_bind_change(Molecule * mol_component, map <string, int> new_states) {
    if(new_states.empty()) return;
    // replace the relevant binding site statuses of reacting molecule with new states
    for(auto & changed_site : new_states) { 
        mol_component->mol_config.molecule_binding_map[changed_site.first] = changed_site.second;
    }
}
