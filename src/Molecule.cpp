//============================================================
// Molecule.cpp
// Donny Chan, Billy Taj

#include "../inc/Molecule.h"
#include <unordered_map>


extern map<string, Species_Attributes *> species_details_map;   // store reactant details. key: reactant ID
extern unordered_map <string, unordered_map <unsigned long, Molecule *> > Mol_Table;


//////////////////////////////////////////////////////////////////////
// Constructors for the Molecule object
//

Molecule::Molecule(string name) {
    set_name(name);
    initialize();
}

Molecule::Molecule(string name, unsigned long old_mol_id) {
    set_name(name);
    initialize(old_mol_id);
}

Molecule::Molecule(Species_Attributes * mol_info) {
    set_name(mol_info->ID);
    initialize(mol_info);
}


//////////////////////////////////////////////////////////////////////
// Destructor for the Molecule object
//
Molecule::~Molecule(){
    if(Mol_Table[molecule_name].find(molecule_id) == Mol_Table[molecule_name].end()) {
        Mol_Table[molecule_name].erase(molecule_id);
    }
}

void Molecule::initialize(Species_Attributes * mol_info, unsigned long mol_id){
    lEnv = LatticeEnvironment :: getLatticeEnv(); // Environment Manager
    diffusion_constant = species_details_map[molecule_name]->diffusion_constant; // Default mobility
    diffusible = true; // Default Diffusibility
    activated = false; // Molecule initially non-active

    for(auto & state : mol_info->default_states) { // insert the default state of the species as the molecule is initialized
        mol_config.molecule_state_map[state.first] = state.second;
    }
    for(auto & site : mol_info->default_bind) { // insert the default state of the species as the molecule is initialized
        mol_config.molecule_binding_map[site.first] = site.second;
    }
    if(mol_id == 0) {
        set_molecule_id();
    } else {
        molecule_id = mol_id;
    }
    if(Mol_Table[molecule_name].find(molecule_id) == Mol_Table[molecule_name].end()) {
        Mol_Table[molecule_name][molecule_id] = this;
    }
    mol_config.molecule_name = molecule_name;
}

void Molecule::initialize(unsigned long mol_id){
    // call the overloaded initializer using default molecule info from species details map, when none was provided
    Species_Attributes * default_mol_info = species_details_map[molecule_name];
    initialize(default_mol_info, mol_id);
}

bool Molecule::mol_point_comp(Molecule * a, Molecule * b) {
    return a->molecule_name < b->molecule_name;
}


//////////////////////////////////////////////////////////////////////
// Set the molecule's name
void Molecule::set_name(string name){
    molecule_name = name;
}

//////////////////////////////////////////////////////////////////////
// Get the molecule's name
string Molecule::get_name(){
    return molecule_name;
}

// set molecule id
void Molecule::set_molecule_id() {
    // using mt19937_64 (high period random generator), make a molecule id.
    molecule_id = RandomNG::randLong();
}

//////////////////////////////////////////////////////////////////////
// Print the current location of the molecule
//
string Molecule::position(){

    // String representation of the current position of the molecule
    // in the Lattice
    return "(" + floatToStr(pos.x) + ", " +
            floatToStr(pos.y) + ", " +
            floatToStr(pos.z) + ")";
}



//////////////////////////////////////////////////////////////////////
// Molecule translates/moves by the specified vector
//
void Molecule::translate(Vector3D & v){

    // Apply translation
    pos.x += v.x; pos.y += v.y; pos.z += v.z;
}



//////////////////////////////////////////////////////////////////////
// Set molecule's Diffusibility
//
void Molecule::setDiffusible(bool diff){
    diffusible = diff;
}



//////////////////////////////////////////////////////////////////////
// Get molecule's Diffusibility
//
bool Molecule::isDiffusible(){
    return diffusible;
}



//////////////////////////////////////////////////////////////////////
// Set molecule to Position
//
void Molecule::set_position(Point & p){
    pos.x = p.x; pos.y = p.y; pos.z = p.z;
}



//////////////////////////////////////////////////////////////////////
// Get the molecule's mobility
//
double Molecule::get_difc() {
    return diffusion_constant;
}

//////////////////////////////////////////////////////////////////////
// Is the molecule activated?
//
bool Molecule::isActivated(){

    // return activation status of the molecule
    return activated;
}



//////////////////////////////////////////////////////////////////////
// Activate the molecule
//
void Molecule::activate(){

    // change status of molecule to active
    activated = true;
}



//////////////////////////////////////////////////////////////////////
// Activate the activatee molecule
//
void Molecule::activates(Molecule * activatee){
    // change status of activatee to active
    activatee->activate();
}



//////////////////////////////////////////////////////////////////////
// De-activate the molecule
//
void Molecule::deactivate(){
    // change status of molecule to non-active
    activated = false;
}



//////////////////////////////////////////////////////////////////////
// Is the molecule activated?
//
bool Molecule::isActivated(Molecule * activator){
    return activated;

}


//////////////////////////////////////////////////////////////////////
// Determine if a molecule is a "neighbor," which is probabilistic based on distance
// and the threshold of the simulation
bool Molecule::isNear(Molecule * m){

    Vector3D diff;
    // Get relative distances of the 2 molecules
    diff.setDiff(pos, m->pos);

    float dist = diff.magnitude();

    float threshold = 1.5 * resolution;

    // Molecules are not within the attraction distance
    if(dist > threshold) return false;

    float pActivate = (pow(threshold,2) - pow(dist,2)) / pow(threshold,2);

    // Return whether their distances is less than a magnitude of 1
    return prob(pActivate);
}

float Molecule::neighbor_reaction_probability(Molecule * m) {
	Vector3D diff;
	// Get relative distances of the 2 molecules
	diff.setDiff(pos, m->pos);
	float dist = diff.magnitude();
	float threshold = 1.5 * resolution;
	// Molecules are not within the attraction distance
	if(dist > threshold) return 0;
	float pActivate = (pow(threshold,2) - pow(dist,2)) / pow(threshold,2);
	return pActivate;
}

//////////////////////////////////////////////////////////////////////
