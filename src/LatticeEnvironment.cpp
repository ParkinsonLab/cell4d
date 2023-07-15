//===============================================================================
//
// LatticeEnvironment.cpp
// Donny Chan

// The simulation works on a superimposed grid, of xdim * ydim * zdim dimensions 
// (these parameters are located in the sbml input file of the model)
// This file constructs the grid.
// By constructing the grid, we mean it will label each grid by some category, 
// and this category represents boundaries, with rules for particle interaction.

// The extern vars originate in Simulation.cpp, and its fill functions most 
// likely reside in ParameterManager.cpp


#include <math.h>
#include <cassert>

#include "../inc/LatticeEnvironment.h"
//these come from Simulation.cpp
extern Model_t *model;
extern map<int, Species_t *> enzyme_map;
extern vector<string> small_molecule_id_list;
//extern map<int, Species_t *> final_metabolite_map;

extern vector<Compart *> compartments;
// 2018 additions
extern map <string, Species_Attributes *> species_details_map;
extern map<string, Compart> compartments_map;   
extern bool sorted_compart;
extern map<int, string> small_molecule_int_to_id_map;
extern map<string, vector<string>> container_components_map; // stores the vector of molecule names in a container

float default_membrane_offset = 0.5 - 1e-6; // the default offset position of a particle in membrane is 1e-5 of a voxel length

//these are special
extern int MODEL_SHAPE;
extern bool NUCLEUS_ENABLED, PLASMA_MEMBRANE_ENABLED;

bool LatticeEnvironment::isInstantiated = false;
//global var
// struct diffuseInto *** diffusionPatterns;
////////////////////////////////////////////////////////
// 2018 new maps
// extern map <string, int> new_metabolite_index_map;

// Pointer to the Lattice Environment
LatticeEnvironment * LatticeEnvironment :: lEnv;

int NUCLEUS_VOLUME_TO_EXPOSE = 0;

ParameterManager *pm = ParameterManager::getParameterManager(); // this is the only thing that isn't extern

// Dimensions of the Lattice
int xDim = pm->get_int(X_DIM);
int yDim = pm->get_int(Y_DIM);
int zDim = pm->get_int(Z_DIM);

/* The length of time that one timestep represents, in seconds */
double timescale;

// in biological terms, spacescale is length of one voxel, so volume of simulation space is spacescale^3
double spacescale; /* 10^-7 sized lattice cells means that 10 of them makes up about the diameter of an E. coli cell */
/* The distance between neighbouring lattice cells, in metres */

double inaccessible_space_percent;
//////////////////////////////////////////////////////////////////////
//  Constructor for the Lattice Environment
//
LatticeEnvironment::LatticeEnvironment(){
// Dimensions of the Lattice
    xDim = pm->get_int(X_DIM);
    yDim = pm->get_int(Y_DIM);
    zDim = pm->get_int(Z_DIM);
    timescale = pm->get_double(TIMESCALE);
    /* The length of time that one timestep represents, in seconds */

    // in biological terms, spacescale is length of one voxel
    spacescale = pm->get_double(SPACESCALE); /* 10^-7 sized lattice cells means that 10 of them makes up about the diameter of an E. coli cell */
    /* The distance between neighbouring lattice cells, in metres */

    inaccessible_space_percent = pm->get_double(INACCESSIBLE_SPACE_PERCENT);

    int numLatticeSites = xDim * yDim * zDim;
    
    stateMatrix = new Matrix(numLatticeSites, small_molecule_id_list.size());
    updateMatrix = new Matrix(numLatticeSites, small_molecule_id_list.size());
    statisticsMatrix = new Matrix(0,small_molecule_id_list.size());

    // Create the 3 Dimensional array of Lattice Sites
    buildLatticeSpace();

    // Build the list of neighbors for each lattice site
    build_voxel_neighbors();
    buildDiffusionRateTable();

    // for compartment-defined location events, set voxels here after compartments have been decided
    ParameterManager::getParameterManager()->set_event_compartment_voxels();
    
}



//////////////////////////////////////////////////////////////////////
//  Get instance of the Lattice Environment
//
LatticeEnvironment * LatticeEnvironment::getLatticeEnv(){

    // If Lattice Environment already instantiated, then return existing Environment
    // Else create new Lattice Environment and return the pointer
    if(!isInstantiated){

        // Create new Lattice Environment
        LatticeEnvironment * new_lEnv = new LatticeEnvironment();

        // Assign pointer
        lEnv = new_lEnv;

        // Lattice Environment exists
        isInstantiated = true;
    }

    // Return the pointer to the Lattice Environment
    return lEnv;
}



//////////////////////////////////////////////////////////////////////
//  Construct the Lattice Space 
//
void LatticeEnvironment :: buildLatticeSpace(){

    // make a 3D array of LatticeSites (Pointers)
    c_voxel_list = new LatticeSite ***[xDim];

    // Step through each block of the array

    // Construct the X Dimension
    for (int x = 0 ; x < xDim ; x ++){
        c_voxel_list[x] = new LatticeSite ** [yDim];

        // Construct the Y Dimension
        for (int y = 0 ; y < yDim ; y ++){
            c_voxel_list[x][y] = new LatticeSite * [zDim];

            // Construct the Z Dimension
            for(int z = 0 ; z < zDim ; z ++){
                // Create a new LatticeSite at location (i,j,k)
                Point pt(x,y,z);

                // Create a Lattice Site at the specified Point
                c_voxel_list[x][y][z] = new LatticeSite(pt);
            }
        }
    }
    // Construct the environment for the Lattice Space
    buildLatticeEnvironment();

}


void LatticeEnvironment :: build_voxel_neighbors(){
// Donny - June 15, 2018 -> This goes through each of the lattice sites and checks all 26 possible directions for neighboring c-voxels 
//                          to add onto the neighbors map of each c-voxel.
//                          Organized by degree (how much contact the neighbor has with the origin c-voxel)
//                          degree 0: neighbor sharing sides, degree 1: neighbor sharing edges, degree 2: neighbour sharing corners
    for (string & small_mol_itr : small_molecule_id_list) {
        for(int x = 0 ; x < xDim ; x ++){
        for(int y = 0 ; y < yDim ; y ++){
        for(int z = 0 ; z < zDim ; z ++){
            LatticeSite * current = c_voxel_list[x][y][z];
            for(int i = -1 ; i <= 1 ; i ++){
            for(int j = -1 ; j <= 1 ; j ++){
            for(int k = -1 ; k <= 1 ; k ++){
                if(i == 0 && j == 0 && k == 0) continue;
                Point p((float)i,(float)j,(float)k);
                float mag = p.distFromOrigin();
                p.add(*(current->getLocation()));
                int degree = -1;

                if(mag == (float)sqrt(1.0)) degree = 0;
                else if (mag == (float)sqrt(2.0)) degree = 1;
                else if (mag == (float)sqrt(3.0)) degree = 2;

                if(degree >= 0 ){
                    if(withinCell(p) && isCompatibleEnvironment(small_mol_itr,get_compart_name(p))){
                        LatticeSite * neighbor = getLatticeSite(p);
                        current->neighbors[small_mol_itr][degree].push_back(neighbor);
                        // current->neighbors[degree].push_back(neighbor);
                        } // end withinCell
                } // end degree
            } // end k
            } // end j
            } // end i
        } // end z
        } // end y
        } // end x
    }
}

//////////////////////////////////////////////////////////////////////
//  Construct the environment for the Lattice space
//  Iterating through each compartment, the voxels specified in the compartment's lattice definition will be assigned the corresponding compartment label
//  Because this is first-come-first-serve, voxels cannot be reassigned if they already have a compartment label. Order of assignment is significant here.
//
void LatticeEnvironment :: buildLatticeEnvironment(){
    // names of compartments in order that compartments should be assigned. by default, alphanumeric order. if sorted compartments, then sorted in that indexed order
    vector<string> compart_read_order;
    compart_read_order.reserve(compartments_map.size());
    if(sorted_compart) {
        // compart name, index
        vector<std::pair<string, string>> compart_vec;
        compart_vec.reserve(compartments_map.size());
        for(auto & compart_pair : compartments_map) {
            compart_vec.push_back(make_pair(compart_pair.first, compart_pair.second.name));
        }
        // sort by the index
        std::sort(compart_vec.begin(),
                  compart_vec.end(),
                  [](const std::pair<string, string> &lhs, const std::pair<string, string> &rhs)
                  {
                      return lhs.second < rhs.second;
                  });
        // push in the sorted compart names, which is .first
        for(auto & compart_pair: compart_vec) {
            compart_read_order.push_back(compart_pair.first);
        }
    } else {
        for (auto & compart_pair : compartments_map) {
            compart_read_order.push_back(compart_pair.first);
        }
    }

    // draw the compartments in order specified by compart_read_order
    for(auto & compart_key : compart_read_order) {
        Compart& cur_compartment = compartments_map[compart_key];
        string compart_name = cur_compartment.id;
        cout << "Drawing compartment:" << compart_name << endl;
        for(auto & shape_definition : cur_compartment.compart_shape_definitions) {
            if(shape_definition["shape"] == SOLID) {
                // need to saturate the c-voxels, so they don't expand out of the simulation
                // Forces each axis to have a minimum of 0, and max of x/y/zDim-1, then makes the first value smaller than the second
                shape_definition["x1"] = (shape_definition["x1"] < 0) ? 0 : (shape_definition["x1"] >= xDim) ? xDim-1 : shape_definition["x1"];
                shape_definition["x2"] = (shape_definition["x2"] < 0) ? 0 : (shape_definition["x2"] >= xDim) ? xDim-1 : shape_definition["x2"];
                if(shape_definition["x1"] > shape_definition["x2"]) swap(shape_definition["x1"], shape_definition["x2"]);

                shape_definition["y1"] = (shape_definition["y1"] < 0) ? 0 : (shape_definition["y1"] >= yDim) ? yDim-1 : shape_definition["y1"];
                shape_definition["y2"] = (shape_definition["y2"] < 0) ? 0 : (shape_definition["y2"] >= yDim) ? yDim-1 : shape_definition["y2"];
                if(shape_definition["y1"] > shape_definition["y2"]) swap(shape_definition["y1"], shape_definition["y2"]);
                
                shape_definition["z1"] = (shape_definition["z1"] < 0) ? 0 : (shape_definition["z1"] >= zDim) ? zDim-1 : shape_definition["z1"];
                shape_definition["z2"] = (shape_definition["z2"] < 0) ? 0 : (shape_definition["z2"] >= zDim) ? zDim-1 : shape_definition["z2"];
                if(shape_definition["z1"] > shape_definition["z2"]) swap(shape_definition["z1"], shape_definition["z2"]);

                cout << "adding solid with dimensions"
                    << " x:" << shape_definition["x1"] << ", " << shape_definition["x2"]
                    << " y:" << shape_definition["y1"] << ", " << shape_definition["y2"]
                    << " z:" << shape_definition["z1"] << ", " << shape_definition["z2"] << endl;
                for (int x = shape_definition["x1"]; x <= shape_definition["x2"]; x++)
                    for (int y = shape_definition["y1"]; y <= shape_definition["y2"]; y++)
                        for (int z = shape_definition["z1"]; z <= shape_definition["z2"]; z++) {
                            if (c_voxel_list[x][y][z]->get_compart_name().compare("default") != 0) {
                                // cout << "C-voxel " << x << "," << y << "," << z << " already labelled as " << c_voxel_list[x][y][z]->get_compart_name() << ", skipped." << endl;
                                continue;
                            }
                            LatticeSite * ls = c_voxel_list[x][y][z];
                            ls->set_compartment(compart_name);
                            compartments_map[compart_name].list_of_voxels.push_back(ls->getLocation());
                        }
            }
            else if(shape_definition["shape"] == POINT) {
                cout << "adding point at"
                    << " x:" << shape_definition["x1"]
                    << " y:" << shape_definition["y1"]
                    << " z:" << shape_definition["z1"] << endl;
                if (c_voxel_list[shape_definition["x1"]][shape_definition["y1"]][shape_definition["z1"]]->get_compart_name().compare("default") != 0) {
                    // cout << "C-voxel " << x << "," << y << "," << z << " already labelled as " << c_voxel_list[x][y][z]->get_compart_name() << ", skipped." << endl;
                    continue;
                }
                LatticeSite * ls = c_voxel_list[shape_definition["x1"]][shape_definition["y1"]][shape_definition["z1"]];
                ls->set_compartment(compart_name);
                compartments_map[compart_name].list_of_voxels.push_back(ls->getLocation());
            } else {
                cerr << "undefined compartment shape, exit." << endl;
                exit(1);
                continue;
            }
        }
        // if this compartment is a membrane, check parameters to see what axis/side its particles should be locked to
        if(cur_compartment.is_membrane) {
            set_membrane_properties(cur_compartment);
        }
    }
}

// if this compartment is a membrane, check parameters to see what axis/side its particles should be locked to
void LatticeEnvironment::set_membrane_properties(Compart& cur_compartment) {
    if(cur_compartment.axis == "n" || cur_compartment.face == "n") {
        cerr << "something went wrong with membrane compartment declaration, exiting." << endl;
        exit(1);
    }
    
    for(auto & loc : cur_compartment.list_of_voxels) {
        LatticeSite * cur_voxel = c_voxel_list[roundToInt(loc->x)][roundToInt(loc->y)][roundToInt(loc->z)];
        cur_voxel->locked_axis.first = cur_compartment.axis;
        if(cur_compartment.face == "front") {
            if(cur_voxel->locked_axis.first == "x") {
                loc->x = roundToInt(loc->x) - default_membrane_offset;
                cur_voxel->locked_axis.second = loc->x;
            } else if(cur_voxel->locked_axis.first == "y") {
                loc->y = roundToInt(loc->y) - default_membrane_offset;
                cur_voxel->locked_axis.second = loc->y;
            } 
            else if(cur_voxel->locked_axis.first == "z") {
                loc->z = roundToInt(loc->z) - default_membrane_offset;
                cur_voxel->locked_axis.second = loc->z;
            } 
        }
        else if(cur_compartment.face == "back") {
            if(cur_voxel->locked_axis.first == "x") {
                loc->x = roundToInt(loc->x) + default_membrane_offset;
                cur_voxel->locked_axis.second = loc->x;
            } else if(cur_voxel->locked_axis.first == "y") {
                loc->y = roundToInt(loc->y) + default_membrane_offset;
                cur_voxel->locked_axis.second = loc->y;
            } 
            else if(cur_voxel->locked_axis.first == "z") {
                loc->z = roundToInt(loc->z) + default_membrane_offset;
                cur_voxel->locked_axis.second = loc->z;
            } 
        }
    }
}

void LatticeEnvironment :: buildInaccessibleSpace(){
    int num_inaccessible_cells = round(inaccessible_space_percent * xDim * yDim * zDim / 100);
    int c = 0;
    while (c < num_inaccessible_cells){
        int x, y, z;
        x = RandomNG::randInt(0, xDim - 1);
        y = RandomNG::randInt(0, yDim - 1);
        z = RandomNG::randInt(0, zDim - 1);
        if (!c_voxel_list[x][y][z]->get_compart_name().compare("inaccessible") == 0){
            // cout << "Building inaccessible space at c-voxel" << x << y << z << " in compartment " << c_voxel_list[x][y][z]->get_compart_name() << endl;
            c_voxel_list[x][y][z]->set_compartment("inaccessible");
            c++;
        }
    }

    if (MODEL_SHAPE == 2) {
        // For two cubes joined, do xDim=zDim=10 and yDim=24
        for (int x = 0; x < xDim; x++) {
            for (int y = xDim; y < yDim-xDim; y++) {
                for (int z = 0; z < zDim; z++) {
                    if ((x < (int)(xDim*0.3) || x > (int)(xDim*0.6)) || (z < (int)(zDim*0.3) || z > (int)(zDim*0.6))) {
                        if (c_voxel_list[x][y][z]->get_compart_name().compare("default") == 0) {
                            c_voxel_list[x][y][z]->set_compartment("inaccessible");
                        }
                    }
                }
            }
        }
    }

    if (MODEL_SHAPE == 1) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                for (int z = 0; z < zDim; z++) {
                    if ((x-xDim/2+0.5) * (x-xDim/2+0.5) + (y-yDim/2+0.5) * (y-yDim/2+0.5)  + (z-zDim/2+0.5) * (z-zDim/2+0.5) > zDim * zDim / 4) {
                        if (c_voxel_list[x][y][z]->get_compart_name().compare("default") == 0) {
                            c_voxel_list[x][y][z]->set_compartment("inaccessible");
                        }
                    }
                }
            }
        }
    }
}

void LatticeEnvironment :: buildMicrotubules(bool *** accessible_nucleus_lattice, int max_nucleus_to_block){
    int num_inaccessible_cells = inaccessible_space_percent * xDim * yDim * zDim / 100;
    int cells_blocked = 0;

#ifdef RADIAL_MICROTUBULES
    while (cells_blocked < num_inaccessible_cells){
        Point cell_membrane_end;
        cell_membrane_end = pickRandomPointOnCellMembrane();

        cells_blocked += buildRadialCylinder(   cell_membrane_end,
                                                NUCLEUS_RADIUS,
                                                1, 
                                                INACCESSIBLE_ENVIRONMENT, 
                                                num_inaccessible_cells - cells_blocked, 
                                                accessible_nucleus_lattice, 
                                                &max_nucleus_to_block
                                            );
    }
#endif /* RADIAL_MICROTUBULES */

#ifdef RANDOM_MICROTUBULES  
    while (cells_blocked < num_inaccessible_cells){
        Point end1;
        end1.x = RandomNG::randInt(1, xDim - 1);
        end1.y = RandomNG::randInt(1, yDim - 1);
        end1.z = RandomNG::randInt(1, zDim - 1);

        Point end2;
        end2.x = RandomNG::randInt(1, xDim - 1);
        end2.y = RandomNG::randInt(1, yDim - 1);
        end2.z = RandomNG::randInt(1, zDim - 1);

        cells_blocked += buildCylinder  (   end1, 
                                            end2,
                                            1, 
                                            INACCESSIBLE_ENVIRONMENT, 
                                            num_inaccessible_cells - cells_blocked, 
                                            accessible_nucleus_lattice, 
                                            &max_nucleus_to_block
                                        );
    }
#endif /* RANDOM_MICROTUBULES */
}


//////////////////////////////////////////////////////////////////////
//  Set the LatticeSite object at the specified point in the Lattice Space
//
LatticeSite * LatticeEnvironment :: getLatticeSite(Point & p){

    if(withinCell(p)){
        // Lattice Site at Point p (Make sure its an INTEGER value)
        int x = roundToInt(p.x);
        int y = roundToInt(p.y);
        int z = roundToInt(p.z);

        // Return the pointer to the Latticesite at (x,y,z)
        return c_voxel_list[x][y][z];
    }

    // LatticeSite is not within the cell, return NULL
    return (LatticeSite *)0;
}



//////////////////////////////////////////////////////////////////////
//  Get the Environment type at location p
//
string LatticeEnvironment :: get_compart_name(Point & p){

    LatticeSite * l = getLatticeSite(p);
    if((LatticeSite *)l != 0) {
        return getLatticeSite(p)->get_compart_name();
    }
    // location not within the cell
    return "inaccessible";
}

//////////////////////////////////////////////////////////////////////
//  Determine whether the specified molecule type is compatible in 
//  the specified environment
//

// find the list (vector of strings) of allowed compartments of this container
vector <string> LatticeEnvironment::get_compatible_env(string container_name) {
    if(container_components_map[container_name].empty()) {
        string species_name = container_name;
        Species_Attributes * reactant_detail = species_details_map[species_name];
        return reactant_detail->allowable_compartments;
    }
    vector<string> container_species_list = container_components_map[container_name];
    set <string> unique_compartments;
    for(auto & species_name : container_species_list) {
        Species_Attributes * reactant_detail = species_details_map[species_name];
        unique_compartments.insert(reactant_detail->allowable_compartments.begin(), reactant_detail->allowable_compartments.end());
    }
    vector<string> allowed_compartments(unique_compartments.begin(), unique_compartments.end());
    return allowed_compartments;
}

// this finds the union of the compatible env's of the container. if any of the species are allowed to be in this space, then it is allowed
bool LatticeEnvironment :: isCompatibleEnvironment(string container_name, string compartment_name) {
    bool match_environment = false;
    if(container_components_map[container_name].empty()) {
        string species_name = container_name;
        Species_Attributes * reactant_detail = species_details_map[species_name];
        if(find(reactant_detail->allowable_compartments.begin(), reactant_detail->allowable_compartments.end(), compartment_name) != reactant_detail->allowable_compartments.end()){
            match_environment = true;
        }
    }
    vector<string> container_species_list = container_components_map[container_name];
    for(auto & species_name : container_species_list) {
        Species_Attributes * reactant_detail = species_details_map[species_name];
        if(find(reactant_detail->allowable_compartments.begin(), reactant_detail->allowable_compartments.end(), compartment_name) != reactant_detail->allowable_compartments.end()){
            match_environment = true;
            break;
        }
    }

    return match_environment;
}

bool LatticeEnvironment :: isCompatibleEnvironment(vector<Species_Attributes *> species_attribute_vector, string compartment_name) {
    bool match_environment = false;
    for(auto & species_attr : species_attribute_vector) {
        if(find(species_attr->allowable_compartments.begin(), species_attr->allowable_compartments.end(), compartment_name) != species_attr->allowable_compartments.end()){
            match_environment = true;
        }
    }
    return match_environment;
}

bool LatticeEnvironment :: is_membrane(string compartment_name) {
    return compartments_map[compartment_name].is_membrane;
}

//////////////////////////////////////////////////////////////////////
//  Determine whether the provided location is within the defined dimensions
//  of the Lattice space
//
bool LatticeEnvironment :: withinCell(Point & p){

    int x = roundToInt(p.x);
    int y = roundToInt(p.y);
    int z = roundToInt(p.z);
    // Is specified Point within dimensions
    return     withinCell(x,y,z);
}

bool LatticeEnvironment :: withinCell(int x, int y, int z){

    // Is specified Point within dimensions
        return  (x < xDim) && (x >= 0) && (y < yDim) && (y >= 0) && (z < zDim) && (z >= 0);
}

//////////////////////////////////////////////////////////////////////
//  Get a list of available locations for the specified environment
//
vector <Point *> & LatticeEnvironment :: getCompartmentSites(string compartment_name){

    // return a list of available locations with the compartment
    return compartments_map[compartment_name].list_of_voxels;
}

//////////////////////////////////////////////////////////////////////
//  Display the Lattice Environments of each LatticeSite in the Lattice Space
//
void LatticeEnvironment :: displayLatticeEnv(){
    for(int x = 0 ; x < xDim ; x ++)
        for(int y = 0 ; y < yDim ; y ++)
            for(int z = 0 ; z < zDim ; z ++)
                // Output Location and its environment
                printf("(%d,%d,%d) -> %d\n",x,y,z,(c_voxel_list[x][y][z])->get_compart_name());
}



// build table of diffusion constants
void LatticeEnvironment :: buildDiffusionRateTable(){
    for(auto & small_mol : small_molecule_id_list) {
        // fill difc table
        diffusion_coefficient_table[small_mol] = species_details_map[small_mol]->diffusion_constant;
        // diffusion per timestep to face-sharing neighbors is D(TS/ss^2) * current_conc
        // this multiplier * conc of a voxel gives amount to diffuse out to each neighbor
        diff_multiplier[small_mol] = (diffusion_coefficient_table[small_mol] * timescale) / (spacescale * spacescale);
    }
}


// (particles / 6.02e23 particles/mol) / (m^3 * 1000 L/m) = mol/L = M
double LatticeEnvironment::particle_to_conc(double particles) {
    double converted_conc = (particles / AVO_CONST) / (spacescale * spacescale * spacescale * 1000);
    return converted_conc;
}

double LatticeEnvironment::particle_to_mole(double particles) {
    double converted_mole = particles / AVO_CONST;
    return converted_mole;
}

// conc is mol/L. mol/L * (m^3 * 1000 L / m^3) * 6.02e23 particles/mol = particles
// particle counts are approximated, not rounded to int
double LatticeEnvironment::conc_to_particle(double conc) {
    // example: 100 uM in 1e-7 m length voxel turns into 60 particles
    double converted_particle_count = conc * (spacescale * spacescale * spacescale * 1000) * AVO_CONST; // multiplies input by about 6e5 for spacescale of 1e-7
    return converted_particle_count;
}

double LatticeEnvironment::mole_to_particle(double mole) {
    double converted_particle_count = mole * AVO_CONST; // multiplies input by avogadro's constant
    return converted_particle_count;
}

double LatticeEnvironment::conc_to_mole(double conc) {
    double converted_mole = conc * (spacescale * spacescale * spacescale * 1000);
    return converted_mole;
}

double LatticeEnvironment::mole_to_conc(double mole) {
    double converted_conc = mole / (spacescale * spacescale * spacescale * 1000);
    return converted_conc;
}

// molecular flux is amount of substance diffused per unit area per time
// amount to diffuse to neighbors per timestep = Difc * ((origin_conc - neighbor_conc) / distance) * area of interface

// Diffusion happens between an origin voxel and its 26 potential neighbors that share faces, edges, and vertices. However, the use of cubic cellulation
// (tessellation of regular cubes that fills all space) combined with using Fick's laws of diffusion, which requires an interface area for flux, means that diffusion will only
// occur between face-sharing voxels, so 6 neighbors out of 26.

//this is fine, just turn off diffusion for other neighbors.


void LatticeEnvironment :: buildDiffusionPatterns() {

    // Number of Dimensions
    int numDegrees = 1;
    // Step through each lattice site in the system
    for(int i = 0 ; i < xDim ; i ++) {
    for(int k = 0 ; k < yDim ; k ++) {
    for(int m = 0 ; m < zDim ; m ++) {
        LatticeSite * origin = c_voxel_list[i][k][m];
        for (auto & bulk_pair : origin->bulk_moles_map) {
            string bulk_name = bulk_pair.first;
            double current_conc = mole_to_conc(origin->bulk_moles_map[bulk_name]);
            if(current_conc <= 0) continue;
            // amount to distribute to each face-sharing neighbor
            double diffuse_out_conc =  diff_multiplier[bulk_name] * current_conc;

            /////////////////////////////////////////////////////////////////////
            // Calculate diffusion rates to neighbors with each degree
            // This goes through each of the lattice sites and checks all 26 possible directions for neighboring c-voxels 
            // to add onto the neighbors map of each c-voxel.
            // Organized by degree (how much contact the neighbor has with the origin c-voxel)
            // degree 0: neighbor sharing sides, degree 1: neighbor sharing edges, degree 2: neighbour sharing corners

            // Oct 20 2019: only allow 0th degree neighbors
            for(int degree = 0 ; degree < numDegrees ; degree++) {
                // Get the neighbors of the source lattice
                vector<LatticeSite *> * neighborLattices = &(origin->neighbors[bulk_name][degree]);
                int face_neighbors = neighborLattices->size();
                // subtract out all the molecules leaving origin
                origin->conc_change_timestep_map[bulk_name] -= (diffuse_out_conc * face_neighbors);

                for (int num = 0; num < face_neighbors; num++) {
                    // Get the neighboring lattice site
                    LatticeSite * neighbor = (*neighborLattices)[num];
                    // move the molecules to the neighbor (add to tally)
                    neighbor->conc_change_timestep_map[bulk_name] += diffuse_out_conc;
                }
            }
        }
    }
    }
    }

}



void LatticeEnvironment :: diffuse() {

    // figure out how the system-wide diffusive flux is
    // creates conc_change_timestep_map for each voxel, indicating how much of each bulk mol is changed during this timestep
    buildDiffusionPatterns();
        
    for(int i = 0 ; i < xDim ; i ++) {
    for(int k = 0 ; k < yDim ; k ++) {
    for(int m = 0 ; m < zDim ; m ++) {
        LatticeSite * origin = c_voxel_list[i][k][m];
        // ID of the source of diffusion
        int origin_id = origin->getID();
        // Get the source environment
        string source_env = origin->get_compart_name();
        if(source_env == "default") continue;

        // now that diffusive flux is determined, move all small molecules simutaneously during this timestep
        // amend the current voxel's molecule accounting as well as the compartment-wide count
        for(auto & small_mol_change_pair : origin->conc_change_timestep_map) {
            string small_mol = small_mol_change_pair.first;
            origin->add_mole_to_voxel(small_mol, conc_to_mole(origin->conc_change_timestep_map[small_mol]));
        }

        // done with temporary tally for this voxel, clear it for next timestep
        origin->conc_change_timestep_map.clear();
    }
    }
    }
}

// Inserts all small molecules into the simulation
void LatticeEnvironment::add_mole_to_voxel(string molecule_name, double input_moles, LatticeSite * c_voxel){
    c_voxel->bulk_moles_map[molecule_name] += input_moles;
    compartments_map[c_voxel->get_compart_name()].current_bulk_count[molecule_name] += mole_to_particle(input_moles);
}


void LatticeEnvironment :: reactWithMetabolite(Point & point, string reactant, string product, float Ka, float Kp, float equilibrium_constant, float maximum_rate){

	// This function performs the actual reaction on the enzyme level.
    // It'll scan the concentration of reactants for that particular c-voxel, and perform the reaction based on the equation

    LatticeSite * lattice_site = getLatticeSite(point);
    int id = lattice_site->getID();
    double amount = 0.0;
    double Vm = maximum_rate * timescale;
    double reactant_conc = mole_to_conc(lattice_site->bulk_moles_map[reactant]);
    double product_conc = mole_to_conc(lattice_site->bulk_moles_map[product]);

    double Keq = equilibrium_constant;

    // equation for reversible michaelis menten is rate = ((Vmaxf * reactant_c / Ka) - (Vmaxb * prod_c / Kp)) / (1 + reactant_c / Ka + prod_c / Kp)
    // converts into rate = ((Vm / Ka) * (reactant_c - (prod_c / Keq)) / (1 + (reactant_c / Ka) + (prod_c / Kp)))

    const static double MINIMUM_CONCENTRATION = 1.0e-10;
    // cout << "before - reactant conc:" << reactant_conc << " prod conc:" << product_conc << " in voxel " << id << endl;
    if (reactant_conc >= MINIMUM_CONCENTRATION) {
        if (product_conc >= MINIMUM_CONCENTRATION) {
            amount = Vm * ((reactant_conc/Ka)*(1-((product_conc/reactant_conc)/Keq))) / (1 + (reactant_conc/Ka) + (product_conc/Kp));
        } else {
            amount = Vm * (reactant_conc/Ka) / (1 + (reactant_conc/Ka));
        }
      //forward_amount = Vm * ((a/Ka)*(1-((MAR)/Keq))) / (1 + (a/Ka) + (p/Kp));
    } else {
        amount = 0.0;
    }
    if(amount == 0)    return;
    // cout << "reaction:" << Vm << "  " << ((reactant_conc/Ka)*(1-((product_conc/reactant_conc)/Keq))) / (1 + (reactant_conc/Ka) + (product_conc/Kp)) << endl;
    double cube_volume = spacescale * spacescale * spacescale * 1000; /* volume in cubic metres then litres */
    amount /= cube_volume; /* divide by VOLUME_OF_CUBE to get concentration change */;

    // cout << "pre-amount:" << amount << " cube volume:" << cube_volume << endl;
    if (amount > reactant_conc) {
        amount = reactant_conc;
    } else if (-amount > product_conc) {
        amount = -product_conc;
    }

    /* This is equation (A1) from page 5326 of Eur. J. Biochem. 267, (2000).
     * It is from the article "Can yeast glycolysis be understood in terms of in vitro
     * kinetics of the constituent enzymes? Testing biochemistry."
     */

    // cout << "amount changed:" << amount << endl;

    lattice_site->add_mole_to_voxel(reactant, -conc_to_mole(amount));
    lattice_site->add_mole_to_voxel(product, conc_to_mole(amount));


}

void LatticeEnvironment :: displayState(){
    stateMatrix->display();
}

Matrix * LatticeEnvironment :: getStateMatrix(){
    return     stateMatrix;
}

Matrix * LatticeEnvironment :: getStatisticsMatrix(){
    return     statisticsMatrix;
}

Tuple * LatticeEnvironment :: getCurStat(){
    int num = statisticsMatrix->getNumRows();
    if(num > 0)
        return statisticsMatrix->matrix[num-1];
    return (Tuple *)0;
}
