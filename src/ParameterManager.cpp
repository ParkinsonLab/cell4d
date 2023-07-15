//==================================================
// ParameterManager.cpp
// Donny Chan, Billy Taj

#include "../inc/ParameterManager.h"
#include "../inc/definition.h"
#include "../inc/Compartment.h"
#include "../inc/Cell_Reaction.h"
#include "../inc/Molecule_Container.h"
#include "../inc/Container_Rules.h"


extern SBMLDocument_t *sbml_document;
extern Model_t *model;
extern map<int, Species_t *> enzyme_map;
extern map<int, Species_t *> metabolite_map;
extern string classInfo;
extern int numMetabolite;
extern Reaction_t **rxns;
extern map<int, Species_t *> reactants_map;
extern map<int, Species_t *> products_map;
extern map<int, KineticLaw_t *> kinetic_law_map;
extern map<int, Species_t *> initial_metabolite_map;
extern map<int, Species_t *> final_metabolite_map;
extern map<int, Species_t *> final_signal_enzyme_map;
extern vector<string> small_molecule_id_list;
extern vector<Compart *> compartments;

// 2018 edits
extern map<string, vector<Species_t*>> full_products_map;  
extern map<string, vector<Species_t*>> full_modifiers_map; // used as a filter on whether or not to consure the reactants
extern map<string, Species_Attributes *> species_details_map;
extern map<string, Compart> compartments_map; // stores all the details of the compartments. key = compartment ID
extern map<int, string> small_molecule_int_to_id_map; // July 17 2018 - store small molecule names as ints, used for tuple
extern map<string, vector<string>> container_components_map; // stores the vector of molecule names in a container
extern map<Container_Config, Container_Rules> container_rules_map;

// TODO: get rid of these extern vars. they're persisting because there's no resources to fundamentally change the
// organization

//
extern map<string, Cell_Reaction *> bulk_reaction_details_map;
// stores details about unimolecular reactions (products, reactants, etc.etc.)
extern map<string, Cell_Reaction *> unimolecular_reaction_details_map;
// stores details about bimolecular reactions (products, reactants, etc.etc.)
extern map<string, Cell_Reaction *> bimolecular_reaction_details_map; 
extern map<string, Cell_Reaction *> enzymatic_reaction_details_map; // protein_id as the key, stores vector of all reaction names it participates in as a modifier
// State of the Parameter Manager
bool ParameterManager::isInstantiated = false;
int numSpecies;
json checkpoint_json;
extern int new_timestep;
extern double new_timescale;
extern bool poslog;
extern vector<string> poslog_comparts;
extern map<string, json> poslog_json;

bool sorted_compart = false;

// Pointer to the Parameter Manager
ParameterManager * ParameterManager::params;

//////////////////////////////////////////////////////////////////////
// Get an instance of the Parameter Manager
//
ParameterManager * ParameterManager::getParameterManager(){

    // Create new Parameter Manager if not yet exist
    if(!isInstantiated){

        // Create Parameter

        ParameterManager * new_param = new ParameterManager();

        // Assign Pointer to existing Parameter Manager
        params = new_param;

        // Parameter Manager created
        isInstantiated = true;
    }

    // Return pointer to the current Parameter Manager
    return params;
}

//////////////////////////////////////////////////////////////////////
// Constructor for the Parameter Manager
// Population the parameter table
//
ParameterManager::ParameterManager(){
    // Load default parameter file
    //load(INPUT_FILE_PARAMS);
    loaded = false;
}


void ParameterManager::loadSBML(const char* fName, bool load_checkpoint, string checkpoint_filename){

    const char *checkMoleculeType = NULL;
    if (loaded){
        return;
    }
    sbml_document = readSBML(fName);
    if (sbml_document->getNumErrors() > 0) {
        sbml_document->printErrors(cout);
        exit(1);
    }

    // if loading from checkpoint, do setup of json loading into memory
    if(load_checkpoint) {
        ifstream check_json(checkpoint_filename);
        if(!json::accept(check_json)) {
            cerr << "Incorrect syntax or non-existent json file. exiting." << endl;
            exit(1);
        }
        ifstream json_input_string(checkpoint_filename);
        checkpoint_json = json::parse(json_input_string);
    }


    model = SBMLDocument_getModel(sbml_document);
    string annot = SBase_getAnnotationString(model);
    rapidxml::xml_document<> annot_level_doc;
    annot_level_doc.parse<0>(&annot[0]); // parse annotations here so it's only done once

    // this version determines if the compartment write order is specified by the name attribute. If missing, verison will = 0.
    int version = extract_int("XMLversion", annot);;
    if(version == 2) {
        sorted_compart = true;
    };
    // set model parameters such as dimensions, scales
    set_global_parameters(annot);

    // numSpecies listing all configurations of molecule containers (and small molecules?)
    numSpecies = Model_getNumSpecies(model);

    // take apart the gigantic cell list of species portion in the model-layer annotation
    set_annot_species(annot_level_doc);

    // grab the compartments info and feed it into the compartments_map
    set_compartments(model);

    set_species(model, load_checkpoint);

    set_reactions(model);

    // try to read any optional sim-time events in the model
    set_events(annot_level_doc);

    
    loaded = true;
    cout << "XML input file loaded." << endl;
}


void ParameterManager::set_global_parameters(string & annot) {

    parameterTable[X_DIM] = extract_int(X_DIM, annot);
    parameterTable[Y_DIM] = extract_int(Y_DIM, annot);
    parameterTable[Z_DIM] = extract_int(Z_DIM, annot);
    parameterTable[MAX_CYCLES] = extract_int(MAX_CYCLES, annot);
    parameterTable[TIMESCALE] = extract_double(TIMESCALE, annot);
    parameterTable[SPACESCALE] = extract_double(SPACESCALE, annot);
    parameterTable[INACCESSIBLE_SPACE_PERCENT] = extract_double(INACCESSIBLE_SPACE_PERCENT, annot);

    // if new_timestep isn't the default value, set to new value from cml input
    if(new_timestep > 0) {
        parameterTable[MAX_CYCLES] = new_timestep;
    }

    if(new_timescale != -1) {
        parameterTable[TIMESCALE] = new_timescale;
    }

    cout << "timescale: " << parameterTable[TIMESCALE] << endl;
    cout << "spacescale: " << parameterTable[SPACESCALE] << endl;
    cout << "numCycles: " << parameterTable[MAX_CYCLES] << endl;
}

void ParameterManager::set_annot_species(const rapidxml::xml_document<> & annotation_xml_doc) {

    // annotation -> listOfAnnotationSpeciesTypes
    for(rapidxml::xml_node<>* lvl_1 =  annotation_xml_doc.first_node()->first_node(); lvl_1; lvl_1 = lvl_1->next_sibling()){
        //
        string lvl_1_name = lvl_1->name();
        if(lvl_1_name.compare("cell4d:listOfAnnotationSpeciesTypes") == 0){
            // loop through each annotation species type
            for(rapidxml::xml_node<>* SpeciesType_node = lvl_1->first_node(); SpeciesType_node; SpeciesType_node = SpeciesType_node->next_sibling()){
                // this is where the species details are stored
                Species_Attributes * new_species_attr = new Species_Attributes;
                new_species_attr->ID = SpeciesType_node->first_attribute("id")->value();
                // commented out unused attribute "speciesTypeName"
                // new_species_attr->name = SpeciesType_node->first_attribute("speciesTypeName")->value();
                new_species_attr->species_type = SpeciesType_node->first_attribute("speciesMoleculeType")->value();
                if(new_species_attr->species_type.compare("SIMPLE_MOLECULE") == 0) {
                    small_molecule_id_list.push_back(new_species_attr->ID);
                    numMetabolite++;
                }
                // SpeciesType -> ValidCompartments/BindingSites
                bool found_difc = false;
                for(rapidxml::xml_node<>* SpeciesType_child = SpeciesType_node->first_node(); SpeciesType_child; SpeciesType_child = SpeciesType_child->next_sibling()) {
                    string SpeciesType_child_name = SpeciesType_child->name();
                    if(SpeciesType_child_name.compare("cell4d:listOfValidCompartments") == 0) {
                        // listOfCompartments -> compartment node
                        for(rapidxml::xml_node<>* compart_node = SpeciesType_child->first_node(); compart_node; compart_node = compart_node->next_sibling()){
                            string compartment_id = compart_node->first_attribute("id")->value();
                            new_species_attr->allowable_compartments.push_back(compartment_id);
                        }
                    } else if(SpeciesType_child_name.compare("cell4d:listOfBindingSites") == 0) {
                        // listOfBindingSites -> bindingSite node
                        for(rapidxml::xml_node<>* bindingSite_node = SpeciesType_child->first_node(); bindingSite_node; bindingSite_node = bindingSite_node->next_sibling()) {
                            string bind_site_name = bindingSite_node->first_attribute("id")->value();
                            if(!bindingSite_node->first_node()) {
                                new_species_attr->binding_sites.push_back(bind_site_name);
                            } else {
                                // bindingSite -> listOfPossibleStates -> state node
                                for(rapidxml::xml_node <> * state_node = bindingSite_node->first_node()->first_node(); state_node; state_node = state_node->next_sibling()) {
                                    new_species_attr->possible_states[bind_site_name].push_back(state_node->first_attribute("value")->value());
                                }
                            }
                        }
                    } else if(SpeciesType_child_name.compare("cell4d:diffusionConstant") == 0) {
                        rapidxml::xml_attribute<>* difc_attr = SpeciesType_child->first_attribute("value");
                        if(difc_attr) {
                            found_difc = true;
                            new_species_attr->diffusion_constant = atof(difc_attr->value());
                        }
                    } else if(SpeciesType_child_name.compare("cell4d:displayProperties") == 0) {
                        new_species_attr->color_red = stoi(SpeciesType_child->first_attribute("redValue")->value());
                        new_species_attr->color_green = stoi(SpeciesType_child->first_attribute("greenValue")->value());
                        new_species_attr->color_blue = stoi(SpeciesType_child->first_attribute("blueValue")->value());
                        // optional field for density gradient display on membranes rather than particles
                        if(SpeciesType_child->first_attribute("membraneDisplay")) {
                            string display_value = SpeciesType_child->first_attribute("membraneDisplay")->value();
                            if(display_value == "true") {
                                new_species_attr->membrane_display = true;
                            }
                        }
                    } else if(SpeciesType_child_name.compare("cell4d:compartmentEmissionRate") == 0) {
                        string origin_compartment = SpeciesType_child->first_attribute("origin")->value();
                        string target_compartment = SpeciesType_child->first_attribute("target")->value();
                        float probability = stof(SpeciesType_child->first_attribute("probability")->value());
                        new_species_attr->compartment_emission_rates[origin_compartment][target_compartment] = probability;
                    }
                }
                if(!found_difc) {
                    new_species_attr->diffusion_constant = 3.0e-11; // default diffusion constant
                }                // put the species attribute in the map
                species_details_map[new_species_attr->ID] = new_species_attr;
            }
        }
    }

    // Creating small mol name to int map
    sort(small_molecule_id_list.begin(),small_molecule_id_list.end());
    for(int i = 0; i < small_molecule_id_list.size(); i++) {
        small_molecule_int_to_id_map[i] = small_molecule_id_list[i];
    }

}

void ParameterManager::set_events(const rapidxml::xml_document<> & annotation_xml_doc) {

    // annotation -> events
    rapidxml::xml_node<>* events_node =  annotation_xml_doc.first_node()->first_node("cell4d:events");
    if(!events_node) return;
    
    // go through every event listed
    for(rapidxml::xml_node<>* event_node = events_node->first_node("event"); event_node; event_node = event_node->next_sibling()) {
        bool read_fail = false;
        if(!event_node->first_attribute("name") || !event_node->first_attribute("type") || !event_node->first_attribute("trigger") || !event_node->first_attribute("probability")) {
            read_fail = true;
            cerr << "event read fail at parent node" << endl;
            exit(1);
        }
        Event_Parameters * cur_event = new Event_Parameters;
        cur_event->event_name = event_node->first_attribute("name")->value();
        cur_event->event_type = event_node->first_attribute("type")->value();
        cur_event->event_trigger = event_node->first_attribute("trigger")->value();
        cur_event->event_probability = stod(event_node->first_attribute("probability")->value());

        // ensure the event type is one of the valid ones: add_mols, remove_mols, transport_mols
        if(cur_event->event_type == "add_mols" || cur_event->event_type == "remove_mols" || cur_event->event_type == "transport_mols") {
            string container_name;
            if(event_node->first_node("add_mols")) {
                container_name = event_node->first_node("add_mols")->first_attribute("id")->value();
                cur_event->event_container_amount_total = stoi(event_node->first_node("add_mols")->first_attribute("amount")->value());
            } else if(event_node->first_node("remove_mols")) {
                container_name = event_node->first_node("remove_mols")->first_attribute("id")->value();
                cur_event->event_container_amount_total = stoi(event_node->first_node("remove_mols")->first_attribute("amount")->value());
            } else if(event_node->first_node("transport_mols")) {
                container_name = event_node->first_node("transport_mols")->first_attribute("id")->value();
                cur_event->event_container_amount_total = stoi(event_node->first_node("transport_mols")->first_attribute("amount")->value());
                if(find(small_molecule_id_list.begin(), small_molecule_id_list.end(), container_name) != small_molecule_id_list.end()) {
                    cerr << "cannot transport bulk molecules" << endl;
                    exit(1);
                }
                rapidxml::xml_node<>* destination_node;
                if(!event_node->first_node("transport_destination")) {
                    cerr << "missing transport destination nodes" << endl;
                    exit(1);
                }
                for(destination_node = event_node->first_node("transport_destination"); destination_node; destination_node = destination_node->next_sibling("transport_destination")) {
                    if(!destination_node->first_attribute("xloc") || !destination_node->first_attribute("yloc") || !destination_node->first_attribute("zloc")) {
                        if(!destination_node->first_attribute("compartment")) {
                            read_fail = true;
                            cerr << "missing transport destination coordinates and compartment" << endl;
                            exit(1);
                        }
                        cur_event->transport_destination_compartments.push_back(destination_node->first_attribute("compartment")->value());
                        // setting compartment name allows the rest of the event location info to be filled by set_event_compartment_voxels() called in lEnv

                    } else {
                        int xloc = stoi(destination_node->first_attribute("xloc")->value());
                        int yloc = stoi(destination_node->first_attribute("yloc")->value());
                        int zloc = stoi(destination_node->first_attribute("zloc")->value());
                        Point new_point = Point(xloc, yloc, zloc);
                        cur_event->transport_destination_voxels.push_back(new_point);
                    }
                }

            }
            auto rule_map_itr = find_if(list_of_species_map.begin(), list_of_species_map.end(), [container_name](pair<string, vector<Species_Attributes *>> attr_pair) {
                return(attr_pair.first == container_name);
            });
            if(rule_map_itr == list_of_species_map.end()) {
                read_fail = true;
                cerr << "event read fail, cannot find event molecule in list of species" << endl;
                exit(1);
            }
            cur_event->event_container_species_name = rule_map_itr->first;

        } else {
            cerr << "event read fail at event type name" << endl;
            read_fail = true;
            exit(1);
        }

        if(cur_event->event_trigger == "time") {
            double interval, initial_time;
            // required attributes: repeat, initial
            // optional: interval (dependent on repeat), end
            if(!event_node->first_node("time_trigger")->first_attribute("repeat") || !event_node->first_node("time_trigger")->first_attribute("initial")) {
                read_fail = true;
                cerr << "event read fail at time trigger node" << endl;
                exit(1);
            }
            string event_repeat = event_node->first_node("time_trigger")->first_attribute("repeat")->value();
            if(event_repeat == "true" || event_repeat == "True" || event_repeat == "TRUE") {
                cur_event->event_repeat=true;
                if(!event_node->first_node("time_trigger")->first_attribute("interval")) {
                    cerr << "event read fail, must have interval for repeating events" << endl;
                    exit(1);
                }
                interval = stod(event_node->first_node("time_trigger")->first_attribute("interval")->value());
                if(interval < 1) interval = interval / parameterTable[TIMESCALE]; // if decimal, treat as real seconds
                cur_event->event_interval = interval;
                // parse optional attribute end time
                if(event_node->first_node("time_trigger")->first_attribute("end")) {
                    cur_event->event_final_time = stoi(event_node->first_node("time_trigger")->first_attribute("end")->value());
                }
            } else if(event_repeat == "false" || event_repeat == "False" || event_repeat == "FALSE") {
                cur_event->event_repeat=false;

            } else {
                read_fail = true;
                cerr << "event read fail, repeat must be true or false" << endl;
                exit(1);
            }
            initial_time = stod(event_node->first_node("time_trigger")->first_attribute("initial")->value());
            if(initial_time < 1) initial_time = initial_time / parameterTable[TIMESCALE]; // if decimal, treat as real seconds
            cur_event->event_initial_time = initial_time;

            if(!event_node->first_node("time_trigger")->first_attribute("initial")) {
                read_fail = true;
                cerr << "event read fail at time trigger node, missing initial ts" << endl;
                exit(1);
            }

        } else if(cur_event->event_trigger == "event") {
            if(!event_node->first_node("event_trigger")->first_attribute("name")) {
                read_fail = true;
                cerr << "event read fail, event_trigger's event must be named" << endl;
                exit(1);
            }
            string event_trigger_name = event_node->first_node("event_trigger")->first_attribute("name")->value();
            cur_event->event_triggering_event = event_trigger_name;
            double delay = 0;
            if(event_node->first_node("event_trigger")->first_attribute("delay")) {
                delay = stod(event_node->first_node("event_trigger")->first_attribute("delay")->value());
                if(delay < 1) delay = delay / parameterTable[TIMESCALE]; // if decimal, treat as real seconds
            }
            cur_event->event_interval = delay;

        } else if(cur_event->event_trigger == "state") { // process state-triggered events
            rapidxml::xml_node<>* state_trigger_node = event_node->first_node("state_trigger");
            if(!state_trigger_node->first_attribute("condition") || !state_trigger_node->first_attribute("id") || !state_trigger_node->first_attribute("amount")) {
                cerr << "missing state trigger node attributes" << endl;
                exit(1);
            }
            string event_condition = state_trigger_node->first_attribute("condition")->value();
            if(event_condition == "greater_than") {
                cur_event->event_triggering_condition = 1;
            } else if(event_condition == "less_than") {
                cur_event->event_triggering_condition = 0;
            } else {
                cerr << "invalid state trigger condition, must be \"greater_than\" or \"less_than\"" << endl;
            }

            string mol_name = state_trigger_node->first_attribute("id")->value();
            cur_event->state_trigger_molecule = mol_name;
            auto rule_map_itr = find_if(list_of_species_map.begin(), list_of_species_map.end(), [mol_name](pair<string, vector<Species_Attributes *>> attr_pair) {
                return(attr_pair.first == mol_name);
            });
            if(rule_map_itr == list_of_species_map.end()) {
                read_fail = true;
                cerr << "state trigger read fail, cannot find trigger molecule in list of species" << endl;
                exit(1);
            }

            string event_repeat = state_trigger_node->first_attribute("repeat")->value();
            if(event_repeat == "true" || event_repeat == "True" || event_repeat == "TRUE") {
                cur_event->event_repeat=true;
            } else if(event_repeat == "false" || event_repeat == "False" || event_repeat == "FALSE") {
                cur_event->event_repeat=false;
            } else {
                read_fail = true;
                cerr << "event state trigger read fail, repeat must be true or false" << endl;
                exit(1);
            }
            
            cur_event->state_trigger_mol_amount = stoi(state_trigger_node->first_attribute("amount")->value());

            if(!state_trigger_node->first_attribute("interval")) {
                cerr << "missing state trigger node check interval" << endl;
                exit(1);
            } 
            cur_event->event_interval = stoi(state_trigger_node->first_attribute("interval")->value());
            if(cur_event->event_interval < 1) {
                cerr << "state check interval must >= 1" << endl;
                exit(1);
            }
            if(state_trigger_node->first_attribute("initial")) {
                cur_event->event_initial_time = stoi(state_trigger_node->first_attribute("initial")->value());
            }

            // parse optional attribute end time
            if(state_trigger_node->first_attribute("end")) {
                cur_event->event_final_time = stoi(state_trigger_node->first_attribute("end")->value());
            }
            
            rapidxml::xml_node<>* state_location_node;
            if(!event_node->first_node("state_trigger_loc")) {
                cerr << "missing state trigger location nodes" << endl;
                exit(1);
            }
            for(state_location_node = event_node->first_node("state_trigger_loc"); state_location_node; state_location_node = state_location_node->next_sibling("state_trigger_loc")) {
                if(!state_location_node->first_attribute("xloc") || !state_location_node->first_attribute("yloc") || !state_location_node->first_attribute("zloc")) {
                    if(state_location_node->first_attribute("compartment")) {
                        cur_event->state_trigger_compartment.push_back(state_location_node->first_attribute("compartment")->value());
                        // setting compartment name allows the rest of the event location info to be filled by set_event_compartment_voxels() called in lEnv
                    } else if(state_location_node->first_attribute("xloc_1") && state_location_node->first_attribute("xloc_2")) {
                        if(!state_location_node->first_attribute("yloc_1") || !state_location_node->first_attribute("yloc_2") || 
                        !state_location_node->first_attribute("zloc_1") || !state_location_node->first_attribute("zloc_2")) {
                            cerr << "missing attributes in state trigger location declaration in event " << cur_event->event_name << endl;
                            exit(1);
                        }
                        int x1 = stoi(state_location_node->first_attribute("xloc_1")->value());
                        int x2 = stoi(state_location_node->first_attribute("xloc_2")->value());
                        int y1 = stoi(state_location_node->first_attribute("yloc_1")->value());
                        int y2 = stoi(state_location_node->first_attribute("yloc_2")->value());
                        int z1 = stoi(state_location_node->first_attribute("zloc_1")->value());
                        int z2 = stoi(state_location_node->first_attribute("zloc_2")->value());
                        for(int x = x1; x <= x2; x++) {
                        for(int y = y1; y <= y2; y++) {
                        for(int z = z1; z <= z2; z++) {
                            Point new_point = Point(x, y, z);
                            cur_event->state_trigger_voxels.push_back(new_point);
                        }
                        }
                        }
                    } else {
                        read_fail = true;
                        cerr << "missing state location coordinates and compartment" << endl;
                        exit(1);
                    }

                } else {
                    if(state_location_node->first_attribute("xloc_1") || state_location_node->first_attribute("xloc_2") || 
                        state_location_node->first_attribute("yloc_1") || state_location_node->first_attribute("yloc_2") ||
                        state_location_node->first_attribute("zloc_1") || state_location_node->first_attribute("zloc_2")) {
                            cerr << "state location definition should have either xloc or xloc_1, not both. exiting" << endl;
                            exit(1);
                    }
                    int xloc = stoi(state_location_node->first_attribute("xloc")->value());
                    int yloc = stoi(state_location_node->first_attribute("yloc")->value());
                    int zloc = stoi(state_location_node->first_attribute("zloc")->value());
                    Point new_point = Point(xloc, yloc, zloc);
                    cur_event->state_trigger_voxels.push_back(new_point);
                }
            }


        } else {
            read_fail = true;
            exit(1);
        }

        // process the location nodes, where the events will occur
        rapidxml::xml_node<>* location_node = event_node->first_node("location");
        if(!location_node) {
            cerr << "missing event location node" << endl;
            exit(1);
        } 
        int voxel_counter = 0;
        for(rapidxml::xml_node<>* location_node = event_node->first_node("location"); location_node; location_node = location_node->next_sibling("location")) {
            if(!location_node->first_attribute("xloc") || !location_node->first_attribute("yloc") || !location_node->first_attribute("zloc")) {
                if(location_node->first_attribute("compartment")) {
                    cur_event->event_compartment = location_node->first_attribute("compartment")->value();
                    // setting compartment name allows the rest of the event location info to be filled by set_event_compartment_voxels() called in lEnv

                } else if(location_node->first_attribute("xloc_1") && location_node->first_attribute("xloc_2")) {
                    if(!location_node->first_attribute("yloc_1") || !location_node->first_attribute("yloc_2") || 
                    !location_node->first_attribute("zloc_1") || !location_node->first_attribute("zloc_2")) {
                        cerr << "missing attributes in location declaration in event " << cur_event->event_name << endl;
                        exit(1);
                    }
                    int x1 = stoi(location_node->first_attribute("xloc_1")->value());
                    int x2 = stoi(location_node->first_attribute("xloc_2")->value());
                    int y1 = stoi(location_node->first_attribute("yloc_1")->value());
                    int y2 = stoi(location_node->first_attribute("yloc_2")->value());
                    int z1 = stoi(location_node->first_attribute("zloc_1")->value());
                    int z2 = stoi(location_node->first_attribute("zloc_2")->value());
                    for(int x = x1; x <= x2; x++) {
                    for(int y = y1; y <= y2; y++) {
                    for(int z = z1; z <= z2; z++) {
                        Point new_point = Point(x, y, z);
                        cur_event->event_molecule_placement_pairs.push_back(make_pair(new_point, 0));
                        cur_event->event_placement_voxels.push_back(new_point);
                        voxel_counter++;
                    }
                    }
                    }
                } else {
                    read_fail = true;
                    cerr << "missing state location coordinates and compartment" << endl;
                    exit(1);
                }

            } else {
                if(location_node->first_attribute("xloc_1") || location_node->first_attribute("xloc_2") || 
                    location_node->first_attribute("yloc_1") || location_node->first_attribute("yloc_2") ||
                    location_node->first_attribute("zloc_1") || location_node->first_attribute("zloc_2")) {
                        cerr << "location definition should have either xloc or xloc_1, not both. exiting" << endl;
                        exit(1);
                }
                int xloc = stoi(location_node->first_attribute("xloc")->value());
                int yloc = stoi(location_node->first_attribute("yloc")->value());
                int zloc = stoi(location_node->first_attribute("zloc")->value());
                Point new_point = Point(xloc, yloc, zloc);
                cur_event->event_molecule_placement_pairs.push_back(make_pair(new_point, 0));
                cur_event->event_placement_voxels.push_back(new_point);
                voxel_counter++;
            }
        }
        
        // deal with event amounts not divisible by number of voxels
        if(voxel_counter > 0) {
            int total_amount = cur_event->event_container_amount_total;
            int remainder_amount = total_amount % voxel_counter;
            int spread_amount = (total_amount - remainder_amount) / voxel_counter;
            for(auto & voxel_pair : cur_event->event_molecule_placement_pairs) {
                voxel_pair.second = spread_amount;
            }
            cur_event->event_molecule_placement_pairs[0].second += remainder_amount;
        }

        event_map[cur_event->event_name] = cur_event;
    }

    // link together event dependencies
    // if B is triggered by A, A knows that it triggers B
    for(auto & event_pair : event_map) {
        Event_Parameters * event = event_pair.second;
        if(event->event_trigger == "event") {
            event_map[event->event_triggering_event]->events_to_trigger.push_back(event);
        }
    }

}

// set the compartment-based event voxel info. this is separated from list of voxel definitions in XML because the compartments have not yet been defined.
// this will be called right after compartments have been decided in latenv.cpp
void ParameterManager::set_event_compartment_voxels() {
    for(auto & event : event_map) {
        if(!event.second->event_compartment.empty() && event.second->event_molecule_placement_pairs.empty()) {
            // split compartment-wide counts to each voxel, put remainders in first one in list 
            int total_amount = event.second->event_container_amount_total;
            int num_points = compartments_map[event.second->event_compartment].list_of_voxels.size();
            int remaining_amount = event.second->event_container_amount_total % num_points;
            if(remaining_amount != 0) {
                total_amount = event.second->event_container_amount_total - remaining_amount;
            }
            int spread_amount = total_amount / num_points;
            for(auto cur_point : compartments_map[event.second->event_compartment].list_of_voxels) {
                // event_molecule_placement_pairs only used for adding mols
                // event_placement_voxels used for removing mols, picks randomly for containers, and for bulk evenly removes from selected
                event.second->event_molecule_placement_pairs.push_back(make_pair(*cur_point, spread_amount));
                event.second->event_placement_voxels.push_back(*cur_point);
            }
            // try to spread remainders evenly across all voxels
            while(remaining_amount != 0) {
                for(auto & point_val_pair : event.second->event_molecule_placement_pairs) {
                    point_val_pair.second += 1;
                    remaining_amount -= 1;
                    if(remaining_amount == 0) break;
                }
            }
            // remove voxels where 0 particles are being placed from the list
            // associative containers like maps return next iterator after an erase
            for(auto it = event.second->event_molecule_placement_pairs.begin(); it != event.second->event_molecule_placement_pairs.end();) {
                if (it->second == 0) {
                    it = event.second->event_molecule_placement_pairs.erase(it);
                } else {
                    it++;
                }
            }
        }

        // build state trigger location compartment info
        // allows multiple compartments to be used
        if(!event.second->state_trigger_compartment.empty() && event.second->state_trigger_voxels.empty()) {
            for(auto & compartment : event.second->state_trigger_compartment)
            for(auto cur_point : compartments_map[compartment].list_of_voxels) {
                event.second->state_trigger_voxels.push_back(*cur_point);
            }
        }

        // do the same for transport destination voxel locations, if only compartment was specified
        if(!event.second->transport_destination_compartments.empty() && event.second->transport_destination_voxels.empty()) {
            for(auto & compartment : event.second->transport_destination_compartments)
            for(auto cur_point : compartments_map[compartment].list_of_voxels) {
                event.second->transport_destination_voxels.push_back(*cur_point);
            }
        }
    }
}


void ParameterManager::set_compartments(Model_t * model) {
    int number_of_compartments = Model_getNumCompartments(model);
    cout << "compartment_count: " << number_of_compartments << endl;
    for(int i = 0; i < number_of_compartments; i++){
        Compart new_compartment;
        Compartment * compartment = Model_getCompartment(model, i);
        new_compartment.id = Compartment_getId(compartment);
        string compartment_id = Compartment_getId(compartment);
        if(compartment->isSetName()){
            new_compartment.name = Compartment_getName(compartment);
        } else{
        }
        if(compartment->isSetOutside()){
            new_compartment.outside = Compartment_getOutside(compartment);
        } else {
        }
        if(compartment->isSetAnnotation()){
			string compartment_annotation = SBase_getAnnotationString(compartment);
            rapidxml::xml_document<>compartment_doc;
            compartment_doc.parse<0>(&compartment_annotation[0]);
            for(rapidxml::xml_node<>* annot_node = compartment_doc.first_node(); annot_node; annot_node = annot_node->next_sibling()) {
                // get all lattice definitions of this compartment, each def will declare multiple distinct regions as the same compartment
                for(auto lat_def_node = annot_node->first_node("cell4d:latticePointDefinition"); lat_def_node; lat_def_node = lat_def_node->next_sibling("cell4d:latticePointDefinition")) {
                    map<string, int> current_lat_point_def;
                    // navigate through each attribute of lattice definition
                    for(rapidxml::xml_attribute<>* compart_attribute = lat_def_node->first_attribute(); compart_attribute; compart_attribute = compart_attribute->next_attribute()) {
                        string attr_label = compart_attribute->name();
                        string attr_value = compart_attribute->value();
                        if(attr_label.compare("type") == 0) {
                            // using int identifiers for the type, 0 is rectangular, 1 is point, 2 is spheroid.
                            if(attr_value.compare("solid") == 0) current_lat_point_def["shape"] = SOLID;
                            else if(attr_value.compare("point") == 0) current_lat_point_def["shape"] = POINT;
                            else if(attr_value.compare("spheroid") == 0) current_lat_point_def["shape"] = SPHEROID;
                        } else if((attr_label.compare("x-fixed") == 0) || (attr_label.compare("y-fixed") == 0) || (attr_label.compare("z-fixed") == 0)){
                            current_lat_point_def["fixed_location"] = stoi(attr_value);
                        } else if(attr_label.compare("radius") == 0) {
                            current_lat_point_def["radius"] = stoi(attr_value);
                        } else if(attr_label.compare("radius_x") == 0) {
                            current_lat_point_def["radius_x"] = stoi(attr_value);
                        } else if(attr_label.compare("radius_y") == 0) {
                            current_lat_point_def["radius_z"] = stoi(attr_value);
                        } else if(attr_label.compare("radius_z") == 0) {
                            current_lat_point_def["radius_z"] = stoi(attr_value);
                        } else if(attr_label.compare("x1") == 0){
                            current_lat_point_def["x1"] = stoi(attr_value);
                        } else if(attr_label.compare("y1") == 0){
                            current_lat_point_def["y1"] = stoi(attr_value);
                        } else if(attr_label.compare("z1") == 0){
                            current_lat_point_def["z1"] = stoi(attr_value);
                        } else if(attr_label.compare("x2") == 0){
                            current_lat_point_def["x2"] = stoi(attr_value);
                        } else if(attr_label.compare("y2") == 0){
                            current_lat_point_def["y2"] = stoi(attr_value);
                        } else if(attr_label.compare("z2") == 0){
                            current_lat_point_def["z2"] = stoi(attr_value);
                        }
                    }
                    new_compartment.compart_shape_definitions.push_back(current_lat_point_def);
                }

                // get, if any, compartment properties such as membrane
                if(annot_node->first_node("cell4d:compartmentProperties")) {
                    auto property_node = annot_node->first_node("cell4d:compartmentProperties");
                    for(rapidxml::xml_attribute<>* property_attr = property_node->first_attribute(); property_attr; property_attr = property_attr->next_attribute()) {
                        string attr_label = property_attr->name();
                        string attr_value = property_attr->value();
                        if(attr_label == "type") {
                            if(attr_value == "membrane") new_compartment.is_membrane = true;
                        }
                        if(new_compartment.is_membrane && attr_label == "axis") {
                            if(attr_value == "x") new_compartment.axis = attr_value;
                            else if(attr_value == "y") new_compartment.axis = attr_value;
                            else if(attr_value == "z") new_compartment.axis = attr_value;
                            // if undefined by user, try to figure out correct orientation using shortest span of the compartment dimensions
                            else {
                                map<string, int> compart_dim = new_compartment.compart_shape_definitions[0];
                                vector<int> dimensions;
                                int x_length = compart_dim["x2"] - compart_dim["x1"];
                                int y_length = compart_dim["y2"] - compart_dim["y1"];
                                int z_length = compart_dim["z2"] - compart_dim["z1"];
                                if(x_length < y_length && x_length < z_length) new_compartment.axis="x";
                                else if(y_length < x_length && y_length < z_length) new_compartment.axis="y";
                                else if(z_length < x_length && y_length < x_length) new_compartment.axis="z";
                                else {
                                    cerr << "axis was undefined, could not be determined. exiting." << endl;
                                    exit(1);
                                }
                                cout << "compartment " << compartment_id << " axis automatically determined to be " << new_compartment.axis << "-locked." << endl;
                            }
                        } else if(new_compartment.is_membrane && attr_label == "face") {
                            if(attr_value == "front") new_compartment.face = attr_value;
                            else if(attr_value == "back") new_compartment.face = attr_value;
                            else {
                                cerr << "compartment " << new_compartment.name << " face property invalid." << endl;
                                exit(1);
                            }
                        } else if(attr_label == "membraneEmissionRate") {
                            new_compartment.membrane_emission_rate = atof(attr_value.c_str());
                        } else if(attr_label == "absorptionRate") {
                            new_compartment.membrane_absorption_rate = atof(attr_value.c_str());
                        }

                    }
                }

            }
        } else {
        }
        compartments_map[compartment_id] = new_compartment;
    }
    // make sure the poslog compartment actually exists
    // if it does, add it to poslog_json, a map of compart:json poslogs
    if(poslog) {
        for(auto & compart : poslog_comparts) {
            if(compartments_map.find(compart) == compartments_map.end()) {
                cerr << "Experimental poslog feature enabled, but specified compartment" << compart << "does not exist. exiting." << endl;
                exit(1);
            }
            json new_compart_json;
            poslog_json[compart] = new_compart_json;
        }

    }
}

// go through listOfSpecies
void ParameterManager::set_species(Model_t * model, bool load_checkpoint) {
    Species_t* species_matrix[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        species_matrix[i] = Model_getSpecies(model, i);
        string species_id = Species_getId(species_matrix[i]);

        // walk through the listOfSpecies, grab the initial counts, and the compartments that they exist in.
        string species_annotation = SBase_getAnnotationString(species_matrix[i]);
        rapidxml::xml_document <> species_doc;
        species_doc.parse<0>(&species_annotation[0]);
        
        map<string, int> compart_initial_map;
        vector <Species_Attributes *> current_container_components;
        bool create_container_rule = false;
        // container_configs are the combination of states that a container is in, and the "rule" is the corresponding special stats
        // that containers of this config will have, such as diffusion rate, color, etc.
        Container_Config cur_config;
        cur_config.config_name = species_id;
        Container_Rules cur_container_rules;
        bool draw_particle = false;
        // iterate through each node in species annotation
        for(rapidxml::xml_node<>* lvl_0 = species_doc.first_node()->first_node(); lvl_0; lvl_0 = lvl_0->next_sibling()) {
            string lvl_0_node_name = lvl_0->name();
            // Grab the annotation fields within the species node by node
            if(lvl_0_node_name.compare("cell4d:listOfValidCompartments") == 0) {
                // listOfValidCompartments-> cell4d:compartment
                for(rapidxml::xml_node<>* compart_node = lvl_0->first_node(); compart_node; compart_node = compart_node->next_sibling()) {
                    string compart_name = compart_node->first_attribute("id")->value();
                    if(compartments_map.find(compart_name) == compartments_map.end()) {
                        cerr << "Invalid compartment specified for initial species placement. exiting" << endl;
                        exit(1);
                    }
                    int compart_amount = stoi(compart_node->first_attribute("initial")->value());
                    compart_initial_map[compart_name] = compart_amount;
                }
            } else if(lvl_0_node_name.compare("cell4d:listOfSpeciesTypes") == 0) { // fills out the "default" (initial) states of the species based on info from listOfSpecies
                // Species_Attributes species_attr;
                for(rapidxml::xml_node<>* species_node = lvl_0->first_node(); species_node; species_node = species_node->next_sibling()) {
                    Species_Attributes * current_species = new Species_Attributes;
                    string species_name = species_node->first_attribute("id")->value();
                    current_species->ID = species_name;
                    if(species_details_map[species_name]->is_simple_molecule()) {
                        current_species->species_type = "SIMPLE_MOLECULE";
                    }
                    // only look for mod and binding sites if it is a protein
                    if(species_details_map[species_name]->is_protein()) {
                        draw_particle = true;
                        current_species->species_type = "PROTEIN";
                        // speciesType -> bindingSite
                        for(rapidxml::xml_node<>* bind_node = species_node->first_node(); bind_node; bind_node = bind_node->next_sibling()) {
                            string site_name = bind_node->first_attribute("id")->value();
                            for(rapidxml::xml_attribute<>* bind_attr = bind_node->first_attribute(); bind_attr; bind_attr = bind_attr->next_attribute()) {
                                string attr_name = bind_attr->name();
                                if(attr_name.compare("state") == 0) {
                                    if(bind_attr->value()) {
                                        string site_state = bind_node->first_attribute("state")->value();
                                        current_species->default_states[site_name] = site_state;
                                    }
                                } else if(attr_name.compare("binding") == 0) {
                                    if(bind_attr->value()) {
                                        string site_state = bind_node->first_attribute("binding")->value();
                                        int bound;
                                        if(site_state.compare("unbound") == 0) { // unbound as 0
                                            bound = 0;
                                        } else if(site_state.compare("*") == 0) { // save wildcards as -1
                                            bound = -1;
                                        } else if(site_state.compare("bind") == 0 || site_state.compare("bound") == 0) { // unbound as 0
                                            bound = 1;
                                        } else {
                                            bound = stoi(site_state);
                                        }
                                        current_species->default_bind[site_name] = bound;
                                    }
                                }
                            }
                        }
                        // ensure that no wildcards are used when initializing containers in the simulation
                        if(current_species->default_bind.size() != species_details_map[species_name]->binding_sites.size()) {
                            cerr << "not all binding sites info filled in" << species_id << ". exiting." << endl;
                            exit(1);
                        } else if(current_species->default_states.size() != species_details_map[species_name]->possible_states.size()) {
                            cerr << "not all state sites info filled in" << species_id << ". exiting." << endl;
                            exit(1);
                        }
                    }
                    current_container_components.push_back(current_species);
                }
            
            }
            // here begins container_rule parsing, rules that override implicitly generated stats that apply to specific container configs
            // only applies to particles
            if(draw_particle) {
                if(lvl_0_node_name.compare("cell4d:diffusionConstant") == 0) {
                    create_container_rule = true;
                    if(!lvl_0->first_attribute("value")) {
                        cerr << "diffusion const value not found. exiting." << endl;
                        exit(1);
                    }
                    cur_container_rules.diffusion_constant = atof(lvl_0->first_attribute("value")->value());
                    cur_container_rules.fields_filled["diffusionConstant"] = true;

                } else if(lvl_0_node_name.compare("cell4d:compartmentEmissionRate") == 0) {
                    create_container_rule = true;
                    if(!lvl_0->first_attribute("origin") || !lvl_0->first_attribute("target") || !lvl_0->first_attribute("probability")) {
                        cerr << "compart emission rate node missing required fields. origin, target, probability are mandatory." << endl;
                        exit(1);
                    }
                    string origin_compartment = lvl_0->first_attribute("origin")->value();
                    string target_compartment = lvl_0->first_attribute("target")->value();
                    float probability = stof(lvl_0->first_attribute("probability")->value());
                    cur_container_rules.compartment_emission_rates[origin_compartment][target_compartment] = probability;
                    cur_container_rules.fields_filled["compartmentEmissionRate"] = true;

                } else if(lvl_0_node_name.compare("cell4d:displayProperties") == 0) {
                    create_container_rule = true;
                    // optional field for special color display
                    if(lvl_0->first_attribute("redValue") && lvl_0->first_attribute("greenValue") && lvl_0->first_attribute("blueValue")) {
                        cur_container_rules.container_color["red"] = stoi(lvl_0->first_attribute("redValue")->value());
                        cur_container_rules.container_color["green"] = stoi(lvl_0->first_attribute("greenValue")->value());
                        cur_container_rules.container_color["blue"] = stoi(lvl_0->first_attribute("blueValue")->value());
                        cur_container_rules.fields_filled["displayProperties"] = true;
                    } else if(lvl_0->first_attribute("redValue") || lvl_0->first_attribute("greenValue") || lvl_0->first_attribute("blueValue")) {
                        cerr << "color attributes red blue green must all be specified, if any are present. exiting." << endl;
                        exit(1);
                    }

                    // optional field for density gradient display on membranes rather than particles
                    if(lvl_0->first_attribute("membraneDisplay")) {
                        string display_value = lvl_0->first_attribute("membraneDisplay")->value();
                        if(display_value == "true" || display_value == "True" || display_value == "TRUE") {
                            cur_container_rules.membrane_display = true;
                            cur_container_rules.fields_filled["membraneDisplay"] = true;
                        } else if(display_value != "false" && display_value != "False" && display_value != "FALSE") {
                            cerr << "please specify either \"true\" or \"false\" in the optional membraneDisplay attribute." << endl;
                            exit(1);
                        }
                    }
                }
            }
        }
        if(draw_particle) {
            string container_name = "";
            vector<string> name_vec;
            for(auto & species : current_container_components) {
                name_vec.push_back(species->ID);
                // if(create_container_rule) {
                    Species_Config cur_species_config;
                    cur_species_config.molecule_name = species->ID;
                    cur_species_config.molecule_state_map = species->default_states;
                    cur_species_config.molecule_binding_map = species->default_bind;
                    cur_config.species_multiset.insert(cur_species_config);
                // }
            }
            sort(name_vec.begin(),name_vec.end());
            for(string & name_itr : name_vec) {
                container_name += name_itr;
                container_name += ";";
            }
            if(!container_name.empty()) {
                container_name.pop_back();
            }
            // if(create_container_rule) {
                cur_config.container_name = container_name;
                container_rules_map[cur_config] = cur_container_rules;
            // }
        }

        if(!load_checkpoint) {
            initial_species_map.push_back(make_pair(compart_initial_map, current_container_components));
        }
        // temp workaround for loading in bulk even if using checkpoint
        if(load_checkpoint && current_container_components.size() == 1 && current_container_components[0]->is_simple_molecule()) {
            initial_species_map.push_back(make_pair(compart_initial_map, current_container_components));
        }
        list_of_species_map[species_id] = current_container_components;
    }

}


void ParameterManager::set_reactions(Model_t * model) {
    int numRxns = Model_getNumReactions(model);
    rxns = new Reaction_t* [numRxns];

    // pull in the reaction info
    for(int i=0; i<numRxns; i++) {
        rxns[i] = Model_getReaction(model, i);
        string reaction_key = Reaction_getId(rxns[i]);
        bool reversible_reaction_flag = Reaction_getReversible(rxns[i]);
        ListOf* list_of_modifiers = Reaction_getListOfModifiers(rxns[i]);
        ListOf* list_of_reactions = Reaction_getListOfReactants(rxns[i]);
        ListOf* list_of_products = Reaction_getListOfProducts(rxns[i]);
        KineticLaw* kinetic_law = Reaction_getKineticLaw(rxns[i]);

        string reaction_annot = SBase_getAnnotationString(rxns[i]); // annotation of reaction
        string modifiers_annotation = SBase_getAnnotationString(list_of_modifiers);
        string reactants_annotation = SBase_getAnnotationString(list_of_reactions);
        string products_annotation = SBase_getAnnotationString(list_of_products);

        
        rapidxml::xml_document<> annotation_doc;
        rapidxml::xml_document<> modifiers_doc;
        rapidxml::xml_document<> reactants_doc;
        rapidxml::xml_document<> products_doc;
        
        annotation_doc.parse<0>(&reaction_annot[0]);
        modifiers_doc.parse<0>(&modifiers_annotation[0]);
        reactants_doc.parse<0>(&reactants_annotation[0]); // will only take in strings like this
        products_doc.parse<0>(&products_annotation[0]);
        
        vector<string> compart_list;
        vector<string> modifier_name_list;
        map<int, vector<string>> products_name_list;
        map<int, vector<string>> reactants_name_list;

        map<int, string> product_destinations;

        map<int, map<string,int>> modifier_binding_map;
        map<int, map<string,string>> modifier_state_map;
        map<int, map<int, map<string,int>>> reactant_binding_map; // species name as key, map of site name: state name
        map<int, map<int, map<string,string>>> reactant_state_map; // species name as key, map of site name: state name
        map<int, map<int, map<string,int>>> product_binding_map; // species name as key, map of site name: state name
        map<int, map<int, map<string,string>>> product_state_map; // species name as key, map of site name: state name
        
        // annotation -> listOfCompartments
        rapidxml::xml_node<>* compartment_list_node = annotation_doc.first_node()->first_node("cell4d:listOfCompartments");
        // if the list of compartments for reaction are defined, add all of list to allowed reaction vector
        if(compartment_list_node) {
            // listOfCompartments -> compartment
            for(rapidxml::xml_node<>* compartment_node = compartment_list_node->first_node("cell4d:compartment"); compartment_node; compartment_node = compartment_node->next_sibling()) {
                string compartment_name = compartment_node->first_attribute("id")->value();
                compart_list.push_back(compartment_name);
            }
        } else { // if undefined, add all compartments to the vector
            for(auto const & compart : compartments_map) {
                compart_list.push_back(compart.first);
            }
        }
        bool spectator_defined = false;
        int spectator_placement = 0;
        rapidxml::xml_node<>* spectator_node = annotation_doc.first_node()->first_node("cell4d:spectatorPlacement");
        if(spectator_node) {
            if(spectator_node->first_attribute("index")) {
                spectator_defined = true;
                spectator_placement = stoi(spectator_node->first_attribute("index")->value());
            }
        }
        cout << "Reaction: " << reaction_key << endl;
        // extracting all of the modifiers from the reaction, if any.
        // Todo: (find the code to detect the species reference flag)
        for(rapidxml::xml_node<>* speciesRef_node = modifiers_doc.first_node()->first_node(); speciesRef_node; speciesRef_node = speciesRef_node->next_sibling()){
			int species_index = 0;
            // speciesReference -> listofspeciestypes -> speciestype
            for(rapidxml::xml_node<>* speciesType = speciesRef_node->first_node()->first_node(); speciesType; speciesType = speciesType->next_sibling()){
				string modifier_name = speciesType->first_attribute("id")->value();
				if(modifier_name.compare("empty") == 0){
					break;
                }
                modifier_name_list.push_back(modifier_name);
                species_details_map[modifier_name]->species_type = "ENZYME"; 
                rapidxml::xml_node<>* bindingSite_node;
                for(bindingSite_node = speciesType->first_node(); bindingSite_node; bindingSite_node = bindingSite_node->next_sibling()) { // iterate through all binding sites
                    // inserts to a map that stores each reaction a reactant is related to
                    string site_name = bindingSite_node->first_attribute("id")->value();
                    if(site_name.compare("None") == 0) break; // move on if no binding sites
                    rapidxml::xml_attribute<>* site = bindingSite_node->first_attribute("id")->next_attribute();
                    string attr_name = site->name();
                    if(attr_name.compare("binding") == 0) { // translate "bound" and "unbound" from reaction rules in xml into int encoding
                        string attr_string = site->value();
                        int bound;
                        if(attr_string.compare("unbound") == 0) { // unbound as 0
                            bound = 0;
                        } else if(attr_string.compare("*") == 0) { // save wildcards as -1
                            bound = -1;
                        } else if(attr_string.compare("bind") == 0 || attr_string.compare("bound") == 0) { // unbound as 0
                            bound = 1;
                        } else {
                            bound = stoi(attr_string);
                        }
                        modifier_binding_map[species_index][site_name] = bound;
                    } else if(attr_name.compare("state") == 0) {
                        modifier_state_map[species_index][site_name] = site->value();
                    }
                }
			}
		}
        
        bool has_solo_bulk_reactant = false;
        bool has_solo_bulk_product = false;
        int bulk_index = -1;

        // extracting all of the reactants from the reaction.
        int reactant_counter = 0;
        // annot -> speciesReference
        for(rapidxml::xml_node<>* speciesRef_node = reactants_doc.first_node()->first_node(); speciesRef_node; speciesRef_node = speciesRef_node->next_sibling()) { // this will go max. 2 times
            // speciesReference -> listofspeciestypes -> speciestype
            int species_index = 0;
            bool has_bulk_reactant = false;
            for(rapidxml::xml_node<>* speciesType_node = speciesRef_node->first_node()->first_node(); speciesType_node; speciesType_node = speciesType_node->next_sibling()) {
                //finds all the species in one reacting container
                string reactant_name = speciesType_node->first_attribute("id")->value();
                cout << reaction_key << " reactant num " << reactant_counter << " - " << speciesType_node->first_attribute("id")->value() << endl;
                reactants_name_list[reactant_counter].push_back(reactant_name);
                // check if any reactants are actually bulk molecules interacting with particles
                if(find(small_molecule_id_list.begin(), small_molecule_id_list.end(), reactant_name) != small_molecule_id_list.end()) {
                    has_bulk_reactant = true;
                }
                // speciestype -> bindingsite
                rapidxml::xml_node<>* bindingSite_node;
                for(bindingSite_node = speciesType_node->first_node(); bindingSite_node; bindingSite_node = bindingSite_node->next_sibling()) { // iterate through all binding sites
                    // inserts to a map that stores each reaction a reactant is related to
                    string site_name = bindingSite_node->first_attribute("id")->value();
                    if(site_name.compare("None") == 0) break; // move on if no binding sites
                    rapidxml::xml_attribute<>* site = bindingSite_node->first_attribute("id")->next_attribute();
                    string attr_name = site->name();
                    if(attr_name.compare("binding") == 0) { // translate "bound" and "unbound" from reaction rules in xml into int encoding
                        string attr_string = site->value();
                        int bound;
                        if(attr_string.compare("unbound") == 0) { // unbound as 0
                            bound = 0;
                        } else if(attr_string.compare("*") == 0) { // save wildcards as -1
                            bound = -1;
                        } else if(attr_string.compare("bind") == 0 || attr_string.compare("bound") == 0) { // unbound as 0
                            bound = 1;
                        } else {
                            bound = stoi(attr_string);
                        }
                        reactant_binding_map[reactant_counter][species_index][site_name] = bound;
                    } else if(attr_name.compare("state") == 0) {
                        reactant_state_map[reactant_counter][species_index][site_name] = site->value();
                    }
                }
                species_index++;
            }
            // if a reactant is by itself and is small molecule, will be modelled as bulk in simulation. particle needs to react with it specially.
            if(species_index == 1 && has_bulk_reactant) {
                has_solo_bulk_reactant = true;
                bulk_index = reactant_counter;
            }
            reactant_counter++;
        }

        int product_counter = 0;
        // extracting all of the products of the reaction  
        // annot -> speciesReference
        for(rapidxml::xml_node<>* speciesRef_node = products_doc.first_node()->first_node(); speciesRef_node; speciesRef_node = speciesRef_node->next_sibling()) { 
            int species_index = 0;
            bool has_bulk_product = false;
            rapidxml::xml_node<>* act_transport_node = speciesRef_node->first_node("cell4d:activeTransport");
            if(act_transport_node) {
                string transport_dest_compart = act_transport_node->first_attribute("destination")->value();
                product_destinations[product_counter] = transport_dest_compart;
            }
            // speciesReference -> listofspeciestypes -> speciestype
            for(rapidxml::xml_node<>* speciesType_node = speciesRef_node->first_node("cell4d:listOfSpeciesTypes")->first_node(); speciesType_node; speciesType_node = speciesType_node->next_sibling()) {
                string product_name = speciesType_node->first_attribute("id")->value();
                cout << reaction_key << " product num " << product_counter << " - " << product_name << endl;
                products_name_list[product_counter].push_back(product_name);
                // check if any products are actually bulk molecules released from complex
                if(find(small_molecule_id_list.begin(), small_molecule_id_list.end(), product_name) != small_molecule_id_list.end()) {
                    has_bulk_product = true;
                }
                // speciestype -> bindingsite
                rapidxml::xml_node<>* bindingSite_node;
                for(bindingSite_node = speciesType_node->first_node(); bindingSite_node; bindingSite_node = bindingSite_node->next_sibling()) { // iterate through all binding sites
                    // inserts to a map that stores each reaction a product is related to
                    string site_name = bindingSite_node->first_attribute("id")->value();
                    if(site_name.compare("None") == 0) break; // move on if no binding sites
                    rapidxml::xml_attribute<>* site = bindingSite_node->first_attribute("id")->next_attribute();
                    string attr_name = site->name();
                    if(attr_name.compare("binding") == 0) { // translate "bound" and "unbound" from reaction rules in xml into int encoding
                        string attr_string = site->value();
                        int bound;
                        if(attr_string.compare("unbound") == 0) { // unbound as 0
                            bound = 0;
                        } else if(attr_string.compare("*") == 0) { // save wildcards as -1
                            bound = -1;
                        } else if(attr_string.compare("bind") == 0 || attr_string.compare("bound") == 0) { // unbound as 0
                            bound = 1;
                        } else {
                            bound = stoi(attr_string);
                        }
                        product_binding_map[product_counter][species_index][site_name] = bound;
                    } else if(attr_name.compare("state") == 0) {
                        product_state_map[product_counter][species_index][site_name] = site->value();
                    }
                }
                species_index++;
            }
            // if a product is by itself and is small molecule, will be modelled as bulk in simulation.
            if(species_index == 1 && has_bulk_product) {
                has_solo_bulk_product = true;
                bulk_index = product_counter;
            }
            product_counter++;
        }

        if(has_solo_bulk_product && has_solo_bulk_reactant) {
            cerr << "reaction " << reaction_key << "has a bulk reactant and a bulk product, this isn't allowed. exiting." << endl;
        }

        // Cell_Reaction current_reaction;
        Cell_Reaction * current_reaction = new Cell_Reaction;
        current_reaction->reaction_name = reaction_key;
        current_reaction->is_reversible = reversible_reaction_flag;
        current_reaction->allowed_compartments = compart_list;

        if(spectator_defined) {
            current_reaction->spectator_index = spectator_placement;
        } 

        current_reaction->reactants_binding_map = reactant_binding_map;
        current_reaction->reactants_state_map = reactant_state_map;
        current_reaction->products_binding_map = product_binding_map;
        current_reaction->products_state_map = product_state_map;
        current_reaction->reactants_list = reactants_name_list;
        current_reaction->products_list = products_name_list;

        current_reaction->product_container_destinations = product_destinations;

        if(modifier_name_list.empty()) { // labels the reaction type based on the reactive species
            bool param_found = false;
            for(int i = 0; i < KineticLaw_getNumParameters(kinetic_law); i++){
                // this picks apart the kinetic law portion and feeds it into the reaction class
                Parameter* k_law_field = KineticLaw_getParameter(kinetic_law, i);
                double k_law_field_value = k_law_field->getValue();
                const string& k_law_field_name = k_law_field->getName(); // to accommodate HPF behavior, trim string
                string law_name_substring = k_law_field_name.substr(0,8);
                cout << "kinetic name: " << law_name_substring << endl;
                cout << "kinetic value: " << k_law_field_value << endl;
                if(law_name_substring.compare("Kforward") == 0) {
                    if(param_found) {
                        cerr << "Reaction " << reaction_key << " already has a set parameter conflicting with forward reaction rate." << endl;
                        exit(1);
                    }
                    current_reaction->forward_reaction_rate = k_law_field_value;
                    if(current_reaction->forward_reaction_rate < 1) {
                        cerr << "Warning: Forward rate of " << reaction_key << " is less than 1. Was this supposed to be radius?" << endl;
                    }
                    param_found = true;
                } else if(law_name_substring.compare("Kreverse") == 0) {
                    current_reaction->reverse_reaction_rate = k_law_field_value;
                } else {
                    string new_substring = law_name_substring.substr(0,6);
                    if(new_substring == "radius") {
                        if(reactants_name_list.size() != 2) {
                            cerr << "Reaction " << reaction_key << " is not bimolecular, should not have a radius parameter." << endl;
                            exit(1);
                        }
                        if(param_found) {
                            cerr << "Reaction " << reaction_key << " already has a set parameter conflicting with forward reaction radius." << endl;
                            exit(1);
                        }
                        current_reaction->binding_radius = k_law_field_value;
                        if(current_reaction->binding_radius > 1) {
                            cerr << "Warning: Binding radius of " << reaction_key << " is greater than 1. Was this supposed to be a rate constant?" << endl;
                        }
                        param_found = true;
                    } else if(new_substring.compare("unbind") == 0) {
                        current_reaction->set_unbinding_radius(k_law_field_value / parameterTable[SPACESCALE]);
                        cout << "Reaction " << reaction_key << " unbinding radius set to " << k_law_field_value << endl;
                        if(k_law_field_value > 1) {
                            cerr << "Warning: minimum unbinding radius parameter is in units meter. Currently set to " << k_law_field_value << " meters." << endl;
                        }
                    } else {
                        cerr << "Kinetic parameter in reaction " << reaction_key << " not recognized, exiting." << endl;
                        exit(1);
                    }
                }
            }
            if(!param_found) {
                cerr << "kineticparameters not found for interaction reaction in XML, name should be Kforward and optionally radius, Kreverse. exiting" << endl;
                exit(EXIT_FAILURE);
            }
            current_reaction->set_reactant_product_pairing();
            // decide which reaction map to put this reaction in
            if(!has_solo_bulk_reactant) {
                if(reactant_counter > 1) {
                    bimolecular_reaction_details_map[reaction_key] = current_reaction;
                } else {
                    unimolecular_reaction_details_map[reaction_key] = current_reaction;
                }
            } else {
                current_reaction->bulk_index = bulk_index;
                // because Cell_Reaction object isn't a pointer, assign to map LAST, after all changes
                bulk_reaction_details_map[reaction_key] = current_reaction;
            }
            

            if(current_reaction->is_reversible) {
                if(current_reaction->reverse_reaction_rate == 0) {
                    cerr << "Warning: reaction " << reaction_key << " reverse rate is 0." << endl;
                } 
                // Cell_Reaction reverse_reaction = current_reaction;
                Cell_Reaction * reverse_reaction = new Cell_Reaction;
                *reverse_reaction = *current_reaction;

                // set up pointers to the paired reverse reaction
                current_reaction->reverse_reaction = reverse_reaction;
                reverse_reaction->reverse_reaction = current_reaction;

                reverse_reaction->reaction_name = reaction_key + "_rev";
                reverse_reaction->forward_reaction_rate = current_reaction->reverse_reaction_rate;
                reverse_reaction->reverse_reaction_rate = 0;
                reverse_reaction->reactants_binding_map = product_binding_map;
                reverse_reaction->reactants_state_map = product_state_map;
                reverse_reaction->products_binding_map = reactant_binding_map;
                reverse_reaction->products_state_map = reactant_state_map;
                reverse_reaction->reactants_list = products_name_list;
                reverse_reaction->products_list = reactants_name_list;

                reverse_reaction->set_reactant_product_pairing();
                // decide which reaction map to put this reaction in
                if(!has_solo_bulk_product)
                    if(product_counter > 1) {
                        bimolecular_reaction_details_map[reverse_reaction->reaction_name] = reverse_reaction;
                    } else {
                        unimolecular_reaction_details_map[reverse_reaction->reaction_name] = reverse_reaction;
                    }
                else {
                    // because Cell_Reaction object isn't a pointer, assign to map LAST, after all changes
                    reverse_reaction->bulk_index = bulk_index;
                    bulk_reaction_details_map[reverse_reaction->reaction_name] = reverse_reaction;
                }

                if(reverse_reaction->reactants_list.size() > 2) {
                    cerr << "Reactions with more than two products cannot be reversible, since the reverse reaction will have more than two reactants." << endl;
                    exit(1);
                }
            }

        } else if(!modifier_name_list.empty()) { // enzymatic reactions using kinetic values
            current_reaction->modifiers_list = modifier_name_list;
            current_reaction->modifiers_binding_map = modifier_binding_map;
            current_reaction->modifiers_state_map = modifier_state_map;
            for(int i = 0; i < KineticLaw_getNumParameters(kinetic_law); i++){
                // this picks apart the kinetic law portion and feeds it into the reaction class
                Parameter* k_law_field = KineticLaw_getParameter(kinetic_law, i);
                string k_law_field_name = k_law_field->getName();
                double k_law_field_value = k_law_field->getValue();
                //cout << k_law_field_name << ": " << k_law_field_value << endl;

                if(k_law_field_name.compare("Ka") == 0) {
                    current_reaction->forward_reaction_rate = k_law_field_value;
                } else if(k_law_field_name.compare("Kp") == 0) {
                    current_reaction->reverse_reaction_rate = k_law_field_value;
                } else if(k_law_field_name.compare("Keq") == 0) {
                    current_reaction->equilibrium_constant = k_law_field_value;
                } else if(k_law_field_name.compare("Vm") == 0) {
                    current_reaction->maximum_rate = k_law_field_value;
                }
            }
            enzymatic_reaction_details_map[reaction_key] = current_reaction;
        }
    }
}


//////////////////////////////////////////////////////////////////////
// Load parameters from data file 'fName'
//
void ParameterManager::load(string fName){
    string key;
    string value;
    string line;
    string delimiter = "=";

    // Open Parameter data file
    ifstream inFile (fName.c_str());

    if(inFile){
        // Populate the parameter table with data in the file
        while (!inFile.eof())
        {
            getline(inFile,line,'\n');
            int keyPos = line.find(delimiter);
            if(keyPos > 0){
                key = line.substr(0,keyPos);
                value = line.substr(keyPos + delimiter.length());
                // Add parameter to the paramater table
                parameterTable[key] = atoi(value.c_str());
            }
        }

        // Close Parameter data file
        inFile.close();
    }
}

int ParameterManager::extract_int(string key, string annotation){
    //fromat of tag is <cell4d:varname>, +1 to go past the >
    try {
        size_t pos = annotation.find(key) + key.length() + 1;
        string  remain = annotation.substr(pos); 
        //format of end tag is </cell4d:....> just need to go before the <
        size_t lastpos = remain.find(key) - 9;
        //get value 
        string temp = remain.substr(0,lastpos);
        return  atoi(temp.c_str());
    } catch(exception& e) {
        return 0;
    }
}

double ParameterManager::extract_double(string key, string annotation){
    //fromat of tag is <cell4d:varname>, +1 to go past the >
    try {
        size_t pos = annotation.find(key) + key.length() + 1;
        string  remain = annotation.substr(pos); 
        //format of end tag is </cell4d:....> just need to go before the <
        size_t lastpos = remain.find(key) - 9;
        //get value 
        string temp = remain.substr(0,lastpos);
        return  atof(temp.c_str());
    } catch(exception& e) {
        return 0;
    }
}


//////////////////////////////////////////////////////////////////////
// Get value of the Parameter by the specified KEY
//
int ParameterManager::get_int(string key){

    // Iterate through the parameter table
    map<string, double >::const_iterator itr = parameterTable.find(key);

    // If the KEY is found in the parameter table
    if(itr != parameterTable.end()){
        int returned_value = int(round(itr->second));
        return (returned_value);
    }

    // Return NIL if KEY not found in the parameter table
    return 0;
}

double ParameterManager::get_double(string key){

    // Iterate through the parameter table
    map<string, double >::const_iterator itr = parameterTable.find(key);

    // If the KEY is found in the parameter table
    if(itr != parameterTable.end()){
        return (itr->second);
    }

    // Return NIL if KEY not found in the parameter table
    return 0;
}


//////////////////////////////////////////////////////////////////////
// Add parameter (key, value) pair
//
void ParameterManager::add(string key, double value){
    parameterTable[key] = value;
}
