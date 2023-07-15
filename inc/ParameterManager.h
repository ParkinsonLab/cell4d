//========================================================
// ParameterManager.h
// Donny Chan, Billy Taj

#ifndef PMANAGER_H
#define PMANAGER_H

class Point;

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <map>
#include <fstream>
#include <iostream>

#include <sbml/common/extern.h>
#include <sbml/SBase.h>
#include <sbml/common/sbmlfwd.h>
#include <sbml/SBMLTypes.h>
#include "../ext/rapidxml/rapidxml.hpp"
#include "../ext/rapidxml/rapidxml_utils.hpp"
#include "../ext/rapidxml/rapidxml_print.hpp"
#include "definition.h"
#include <sbml/ModifierSpeciesReference.h>
#include "Species_Attributes.h"

// for json parsing
#include "../ext/json.hpp"
// for convenience
using json = nlohmann::json;

using namespace std;

struct Event_Parameters {
    string event_name; // name of event from xml model
    string event_type; // type of event: add_mols, delete_mols, transport_mols
    string event_trigger; // trigger of event: time, other events, simulation state
    double event_probability; // probability of event going when triggered

    bool event_repeat = true;
    string event_triggering_event;
    int event_triggering_condition = -1;
    vector<Event_Parameters *> events_to_trigger;

    string state_trigger_molecule;
    int state_trigger_mol_amount;
    vector<Point> state_trigger_voxels;
    vector<string> state_trigger_compartment;

    int event_interval = 0;
    int event_initial_time = 0;
    int event_final_time = -1;

    int event_container_amount_total;
    string event_compartment = "";

    string event_container_species_name;

    vector<Point> event_placement_voxels;
    vector<pair<Point, int>> event_molecule_placement_pairs;

    vector<Point> transport_destination_voxels;
    vector<string> transport_destination_compartments;

    bool event_done = false; // check off non-repeating events
};

class ParameterManager{

    public:

        // Get an instance of the Parameter Manager
        static ParameterManager * getParameterManager();

        // Get value of the Parameter by the specified KEY
        int get_int(string);

        // Get value of the Parameter by the specified KEY
        double get_double(string);

        // Add parameter (key, value) pair
        void add(string, double);

        // load parameters from data file
        void load(string);

        // load SBML data file 
        void loadSBML(const char*, bool load_checkpoint = false, string checkpoint_filename = NULL);

        // set global simulation parameters
        void set_global_parameters(string & annot);
        void set_annot_species(const rapidxml::xml_document<> & annotation_xml_doc);
        void set_compartments(Model_t * model);
        void set_species(Model_t * model, bool load_checkpoint = false);
        void set_reactions(Model_t * model);
        void set_events(const rapidxml::xml_document<> & annotation_xml_doc);
        
        void set_event_compartment_voxels();

        //save model
        void ParameterManagerModelSave(SBMLDocument * doc, const char *  filename );

        vector <pair<map<string, int>, vector <Species_Attributes *>>> initial_species_map; // pair<<compartment:amount>, container_contents>
        map<string, vector <Species_Attributes *>> list_of_species_map; // map of species id from xml to the vector of species attributes

        map<string, Event_Parameters *> event_map;
        map<int, vector<Event_Parameters *>> event_queue;
        
    protected:

        // Constructor for the Parameter Manager
        ParameterManager();

    private:

        string formatWriteInt(int var, string varName);
        string formatWriteDouble(double var, string varName);

        int extract_int(string key, string annotation);

        double extract_double(string key, string annotation);
        // Does the Parameter Manager already exist (Status)
        static bool isInstantiated;

        // Pointer to the Parameter Manager
        static ParameterManager * params;

        // Parameter Table    (key <-> value)
        map <string, double> parameterTable;

        bool loaded;
}
;

#endif
