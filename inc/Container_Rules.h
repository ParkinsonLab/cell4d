//=================================================
// Container_Rules.h
// Donny Chan, Billy Taj
//
// This class is used for the application of rules or specifications that are container-combination specific, rather than species-specific.
// That is, when certain rules such as color display or diffusion rates are described in the listOfSpecies attributes of the input
// file, rather than in the annotationspeciestypes. In these cases, the XML specified rules will be applied to the matching containers,
// overriding the inferred stats from the combination of species.
//
// Example: species A is displayed as red, so AA containers will also show up as red. the combination AA can be described in listOfSpecies,
// with an alternate color property. Now all containers that match AA will show up as the alternate color instead.
//
// Wildcard rules do not apply to container rules. Container rules are only applicable to containers that completely match the rule
// configuration, that is, the molecule components and all their states match.



#ifndef CONTAINER_RULES
#define CONTAINER_RULES

#include <map>

using namespace std;

class Container_Rules{
    public:
        // string container_name = "";
        Container_Rules();
        ~Container_Rules();

        double diffusion_constant;
        map<string, map<string, float>> compartment_emission_rates; 

        // red, blue, green keys for float RGB values of this container
        map<string, float> container_color;
        bool membrane_display = false;

        // map checklist of which rules are modified within this config container rule
        // valid keys right now are "diffusionConstant", "compartmentEmissionRate", "displayProperties", and "membraneDisplay"
        map<string, bool> fields_filled;
};

#endif