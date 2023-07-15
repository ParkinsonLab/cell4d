//==========================================================
// Molecule_Container_Reactions.h
// Donny Chan, Billy Taj

#ifndef MOL_CONTAINER_REACTIONS
#define MOL_CONTAINER_REACTIONS

#include "Molecule_Container.h"
#include <iostream>

/*
    This class contains code that help make a reaction happen. 
    It will dictate what sort of manipulations that will occur when 2 or more Molecule Containers come into contact with eachother.
*/

class Molecule_Container_Reactions{
    public:
        Molecule_Container_Reactions();
        ~Molecule_Container_Reactions();

        static vector<Molecule_Container *> create_new_empty_containers(int n);

        static Molecule_Container * fusion(Molecule_Container * head, Molecule_Container * tail);
        // static vector<Molecule_Container *> fission(map<Molecule *, pair<int, int>>, map<int, vector <string>>, int spectator_placement);
        static vector<Molecule_Container *> fission(map<Molecule *, pair<int, int>>, Cell_Reaction * reaction_details);

        static Molecule_Container * bulk_fusion(Molecule_Container * head, Molecule * bulk);

        static void molecule_state_change   (Molecule *, map <string, string>);
        static void molecule_bind_change    (Molecule *, map <string, int>);

        static void bind_molecules      (Molecule * mol_0, Molecule * mol_1, string site_0, string site_1);
        static void unbind_molecules    (Molecule * mol_0, Molecule * mol_1, string site_0, string site_1);

        static void add_molecule_to_container      (Molecule_Container *, Molecule *);
        static void remove_molecule_from_container (Molecule_Container *, Molecule *);

        static bool molecule_fills_rule(Molecule * mol, const map<string, string> & state_rules, const map<string, int> & bind_rules);

        static void createMoleculeContainers(int num, string molecule_name, Point *loc);

    private:

        // Environment Manager
        LatticeEnvironment * lEnv;

}
;

#endif
