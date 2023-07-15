//======================================================
// SimEngine.h
// Donny Chan, Billy Taj

#ifndef SIM_ENGINE_H
#define SIM_ENGINE_H

// include OpenGL
#ifdef MAC_OS
    #include <OpenGL/glu.h>
    #include <OpenGL/gl.h>
    #include <GLUT/glut.h>
#else
#if 1
    #include <GL/glu.h>
    #include <GL/gl.h>
    #include <GL/glut.h>
#endif
#endif

#include "Molecule.h"
#include "Matrix.h"
#include "GenSimulation.h"
#include "functionLib.h"
#include "ParameterManager.h"
#include "LatticeEnvironment.h"
#include <string.h>
#include <sstream>
#include <algorithm>
#include <map>

extern int xDim, yDim, zDim;

class SimEngine{

    public:
        SimEngine();
        ~SimEngine();
        bool simulate();
        string statsToString();
        bool continueSimulation;
        bool showMetaConc;
        bool showMetaEnz;
        bool showEnv;
        int actionType;
        bool showEnvCyt;
        bool showEnv1;
        bool showEnv2;
        bool showEnv3;
        bool showEnv4;
        Simulation * cell;
        map<string, bool> disp_compart_toggle;
        map<string, bool> disp_metabolite_toggle;
        map<string, bool> disp_protein_toggle;
        map<string, vector<float>> graphic_display_colours;
        int object_display_count = 1;


#if 1
        void display();
        void displayStats();
        void displayGraph();
        
#endif

       private:

        void displayTotalConcentrations();

        float *r, *g, *b;

//        Simulation * cell;
        ParameterManager * pm;
        LatticeEnvironment * lEnv;
        // stores the list of species that will light up on the membrane
        vector<string> membrane_species_vector;
        // using container name as key, store the multiplier of brightness to display for that container
        // for example if there are 2+ species in that container 
        map<string, int> mem_container_bright_map;

#if 1
        void setDisplayEnvironment();
        void get_and_set_color_display(Molecule_Container *);
        vector<float> get_default_container_color(string container_name);
        void draw_Molecule_Container(Molecule_Container *);
        void displayMetaboliteConcentration();
        void displayEnvironment();
        void displayNumCycles();
        void draw_particles();
        void draw_sim_box();
        void master_highlight();
        void constructLattice();
        void display_membrane_gradient();
        void highlight_compartment();
        void highlight_metabolites();
        vector<float> set_display_colour(int object_count);
        // Background color, array in order or R, G, B values
        // float background_color [3] = {0.85, 0.85, 0.85}; 
        float background_color [3] = {0, 0, 0}; 

        //Colors and Materials
        GLfloat* COL_MOL_ENZYME_A;
        GLfloat* COL_MOL_ENZYME_B;
        GLfloat* COL_MOL_ENZYME_C;
        GLfloat* COL_MOL_ENZYME_D;
        GLfloat* COL_MOL_ENZYME_E;

        GLfloat* LATTICE;

        GLUquadricObj *qobj;
#endif
};

#endif
