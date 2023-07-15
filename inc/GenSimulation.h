//===========================================================
// GenSimulation.h
// Donny Chan, Billy Taj

#ifndef GENSimulation_H
#define GENSimulation_H

#include "Simulation.h"

#include <string>
#include <memory>

using namespace std;

class GenSimulation : public Simulation
{
    public:

        static GenSimulation * getSimulation();

        bool isSimulationCompleted();

    protected:

        GenSimulation();
        virtual ~GenSimulation();

    private:

        //////////////////////////////////////////////////////////////////////
        // Member variables

        //////////////////////////////////////////////////////////////////////
        // Private Operations

        // Input Data
        void initialDataInput();

        // Create new molecule of specified name
        Molecule_Container * createNewMoleculeContainer(string molecule_name, Point loc);
        Molecule_Container * createNewMoleculeContainer(vector<Species_Attributes *>, Point loc);

        // for recreating json molecules with pre-determined states and position
        Molecule_Container * createNewMoleculeContainer(vector<Molecule *>, Point loc, unsigned long old_container_id);

        // create molecules from checkpoint file
        void create_checkpoint_molecules();

        // create small molecules
        void create_bulk_molecules();

        // load bulk molecules from checkpoint
        void load_bulk_molecules();

        // Periodically update Behaviour/Properties of the Cell TNF_Simulation
        void updateSimulation();

        // Display data from simulation
        void displayResults();
}
;

#endif
