//======================================================
// Compartment.h
// Donny Chan, Billy Taj

#ifndef COMPART_H
#define COMPART_H
#include <string>
#include <vector>
#include <map>
using namespace std;
// June 14, 2018: This class has been repurposed to be a datastructure class.  Each class instance will store details of one compartment
// in the simulation

class Point;
class Compart {
	public:
		Compart();
		~Compart();
		string id; // Donny June 25 2018 - Added id, which should be used as the env identifier
		string name;
		string outside;

		// shape only defines how the given coordinates define how the compartment segment is drawn, and nothing else
		string shape;
		// flag for whether or not it is a membrane compartment
		bool is_membrane = false;
		// axis declares what plane the membrane is along, x/y/z
		string axis = "n";
		// face describes which side of the membrane particles will be close to, front means its close to the 0 side of that axis, back means the opposite side
		string face = "n";
		// does this compartment move or not, second variable is the movement speed of this compartment
		std::pair<bool, double> is_mobile = make_pair(false, 0);
		double membrane_emission_rate = 0;
		double membrane_absorption_rate = 0;
		// vector of all shapes to draw that defines all coordinates of this compartment, includes x/y/z1,2 for rectangles, center point and 1-3 radii for spheroids 
		vector<map<string, int>> compart_shape_definitions;

		int x1 = -1;
		int x2 = -1; 
		int y1 = -1; 
		int y2 = -1; 
		int z1 = -1; 
		int z2 = -1;
		int radius = -1;

		// list of points that represent all lattice sites that have this compartment label
		vector <Point *> list_of_voxels;
		
		// keep track of how many small molecules are currently in this compartment
		map<string, double> current_bulk_count;
		
};

#endif

/*
		std::string getName();
		std::string setName(std::string);
		std::string getType();
		std::string setType(std::string);
		int getX1();
		int setX1(int X1);
		int getX2();
		int setX2(int X2);
		int getY1();
		int setY1(int Y1);
		int getY2();
		int setY2(int Y2);
		int getZ1();
		int setZ1(int Z1);
		int getZ2();
		int setZ2(int Z2);
		int getFixedLocation();
		int setFixedLocation(int fixedLocation);
		int getEnv();
		int setEnv(int env);
*/
