
//===============================================================================
// SimEngine.cpp
// Donny Chan
//
// This code controls much of the visual aspects of Cell4D.  
// In particular, it is responsible for the glow effects (highlight)
// and motions of the window.  The main class is the SimEngine class, 
// instantiated with a poorly-named object called "simulation" in SimWindow.


#include "../inc/SimEngine.h"
#include "../inc/Compartment.h"
#include <ctime>
#include <cstdlib>
#include <queue>

//extern map<int, Species_t *> metabolite_map;
//extern map<int, Species_t *> final_metabolite_map;
// 2018 additions
extern map<string, Species_Attributes *> species_details_map;
// extern map<string, int> new_metabolite_index_map;
extern map<string, Compart> compartments_map;   
extern int numMetabolite;
extern vector<string> small_molecule_id_list;
extern map<int, string> small_molecule_int_to_id_map;
extern map<string, vector<string>> container_components_map;
extern map<string, vector<float>> container_color_map; // stores the colors of each container type
extern map<Container_Config, Container_Rules> container_rules_map;

//extern map<string, vector<double>> bulk_counts_map;
//extern map<string, map<string, vector<double>>> Compart_bulk_counts_map;
//extern map<string, vector<int>> Container_counts_map;
//extern map<string, map<string, vector<int>>> Compart_container_counts_map;


int concType = -1;
//string metabolite_selection;
int showEnvType = -1;
float resolution = 0.25;

map<string, float> small_mol_r;
map<string, float> small_mol_g;
map<string, float> small_mol_b;

SimEngine :: SimEngine(){

    cell = GenSimulation::getSimulation();

    lEnv = LatticeEnvironment::getLatticeEnv();

    continueSimulation = false;
    showMetaConc = true;
    showMetaEnz = true;
    showEnv = false;
    showEnvCyt = true;
    showEnv1 = true;
    showEnv2 = true;
    showEnv3 = true;
    showEnv4 = true;

    // fill the display-compartments bool map
    for (map<string, Compart>::iterator compart_it = compartments_map.begin(); compart_it != compartments_map.end(); compart_it++){
        disp_compart_toggle[compart_it->first] = false;
    }
    // fill the metabolites bool map
    for(map<string, Species_Attributes *> ::iterator species_itr = species_details_map.begin(); species_itr!= species_details_map.end(); species_itr++) {
        if (species_itr->second->is_simple_molecule()){
            disp_metabolite_toggle[species_itr->first] = false;
        }
    }
    // fill the proteins bool map
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >  mTable = cell->getMoleculeTable();//simulation->cell->getMoleculeTable();
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >::iterator container_table_itr = mTable.begin();
    for (container_table_itr = mTable.begin(); container_table_itr != mTable.end(); container_table_itr ++){
        disp_protein_toggle[container_table_itr->first] = true;
    }
    
    
    for (string & metabolites_itr : small_molecule_id_list) {
        float small_mol_r = species_details_map[metabolites_itr]->color_red / 255; //small_mol_r[small_mol_itr] = species_details_map[small_mol_itr]->color_red / 255;
        float small_mol_g = species_details_map[metabolites_itr]->color_green / 255; //small_mol_g[small_mol_itr] = species_details_map[small_mol_itr]->color_green / 255;
        float small_mol_b = species_details_map[metabolites_itr]->color_blue / 255; //small_mol_b[small_mol_itr] = species_details_map[small_mol_itr]->color_blue / 255;
        vector<float> small_mol_colours;
        small_mol_colours.push_back(small_mol_r);
        small_mol_colours.push_back(small_mol_g);
        small_mol_colours.push_back(small_mol_b);
        
        graphic_display_colours[metabolites_itr] = small_mol_colours;
    
        cout << metabolites_itr << " red:" << small_mol_r << " green:" << small_mol_g << " blue:" << small_mol_b << endl;
    }

    // sort out which species will have a special display on the membrane
    for(auto & species : species_details_map) {
        if(species.second->membrane_display) {
            membrane_species_vector.push_back(species.first);
        }
    }
}

SimEngine :: ~SimEngine(){
}

// executes 1 loop of the continue sim timestep only. relies on external endless loops to keep calling this
// for graphics, see if the run simulation checkmark is on
// for no-graphics, just have engine->simulate while true.
bool SimEngine :: simulate(){
    if(continueSimulation) {
        if(!(cell->isSimulationCompleted())){
            cell->continueSimulation();
            // displayTotalConcentrations(); // Displays small molecule concentrations in stdout
        }
    }
    if(cell->isSimulationCompleted()) {
        return true;
    }
}

void SimEngine :: displayTotalConcentrations(){
    Matrix * state = lEnv->getStateMatrix();
    Tuple * t = state->collapse();
    t->display();
    delete t;
}

string SimEngine :: statsToString(){
    Matrix * stats = lEnv->getStatisticsMatrix();
    stringstream stat("");
    for(auto & species_attr_pair : species_details_map) {
        Species_Attributes * reactant_details = species_attr_pair.second;
        if(reactant_details->is_simple_molecule()){
            //stat << Species_getName(metItr->second) << '\t';
            stat << species_attr_pair.first << '\t';
        }
    }
    stat << endl;    
    stat << stats->displayToString();
    return stat.str();
}

#if 1

vector<float> SimEngine::set_display_colour(int object_count){
    // function to set the colour based on the ID supplied.  
    // ID is likely the count of each iterative item.
    srand(time(NULL) + object_count);
    float r_colour = rand()/float(RAND_MAX); // % 256 + 1;
    float g_colour = rand()/float(RAND_MAX); // % 256 + 1;
    float b_colour = rand()/float(RAND_MAX); // % 256 + 1;

    cout << "new colour generated r: " << (float) r_colour << " g: " << (float) g_colour << " b: " << (float) b_colour << endl;
    vector<float> return_me {r_colour, g_colour, b_colour};
    return return_me;

}


void SimEngine :: setDisplayEnvironment(){
    // Set properties of the surface material
    GLfloat mat_ambient[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat mat_diffuse[] = {0.6f, 0.6f, 0.6f, 1.0f};
    GLfloat mat_specular[] = {0.0f, 0.0f, 0.0f, 0.5f};
    GLfloat mat_shininess[] = {50.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel (GL_SMOOTH);

}

void SimEngine :: get_and_set_color_display(Molecule_Container * container){
    float r = 0.0f, g = 0.0f, b = 0.0f;
    r = container->container_color["red"] / 255;
    g = container->container_color["green"] / 255;
    b = container->container_color["blue"] / 255;
    // set molecule color to background color if invisible (rgb all at 0)
    if(r == 0 && g == 0 && b == 0) {
        r = background_color[0];
        g = background_color[1];
        b = background_color[2];
    }
    glColor3f(r, g, b);
}

vector<float> SimEngine :: get_default_container_color(string container_name){
    if(container_color_map.find(container_name) == container_color_map.end()) {
        cerr << "can't find colors for " << container_name << endl;
        return vector<float> {0,0,0};
    } 
    return container_color_map[container_name];
}

void SimEngine :: draw_Molecule_Container(Molecule_Container * mol){

    float xScale = (float)xDim/2;
    float yScale = (float)yDim/2;
    float zScale = (float)zDim/2;

    glPushMatrix();
    glTranslated(mol->pos.x-xScale,mol->pos.y-yScale,mol->pos.z-zScale);
    glutSolidSphere(0.02*resolution*zDim,4,4);
    glPopMatrix();
}


// if specified in input file, display membrane proteins as density gradients rather than distinct particles
void SimEngine::display_membrane_gradient() {
    // update the container brightness map, in case new containers were introduced
    for(auto & container : container_components_map) {
        vector<string> & container_vector = container.second;
        int counter = 0;
        for(auto & mem_species : membrane_species_vector) {
            counter += count(container_vector.begin(), container_vector.end(), mem_species);
        }
        if(counter > 0) {
            mem_container_bright_map[container.first] = counter;
        }
    }
    float xScale = ((float)xDim)/2;
    float yScale = ((float)yDim)/2;
    float zScale = ((float)zDim)/2;
    // iterating through each membrane voxel, illuminate it based on the concentration of membrane proteins
    for(auto & compart : compartments_map) {
        if(compart.second.is_membrane) {
            const Compart & cur_membrane = compart.second;
            for(auto & point : cur_membrane.list_of_voxels) {
                LatticeSite * cur_voxel = lEnv->getLatticeSite(*point);
                float intensity[3]; // Used for color display, array elements correspond to RGB values of the current c-voxel.
                intensity[0] = 0.0f; // RGB values
                intensity[1] = 0.0f;
                intensity[2] = 0.0f;
                for(auto & container_type : cur_voxel->lattice_container_map) {
                    auto const & container_type_map = container_type.second;
                    int alt_counter = 0;
                    // make sure the alternate color display rules are followed on a per-container basis
                    for(auto & container_pair : container_type_map) {
                        Molecule_Container * cur_container = container_pair.second;
                        if(cur_container->get_config_status()) {
                            alt_counter++; // track number of containers that follow special rules
                            // perform alternate display rules when appropriate
                            if(cur_container->show_membrane_conc_colors()) {
                                intensity[0] += cur_container->container_color["red"] / 10 / 255;
                                intensity[1] += cur_container->container_color["green"] / 10 / 255;
                                intensity[2] += cur_container->container_color["blue"] / 10 / 255;
                            }
                        }
                    }
                    if(mem_container_bright_map.find(container_type.first) == mem_container_bright_map.end()) continue;
                    // add the colors of membrane particles to the voxel display
                    // moles of particle per area of voxel
                    double container_density = container_type_map.size() - alt_counter; // / 6.02e23) / (spacescale * spacescale);
                    if(container_density == 0) continue;
                    vector<float> container_colors = get_default_container_color(container_type.first);
                    intensity[0] += container_density * container_colors[0] / 10;
                    intensity[1] += container_density * container_colors[1] / 10;
                    intensity[2] += container_density * container_colors[2] / 10;
                }
                if (intensity[0] > 0.02 || intensity[1] > 0.02 || intensity[2] > 0.02) {
                    GLfloat METABOLITECONC[] = {intensity[0], intensity[1], intensity[2], 0.0};
                    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, METABOLITECONC);

                    glPushMatrix();
                    glTranslated(point->x - xScale, point->y - yScale, point->z - zScale);
                    glutSolidCube(1.0);
                    glPopMatrix();
                    
                }
            }
        }
    }
}

// draw particles into the sim
void SimEngine :: draw_particles(){
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >  mTable = cell->getMoleculeTable();

    // Iterate the table of molecules that are in the lattice
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >::iterator container_table_itr = mTable.begin();

    // For each type of molecules in the lattice
    for (container_table_itr = mTable.begin(); container_table_itr != mTable.end(); container_table_itr ++)
    {
        if(disp_protein_toggle[container_table_itr->first]){
            unordered_map<unsigned long, Molecule_Container *>::iterator container_type_itr;
            // For each molecule of that type
            for(container_type_itr = container_table_itr->second.begin(); container_type_itr !=  container_table_itr->second.end() ; container_type_itr++){
                Molecule_Container * m_container = container_type_itr->second;

                if(m_container->show_membrane_conc_colors() && m_container->is_in_membrane()) continue;
                glDisable(GL_LIGHTING);
                get_and_set_color_display(m_container);
                // don't know what colour to supply this.  going to pick the first submmember
                //Molecule * m = m_container->container_vec[0];
                draw_Molecule_Container(m_container);
                glEnable(GL_LIGHTING);
            }
        }
        
    }
}

void SimEngine::master_highlight(){
    //merges all highlighting loops into 1
    time_t now = time(0);
    char* dt = ctime(&now);

    const static float max_mole_display_threshold = lEnv->conc_to_mole(2e-5);
    

    float xScale = ((float)xDim)/2;
    float yScale = ((float)yDim)/2;
    float zScale = ((float)zDim)/2;


    //for the entire sim space... cycle through each voxel, decide what colour to highlight, then check if we're at a voxel to highlight.
    for(int i = 0 ; i < xDim ; i ++){
    for(int j = 0 ; j < yDim ; j ++){
    for(int k = 0 ; k < zDim ; k ++){

        bool display = false;
        float r = 0;
        float g = 0;
        float b = 0;

        //from highlight_metabolites
        Point p(i,j,k);
        LatticeSite * ls = lEnv->getLatticeSite(p);
        int id = ls->getID();
        auto & ls_small_mol_counts = ls->bulk_moles_map;
        
        double mole;
        bool display_flag = false;
        float intensity[3]; // Used for color display, array elements correspond to RGB values of the current c-voxel.
        intensity[0] = 0.0f; // RGB values
        intensity[1] = 0.0f;
        intensity[2] = 0.0f;
        //double display_amount = 0;
        int rgb_division = 0;
        //if(concType != -1) // concType tells which show metabolite button is currently selected, -1 means show all.
        for (map<string, bool> ::iterator metabolite_it = disp_metabolite_toggle.begin(); metabolite_it != disp_metabolite_toggle.end(); metabolite_it++){
            string metabolite_name = metabolite_it->first;
            //if(metabolite_it->first != "All"){
            if(metabolite_it->second){
                mole = ls_small_mol_counts[metabolite_name];//[small_molecule_int_to_id_map[concType-1]];
                if(lEnv->mole_to_particle(mole) < 1) continue;
                double display_amount = mole / max_mole_display_threshold;

                if(display_amount > 1) display_amount = 1;
                // add something to scale colour by concentration (normalization)  r, g, b is max of 1
                intensity[0] += display_amount * graphic_display_colours[metabolite_name][0];//small_mol_r[metabolite_name];
                intensity[1] += display_amount * graphic_display_colours[metabolite_name][1];//small_mol_g[metabolite_name];
                intensity[2] += display_amount * graphic_display_colours[metabolite_name][2];//small_mol_b[metabolite_name];
                if(display_amount != 0){
                    rgb_division++;
                }

            }
  
            if(metabolite_it->second){
                display_flag = true;
                
            }
                
            
            if(rgb_division > 0){
                intensity[0] = intensity[0] / rgb_division; 
                intensity[1] = intensity[1] / rgb_division; 
                intensity[2] = intensity[2] / rgb_division;
            }
        }
        
        if(display_flag){
            if (intensity[0] > 0.01 || intensity[1] > 0.01 || intensity[2] > 0.01) {
                GLfloat METABOLITECONC[] = {intensity[0], intensity[1], intensity[2], 0.0};
                glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, METABOLITECONC);

                glPushMatrix();
                glTranslated(i-xScale,j-yScale,k-zScale);
                glutSolidCube(1.0);
                glPopMatrix();
                
            }
        }

        // highlight_compartment

        vector<float> colour_map;
        
        int compart_intensity = 0;
        for(auto & show_compart: disp_compart_toggle) {
            string compart_name = show_compart.first;
            time_t now = time(0);
            char* dt = ctime(&now);
            if(show_compart.second){
                Compart selected_compartment = compartments_map[show_compart.first];
                if(show_compart.first == "off"){
                    display = false;
                }   
                else{
                    map<string, vector<float>>::iterator colours_itr;
                    colours_itr = graphic_display_colours.find(compart_name);
                    
                    if(colours_itr != graphic_display_colours.end()){
                        colour_map = graphic_display_colours[compart_name];
                    }
                    else{
                        colour_map = set_display_colour(object_display_count);
                        cout << "new display clr set for compart: " << show_compart.first << " R[" << colour_map[0] << "]G[" << colour_map[1] << "]B[" << colour_map[2] << "]";
                        graphic_display_colours[show_compart.first] = colour_map;
                    }
                    
                    compart_intensity = 2;// / 8;
                    //cout << "intensity: " << intensity << endl;
                    r = (r == 0) ? ((colour_map[0] + compart_intensity) * 0.1) : (r + (colour_map[0] + compart_intensity) * 0.1);
                    g = (g == 0) ? ((colour_map[1] + compart_intensity) * 0.1) : (g + (colour_map[1] + compart_intensity) * 0.1);
                    b = (b == 0) ? ((colour_map[2] + compart_intensity) * 0.1) : (b + (colour_map[2] + compart_intensity) * 0.1);
                    
                    //check to see if this voxel is to be displayed
                    for (auto & disp_voxel: selected_compartment.list_of_voxels){
                        if((i == disp_voxel->x) && (j == disp_voxel->y) && (k == disp_voxel->z)){
                            //cout << show_compart.first << " Highlight:" + to_string(disp_voxel->x) + "," + to_string(disp_voxel->y) +  "," + to_string(disp_voxel->z) << endl;
                            display = true;
                            break;
                        }   
                    }
                }
            }
            object_display_count++;
        }
        

        if(display){
            GLfloat ENV[] = {r, g, b, 0.0};
            glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, ENV);

            float xScale = (float)xDim/2;
            float yScale = (float)yDim/2;
            float zScale = (float)zDim/2;

            glPushMatrix();
            glTranslated(i-xScale,j-yScale,k-zScale);
            glutSolidSphere(0.4,5,5);
            glPopMatrix();
        }

    //    }
    }
    }
    }

}

void SimEngine :: draw_sim_box(){
    // Construct the Cell Volume
    glPushMatrix();
    int maxDim = max(xDim,max(yDim,zDim));
    glScaled(((float)xDim)/maxDim,((float)yDim)/maxDim,((float)zDim)/maxDim);
    glTranslated(-0.5, -0.5, -0.5);
    glDisable(GL_LIGHTING);
    glColor3f(0.5,0.5,0.5);
    glutWireCube(maxDim);
    glEnable(GL_LIGHTING);
    glPopMatrix();

    GLfloat outwardFace[] = {0.0, 0.25, 0.25, 0.0};
    GLfloat inwardFace[]  = {0.0, 0.10, 0.10, 0.0};
}

void SimEngine :: constructLattice(){
//fires at every screen
    glEnable (GL_BLEND);
    glBlendFunc (GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc (GL_ONE_MINUS_SRC_ALPHA, GL_ONE);

    glPushMatrix();

    draw_particles();
    draw_sim_box();
    master_highlight();

    if(!membrane_species_vector.empty()){
        display_membrane_gradient();
    }

    

    //the metabolites
    //highlight_metabolites();

    //the compartments
    //highlight_compartment();

    glPopMatrix();

    glDisable (GL_BLEND);
}

void SimEngine :: displayStats(){
    
    Matrix * stats = lEnv->getStatisticsMatrix();
    int num = stats->getNumRows();

    Tuple * t = stats->getRow(num-1);

    float max = 800;
    float one = 1;

    float length = 50.0;
    float factor = 1;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -20.0);
    glTranslated(0.0,-7.0,0.0);

    glDisable(GL_LIGHTING);

    for(int i = 0 ; i < numMetabolite ; i ++){
        glPushMatrix();
        glColor3f(r[i],g[i],b[i]);
        factor = 1.0 * min(t->tuple[i]/max,one) * length;
        glScaled(factor,1.0,1.0);
        glutSolidCube(0.2);
        glPopMatrix();

        glTranslated(0.0,-0.5,0.0);
    }
    glEnable(GL_LIGHTING);
}

void SimEngine :: displayGraph(){
    //this function graphs the data.
    // june 28, 2021
    // need to remove all connections to lEnv stats matrix.  it's no longer used

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
     gluOrtho2D( 0.0, 6500.0, 0.0, 6500.0 );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef( 25.0, 50.0, 0.0 );
    glDisable(GL_LIGHTING);
    //also need to change this to display protein concentrations too, if selected
    string metName;
    string graph_display_name;
    int colour_count = 1; //must be larger than 0

    vector<string> selected_proteins; //stores the names of things to be graphed
    vector<string> selected_metabolites;
    vector<string> selected_comparts;

    map<string, vector<int>> graph_display_proteins;  //stores the graphing data
    map<string, vector<double>> graph_display_metabolites;

    
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >  mTable = cell->getMoleculeTable();//simulation->cell->getMoleculeTable();
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >::iterator container_table_itr;// = mTable.begin();
    for (container_table_itr = mTable.begin(); container_table_itr != mTable.end(); container_table_itr ++){
        string protein_name = container_table_itr->first;
        if(disp_protein_toggle[protein_name]){
            //cout << "displaying protein in graph: " << protein_name << endl;
            graph_display_name = protein_name;
            selected_proteins.push_back(protein_name);
            vector<float> colour_vec;
            std::map<string, vector<float>>::iterator disp_colour_itr = graphic_display_colours.find(protein_name);
            if(disp_colour_itr != graphic_display_colours.end()){
                colour_vec = graphic_display_colours[protein_name];
            } else{
                colour_vec = set_display_colour(colour_count);
                graphic_display_colours[protein_name] = colour_vec;
            }
            cout << "colour for[" << protein_name << "] R[" << colour_vec[0] << "]G[" << colour_vec[1] << "]B[" << colour_vec[2] << "]" << endl;
            glColor3f(colour_vec[0], colour_vec[1], colour_vec[2]);
            colour_count++;
            for (int letter_count = 0; letter_count < graph_display_name.length(); letter_count++){
                glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, graph_display_name[letter_count]);
            }
            glColor3f(1.0, 1.0, 1.0);
            glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, '|');
            
        }
    }
    
    for(map<string, Species_Attributes *> ::iterator species_itr = species_details_map.begin(); species_itr!= species_details_map.end(); species_itr++) {
        if (species_itr->second->is_simple_molecule()){
            string metabolite_name = species_itr->first;
            if(disp_metabolite_toggle[metabolite_name]){
                //cout << "displaying metabolite in graph: " << metabolite_name << endl;
                graph_display_name = metabolite_name;
                selected_metabolites.push_back(metabolite_name);
                std::map<string, vector<float>>::iterator disp_colour_itr = graphic_display_colours.find(metabolite_name);
                vector<float> colour_vec;
                if(disp_colour_itr != graphic_display_colours.end()){
                    colour_vec = graphic_display_colours[metabolite_name];
                } else{
                    colour_vec = set_display_colour(colour_count);
                    
                    graphic_display_colours[metabolite_name] = colour_vec;
                }

                glColor3f(colour_vec[0], colour_vec[1], colour_vec[2]);
                colour_count++;
                for (int letter_count = 0; letter_count < graph_display_name.length(); letter_count++){
                    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, graph_display_name[letter_count]);
                }
                glColor3f(1.0, 1.0, 1.0);
                glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, '|');
            }
        }
    }

    
    for (map<string, Compart>::iterator compart_it = compartments_map.begin(); compart_it != compartments_map.end(); compart_it++){
        if(disp_compart_toggle[compart_it->first]){
            selected_comparts.push_back(compart_it->first);
            //cout << "graph display: compartment selected: " << compart_it->first << endl;
        }
    }
    //sanity check on no-graph
    if(selected_metabolites.empty()){
        if(selected_proteins.empty()){
            //cout << "nothing to be graphed" << endl;
            return;
        }
    }

    if(selected_comparts.empty()){
        vector<string>::iterator protein_itr;
        vector<string>::iterator metabolite_itr;
        for(protein_itr = selected_proteins.begin(); protein_itr != selected_proteins.end(); protein_itr++){
            int orig_size = cell->Container_counts_map[*protein_itr].size();
            graph_display_proteins[*protein_itr] = cell->Container_counts_map[*protein_itr];
            //cout << "no compart: " << *protein_itr << " orig size: " << orig_size << endl; 
        }   
        for(metabolite_itr = selected_metabolites.begin(); metabolite_itr != selected_metabolites.end(); metabolite_itr++){
            graph_display_metabolites[*metabolite_itr] = cell->bulk_counts_map[*metabolite_itr];
            int orig_size = cell->bulk_counts_map[*metabolite_itr].size();
            //cout << "no compart: " << *metabolite_itr << " orig size: " << orig_size << endl; 
            
        } 
    }else{
        //cout << "compartment selection not empty" << endl;
        vector<string>::iterator compart_itr;
        for (compart_itr = selected_comparts.begin(); compart_itr != selected_comparts.end(); compart_itr++){
            
            map<string, map<string, vector<int>>>::iterator search_compart_protein_itr = cell->Compart_container_counts_map.find(*compart_itr);
            if(search_compart_protein_itr != cell->Compart_container_counts_map.end()){
                //cout << "found " << *compart_itr << " in compart protein map" << endl;  
                map<string, vector<int>> selected_inner_protein_map = cell->Compart_container_counts_map[*compart_itr];  
                vector<string>::iterator protein_itr;
                for(protein_itr = selected_proteins.begin(); protein_itr != selected_proteins.end(); protein_itr++){
                    map<string, vector<int>>::iterator search_protein_itr = selected_inner_protein_map.find(*protein_itr);
                    if(search_protein_itr != selected_inner_protein_map.end()){
                        //cout << "compartment[" << *compart_itr << "]: has data for protein[" << *protein_itr << "]" << endl;
                        map<string, vector<int>>::iterator search_graph_display_protein_itr = graph_display_proteins.find(*protein_itr);
                        if(search_graph_display_protein_itr != graph_display_proteins.end()){
                            //cout << "protein[" << *protein_itr << "] already found in graphable proteins map in another compartment.  combining from [" << *compart_itr << "]" << endl;
                            std::transform(
                                selected_inner_protein_map[*protein_itr].begin(), 
                                selected_inner_protein_map[*protein_itr].end(), 
                                graph_display_proteins[*protein_itr].begin(),  
                                graph_display_proteins[*protein_itr].begin(),
                                std::plus<int>()
                            );
                        } else{
                            //cout << "protein[" << *protein_itr << "freshly added to graphable protein map, from [" << *compart_itr << "]" << endl;
                            graph_display_proteins[*protein_itr] = selected_inner_protein_map[*protein_itr];
                        }

                        
                    } //else{
                        //cout << "compartment[" << *compart_itr << "]: has NO data for protein[" << *protein_itr << "]" <<  endl;
                    //}
                }
            }

            map<string, map<string, vector<double>>>::iterator search_compart_metabolite_itr = cell->Compart_bulk_counts_map.find(*compart_itr);
            if(search_compart_metabolite_itr != cell->Compart_bulk_counts_map.end()){
                //cout << "found " << *compart_itr << " in compart metabolite map" << endl;
                map<string, vector<double>> selected_inner_metabolite_map = cell->Compart_bulk_counts_map[*compart_itr];
                vector<string>::iterator metabolite_itr;
                for(metabolite_itr = selected_metabolites.begin(); metabolite_itr != selected_metabolites.end(); metabolite_itr++){
                    map<string, vector<double>>::iterator search_metabolite_itr = selected_inner_metabolite_map.find(*metabolite_itr);
                    if(search_metabolite_itr != selected_inner_metabolite_map.end()){
                        //cout << "compartment[" << *compart_itr << "]: has data for metabolite[" << *metabolite_itr << "]" << endl;
                        map<string, vector<double>>::iterator search_graph_display_metabolite_itr = graph_display_metabolites.find(*metabolite_itr);
                        if(search_graph_display_metabolite_itr != graph_display_metabolites.end()){
                            //cout << "metabolite[" << *compart_itr << "] already found in graphable metabolite map in another compartment. combining from [" << *compart_itr << "]" << endl;
                            std::transform(
                                selected_inner_metabolite_map[*metabolite_itr].begin(),
                                selected_inner_metabolite_map[*metabolite_itr].end(),
                                graph_display_metabolites[*metabolite_itr].begin(),
                                graph_display_metabolites[*metabolite_itr].begin(), 
                                std::plus<double>()
                            );
                        } else{
                            //cout << "metabolite[" << *metabolite_itr << " freshly added to graphable metabolite map. from[" << *compart_itr << "]" << endl;
                            graph_display_metabolites[*metabolite_itr] = selected_inner_metabolite_map[*metabolite_itr];
                        }
                    } //else {
                        //cout << "compartment[" << *compart_itr << "]: has NO data for metabolite[" << *metabolite_itr << "]" << endl;
                    //}
                }

            }
            
        }

    }
    
    

    glEnable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-0.5, 0.5, -0.5, 0.5, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -40.0);

    float x = -19;
    float y = -16;
    float height = 5.0;
    float scale = 0.04;

    glDisable(GL_LIGHTING);

    glPointSize(2);

    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINES);
        glVertex3f(x,y+6.5*height,0.0);
        glVertex3f(x,y,0.0);
    glEnd();
    glBegin(GL_LINES);
        glVertex3f(x,y,0.0);
        glVertex3f(-x,y,0.0);
    glEnd();

    float alter = 0;
    int time_steps = 0;
    if(!graph_display_proteins.empty()){
        time_steps = graph_display_proteins.begin()->second.size();
        //cout << "timestep first read: " << graph_display_proteins.begin()->first << ": " << time_steps << endl;
    } else {
        time_steps = graph_display_metabolites.begin()->second.size();
        //cout << "timestep first read: " << graph_display_metabolites.begin()->first << ": " << time_steps << endl;
    } 

    
    

    for (int i = 0; i < time_steps; i++){
        glBegin(GL_POINTS);
        map<string, vector<int>>::iterator graph_disp_protein_itr;
        map<string, vector<double>>::iterator graph_disp_metabolite_itr;
        //graphing needs normalization.  y+10 is already huge

        for (graph_disp_protein_itr = graph_display_proteins.begin(); graph_disp_protein_itr != graph_display_proteins.end(); graph_disp_protein_itr++){
            vector<float> colour_vec = graphic_display_colours[graph_disp_protein_itr->first];
            glColor3f(colour_vec[0], colour_vec[1], colour_vec[2]);
            //glColor3f(1.0,1.0,1.0);
            vector<int> display_protein_vec = graph_disp_protein_itr->second;
            glVertex3f(x+alter, y+(display_protein_vec[i] * scale), 0.0f);
            //glVertex3f(x+alter, y+(display_protein_vec[i]), 0.0f);
            //glVertex3f(x+alter, y+(5*inner_count * height), 0.0f);
            
            //cout << graph_disp_protein_itr->first <<  " x: " << x+alter << " y: " << y+(display_protein_vec[i]) << endl;

        }
        if(!graph_display_metabolites.empty()){

            for (graph_disp_metabolite_itr = graph_display_metabolites.begin(); graph_disp_metabolite_itr != graph_display_metabolites.end(); graph_disp_metabolite_itr++){
                vector<float> colour_vec = graphic_display_colours[graph_disp_metabolite_itr->first];
                glColor3f(colour_vec[0], colour_vec[1], colour_vec[2]);
                vector<double> display_metabolite_vec = graph_disp_metabolite_itr->second;
                glVertex3f(x+alter, y+(display_metabolite_vec[i] * scale), 0.0f);

            }
        }
        
        alter +=0.02;
        glEnd();
    }
    glEnable(GL_LIGHTING);

    
}

void SimEngine :: displayNumCycles() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D( 0.0, 6500.0, 0.0, 6500.0 );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef( 25.0, 50.0, 0.0 );
    glDisable(GL_LIGHTING);
    glColor3f( 0.0, 0.8, 1.0 );

    char buf[80];
    sprintf( buf, "%d", cell->getCycles() );

    int i = 0;
    while (buf[i] != '\0') {
        glutStrokeCharacter( GLUT_STROKE_MONO_ROMAN, buf[i] );
        i++;
    }
}


void SimEngine :: display(){
    // background display color
    glClearColor(background_color[0], background_color[1], background_color[2], 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);

    glPushMatrix();
    glTranslated(1000.0,1000.0,1000.0);
    glutSolidSphere(1.0,8,8);
    glPopMatrix();

    glPushMatrix();
    constructLattice();
    glPopMatrix();

//    glPushMatrix();
//    displayStats();
//    glPopMatrix();

    displayNumCycles();
}

#endif
