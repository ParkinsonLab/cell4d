//------------------------------------------------------------
// CellApp.cpp
// Donny Chan, Billy Taj
// This is the top-level of the cell4d simulator.  It's also the splash page.


#include "../inc/CellApp.h"
#include <ctime>
extern ParameterManager *pm;
//extern multimap<int, string> molEnzyme;
//extern multimap<int, string> molSubstrate;

extern int seed, max_cycles, inaccessible_space_percent;
extern int numSpecies;
bool from_checkpoint;
extern vector<Compart *> compartments;
extern map<string, Species_Attributes *> species_details_map;

int MODEL_SHAPE = 0; 
bool NUCLEUS_ENABLED = true, PLASMA_MEMBRANE_ENABLED = true; 

// global input parameters, can refactor in future
string model_filename;
string checkpoint_filename;
int new_timestep = 0;
double new_timescale = -1;
long new_seed = 0;
string whole_bulk_out_name = "total_bulk_stats.tsv";
string compart_bulk_out_name = "compart_bulk_stats.tsv";
string whole_particle_out_name = "total_container_stats.tsv";
string compart_particle_out_name = "compart_container_stats.tsv";
string checkpoint_out_name = "checkpoint";
string config_path = "none";
int checkpoint_every_timestep = 25000;
string model_editor_command = "firefox https://compsysbio.org/cell4d";
bool graphics = true;
bool multi_out = false;

bool poslog = false;
int poslog_interval = 500;
string poslog_comparts_in;
vector<string> poslog_comparts;
string poslog_out = "poslog.json";

#if 1

SimCanvas *simCanvas;
SimSidebar *simSidebar;

//make the object
IMPLEMENT_APP(CellApp)

// CellApp init, and effectively the entire program
bool CellApp::OnInit() {
    #ifndef MAC_OS
        int i = wxApp::argc;
        glutInit(&i, wxApp::argv);
    #endif

    // parser for command line arguments
    // the only mandatory argument is xml model
    auto cmdl = parameter_parsing(argc, argv);
    RandomNG::setSeed(new_seed);
    if(!graphics) {

        pm->loadSBML(model_filename.c_str(), from_checkpoint, checkpoint_filename);
        SimEngine *simulation = new SimEngine(); // starts the engine
        simulation->continueSimulation = true;
        // start point
        while (1) { 
            bool completed_sim = simulation->simulate();
            if(completed_sim) break;
        } 
        return 0;
    }
    else if (graphics) {
        if(!model_filename.empty()) {
            ShowSimulation(model_filename.c_str(), model_filename);
        } else {
            ShowMenu();
        }
    }
    return true;
}

void CellApp::ShowMenu() {
    menuFrame = new wxFrame((wxFrame *)NULL, wxID_ANY,  wxT("Cell4D Diffusion Simulation"), wxPoint(200, 200), MENU_SIZE);
        wxPanel *panel = new wxPanel(menuFrame, wxID_ANY);
            main_logo.LoadFile("assets/cell4d_sign.jpg", wxBITMAP_TYPE_JPEG);
            
            main_image = new wxStaticBitmap(panel, wxID_ANY, main_logo, wxPoint(10, 10), wxSize(250, 250));
            int text_x = 300;

            heading = new wxStaticText(panel, wxID_ANY, wxT("Cell4D"), wxPoint(text_x, 50), wxSize(200, 80), wxALIGN_LEFT);
            wxFont heading_font = heading->GetFont();
            heading_font.SetPointSize(40);
            //heading_font.SetPointSize(32);
            heading_font.SetWeight(wxFONTWEIGHT_BOLD);
            heading->SetFont(heading_font);

            //subheading = new wxStaticText(panel, wxID_ANY, wxT("Spatio-Temporal\nModeling Platform"), wxPoint(text_x, 110), wxSize(200, 40), wxALIGN_LEFT);
            subheading = new wxStaticText(panel, wxID_ANY, wxT("Spatio-Temporal\nModeling Platform"), wxPoint(text_x, 110), wxSize(200, 80), wxALIGN_LEFT);
            wxFont subheading_font = subheading ->GetFont();
            subheading_font.SetPointSize(18);
            //subheading_font.SetPointSize(13);
            subheading_font.SetWeight(wxFONTWEIGHT_NORMAL);
            subheading->SetFont(subheading_font);

            version_heading = new wxStaticText(panel, wxID_ANY, wxT("version 1.0.0"), wxPoint(text_x, 200), wxSize(200, 40), wxALIGN_LEFT);
            wxFont version_heading_font = version_heading->GetFont();
            version_heading_font.SetPointSize(10);
            version_heading_font.SetWeight(wxFONTWEIGHT_LIGHT);
            version_heading -> SetFont(version_heading_font);

            

            MACRO_INIT(-1);
            int menu_y = 275;
            int menu_header_x = 100;
            demo_label = new wxStaticText(panel, wxID_ANY, wxT("Demos"), wxPoint(menu_header_x, menu_y), wxSize(200, 40), wxALIGN_LEFT);
            wxFont label_font = demo_label->GetFont();
            label_font.SetPointSize(15);
            label_font.SetWeight(wxFONTWEIGHT_NORMAL);
            demo_label ->SetFont(label_font);

            model_label = new wxStaticText(panel, wxID_ANY, wxT("Models"), wxPoint(menu_header_x + 280, menu_y), wxSize(200, 40), wxALIGN_LEFT);
            model_label ->SetFont(label_font);

            menu_y += 40;
            int left_column_x = 35;
            int right_column_x = 310;
            BUTTON(panel, "Glycolysis",                 left_column_x, menu_y, 200, 25, CellApp::DemoGlycolysis);
            BUTTON(panel, "Workbench - Simulation",     right_column_x, menu_y, 200, 25, CellApp::SimulationFromFile);
            menu_y+= 30;

            BUTTON(panel, "Cascade",                    left_column_x, menu_y, 200, 25, CellApp::DemoCascade);
            BUTTON(panel, "Open Model Editor",          right_column_x, menu_y, 200, 25, CellApp::OpenEditor);
            menu_y+=30;

            BUTTON(panel, "CEACAM",                     left_column_x, menu_y, 200, 25, CellApp::DemoCeacam);
            menu_y+=30;
            
            BUTTON(panel, "AB Fusion",                  left_column_x, menu_y, 200, 25, CellApp::DemoABFusion);
            menu_y+=60;
            
            
            //BUTTON(panel, "Close simulation",         100, menu_y, 200, 25, CellApp::CloseSimulation);menu_y+=30;
            BUTTON(panel, "About",                      10, menu_y, 100, 25, CellApp::AboutButton);//menu_y+=30;
            BUTTON(panel, "License",                    100, menu_y, 100, 25, CellApp::LicenseButton);menu_y+=30;

    menuFrame->Show();
}
void CellApp::ShowAbout() {

    if (aboutFrame == NULL){
      
       aboutFrame = new wxFrame((wxFrame *)NULL, -1,  wxT("Cell4D Diffusion Simulation"), wxPoint(500, 100), ABOUT_SIZE);
    
       wxHtmlWindow *html = new wxHtmlWindow(aboutFrame);
       html->LoadPage("about.html");
       aboutFrame->Bind( wxEVT_CLOSE_WINDOW , &CellApp::AboutPageClosed, this);       
    }

    aboutFrame->Show(true);
  
   
}

void CellApp::AboutPageClosed(wxCloseEvent& evt){
    aboutFrame = NULL;
    evt.Skip();
}

void CellApp::ShowLicense(){
    if (licenseFrame == NULL){
      
        licenseFrame = new wxFrame((wxFrame *)NULL, -1,  wxT("Cell4D Diffusion Simulation"), wxPoint(500, 100), LICENSE_SIZE);
        ifstream licenseFile("LICENSE");
        stringstream buffer;
        buffer << licenseFile.rdbuf();
        wxTextCtrl *licenseInformation = new wxTextCtrl(licenseFrame, -1, buffer.str(), wxPoint(0,0), LICENSE_SIZE, wxTE_MULTILINE | wxHSCROLL );   
        licenseFrame->Bind( wxEVT_CLOSE_WINDOW , &CellApp::LicensePageClosed, this);       
    }

    licenseFrame->Show(true);
}

void CellApp::OpenEditor(wxCommandEvent& event){
    cout << "open editor button pushed. calling: " << model_editor_command << endl;
    system(model_editor_command.c_str());
}

void CellApp::LicensePageClosed(wxCloseEvent& evt){
    licenseFrame = NULL;
    evt.Skip();
}

////////////////////////////////////////////////////////////////////

void CellApp::ShowSimulation(const char *SBMLFilePath, wxString title) {
    time_t now = time(0);
    char* dt = ctime(&now);
    
    if(simCanvas){
        cout << "instance of simulation already run before " << dt;
        delete simFrame;
        delete simCanvas;
        delete simSidebar;

        
    }
    else{
        pm->loadSBML(SBMLFilePath, from_checkpoint, checkpoint_filename);;
        title = "Cell4D Diffusion Simulation - " + title;
        simFrame = new wxFrame((wxFrame *)NULL, -1,  title.c_str(), wxPoint(100, 100), SIM_SIZE);
        int args[] = {WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16, 0};
        simCanvas = new SimCanvas((wxFrame*) simFrame, args, SBMLFilePath);
        simSidebar = new SimSidebar((wxFrame*) simFrame);

        wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
        sizer->Add(simCanvas, 1, wxEXPAND);
        sizer->Add(simSidebar, 0, wxEXPAND);
        simFrame->SetSizer(sizer);
        simFrame->SetAutoLayout(true);
        simFrame->Show();
        
        
        if (menuFrame != NULL){
            menuFrame->Close();
        }   
        
        if (aboutFrame != NULL) {
            aboutFrame->Close();
        }
    
        if (pFrame != NULL) {
            pFrame->Close();
        }
        
        
        cout << "control returns to ShowSimulation " << dt; 
        //exit(0);
    }
}

void CellApp::CloseSimulation(wxCommandEvent& event){
    if(simFrame){
        if(simFrame->IsShownOnScreen()){
            cout << "close triggered" << endl;
            simFrame->Hide();            
        } else{
            cout << "simFrame exists but not shown onscreen" << endl;
        }
    } else{
        cout << "Close Simulation pressed but we're not reacting" << endl;
    }
}

void CellApp::SimulationFromFile(wxCommandEvent& event) {
    //opens up a new UI for file-making
    filename = wxFileSelector(wxT("Choose SBML XML Document"), wxT(""), wxT(""), wxT(""), 
        wxT("SBML XML Document (*.xml)|*.xml"), wxFD_OPEN);
    if (filename == "") return ;
    SBMLDocument_t *SBMLDoc = readSBML(filename);


    wxFrame * workbenchFrame = new wxFrame((wxFrame *)NULL, -1,  wxT("Workbench"), wxPoint(100, 100), wxSize(600, 500)); 
    wxPanel *panel = new wxPanel(workbenchFrame, wxID_ANY);

    new wxButton(panel, 207, "Launch Simulation", wxPoint(400, 400), wxSize(150, 20));
    Connect(205, wxEVT_BUTTON, wxCommandEventHandler(CellApp::ShowParamFrame));
    Connect(207, wxEVT_BUTTON, wxCommandEventHandler(CellApp::LaunchSimulation));

    workbenchFrame->Show();


}

void CellApp::LaunchSimulation(wxCommandEvent& event) {
    //This is the part where the file-maker hooks back to the regular flow of things
    ShowSimulation(filename.c_str(), filename);
}

void CellApp::ShowParamFrame(wxCommandEvent& event) {
    pFrame = new wxFrame((wxFrame *)NULL, -1,  wxT("One Time Parameters"), wxPoint(700, 100), wxSize(250, 600)); 
    paramFrame = new wxScrolledWindow(pFrame, -1, wxPoint(-1,-1), wxSize(250, 20000), wxScrolledWindowStyle, "One time Parameters");
    pm->loadSBML(filename);

    int y = 20;
    new wxStaticText(paramFrame, -1, "Seed (numbers only)", wxPoint(20, y)); y += 20;
    seedText = new wxTextCtrl(paramFrame, wxID_ANY, "RANDOM", wxPoint(20, y), wxSize(210, 20)); y += 40;

    new wxStaticText(paramFrame, -1, "Maximum cycles", wxPoint(20, y)); y += 20;
    maxCycleText = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(pm->get_int(MAX_CYCLES)), wxPoint(20, y), wxSize(210, 20)); y += 40;

    new wxStaticText(paramFrame, -1, "Time per step (in microseconds)", wxPoint(20, y)); y += 20;
    timeScaleText = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(int(pm->get_int(TIMESCALE))), wxPoint(20, y), wxSize(210, 20)); y += 40;

    new wxStaticText(paramFrame, -1, "Cell size (in nanometers)", wxPoint(20, y)); y += 20;
    spaceScaleText = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(int(pm->get_int(SPACESCALE))), wxPoint(20, y), wxSize(210, 20)); y += 40;

    new wxStaticText(paramFrame, -1, "Lattice Resolution (even number)", wxPoint(20, y)); y += 20;
    letticeDimText = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(pm->get_int(X_DIM)), wxPoint(20, y), wxSize(210, 20)); y += 40;

    new wxStaticText(paramFrame, -1, "Inaccessible space (percent)", wxPoint(20, y)); y += 20;
    inaccessibleSpaceText = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(pm->get_int(INACCESSIBLE_SPACE_PERCENT)), wxPoint(20, y), wxSize(210, 20)); y += 40;

    speciesInput = new wxTextCtrl * [numSpecies];

    int species_count = 0;
    for(map<string, Species_Attributes *> :: const_iterator reactant_itr = species_details_map.begin(); reactant_itr != species_details_map.end(); reactant_itr++){
        
        string species_name = reactant_itr->first + " initial amount";
        new wxStaticText(paramFrame, -1, species_name, wxPoint(20, y)); y += 20;
        speciesInput[species_count] = new wxTextCtrl(paramFrame, wxID_ANY, std::to_string(reactant_itr->second->current_total_amount), wxPoint(20, y), wxSize(210, 20)); y += 40;
        species_count++;
    }
    new wxStaticText(paramFrame, -1, "Model:", wxPoint(20, y)); y += 20;
    cubeButton = new wxRadioButton(paramFrame, wxID_ANY, "Cube", wxPoint(20, y), wxDefaultSize, wxRB_GROUP); y += 20;
    cubeButton->SetValue(true);
    oneLayerButton = new wxRadioButton(paramFrame, wxID_ANY, "Flat, one layer", wxPoint(20, y)); y += 20;
    threeLayerButton = new wxRadioButton(paramFrame, wxID_ANY, "Flat, three layer", wxPoint(20, y)); y += 40;

    (new wxButton(paramFrame, 99, "Save Parameters", wxPoint(20, y), wxSize(210, 20)))->SetDefault();
    Connect(99, wxEVT_BUTTON, wxCommandEventHandler(CellApp::SimParam));
    paramFrame->SetScrollbars(20,20,50,50); 
    paramFrame->SetVirtualSize(250, y+40);
    pFrame->Show();
}

void CellApp::DemoCascade(wxCommandEvent& event) {
    ShowSimulation("demo_files/Cascade.xml", "Cascade Ver.5"); 
}
void CellApp::DemoABFusion(wxCommandEvent& event){
    ShowSimulation("demo_files/AB_fusion_demo.xml", "AB fusion Ver.1");
}
void CellApp::DemoCeacam(wxCommandEvent& event){
    ShowSimulation("demo_files/CEACAM_transport.xml", "CEACAM Ver.1");
}
void CellApp::DemoGlycolysis(wxCommandEvent& event) {
    ShowSimulation("demo_files/Glycolysis_V5.xml", "Glycolysis Ver.5"); 
}
void CellApp::SimParam(wxCommandEvent& event) {
    if (seedText->GetValue() != "RANDOM" && seedText->GetValue() != "Random" && seedText->GetValue() != "random"){
        try {
            seed = stoi((string)seedText->GetValue());
        } catch (const std::invalid_argument& ia) {
            wxMessageBox( "Error found in seed \"" + seedText->GetValue() + "\"!", 
                      "One Time Parameters", wxOK | wxICON_ERROR );
            return;
        }
    }

    try {
        pm->add(MAX_CYCLES, stoi((string)maxCycleText->GetValue()));
    } catch (const std::invalid_argument& ia) {
        wxMessageBox( "Error found in maximum cycles \"" + maxCycleText->GetValue() + "\"!", 
                  "One Time Parameters", wxOK | wxICON_ERROR );
        return;
    }

    try {
        pm->add(X_DIM,stoi((string)letticeDimText->GetValue()));
        pm->add(Y_DIM,stoi((string)letticeDimText->GetValue()));
        pm->add(Z_DIM,stoi((string)letticeDimText->GetValue()));
    } catch (const std::invalid_argument& ia) {
        wxMessageBox( "Error found in lattice resolution \"" + letticeDimText->GetValue() + "\"!", 
                  "One Time Parameters", wxOK | wxICON_ERROR );
        return;
    }

    try {
        pm->add(TIMESCALE, stoi((string)timeScaleText->GetValue()));
    } catch (const std::invalid_argument& ia) {
        wxMessageBox( "Error found in time per step \"" + timeScaleText->GetValue() + "\"!", 
                  "One Time Parameters", wxOK | wxICON_ERROR );
        return;
    }

    try {
        pm->add(SPACESCALE, stoi((string)spaceScaleText->GetValue()));
    } catch (const std::invalid_argument& ia) {
        wxMessageBox( "Error found in cell size \"" + spaceScaleText->GetValue() + "\"!", 
                  "One Time Parameters", wxOK | wxICON_ERROR );
        return;
    }

    try {
        pm->add(INACCESSIBLE_SPACE_PERCENT, stoi((string)inaccessibleSpaceText->GetValue()));
    } catch (const std::invalid_argument& ia) {
        wxMessageBox( "Error found in inaccessible space \"" + inaccessibleSpaceText->GetValue() + "\"!", 
                  "One Time Parameters", wxOK | wxICON_ERROR );
        return;
    }
    int i = 0;

    if (oneLayerButton->GetValue()) {
        pm->add(Z_DIM,((int) pm->get_int(X_DIM) * .1) * 2);
       
    }
    if (threeLayerButton->GetValue()) {
         pm->add(Z_DIM,((int) pm->get_int(X_DIM) * .33) * 2);
       
    }
    std::ofstream out("input/params.in");

    pFrame->Close();
    pFrame = NULL;
}
void CellApp::AboutButton(wxCommandEvent& event) {
    ShowAbout();
}

void CellApp::LicenseButton(wxCommandEvent& event){
    ShowLicense();
}



#else

int main(int argc, char * argv[]) {
    graphics = false;
    // parser for command line arguments
    // the only mandatory argument is xml model
    auto cmdl = parameter_parsing(argc, argv);

    pm->loadSBML(model_filename.c_str(), from_checkpoint, checkpoint_filename);
    RandomNG::setSeed(new_seed);
    SimEngine *simulation = new SimEngine(); // starts the engine
    simulation->continueSimulation = true;
    // start point
    while (1) { 
        bool completed_sim = simulation->simulate();
        if(completed_sim) break;
    } 
    return 0;
    auto mol_Table = GenSimulation::getSimulation()->getMoleculeTable();
}

#endif

// split algorithm from tbundle
vector<string> string_split(const std::string& str, const std::string& delims = ",", bool skip_empty = true) {
    std::vector<std::string> output;
    auto first = str.cbegin();

    while (first != str.cend()) {
        const auto second = std::find_first_of(first, str.cend(),
                                               delims.cbegin(), delims.cend());
        if (first != second || !skip_empty) {
            output.emplace_back(first, second);
        }
        if (second == str.cend()) break;
        first =  std::next(second);
    }

    return output;
}



argh::parser parameter_parsing(int argc, char * argv[]) {

    // batch pre-register multiple params: name + value
    argh::parser cmdl({ "--whole-particle-output", "--config", "--compartment-particle-output", "--whole-bulk-output", "--compartment-bulk-output", 
    "--checkpoint", "--check-in", "--checkpoint-output", "--check-out", "--checkpoint-timesteps", "--timestep", "--timescale", "--seed" }); 
    cmdl.parse(argc, argv);


    if(cmdl["-h"] || cmdl["--help"] || cmdl["--usage"]) {
        usage(-1);
        exit(0);
    }

    // split outputs based on compartment
    if(cmdl["multi"]) {
        multi_out = true;
        cout << "compartment output files are split." << endl;
    } 

    // graphical option
    if(cmdl["no-graphics"]) {
        graphics = false;
        cout << "graphics are off." << endl;
    } 

    // if there are more than 1 positional arguments (the first being simulation executable), parse in order of model_input, checkpoint_input
    if(cmdl.pos_args().size() > 1) {
        if(cmdl[1].find(".xml") != string::npos) {
            model_filename = cmdl[1];
            cout << "loading " << model_filename << " as model file." << endl;
        } else {
            if(!graphics) {
                // usage help for missing model
                usage(0);
                exit(1);
            }
        }

        // look for checkpoint file to load in
        if(cmdl[2].find(".json") != string::npos) {
            checkpoint_filename = cmdl[2];
        }
    }
    //the config parser
    if(cmdl({ "--config" })) {
        config_path = cmdl({ "--config" }).str();
        cout << "Using " << config_path << " as the config file." << endl;
        ifstream config_file(config_path);
        string config_line;

        //flag info.  need both the length and the flag itself.  don't forget the lagging space
        string no_graphics_header = "--no-graphics ";
        string whole_out_header = "--whole-particle-output ";
        string compart_out_header = "--compartment-particle-output ";
        string whole_bulk_out_header = "--whole-bulk-output ";
        string compart_bulk_out_header = "--compartment-bulk-output ";
        string checkpt_in_header = "--checkpoint-input ";
        string checkpt_out_header = "--checkpoint-output ";
        string checkpt_timestep_header = "--checkpoint-timesteps ";
        string timestep_header = "--timestep ";
        string timescale_header = "--timescale ";
        string seed_header = "--seed ";
        string poslog_out_header = "--poslog-out ";
        string poslog_compart_header = "--poslog-comparts ";
        string poslog_interval_header = "--poslog-interval ";
        string model_editor_header = "--model-editor-command ";

        size_t no_graphics_length = no_graphics_header.length();
        size_t whole_out_length = whole_out_header.length();
        size_t compart_out_length = compart_out_header.length();
        size_t whole_bulk_out_length = whole_bulk_out_header.length();
        size_t compart_bulk_out_length = compart_bulk_out_header.length();
        size_t checkpt_in_length = checkpt_in_header.length();
        size_t checkpt_out_length = checkpt_out_header.length();
        size_t checkpt_timestep_length = checkpt_timestep_header.length();
        size_t timestep_length = timestep_header.length();
        size_t timescale_length = timescale_header.length();
        size_t seed_length = seed_header.length();
        size_t poslog_out_length = poslog_out_header.length();
        size_t poslog_compart_length = poslog_compart_header.length();
        size_t poslog_interval_length = poslog_interval_header.length();
        size_t model_editor_length = model_editor_header.length();

        bool poslog_out_flag = false;
        bool poslog_editor_flag = false;
        bool poslog_compart_flag = false;
        bool poslog_interval_flag = false;

        //walk through the file, search for the flags
        if(config_file.is_open()){
            while(getline(config_file, config_line)){
                size_t no_graphics_flag = config_line.find(no_graphics_header);
                size_t whole_out_flag = config_line.find(whole_out_header);
                size_t compart_out_flag = config_line.find(compart_out_header);
                size_t whole_bulk_out_flag = config_line.find(whole_bulk_out_header);
                size_t compart_bulk_out_flag = config_line.find(compart_bulk_out_header);
                size_t checkpt_in_flag = config_line.find(checkpt_in_header);
                size_t checkpt_out_flag = config_line.find(checkpt_out_header);
                size_t checkpt_timestep_flag = config_line.find(checkpt_timestep_header);
                size_t timestep_flag = config_line.find(timestep_header);
                size_t timescale_flag = config_line.find(timescale_header);
                size_t seed_flag = config_line.find(seed_header);
                size_t poslog_out_flag = config_line.find(poslog_out_header);
                size_t poslog_compart_flag = config_line.find(poslog_compart_header);
                size_t poslog_interval_flag = config_line.find(poslog_interval_header);
                size_t model_editor_flag = config_line.find(model_editor_header);


                if(no_graphics_flag != string::npos){
                    string sub = config_line.substr(no_graphics_flag + no_graphics_length);
                    cout << " no-graphics flag detected " << "substring: " << sub << endl;
                    if(sub == "true" || sub == "True"){
                        graphics = false;
                    } else{
                        graphics  = true;
                    }

                } else if(whole_out_flag != string::npos){
                    whole_particle_out_name = config_line.substr(whole_out_flag + whole_out_length);

                } else if(compart_out_flag != string::npos){
                    compart_particle_out_name = config_line.substr(compart_out_flag + compart_out_length);

                } else if(whole_bulk_out_flag != string::npos){
                    whole_bulk_out_name = config_line.substr(whole_bulk_out_flag + whole_bulk_out_length);

                } else if(compart_bulk_out_flag != string::npos){
                    compart_bulk_out_name = config_line.substr(compart_bulk_out_flag + compart_bulk_out_length);

                } else if(checkpt_in_flag != string::npos){
                    checkpoint_filename = config_line.substr(checkpt_in_flag + checkpt_in_length);
                    //the file name must end in .json
                    size_t find_ext = checkpoint_filename.find(".json");
                    if(find_ext != string::npos){
                        from_checkpoint = true;
                    } else{
                        std::cout << "invalid checkpoint file.  Must be a JSON. exiting" << endl;
                        exit(1);
                    }

                } else if(checkpt_out_flag != string::npos){
                    checkpoint_out_name = config_line.substr(checkpt_out_flag + checkpt_out_length);
                    size_t ext_found_at = checkpoint_out_name.find(".");
                    //clean up the filename
                    if(ext_found_at != string::npos) {
                        checkpoint_out_name.replace(checkpoint_out_name.find("."),checkpoint_out_name.length(),"");
                    }

                } else if(checkpt_timestep_flag != string::npos){
                    checkpoint_every_timestep = stoi(config_line.substr(checkpt_timestep_flag + checkpt_timestep_length));
                    if(checkpoint_every_timestep > 0) {
                        std::cout << "simulation will be checkpointed every " << checkpoint_every_timestep << " timesteps." << endl;
                    }

                } else if(timestep_flag != string::npos){
                    new_timestep = stoi(config_line.substr(timestep_flag + timestep_length));
                    if(new_timescale > 1e-1) {
                        std::cerr << "a timescale greater than 0.1 seconds is not recommended. exiting." << endl;
                        exit(1);
                    }

                } else if(timescale_flag != string::npos){
                    new_timescale = stod(config_line.substr(timescale_flag + timescale_length));

                } else if(seed_flag != string::npos){
                    new_seed = stol(config_line.substr(seed_flag + seed_length));

                } else if(poslog_out_flag != string::npos){
                    poslog = true;
                    poslog_out_flag = true;
                    poslog_out = config_line.substr(poslog_out_flag + poslog_out_length);
                    size_t ext_found_at = poslog_out.find(".");
                    if(ext_found_at != string::npos) {
                        poslog_out.replace(poslog_out.find("."),poslog_out.length(),"");
                    }
                    cout << "Position log file name set to " << poslog_out + ".json" << endl;


                } else if(poslog_compart_flag != string::npos){
                    poslog_comparts_in = config_line.substr(poslog_compart_flag + poslog_compart_length);
                    poslog_comparts = string_split(poslog_comparts_in);
                    for(auto & compart : poslog_comparts) {
                        if(compart == poslog_comparts.front()) {
                            cout << compart;
                        } else {
                            cout << ", " << compart;
                        }
                    }
                    poslog_compart_flag = true;

                } else if(poslog_interval_flag != string::npos){
                    poslog_interval = stoi(config_line.substr(poslog_interval_flag + poslog_interval_length));
                    
                    poslog_interval_flag = true;

                } else if(model_editor_flag != string::npos){
                    model_editor_command = config_line.substr(model_editor_flag + model_editor_length);
                    cout << "model editor is: " << model_editor_command << endl;

                }


            }
            // check the integrity of the poslog arguments
            if(!poslog_out_flag || !poslog_compart_flag || !poslog_interval_flag){
                cout << "All poslog options must be set: poslog-out, poslog-comparts, poslog-interval" << endl;
            }

        }
        else{
            cout << "cannot open config file at: " << config_path << " skipping" << endl;
        }
    }

    cout << endl;
    return(cmdl);
}

void usage(int usage_in) {
    // usage info
    // default: all help options

    // -m or --model for model xml
    // -s or --checkpoint for checkpoint file input
    // --whole-particle-output for name of total particle counts
    // --compartment-particle-output for name of compartmentalized particle counts
    // --no-graphics for graphics toggle, should be last option
    // -t or --timesteps or --timestep for model timestep override, specify number of timecycles for simulation
    switch (usage_in) {
    
    // missing model file
    case 0:
        cerr << "model file input not found, must have xml extension." << endl;
        break;
    
    // incorrect checkpoint json
    case 1:
        cerr << "checkpoint file must have json extension." << endl;
        break;

    // help
    default:
        cout << "Cell4D Stochastic Spatial Biological Simulator" << endl;
        cout << "Usage: simulation <model> [<sim-parameters>] [<flags>]" << endl;
        cout << "where <model> is a cell4D formatted XML file, and <sim-parameters> include:" << "\n\t";
        cout << "--check-in or --checkpoint: checkpointed molecule file from previous cell4D runs, formatted as json." << "\n\t";
        cout << "--whole-particle-output: name of particle count .tsv at each time step. default: " << whole_particle_out_name << "\n\t";
        cout << "--compartment-particle-output: name of particle count file for each compartment. default: " << compart_particle_out_name << "\n\t";
        cout << "--whole-bulk-output: name of bulk count .tsv at each time step. default: " << whole_bulk_out_name << "\n\t";
        cout << "--compartment-bulk-output: name of bulk count file for each compartment. default: " << compart_bulk_out_name << "\n\t";
        cout << "--check-out or --checkpoint-output: name of checkpoint file created. default: " << checkpoint_out_name << "\n\t";
        cout << "--checkpoint-timesteps: timestep interval for regular checkpointing. default: " << checkpoint_every_timestep << "\n\t";
        

        cout << "--timestep: number of time cycles to run this simulation, overrides MAX_CYCLES field in model." << "\n\t";
        cout << "--timescale: timestep length (seconds) to use in this simulation, overrides TIMESCALE field in model." << "\n\t";
        cout << "--seed: simulation seed to use for RNG." << endl;

        cout << "<flags> are:" << "\n\t";
        cout << "--no-graphics: turns off graphical display of simulation." << "\n\t";
        cout << "--multi: compartmental molecule output will be separate files by compartment rather than long concatenated file" << endl;

        cout << "<sim-parameters> can be called as \"--parameter_name\" value OR \"--parameter_name=value\"" << "\n\t";
        cout << endl;

        cout << "Experimental features:" << "\n\t";
        cout << "--poslog-out: name of molecule position log output" << "\n\t";
        cout << "--poslog-comparts: name of compartment of logged molecule positions" << "\n\t";
        cout << "--poslog-interval: interval of molecule log recordings, default at " << poslog_interval << endl;
        

        break;

    }
    

}
