//====================================================================
// SimSidebar.cpp
// Donny Chan, Billy Taj
//This code controls/specifies the sidebar of the simulation window.

#include "../inc/_wxHeader.h"
#include <ctime>
extern map<string, Species_Attributes *> species_details_map;
extern map<string, Compart>compartments_map;

extern SimCanvas *simCanvas;
extern int concType;
StatsCanvas *statsCanvas = NULL;

wxFrame *statsFrame;
extern SBMLDocument_t *sbml_document;
extern SimEngine *simulation;



//checklistbox position counters
int clb_proteins_position = 0;  
int clb_metabolites_position = 0;
int clb_compart_position = 0;


SimSidebar::SimSidebar(wxFrame *parent)
       : wxPanel(parent, wxID_ANY, wxDefaultPosition, 
         wxSize(200, 540), wxBORDER_SUNKEN) {

    
    MACRO_INIT(10);

    STATICTEXT( "Control Panel", 10)

   
    play_logo.LoadFile("assets/play.jpg", wxBITMAP_TYPE_JPEG);
    pause_logo.LoadFile("assets/pause.jpg", wxBITMAP_TYPE_JPEG);
    stop_logo.LoadFile("assets/stop.jpg", wxBITMAP_TYPE_JPEG);
    all_logo.LoadFile("assets/all.jpg", wxBITMAP_TYPE_JPEG);
    none_logo.LoadFile("assets/none.jpg", wxBITMAP_TYPE_JPEG);


    pbb_play = BITMAPBUTTON(this, play_logo, 10, SimSidebar::ContinueSimulation);
    pbb_pause = BITMAPBUTTON(this, pause_logo, 50, SimSidebar::PauseSimulation);
    pbb_stop = BITMAPBUTTON(this, stop_logo, 90, SimSidebar::StopSimulation);
    y+= 50;
    
    STATICTEXT("Simulation Speed", 10)
    spc_time_interval = new wxSpinCtrl(this, obj_id, "time interval", wxPoint(10, y), wxSize(160, 20), wxSP_ARROW_KEYS, 0, 160, 160);
    Connect(obj_id, wxEVT_SPIN, wxCommandEventHandler(SimSidebar::spc_simulation_speed));
    
    y+=20;
                            

    STATICLINE()
    STATICTEXT("Display Proteins", 10)
    // need to have custom names because reacted-particle names will be complicated
    pbb_all_proteins = BITMAPBUTTON(this, all_logo, 10, SimSidebar::display_all_proteins);
    pbb_no_proteins = BITMAPBUTTON(this, none_logo, 45, SimSidebar::display_no_proteins);
    y+=30;
    clb_proteins = CHECKLISTBOX(this);
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >  mTable = simulation->cell->getMoleculeTable();//simulation->cell->getMoleculeTable();
    unordered_map<string, unordered_map <unsigned long, Molecule_Container *> >::iterator container_table_itr = mTable.begin();

    for (container_table_itr = mTable.begin(); container_table_itr != mTable.end(); container_table_itr ++){
        clb_proteins->Insert(container_table_itr->first, clb_proteins_position);
        clb_proteins->Check(clb_proteins_position, true);
        clb_proteins_position++;
    }
    Connect(obj_id, wxEVT_CHECKLISTBOX, wxCommandEventHandler(SimSidebar::GetProteinSelection));
    obj_id++;

    STATICLINE()
    //Controls to display metabolites
    STATICTEXT("Display Metabolites", 10)
    pbb_all_metabolites = BITMAPBUTTON(this, all_logo, 10, SimSidebar::display_all_metabolites);
    pbb_no_metabolites = BITMAPBUTTON(this, none_logo, 45, SimSidebar::display_no_metabolites);
    y+= 30;
    clb_metabolites = CHECKLISTBOX(this);

    for(map<string, Species_Attributes *> ::iterator species_itr = species_details_map.begin(); species_itr!= species_details_map.end(); species_itr++) {
        if (species_itr->second->is_simple_molecule()){
            clb_metabolites->Insert(species_itr->first, clb_metabolites_position);
            
            clb_metabolites_position++;        
        }
    }

    Connect(obj_id, wxEVT_CHECKLISTBOX, wxCommandEventHandler(SimSidebar::GetMetaboliteSelection));
    obj_id++;
    STATICLINE()
    STATICTEXT("Display Compartments", 10)
    //make the combobox. populate it with the Compartments to highlight
    pbb_all_comparts = BITMAPBUTTON(this, all_logo, 10, SimSidebar::display_all_comparts);
    pbb_no_comparts = BITMAPBUTTON(this, none_logo, 45, SimSidebar::display_no_comparts);
    y+=30;
    clb_compart = CHECKLISTBOX(this);
    
    
    for (map<string, bool>::iterator it = simulation->disp_compart_toggle.begin(); it != simulation->disp_compart_toggle.end(); it++){
        clb_compart->Insert(it->first, clb_compart_position);
        clb_compart_position++;
        
     
    }
    Connect(obj_id, wxEVT_CHECKLISTBOX, wxCommandEventHandler(SimSidebar::GetCompartSelection));
    obj_id++;
    


                            STATICLINE()

    cb_rotate =             CHECKBOX_E( this, "Rotate",                    SimSidebar::Rotate);            cb_rotate->SetValue(true);
    // cb_translate =          CHECKBOX_E( this, "Translate",                 SimSidebar::Translate);         cb_translate->Disable();
    cb_scale =              CHECKBOX_E( this, "Scale",                     SimSidebar::Scale);             cb_scale->SetValue(true);
    
    cb_displayStats =       CHECKBOX_E( this, "Display Statistics",        SimSidebar::DisplayStats);
    //pb_endEarly =           BUTTON(     this, "End Early", 10, y, 175, 25, SimSidebar::EndEarly); y+=30;
    STATICLINE()
    pb_close =              BUTTON(     this, "Close", 10, y, 175, 25, SimSidebar::Close);
}

void SimSidebar::ContinueSimulation(wxCommandEvent& event){
    simulation->continueSimulation = true;
    simCanvas->Refresh();
}

void SimSidebar::PauseSimulation(wxCommandEvent& event){
    simulation->continueSimulation = false;
    simCanvas->Refresh();
}

void SimSidebar::StopSimulation(wxCommandEvent& event){
    simulation->continueSimulation = false;
    wxMessageDialog * exitEarly = new wxMessageDialog(this, "You have stopped the simulation.  Quit the program?", "Stop simulation", wxYES_NO);
    int option = exitEarly->ShowModal();
    if (option == wxID_YES){
        string new_command = "/home/billy/cell4d/cell-wxwidgets/CellApp --seed " + to_string(RandomNG::getSeed()); 
        
        exit(0);

    } else if (option == wxID_NO){
        return;
    }

}

void SimSidebar::SimulationSpeed(wxCommandEvent& event){
    simCanvas->simulation_timeout = (int) pow(2.0, (160-sl_timeInterval->GetValue())/16) - 1;
}

void SimSidebar::spc_simulation_speed(wxCommandEvent &event){
    simCanvas->simulation_timeout = (int) pow(2.0, (160-spc_time_interval->GetValue())/16) - 1;
    cout << "timestep orig[" << spc_time_interval->GetValue() << "]| real num[" << simCanvas->simulation_timeout << "]" << endl;
}

void SimSidebar::ShowMetaConc(wxCommandEvent& event){
    simulation->showMetaConc = cb_showMetaConc->IsChecked();
    simCanvas->Refresh();
}
void SimSidebar::ShowMetaEnz(wxCommandEvent& event){
    simulation->showMetaEnz = cb_showMetaEnz->IsChecked();
    simCanvas->Refresh();
}
void SimSidebar::ShowEnv(wxCommandEvent& event){
    time_t now = time(0);
    char* dt = ctime(&now);
    cout << "env shown: " << dt << endl;
    simulation->showEnv = cb_showEnv->IsChecked();
    simCanvas->Refresh();
}
void SimSidebar::ShowConcType(wxCommandEvent& event) {
    time_t now = time(0);
    char* dt = ctime(&now);
    concType = event.GetId()-100;
    
    simCanvas->Refresh();
}
void SimSidebar::GetComboboxSelection(wxCommandEvent& event){
    time_t now = time(0);
    char* dt = ctime(&now);
    string compart_reader = string(cb_env_combobox->GetValue().mb_str());
    for (map<string, bool>::iterator it = simulation->disp_compart_toggle.begin(); it != simulation->disp_compart_toggle.end(); it++){
    }


    for (map<string, bool>::iterator it = simulation->disp_compart_toggle.begin(); it != simulation->disp_compart_toggle.end(); it++){
        if(compart_reader == "off"){
            simulation->disp_compart_toggle[it->first] = false;
            
        }
        else
            if (compart_reader == it->first){
                simulation->disp_compart_toggle[it->first] = !simulation->disp_compart_toggle[it->first];//true;
                simulation->disp_compart_toggle["off"] = false;
                break;
            }
    }
    
    simCanvas->Refresh();
}
void SimSidebar::GetProteinSelection(wxCommandEvent& event){
    time_t now = time(0);
    char* dt = ctime(&now);
 
    for (int i = 0; i < clb_proteins_position; i++){
        wxString selected_protein_string = clb_proteins->GetString(i);
        string protein_selection = string(selected_protein_string.mb_str());
        simulation->disp_protein_toggle[protein_selection] = clb_proteins->IsChecked(i);
           
    }
    simCanvas->Refresh();
    for (map<string, bool>::iterator it = simulation->disp_protein_toggle.begin(); it!= simulation->disp_protein_toggle.end(); it++){
        cout << it->first << " protein state: " << to_string(simulation->disp_protein_toggle[it->first]) << " " << dt;
    }
}

void SimSidebar::display_all_proteins(wxCommandEvent& event){
    for (int i = 0; i < clb_proteins_position; i++){
        wxString selected_protein_string = clb_proteins->GetString(i);
        string protein_selection = string(selected_protein_string.mb_str());
        clb_proteins->Check(i, true);
        simulation->disp_protein_toggle[protein_selection] = true;
    }
    simCanvas->Refresh();
}

void SimSidebar::display_no_proteins(wxCommandEvent& event){
    for (int i = 0; i < clb_proteins_position; i++){
        wxString selected_protein_string = clb_proteins->GetString(i);
        string protein_selection = string(selected_protein_string.mb_str());
        clb_proteins->Check(i, false);
        simulation->disp_protein_toggle[protein_selection] = false;
        
    }
    simCanvas->Refresh();
}

void SimSidebar::GetMetaboliteSelection(wxCommandEvent& event){
    time_t now = time(0);
    char* dt = ctime(&now);
    
    for (int i = 0; i < clb_metabolites_position; i++){
        wxString selected_metabolite_string = clb_metabolites->GetString(i);
        string metabolite_selection = string(selected_metabolite_string.mb_str());
        simulation->disp_metabolite_toggle[metabolite_selection] = clb_metabolites->IsChecked(i);
           
    }

    for (map<string, bool>::iterator it = simulation->disp_metabolite_toggle.begin(); it != simulation->disp_metabolite_toggle.end(); it++){
        cout << it->first << " metabolite state: " << to_string(simulation->disp_metabolite_toggle[it->first]) << endl;
    }
    simCanvas->Refresh();
}

void SimSidebar::display_all_metabolites(wxCommandEvent& event){

    for (int i = 0; i < clb_metabolites_position; i++){
        wxString selected_metabolite_string = clb_metabolites->GetString(i);
        string metabolite_selection = string(selected_metabolite_string.mb_str());
        clb_metabolites->Check(i, true);
        simulation->disp_metabolite_toggle[metabolite_selection] = true;
    }
    simCanvas->Refresh();
}

void SimSidebar::display_no_metabolites(wxCommandEvent& event){

    for (int i = 0; i < clb_metabolites_position; i++){
        wxString selected_metabolite_string = clb_metabolites->GetString(i);
        string metabolite_selection = string(selected_metabolite_string.mb_str());
        clb_metabolites->Check(i, false);
        simulation->disp_metabolite_toggle[metabolite_selection] = false;
    }
    simCanvas->Refresh();
}

void SimSidebar::GetCompartSelection(wxCommandEvent& event){
    time_t now = time(0);
    char* dt = ctime(&now);
    for (int i = 0; i < clb_compart_position; i++){
        wxString selected_compart_string = clb_compart->GetString(i);
        string compart_selection = string(selected_compart_string.mb_str());
        simulation->disp_compart_toggle[compart_selection] = clb_compart->IsChecked(i);
           
    }
    simCanvas->Refresh();
    for (map<string, bool>::iterator it = simulation->disp_compart_toggle.begin(); it != simulation->disp_compart_toggle.end(); it++){
        cout << it->first << " compart state: " << to_string(simulation->disp_compart_toggle[it->first]) << endl;
    }
}

void SimSidebar::display_all_comparts(wxCommandEvent& event){
    for (int i = 0; i < clb_compart_position; i++){
        wxString selected_compart_string = clb_compart->GetString(i);
        string compart_selection = string(selected_compart_string.mb_str());
        clb_compart->Check(i, true);
        simulation->disp_compart_toggle[compart_selection] = true;
    }
    simCanvas->Refresh();
}

void SimSidebar::display_no_comparts(wxCommandEvent& event){
    for (int i = 0; i < clb_compart_position; i++){
        wxString selected_compart_string = clb_compart->GetString(i);
        string compart_selection = string(selected_compart_string.mb_str());
        clb_compart->Check(i, false);
        simulation->disp_compart_toggle[compart_selection] = false;
    }
    simCanvas->Refresh();
}

void SimSidebar::Rotate(wxCommandEvent& event){
    simCanvas->rotate_enabled = cb_rotate->IsChecked();
}
void SimSidebar::Translate(wxCommandEvent& event){}
void SimSidebar::Scale(wxCommandEvent& event){
    simCanvas->scale_enabled = cb_scale->IsChecked();
}
void SimSidebar::DisplayStats(wxCommandEvent& event){
    if (cb_displayStats->IsChecked()) {
        statsFrame = new wxFrame((wxFrame *)NULL, -1,  wxT("Statistics"), wxPoint(800, 100), STATS_SIZE);
        statsCanvas = new StatsCanvas(statsFrame, simCanvas->gl_args);
        statsFrame->Show();
    } else {
        statsFrame->Close();
        delete statsFrame;
    }
}

void SimSidebar::EndEarly(wxCommandEvent& event){
    simulation->continueSimulation = false;
    wxMessageDialog * exitEarly = new wxMessageDialog(this, "You are exiting the program early. Would you like to save the data to a file", "Exiting", wxCANCEL|wxYES_NO);
    int option = exitEarly->ShowModal();
    if (option == wxID_CANCEL){
        return;
    } else if (option == wxID_YES){
       string filename = sbml_document->getLocationURI() + "_results";
       string results = simulation->statsToString();
       ofstream f;
       f.open(filename, fstream::out | ofstream::trunc);
       f << results;
       f.close();

    } 
    exit(0);
}

void SimSidebar::Close(wxCommandEvent& event){
    exit(0);
}
