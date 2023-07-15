//=======================================================
// SimSidebar.h
// Donny Chan, Billy Taj

#include <wx/wx.h>
#include <wx/statline.h>
#include <wx/spinctrl.h>
#include "_wxClass.h"
#include <sstream>
#include <fstream>

class SimSidebar : public wxPanel {
    wxCheckBox      *cb_continueSimulation;
    wxSlider        *sl_timeInterval;
    wxCheckBox      *cb_showMetaConc, *cb_showMetaEnz, *cb_showEnzyme; //
    wxCheckBox      *cb_showEnv;

    wxCheckBox      *cb_dispConc;

    wxRadioButton   **rb_showMetabolites;

    wxCheckBox      *cb_showCytosol, *cb_showEnv1, *cb_showEnv2, *cb_showEnv3, *cb_showEnv4;

    wxCheckBox      *cb_rotate, *cb_translate, *cb_scale;
    wxBitmapButton  *pbb_all_proteins, *pbb_no_proteins, *pbb_all_metabolites, *pbb_no_metabolites, *pbb_all_comparts, *pbb_no_comparts;
    wxButton        *pb_close, *pb_endEarly;
    //wxButton        *pb_play, *pb_pause, *pb_stop;
    wxComboBox      *cb_env_combobox;
    wxBitmapButton  *pbb_play, *pbb_pause, *pbb_stop;

    wxCheckListBox  *clb_compart, *clb_metabolites, *clb_proteins;

    wxSpinCtrl      *spc_time_interval;
    
    wxBitmap play_logo, pause_logo, stop_logo;
    wxBitmap all_logo, none_logo;

public:
    SimSidebar(wxFrame *parent);
    wxCheckBox      *cb_displayStats;

    // events
    void ContinueSimulation(wxCommandEvent& event);
    void SimulationSpeed(wxCommandEvent& event);
    void ShowMetaConc(wxCommandEvent& event);
    void ShowMetaEnz(wxCommandEvent& event);
    void ShowEnv(wxCommandEvent& event);
    void ShowConcType(wxCommandEvent& event);
    void Rotate(wxCommandEvent& event);
    void Translate(wxCommandEvent& event);
    void Scale(wxCommandEvent& event);
    void DisplayStats(wxCommandEvent& event);
    void EndEarly(wxCommandEvent& event);
    void GetComboboxSelection(wxCommandEvent& event);
    void GetCompartSelection(wxCommandEvent& event);
    void GetMetaboliteSelection(wxCommandEvent& event);
    void GetProteinSelection(wxCommandEvent& event);
    void Close(wxCommandEvent& event);
    void display_all_proteins(wxCommandEvent& event);
    void display_no_proteins(wxCommandEvent& event);
    void display_all_metabolites(wxCommandEvent& event);
    void display_no_metabolites(wxCommandEvent& event);
    void display_all_comparts(wxCommandEvent& event);
    void display_no_comparts(wxCommandEvent& event);
    void PauseSimulation(wxCommandEvent& event);
    void StopSimulation(wxCommandEvent& event);
    void spc_simulation_speed(wxCommandEvent &event);


};  
