//=====================================================
// CellApp.h
// Donny Chan, Billy Taj

// include OpenGL
#ifdef MAC_OS
    #include <OpenGL/glu.h>
    #include <OpenGL/gl.h>
    #include <GLUT/glut.h>
#else
    #if 1
    // #include <GL/glu.h>
    // #include <GL/gl.h>
    // #include <GL/glut.h>
    #endif
#endif

#include <algorithm>
#include <iterator>
#include "ParameterManager.h"
#include "SimEngine.h"
#include "../ext/argh.h"

#if 1
#include <wx/wx.h>
#include <wx/sizer.h>
#include <wx/glcanvas.h>
#include <wx/wxhtml.h>
#include <wx/event.h>
#include "_wxHeader.h"
#include <iostream>
#include <fstream>

class CellApp: public wxApp {
    wxFrame *menuFrame, *simFrame;
    wxFrame *aboutFrame = NULL;
    wxFrame *pFrame = NULL;
    wxFrame *gFrame = NULL;
    wxScrolledWindow *paramFrame = NULL;
    wxFrame *licenseFrame = NULL;
    wxString filename;
    wxTextCtrl *seedText, *maxCycleText, *timeScaleText, *spaceScaleText, *letticeDimText, *nucleusSizeText, *inaccessibleSpaceText;
    wxRadioButton *originalButton, *cubeButton, *oneLayerButton, *threeLayerButton, *sphereButton, *twoCubesButton;
    wxCheckBox *nucleusEnabled, *plasmaMembraneEnabled;
    wxTextCtrl ** speciesInput;

    wxBitmap main_logo;
    wxStaticBitmap *main_image;

    wxStaticText *heading;
    wxStaticText *subheading;
    wxStaticText *version_heading;

    wxStaticText * demo_label;
    wxStaticText * model_label;

    virtual bool OnInit();
    void ShowMenu();
    void ShowAbout();
    void ShowLicense();
    void ShowSimulation(const char*SBMLFile, wxString title);
    void ShowParamFrame(wxCommandEvent& event);
    void ShowGeometryFrame(wxCommandEvent& event);
    void LaunchSimulation(wxCommandEvent& event);
    void SaveGeometry(wxCommandEvent& event);

    // events
    void SimulationFromFile(wxCommandEvent& event);
    void DemoCascade(wxCommandEvent& event);
    void DemoCeacam(wxCommandEvent& event);
    void DemoABFusion(wxCommandEvent& event);
    void DemoGlycolysis(wxCommandEvent& event);
    void SimParam(wxCommandEvent& event);
    void AboutButton(wxCommandEvent& event);
    void AboutPageClosed(wxCloseEvent& evt);
    void LicenseButton(wxCommandEvent& event);
    void LicensePageClosed(wxCloseEvent& evt);
    void CloseSimulation(wxCommandEvent& event);
    void OpenEditor(wxCommandEvent& event);

};
#endif

argh::parser parameter_parsing(int argc, char * argv[]);
void usage(int usage_in);
