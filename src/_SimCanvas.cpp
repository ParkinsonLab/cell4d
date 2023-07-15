//========================================================
// SimCanvas.cpp
// Donny Chan, Billy Taj

//Simulation Canvas: controls simulation space. 

#include "../inc/_wxHeader.h"

extern ParameterManager *pm;

extern map<string, Species_Attributes *> species_details_map; // replaces metabolite map.  now we just check what type it is


extern SimSidebar *simSidebar;
extern StatsCanvas *statsCanvas;

SimEngine *simulation;

SimCanvas::SimCanvas(wxFrame *parent, int *args, const char *SBMLFile) :
    wxGLCanvas(parent, wxID_ANY, args, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE) {
    gl_args = args;
    r = 50, theta = .87, phi = 1.;
    rotate_enabled = true;
    scale_enabled = true;
    m_context = new wxGLContext(this);
    SetBackgroundStyle(wxBG_STYLE_CUSTOM);

    t0 = T::now();
    simulation_timeout = 1;

    simulation = new SimEngine();
}

SimCanvas::~SimCanvas() {
    delete m_context;
}

BEGIN_EVENT_TABLE(SimCanvas, wxGLCanvas)
    EVT_PAINT(SimCanvas::OnPaint)
    EVT_SIZE(SimCanvas::OnSize)
    EVT_MOUSE_EVENTS(SimCanvas::OnMouse)
    EVT_IDLE(SimCanvas::OnIdle)
END_EVENT_TABLE()

void SimCanvas::OnPaint(wxPaintEvent& evt) {
    wxGLCanvas::SetCurrent(*m_context);
    wxPaintDC(this); // only to be used in paint events. use wxClientDC to paint outside the paint event

    glViewport(0, 0, GetSize().x, GetSize().y);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float ratio_w_h = (float)(GetSize().x)/(float)(GetSize().y);
    gluPerspective(15 /*view angle*/, ratio_w_h, 1 /*clip close*/, 20000 /*clip far*/);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -25.0);

    float maxDim = max(xDim, max(yDim, zDim));

    double x, y, z;
    x = r*sin(phi)*sin(theta);
    y = r*sin(phi)*cos(theta);
    z = r*cos(phi);

    gluLookAt(x * maxDim/10.0 - 0.5, y * maxDim/10.0 - 0.5, z * maxDim/10.0 - 0.5, -0.5, -0.5, -0.5, 0, 0, 1);

    simulation->display();

    glFlush();
    SwapBuffers();
}
void SimCanvas::OnSize(wxSizeEvent& evt) {
    wxGLCanvas::OnSize(evt);
    Refresh();
}
void SimCanvas::OnMouse(wxMouseEvent& event) {
    if (event.Dragging() && rotate_enabled) {
        theta += (event.GetX() - mouseX)*5./GetSize().y;
        phi -= (event.GetY() - mouseY)*5./GetSize().y;
        phi = min(max(phi, 0.00001), 3.14159);
        Refresh();
    }
    if (event.GetWheelRotation() != 0 && scale_enabled){
        r -= event.GetWheelRotation()*1./event.GetWheelDelta();
        r = max(0.1, r);
        Refresh();
    }
    mouseX = event.GetX();
    mouseY = event.GetY();
}
void SimCanvas::OnIdle(wxIdleEvent& event) {
    // Refresh();
    if (simulation->continueSimulation && chrono::duration_cast<ms>(T::now() - t0).count() > simulation_timeout) {
        simulation->simulate();
        if (statsCanvas != NULL){
            statsCanvas->Refresh();
        }
        Refresh();
        t0 = T::now();
    }
    event.RequestMore();
}
