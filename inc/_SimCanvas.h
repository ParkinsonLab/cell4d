//====================================================
// SimCanvas.h
// Donny Chan, Billy Taj

#include <wx/wx.h>
#include <wx/glcanvas.h>

#include <chrono>

#include "ParameterManager.h"
#include "SimEngine.h"

#include "_wxClass.h"

typedef std::chrono::high_resolution_clock T;
typedef std::chrono::milliseconds ms;

class SimCanvas : public wxGLCanvas {
    wxGLContext*    m_context;
    wxCoord         mouseX, mouseY;
    double          r, theta, phi;
    T::time_point   t0;

public:
    SimCanvas(wxFrame *parent, int *args, const char *SBMLFile);
    virtual ~SimCanvas();

    bool            rotate_enabled, scale_enabled;
    int             simulation_timeout;
    int             *gl_args;

    // events
    void            OnPaint(wxPaintEvent& evt);
    void            OnSize(wxSizeEvent& evt);
    void            OnMouse(wxMouseEvent& event);
    void            OnIdle(wxIdleEvent& event);
    DECLARE_EVENT_TABLE()
};
