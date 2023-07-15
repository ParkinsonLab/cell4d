//==========================================
// StatsCanvas.h
// Donny Chan, Billy Taj

// #ifndef _glpane_
// #define _glpane_

#include <wx/wx.h>
#include <wx/glcanvas.h>

class StatsCanvas : public wxGLCanvas {
    wxGLContext*    m_context;

public:
    StatsCanvas(wxFrame *parent, int *args);
    virtual ~StatsCanvas();

    // events
    void            OnPaint(wxPaintEvent& evt);
    void            OnSize(wxSizeEvent& evt);

    DECLARE_EVENT_TABLE()
};
