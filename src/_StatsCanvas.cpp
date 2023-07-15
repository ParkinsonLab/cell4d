//=================================================
// StatsCanvas
// Donny Chan, Billy Taj

#include "../inc/_wxHeader.h"

extern SimSidebar *simSidebar;
extern StatsCanvas *statsCanvas;

extern SimEngine *simulation;

StatsCanvas::StatsCanvas(wxFrame *parent, int *args) :
        wxGLCanvas(parent, wxID_ANY, args, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE) {
    m_context = new wxGLContext(this);
    SetBackgroundStyle(wxBG_STYLE_CUSTOM);
}

BEGIN_EVENT_TABLE(StatsCanvas, wxGLCanvas)
    EVT_PAINT(StatsCanvas::OnPaint)
    EVT_SIZE(StatsCanvas::OnSize)
END_EVENT_TABLE()

StatsCanvas::~StatsCanvas() {
    statsCanvas = NULL;
    simSidebar->cb_displayStats->SetValue(false);
    delete m_context;
}

void StatsCanvas::OnPaint(wxPaintEvent& evt) {
    wxGLCanvas::SetCurrent(*m_context);
    wxPaintDC(this); // only to be used in paint events. use wxClientDC to paint outside the paint event

    glViewport(0, 0, GetSize().x, GetSize().y);
    simulation->displayGraph();

    glFlush();
    SwapBuffers();
}
void StatsCanvas::OnSize(wxSizeEvent& evt) {
    wxGLCanvas::OnSize(evt);
    Refresh();
}
