//=================================================
// SimControls.cpp
// Donny Chan, Billy Taj


#include "../inc/SimControls.h"

extern int shoulder1, shoulder2, shoulder3, shoulder4, lat1, lat2,
  elbow1, elbow2, pivot, tilt, ankle1, ankle2, heel1,
  heel2, hip11, hip12, hip21, hip22, solid_part,
  turn, turn1, lightturn, lightturn1;


SimControls::SimControls(
    int mainWindowID)
{
    _mainWindowID = mainWindowID;

    // Create a control panel window
    _controlWindow = new GlowQuickPaletteWindow("Controls");

    // Add controls
    // First, we'll put a little blurb at the top of the window
    _controlWindow->AddLabel("Metabolite Diffusion Model");

    // The rest of the window appears in three panels. The panels are arranged
    // horizontally, but widgets are arranged vertically within each panel
    GlowQuickPanelWidget* hpanel = _controlWindow->AddArrangingPanel(
        GlowQuickPanelWidget::horizontal);

    // First we have a few general controls
    GlowQuickPanelWidget* panel = hpanel->AddPanel(
        GlowQuickPanelWidget::loweredStyle, "Main",
        GlowQuickPanelWidget::vertical);

    // Note how we use arranging panels to get the arrangement we want.
    GlowQuickPanelWidget* panel2 = panel->AddArrangingPanel(
        GlowQuickPanelWidget::horizontal);
    GlowQuickPanelWidget* panel3 = panel2->AddArrangingPanel(
        GlowQuickPanelWidget::vertical);
    _viewBall = new BallWidget(panel3);
    _viewBall->Notifier().Bind(this);
    panel3->AddLabel("View");
    panel3 = panel2->AddArrangingPanel(GlowQuickPanelWidget::vertical);
    _lightBall = new BallWidget(panel3);
    _lightBall->Notifier().Bind(this);
    panel3->AddLabel("Light");

    // A separator
    panel->AddSeparator();
    // Quit button
    _quitButton = panel->AddPushButton("Quit", this);

    // Arrange controls and show the control panel window
    _controlWindow->Pack();

}


// Destructor

SimControls::~SimControls()
{
    delete _controlWindow;
}


// Receive pushbutton events
void SimControls::OnMessage(
    const GlowPushButtonMessage& message)
{
//    GLOW_DEBUGSCOPE("SimControls::OnMessage(pushbutton)");

    // Was it the quit button?
    if (message.widget == _quitButton)
    {
        exit(0);
    }
}


// Receive checkbox events
void SimControls::OnMessage(
    const GlowCheckBoxMessage& message)
{

}


// Receive slider events
void SimControls::OnMessage(
    const GlowSliderMessage& message)
{

}


// Receive ball events (lesson 8)
void SimControls::OnMessage(
    const BallMessage& message)
{
    GLOW_DEBUGSCOPE("SimControls::OnMessage(ball)");

    // Project the quaternion onto the latitude/longitude specification
    // used by glutmech by applying the quaternion to a sample vector,
    // the forward view (0,0,1), and examining the result
    Vec3f vec = message.rotation * Vec3f(0, 0, 1);
    if (message.widget == _viewBall)
    {
        GLfloat val = vec.Y();
        if (val < -1) val = -1;
        if (val > 1) val = 1;
        GLfloat theta = asin(val);
        val = vec.X()/cos(theta);
        if (val < -1) val = -1;
        if (val > 1) val = 1;
        turn1 = -Math::radiansToDegrees * theta;
        if (vec.Z() > 0)
        {
            turn = Math::radiansToDegrees * asin(val);
        }
        else
        {
            turn = 180 - Math::radiansToDegrees * asin(val);
        }
    }
    else if (message.widget == _lightBall)
    {
        GLfloat val = vec.X();
        if (val < -1) val = -1;
        if (val > 1) val = 1;
        GLfloat theta = asin(val);
        val = vec.Y()/cos(theta);
        if (val < -1) val = -1;
        if (val > 1) val = 1;
        lightturn = Math::radiansToDegrees * theta;
        if (vec.Z() > 0)
        {
            lightturn1 = -Math::radiansToDegrees * asin(val);
        }
        else
        {
            lightturn1 = 180 + Math::radiansToDegrees * asin(val);
        }
    }

    // Redraw scene
    Glow::RefreshGlutWindow(_mainWindowID);
}

