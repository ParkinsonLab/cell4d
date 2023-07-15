//=================================================================
// SimControls.h
// Donny Chan, Billy Taj

#ifndef SIMCONTROLS__H
#define SIMCONTROLS__H


/*
===============================================================================
    Headers
===============================================================================
*/

#include "glow.h"
#include "glowQuickPalette.h"

GLOW_NAMESPACE_USING

#include "BallWidget.h"


/*
===============================================================================
    Controls class
===============================================================================
*/

// SimControls is an object that receives a number of events.
// The idle event receiver is used for background calculating.
// The TextFieldWindow event receiver is used for the save image dialog.
// The other event receivers are widgets in the control palette.

class SimControls :
    public GlowPushButtonReceiver,
    public GlowSliderReceiver,
    public GlowCheckBoxReceiver,
    public BallReceiver
{
    public:

        SimControls(
            int mainWindowID);
        virtual ~SimControls();

    protected:

        // Widget events

        // Implement this method to receive pushbutton events
        virtual void OnMessage(
            const GlowPushButtonMessage& message);

        // Implement this method to receive slider events
        virtual void OnMessage(
            const GlowSliderMessage& message);

        // Implement this method to receive check box events
        virtual void OnMessage(
            const GlowCheckBoxMessage& message);

        // Implement this method to receive ball events (new with lesson 8)
        virtual void OnMessage(
            const BallMessage& message);

    private:

        // Main drawing window id
        int _mainWindowID;

        // Pointer to the control window itself
        GlowQuickPaletteWindow* _controlWindow;

        // Pointers to selected widgets. We store these so that we can tell
        // what widget was hit when we receive events.
        GlowPushButtonWidget* _quitButton;
        BallWidget* _viewBall;
        BallWidget* _lightBall;
};


#endif
