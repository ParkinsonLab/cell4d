//==========================================================
// SimWindow.h
// Donny Chan, Billy Taj

#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include "glow.h"
#include "glowLabelWidget.h"
#include "glowPushButtonWidget.h"
#include "glowRadioButtonWidget.h"
#include "glowCheckBoxWidget.h"
#include "glowTextFieldWindow.h"
#include "glowViewTransform.h"
#include "glowScrollBarWidget.h"
#include "glowTextFieldWidget.h"
#include "glowSeparatorWidget.h"
#include <GL/glut.h>
#include <GL/gl.h>

#include "SimEngine.h"

GLOW_NAMESPACE_USING

#define CLOSE -1
#define ROTATE 0
#define TRANSLATE 1
#define ZOOM 2
#define CONTINUE 3

#define SHOWMETACONC 4
#define SHOWMETAENZ 5

#define SAVE 6

#define SHOWALL 7
//#define SHOWMETAA 8
//#define SHOWMETAB 9
//#define SHOWMETAC 10
//#define SHOWMETAD 11
//#define SHOWMETAE 12


#define SCROLL 8

#define SBMETA 9
#define SBENV 10
#define SBVIEW 11
#define SBSTAT 12

#define SHOWSTATS 13
#define SHOWENV 14

#define SHOWALLENV 15
#define SHOWCYTOSOL    16
#define SHOWENV1 17
#define SHOWENV2 18
#define SHOWENV3 19
#define SHOWENV4 20
#define SHOWENZYME 100

/*
===============================================================================
SimulationViewer
===============================================================================
*/

class SimulationViewer : public GlowComponent{
    public:

        inline SimulationViewer(GlowComponent* parent);
        virtual void OnEndPaint();
};

/*
===============================================================================
StatsWindow
===============================================================================
*/

class StatsWindow : public GlowWindow{

    public:
        inline StatsWindow();

    private:
        void OnEndPaint();
};

/*
===============================================================================
EnzymeWindow
===============================================================================
*/

class EnzymeWindow : public GlowWindow{

    public:
        inline EnzymeWindow();

        GlowCheckBoxWidget ** cbEnzyme;
    protected:

        void displayEnzyme();

    private:

        GlowWindow * enzymeTop;
        GlowWidgetSubwindow * subEnzyme;

};

/*
===============================================================================
SimDisplayWindow
===============================================================================
*/

class SimDisplayWindow : public GlowSubwindow{

    public:

        inline SimDisplayWindow(GlowWindow* parent);

        virtual bool OnBeginPaint();
        virtual void OnMouseDown(Glow::MouseButton button, int x, int y, Glow::Modifiers modifiers);
        virtual void OnMouseDrag(int x, int y);
        virtual void OnMouseUp(Glow::MouseButton button, int x, int y, Glow::Modifiers modifiers);

    private:

        GlowViewManipulator* _manip;
};

/*
===============================================================================
SimUIReceiver
===============================================================================
*/
/*
class SimUIReceiver :
    public GlowPushButtonReceiver,
    public GlowRadioButtonReceiver,
    public GlowCheckBoxReceiver
{

    public:

        virtual void OnMessage(const GlowPushButtonMessage & message);
        virtual void OnMessage(const GlowCheckBoxMessage & message);
        virtual void OnMessage(const GlowRadioButtonMessage & message);

    private:
        StatsWindow * statsWin;
        SimDisplayWindow * displayWin;
};
*/
/*
===============================================================================
SimWindow
===============================================================================
*/

class SimWindow :
    public GlowWindow,
    public GlowPushButtonReceiver,
    public GlowCheckBoxReceiver,
    public GlowRadioButtonReceiver,
    public GlowScrollBarReceiver
{
    public:
        SimWindow();
        ~SimWindow();
        void setup();
        GlowCheckBoxWidget * cb_continueSimulation;
        //GlowScrollBarWidget * sbMeta;
        //GlowScrollBarWidget * sbEnv;

    protected:

        //virtual void OnEndPaint();
        virtual void OnMessage(const GlowPushButtonMessage & message);
        virtual void OnMessage(const GlowCheckBoxMessage & message);
        virtual void OnMessage(const GlowRadioButtonMessage & message);
        virtual void OnMessage(const GlowScrollBarMessage & message);
        virtual void OnReshape(int width, int height);

    private:
        GlowWindow * top;
        StatsWindow * statsWin;
        SimDisplayWindow * displayWin;
        GlowWidgetSubwindow * controls;
        GlowWidgetSubwindow * swMeta;
        GlowWindow * enzyme;
        GlowRadioGroupWidget * showRadioGroup;

        GlowScrollBarWidget * sbMeta;
        GlowScrollBarWidget * sbEnv;

        GlowSeparatorWidget * spSimulation;
        GlowSeparatorWidget * spMeta;
        GlowSeparatorWidget * spEnv;
        //GlowLabelParams lparams;
        //GlowScrollBarParams sbParams;
        //GlowCheckBoxParams cbContinueSim;
        //GlowCheckBoxParams cbShowMetaConc;
        //GlowCheckBoxParams cbShowMetaEnz;
        //GlowCheckBoxParams cbShowEnv;

        //GlowRadioGroupParams viewGrpParams;
        //GlowRadioButtonParams viewRotateParams;
        //GlowRadioButtonParams viewTranslateParams;
        //GlowRadioButtonParams viewZoomParams;

        //GlowRadioGroupParams showGrpParams;
        GlowRadioButtonParams showParamsAll;
        GlowRadioButtonParams *showParams;

       // GlowCheckBoxParams cbShowParams;
        //GlowCheckBoxParams cbShowCytosolParams;
        //GlowCheckBoxParams cbShowMembraneParams;
        //GlowCheckBoxParams cbShowEnv2Params;
        //GlowCheckBoxParams cbShowEnv3Params;
        //GlowCheckBoxParams cbShowEnv4Params;

        //GlowCheckBoxParams cbShowStats;
        //GlowPushButtonParams pbClose;

        GlowCheckBoxWidget * cb_showMetaConc;
        GlowCheckBoxWidget * cb_showMetaEnz;
        GlowCheckBoxWidget * cb_showEnv;
        GlowCheckBoxWidget * cb_showEnzyme;

        GlowCheckBoxWidget * cb_displayStats;

        GlowRadioButtonWidget * rb_showAll;
        GlowRadioButtonWidget ** rb_showMeta;

        GlowCheckBoxWidget * cb_showCytosol;
        GlowCheckBoxWidget * cb_showEnv1;
        GlowCheckBoxWidget * cb_showEnv2;
        GlowCheckBoxWidget * cb_showEnv3;
        GlowCheckBoxWidget * cb_showEnv4;

        GlowRadioButtonWidget * rb_rotate;
        GlowRadioButtonWidget * rb_translate;
        GlowRadioButtonWidget * rb_scale;

        GlowPushButtonWidget * pb_quit;

};

#endif
