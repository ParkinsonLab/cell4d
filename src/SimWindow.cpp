//===========================================================
// SimWindow.cpp
// Donny Chan, Billy Taj

#include "../inc/SimWindow.h"

int winMainWidth = 700;
int winMainHeight = 500;
int winControlWidth = 200;
int winControlHeight = 500;
int winDisplayWidth = 500;
int winDisplayHeight = 500;
int winStatWidth = 500;
int winStatHeight = 200;
int incMainWidth = 0;
int incMainHeight = 0;
int actionType = 0;
bool showStats = false;
bool showEnzyme = false;
bool wasTop = true;
bool wasBottom = false;


/*
===============================================================================
Static Declarations
===============================================================================
*/

extern void simulate();
extern SimDisplayWindow * sub;
extern StatsWindow * stats;
extern SimEngine * simulation;



/*
===============================================================================
Sim Window
===============================================================================
*/

SimWindow :: SimWindow():GlowWindow("Cell4D Diffusion Simulation", 
        GlowWindow::autoPosition, GlowWindow::autoPosition, 
        winMainWidth, winMainHeight, Glow::rgbaBuffer | Glow::doubleBuffer, Glow::noEvents){

    top = this;

    displayWin = new SimDisplayWindow(top);

    statsWin = new StatsWindow();

    controls = new GlowWidgetSubwindow(
            top, // Parent
            winDisplayWidth, 0, // Starting Coord
            winControlWidth, winControlHeight); // Dimension

    //simReceiver = new SimUIReceiver();

    stats = statsWin;
    sub = displayWin;
}

SimWindow :: ~SimWindow(){
    delete top;
    delete displayWin;
    delete statsWin;
    delete controls;
    delete enzyme;
//    delete simReceiver;
}

void 
SimWindow :: setup(){

    if(!showStats)    stats->Hide();

    enzyme = new EnzymeWindow();
    if(!showEnzyme) enzyme->Hide();

    ///////////////////////////////////////////////////////////////////////

    GlowLabelParams lparams; //

    // Label
    lparams.x = 10;
    lparams.y = 10;
    lparams.width = 300;
    lparams.text = "Controls";
    new GlowLabelWidget(controls, lparams);

    // Checkboxes

    GlowCheckBoxParams cbContinueSim; //

    // Run Simulation
    cbContinueSim.x = 10;
    cbContinueSim.y = 30;
    cbContinueSim.height = 15;
    cbContinueSim.width = 160;
    cbContinueSim.refcon = CONTINUE;
    cbContinueSim.text = "Run Simulation";
    cbContinueSim.state = GlowCheckBoxWidget::off;
    cbContinueSim.receiver = this;
    cb_continueSimulation = new GlowCheckBoxWidget(controls, cbContinueSim);

    GlowCheckBoxParams cbShowMetaConc; //

    // Show Metabolic Concentrations
    cbShowMetaConc.x = 10;
    cbShowMetaConc.y = 50;
    cbShowMetaConc.height = 15;
    cbShowMetaConc.width = 160;
    cbShowMetaConc.refcon = SHOWMETACONC;
    cbShowMetaConc.text = "Display Concentrations";
    cbShowMetaConc.state = GlowCheckBoxWidget::on; //for the checkmark to appear.
    cbShowMetaConc.receiver = this;
    cb_showMetaConc = new GlowCheckBoxWidget(controls, cbShowMetaConc);

    GlowCheckBoxParams cbShowMetaEnz; //

    // Show Metabolic Enzymes
    cbShowMetaEnz.x = 10;
    cbShowMetaEnz.y = 70;
    cbShowMetaEnz.height = 15;
    cbShowMetaEnz.width = 160;
    cbShowMetaEnz.refcon = SHOWMETAENZ;
    cbShowMetaEnz.text = "Display Proteins";
    cbShowMetaEnz.state = GlowCheckBoxWidget::on;
    cbShowMetaEnz.receiver = this;
    cb_showMetaEnz = new GlowCheckBoxWidget(controls, cbShowMetaEnz);

    GlowCheckBoxParams cbShowEnv; //

    // Show Environment
    cbShowEnv.x = 10;
    cbShowEnv.y = 90;
    cbShowEnv.height = 15;
    cbShowEnv.width = 160;
    cbShowEnv.refcon = SHOWENV;
    cbShowEnv.text = "Display Environment";
    cbShowEnv.state = GlowCheckBoxWidget::off;
    cbShowEnv.receiver = this;
    cb_showEnv = new GlowCheckBoxWidget(controls, cbShowEnv);

    // Show Enzymes
    GlowCheckBoxParams cbShowEnzyme; //

    cbShowEnzyme.x = 10;
    cbShowEnzyme.y = 110;
    cbShowEnzyme.height = 15;
    cbShowEnzyme.width = 160;
    cbShowEnzyme.refcon = SHOWENZYME;
    cbShowEnzyme.text = "Display Enzymes";
    cbShowEnzyme.state = GlowCheckBoxWidget::off;
    cbShowEnzyme.receiver = this;
    cb_showEnzyme = new GlowCheckBoxWidget(controls, cbShowEnzyme);


    // Separator

    GlowSeparatorParams spParams;

    spParams.x = 5;
    spParams.y = 130;
    spParams.width = controls->Width() - 10;
    spSimulation = new GlowSeparatorWidget(controls, spParams);
    /////////////////////////////////////////////////////////

    // Metabolite Panel

    //GlowPanelParams pMetaParams;
    //pMetaParams.x = 5;
    //pMetaParams.y = 135;
    //pMetaParams.width = controls->Width() - 10;
    //pMetaParams.height = 80;
    //pMetaParams.style = GlowPanelWidget::etchedStyle;

    //GlowPanelWidget * pMeta = new GlowPanelWidget(controls, pMetaParams);
    // Display Metabolites (Radio Buttons)

    GlowRadioGroupParams showGrpParams; //

    showGrpParams.height = 60;
    showGrpParams.width = 200;
    showGrpParams.y = 150;
    showGrpParams.x = 10;
    //showGrpParams.y = 0;
    //showGrpParams.x = 10;
    showGrpParams.receiver = this;
    showRadioGroup = new GlowRadioGroupWidget(controls, showGrpParams);

    //GlowRadioButtonParams* showParams; 

    //GlowRadioButtonParams showParamsAll; // 

    showParamsAll.height = 15;
    showParamsAll.width = 100;
    showParamsAll.x = 0;
    showParamsAll.y = 0;
    showParamsAll.refcon = SHOWALL;
    showParamsAll.text = "Show All Metabolites";
    rb_showAll = new GlowRadioButtonWidget(showRadioGroup, showParamsAll);

    rb_showMeta = (GlowRadioButtonWidget**) calloc(numMetabolite, sizeof(GlowRadioButtonWidget*));
    showParams = new GlowRadioButtonParams [numMetabolite]; 

    int i = 1;
    for(map<int, Species_t*>::const_iterator mItr = mets.begin(); mItr!= mets.end(); mItr++)
    {
        showParams[i-1].x = 10;
        showParams[i-1].y = i * 20;
        showParams[i-1].refcon = i + 20;
        showParams[i-1].text = Species_getName(mItr->second);
        rb_showMeta[i-1] = new GlowRadioButtonWidget(showRadioGroup, showParams[i-1]);
        if(showParams[i-1].y > 60)
            rb_showMeta[i-1]->Hide();
        i++;
    }

    GlowScrollBarParams sbParams;

    sbParams.width = 20;
    sbParams.height = 90;
    sbParams.x = controls->Width() - 25;
    sbParams.y = 150;
    sbParams.span = 1;
    sbParams.max = numMetabolite / 2;
    sbParams.refcon = SBMETA;

    
    sbMeta = new GlowScrollBarWidget(controls, sbParams);

    sbMeta->Notifier().Bind(this);

    cout << "first delay: " << sbMeta->GetFirstDelay() << endl;
    cout << "second delay: " << sbMeta->GetSecondDelay() << endl;
    cout << "arrow step: " << sbMeta->GetArrowStep() << endl;
    cout << "page step: " << sbMeta->GetPageStep() << endl;

    sbMeta->SetFirstDelay(1000);
    sbMeta->SetSecondDelay(1000);

    // Separator metabolite to environment

    spParams.x = 5;
    spParams.y = 245;
    spParams.width = controls->Width() - 10;
    spMeta = new GlowSeparatorWidget(controls, spParams);
    ///////////////////////////////////////////////////

    GlowCheckBoxParams cbShowParams; //

    cbShowParams.x = 10;
    cbShowParams.y = 250;

    cbShowParams.height = 15;
    cbShowParams.width = 100;
    cbShowParams.x += 0;
    cbShowParams.y += 0;

    GlowCheckBoxParams cbShowCytosolParams; //

    cbShowCytosolParams.x = cbShowParams.x;
    cbShowCytosolParams.y = cbShowParams.y;

    cbShowCytosolParams.height = cbShowParams.height;
    cbShowCytosolParams.width = cbShowParams.width;
    cbShowCytosolParams.x += 0;
    cbShowCytosolParams.y += 10;
    cbShowCytosolParams.refcon = SHOWCYTOSOL;
    cbShowCytosolParams.text = "Cytosol";
    cbShowCytosolParams.state = GlowCheckBoxWidget::on;
    cbShowCytosolParams.receiver = this;
    cb_showCytosol = new GlowCheckBoxWidget(controls, cbShowCytosolParams);

    GlowCheckBoxParams cbShowMembraneParams; //

    cbShowMembraneParams.x = cbShowCytosolParams.x;
    cbShowMembraneParams.y = cbShowCytosolParams.y;

    cbShowMembraneParams.height = cbShowCytosolParams.height;
    cbShowMembraneParams.width = cbShowCytosolParams.width;
    cbShowMembraneParams.x += 0;
    cbShowMembraneParams.y += 20;
    cbShowMembraneParams.refcon = SHOWENV1;
    cbShowMembraneParams.text = "Membrane";
    cbShowMembraneParams.state = GlowCheckBoxWidget::on;
    cbShowMembraneParams.receiver = this;
    cb_showEnv1 = new GlowCheckBoxWidget(controls, cbShowMembraneParams);

    GlowCheckBoxParams cbShowEnv2Params; //

    cbShowEnv2Params.x = cbShowMembraneParams.x;
    cbShowEnv2Params.y = cbShowMembraneParams.y;

    cbShowEnv2Params.height = cbShowMembraneParams.height;
    cbShowEnv2Params.width = cbShowMembraneParams.width;
    cbShowEnv2Params.x += 0;
    cbShowEnv2Params.y += 20;
    cbShowEnv2Params.refcon = SHOWENV2;
    cbShowEnv2Params.text = "Env 2";
    cbShowEnv2Params.state = GlowCheckBoxWidget::on;
    cbShowEnv2Params.receiver = this;
    cb_showEnv2 = new GlowCheckBoxWidget(controls, cbShowEnv2Params);

    GlowCheckBoxParams cbShowEnv3Params; //

    cbShowEnv3Params.x = cbShowEnv2Params.x;
    cbShowEnv3Params.y = cbShowEnv2Params.y;

    cbShowEnv3Params.height = cbShowEnv2Params.height;
    cbShowEnv3Params.width = cbShowEnv2Params.width;
    cbShowEnv3Params.x += 0;
    cbShowEnv3Params.y += 20;
    cbShowEnv3Params.refcon = SHOWENV3;
    cbShowEnv3Params.text = "Env 3";
    cbShowEnv3Params.state = GlowCheckBoxWidget::on;
    cbShowEnv3Params.receiver = this;
    cb_showEnv3 = new GlowCheckBoxWidget(controls, cbShowEnv3Params);

    GlowCheckBoxParams cbShowEnv4Params; //

    cbShowEnv4Params.x = cbShowEnv3Params.x;
    cbShowEnv4Params.y = cbShowEnv3Params.y;

    cbShowEnv4Params.height = cbShowEnv3Params.height;
    cbShowEnv4Params.width = cbShowEnv3Params.width;
    cbShowEnv4Params.x += 0;
    cbShowEnv4Params.y += 20;
    cbShowEnv4Params.refcon = SHOWENV4;
    cbShowEnv4Params.text = "Env 4";
    cbShowEnv4Params.state = GlowCheckBoxWidget::on;
    cbShowEnv4Params.receiver = this;
    cb_showEnv4 = new GlowCheckBoxWidget(controls, cbShowEnv4Params);

    sbParams.width = 20;
    sbParams.height = 90;
    sbParams.x = controls->Width() - 25;
    sbParams.y = 260;
    sbParams.span = 1;
    sbParams.max = 4;
    sbParams.refcon = SBENV;

    sbEnv = new GlowScrollBarWidget(controls, sbParams);
    sbEnv->Notifier().Bind(this);

    spParams.x = 5;
    spParams.y = 365;
    spParams.width = controls->Width() - 10;

    spEnv = new GlowSeparatorWidget(controls, spParams);
    //////////////////////////////////////////////////////

    GlowRadioGroupParams viewGrpParams; //

    viewGrpParams.height = 100;
    viewGrpParams.width = 100;
    viewGrpParams.y = spEnv->PositionY() + 10;
    viewGrpParams.x = 10;
    viewGrpParams.receiver = this;
    GlowRadioGroupWidget* viewRadioGroup = new GlowRadioGroupWidget(controls, viewGrpParams);

    GlowRadioButtonParams viewRotateParams; //

    viewRotateParams.height = 15;
    viewRotateParams.width = 100;
    viewRotateParams.x = 0;
    viewRotateParams.y = 0;
    viewRotateParams.refcon = ROTATE;
    viewRotateParams.text = "Rotate";
    rb_rotate = new GlowRadioButtonWidget(viewRadioGroup, viewRotateParams);



    GlowRadioButtonParams viewZoomParams; //

    viewZoomParams.height = 15;
    viewZoomParams.width = 100;
    viewZoomParams.x = 0;
    viewZoomParams.y = 40;
    viewZoomParams.refcon = ZOOM;
    viewZoomParams.text = "Zoom";
    rb_scale = new GlowRadioButtonWidget(viewRadioGroup, viewZoomParams);

    /////////////////////////////////////////////////////////

   
    GlowPushButtonParams pbClose; //
    // Quit Button
    pbClose.x = 10;
    pbClose.y = 470;
    pbClose.width = 180;
    pbClose.text = "Close";
    pbClose.refcon = CLOSE;
    pbClose.receiver = this;
    pb_quit = new GlowPushButtonWidget(controls, pbClose);

}

/*
===============================================================================
SimulationViewer
===============================================================================
*/

inline SimulationViewer::SimulationViewer(GlowComponent* parent) : GlowComponent(parent){}

void SimulationViewer::OnEndPaint(){
    simulation->display();
}

/*
===============================================================================
StatsWindow
===============================================================================
*/

inline StatsWindow::StatsWindow():
    GlowWindow("Statistics", GlowWindow::autoPosition, GlowWindow::autoPosition, 500, 200,
    Glow::rgbBuffer | Glow::depthBuffer | Glow::doubleBuffer,
    Glow::noEvents)
    {OnEndPaint();}

void StatsWindow::OnEndPaint(){
    if(showStats)
        simulation->displayGraph();
}

/*
===============================================================================
EnzymeWindow
===============================================================================
*/

inline EnzymeWindow::EnzymeWindow():
    GlowWindow("Enzymes", GlowWindow::autoPosition, GlowWindow::autoPosition, 300, 500,
            Glow::rgbBuffer | Glow::depthBuffer | Glow::doubleBuffer,
            Glow::noEvents)
{
    enzymeTop = this;
    displayEnzyme();}

void
EnzymeWindow::displayEnzyme()
{
    subEnzyme = new GlowWidgetSubwindow(
            enzymeTop, 0, 0, 300, 500);

    GlowCheckBoxParams cbEnzymeParams;

    cbEnzyme = (GlowCheckBoxWidget**) calloc (prot.size(), sizeof(GlowCheckBoxWidget*));
    int i=0;
    for(map<int, Species_t*>::const_iterator pItr = prot.begin(); pItr!= prot.end(); pItr++)
    {
        cbEnzymeParams.x = 10;
        cbEnzymeParams.y = 10 + 20 * i;
        cbEnzymeParams.width = 200;
        cbEnzymeParams.height = 15;
        cbEnzymeParams.text = Species_getName(pItr->second);
        cbEnzymeParams.state = GlowCheckBoxWidget::on;
        cbEnzyme[i] = new GlowCheckBoxWidget(subEnzyme, cbEnzymeParams);
        i++;
    }
  
}

/*
===============================================================================
SimDisplayWindow
===============================================================================
*/

inline SimDisplayWindow::SimDisplayWindow(GlowWindow* parent) :
    GlowSubwindow(parent, 0, 0, winDisplayWidth, winDisplayHeight,
    Glow::rgbBuffer | Glow::depthBuffer | Glow::doubleBuffer,
    Glow::mouseEvents | Glow::dragEvents)
{
    _manip = new GlowViewManipulator(this,
        GlowViewManipulatorParams::defaults);
    _manip->SetSpinnable(true);
    new SimulationViewer(_manip);
    glEnable(GL_DEPTH_TEST);
}

bool SimDisplayWindow::OnBeginPaint(){

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-0.5, 0.5, -0.5, 0.5, 1.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 3.0, -25.0);

    float maxDim = max( xDim, max( yDim, zDim ));

    gluLookAt(4 * maxDim/10.0, 5 * maxDim/10.0, 6 * maxDim/10.0, 0,1,0,0,1,0);
    return true;
}

void SimDisplayWindow::OnMouseDown(Glow::MouseButton button, int x, int y, Glow::Modifiers modifiers){

    if (!_manip->IsDragging()){
        float xn, yn;
        NormalizeCoordinates(x, y, xn, yn);
        if(actionType == ROTATE)
            _manip->BeginRotationDrag(xn, yn);
        else if(actionType == TRANSLATE)
            _manip->BeginTranslationDrag(xn, yn);
        else if(actionType == ZOOM)
            _manip->BeginScaleDrag(xn, yn);
    }
}

void SimDisplayWindow::OnMouseDrag(int x, int y){

    if (_manip->IsDragging()){
        float xn, yn;
        NormalizeCoordinates(x, y, xn, yn);
        _manip->InDrag(xn, yn);
    }
}

void SimDisplayWindow::OnMouseUp(Glow::MouseButton button, int x, int y, Glow::Modifiers modifiers){

    if (_manip->IsDragging()){
        float xn, yn;
        NormalizeCoordinates(x, y, xn, yn);
        _manip->EndDrag(xn, yn);
    }
}

/*
===============================================================================
SimUIReceiver
===============================================================================
*/

//void SimUIReceiver::OnMessage(const GlowPushButtonMessage& message){
void SimWindow::OnMessage(const GlowPushButtonMessage& message){

    int evt = message.widget->GetRefCon();
    if(evt == CLOSE){
        exit(0);
    }
}

void SimWindow::OnMessage(const GlowCheckBoxMessage& message){

    int evt = message.widget->GetRefCon();
    bool checked = (message.state == GlowCheckBoxWidget::on);

    if(evt == CONTINUE){

        simulation->continueSimulation = checked;

        if(simulation->continueSimulation)
            Glow::SetIdleFunc(simulate);
        else
            Glow::SetIdleFunc(NULL);
    }
    else if(evt == SHOWMETACONC){
        simulation->showMetaConc = checked;
    }
    else if(evt == SHOWMETAENZ){
        simulation->showMetaEnz = checked;
    }
    else if(evt == SHOWENV){
        simulation->showEnv = checked;
    }
    else if(evt == SHOWSTATS){
        showStats = checked;
        if(checked)
            stats->Show();
        else
            stats->Hide();
        stats->Refresh();
    }
    else if(evt == SHOWCYTOSOL){
        simulation->showEnvCyt = checked;
    }
    else if(evt == SHOWENV1){
        simulation->showEnv1 = checked;
    }
    else if(evt == SHOWENV2){
        simulation->showEnv2 = checked;
    }
    else if(evt == SHOWENV3){
        simulation->showEnv3 = checked;
    }
    else if(evt == SHOWENV4){
        simulation->showEnv4 = checked;
    }
    else if(evt == SHOWENZYME){
        showEnzyme = checked;
        if(checked)
            enzyme->Show();
        else
            enzyme->Hide();
        enzyme->Refresh();
    }

    message.widget->Refresh();
    sub->Refresh();
}

void SimWindow::OnMessage(const GlowRadioButtonMessage& message){
    int evt = message.buttonWidget->GetRefCon();

    if(evt == ROTATE){
        actionType = ROTATE;
    }
    else if(evt == TRANSLATE){
        actionType = TRANSLATE;
    }
    else if(evt == ZOOM){
        actionType = ZOOM;
    }
    else if(evt == SHOWALL){
        concType = -1;
    }
    
    else if(evt > SHOWENV4){    // larger than constant 20
        concType = evt-21;
    }
    sub->Refresh();
}


void SimWindow::OnMessage(const GlowScrollBarMessage & message)
{
    int evt = message.widget->GetRefCon();

    switch(evt)
    {
        case SBMETA:    // if scrollbar for metabolites is selected

            if(message.mouseButton == Glow::leftButton) // and left mouse button clicked
            {
                switch(message.part)
                {
                    case GlowScrollBarWidget::upButtonPart: // when up button is pressed

                    ///////////////////////////////////////////////////////////////////////////
                        // Check if the bar is already at the very top before the button was pressed
                        // if true then do nothing
                        if((sbMeta->GetTopValue() == sbMeta->GetMinimum()) && (wasTop == true))
                            break;
                    ///////////////////////////////////////////////////////////////////////////

                        for(int i=1; i<=numMetabolite; i++)
                        {
                            // Scroll down all radio buttons by 10 unit 
                            showParams[i-1].y += 20;

                            // if any radio button (metabolites) are not
                            // between y range of 20 to 60 then hide them
                            // or else show them
                            if(showParams[i-1].y > 60 || showParams[i-1].y < 20)
                            {
                                rb_showMeta[i-1]->Hide();
                                continue;
                            }
                            else
                                rb_showMeta[i-1]->Show();

                            // Move all radio buttons by 10 unit as specified
                            // above
                            rb_showMeta[i-1]->Move(showRadioGroup->PositionX(),showParams[i-1].y);

                            // Set flags for when the bar is at the very top or the bottom
                            // of the scrollbar
                            if(sbMeta->GetBottomValue() == sbMeta->GetMaximum())
                                wasBottom = true;
                            else if(sbMeta->GetTopValue() == sbMeta->GetMinimum())
                                wasTop = true;
                            else
                            {
                                wasBottom = false;
                                wasTop = false;
                            }

                            //cout << rb_showMeta[i-1]->PositionY() << endl;
                        }

                        break;

                    case GlowScrollBarWidget::downButtonPart:

                        // If the bar was at the very bottom of scrollbar before
                        // the down button is clicked then
                        // do nothing
                        if((sbMeta->GetBottomValue() == sbMeta->GetMaximum()) && (wasBottom == true))
                            break;

                        /////////////////////////////////////////////////////////////////////////////

                        for(int i=1; i<=numMetabolite; i++)
                        {
                            // Scroll up all radio buttons by 10 units
                            showParams[i-1].y -= 20;                  

                            // If any radio buttons are not between y coord
                            // of 20 to 60, then hide the buttons
                            // or else show the buttons
                            if(showParams[i-1].y > 60 || showParams[i-1].y < 20)
                            {
                                rb_showMeta[i-1]->Hide();
                                continue;
                            }
                            else
                                rb_showMeta[i-1]->Show();
                            // Move all radio buttons are specified above
                            rb_showMeta[i-1]->Move(showRadioGroup->PositionX(),showParams[i-1].y);

                            // Set flags for when the bar is at the very top or the bottom
                            // of the scrollbar
                            if(sbMeta->GetBottomValue() == sbMeta->GetMaximum())
                                wasBottom = true;
                            else if(sbMeta->GetTopValue() == sbMeta->GetMinimum())
                                wasTop = true;
                            else
                            {
                                wasBottom = false;
                                wasTop = false;

                            }

                            //cout << rb_showMeta[i-1]->PositionY() << endl;
                        }

                        break;

                    case GlowScrollBarWidget::upPagePart:

                        cout << "UpPagePart !!\n";

                        for(int i=1; i<=numMetabolite; i++)
                        {
                            showParams[i-1].y += 15;
                            if(showParams[i-1].y > 60 || showParams[i-1].y < 20)
                            {
                                rb_showMeta[i-1]->Hide();
                                continue;
                            }
                            else
                                rb_showMeta[i-1]->Show();
                            rb_showMeta[i-1]->Move(showRadioGroup->PositionX(),showParams[i-1].y);                            
                        }

                        break;

                    case GlowScrollBarWidget::downPagePart:

                        cout << "DownPagePart!!\n";

                        for(int i=1; i<=numMetabolite; i++)
                        {
                            showParams[i-1].y -= 15;
                            if(showParams[i-1].y > 60 || showParams[i-1].y < 20)
                            {
                                rb_showMeta[i-1]->Hide();
                                continue;
                            }
                            else
                                rb_showMeta[i-1]->Show();
                            rb_showMeta[i-1]->Move(showRadioGroup->PositionX(),showParams[i-1].y);
                        }

                        break;
                }
            }

            break;
        case SBENV:
            break;

    }
}

void SimWindow::OnReshape(int width, int height)
{

    incMainWidth = width - winMainWidth;
    incMainHeight = height - winMainHeight;
    winMainWidth = width;
    winMainHeight = height;

    winDisplayWidth = winMainWidth - winControlWidth;
    winControlHeight = winMainHeight;
    winDisplayHeight = winMainHeight;

    sub->Reshape(winDisplayWidth, winDisplayHeight);
    controls->Reshape(winControlWidth, winControlHeight);
    sub->Refresh();
    controls->Move(winMainWidth - winControlWidth, 0);
    sub->Refresh();
    controls->Refresh();
}
