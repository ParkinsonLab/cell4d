//=========================================================
// wxHeader.h
// Donny Chan, Billy Taj


#include "_SimCanvas.h"
#include "_SimSidebar.h"
#include "_StatsCanvas.h"

#define MACRO_INIT(init_y) \
	int obj_id = 10000; \
	int y = init_y;

#define BUTTON(parent, text, x, y, width, height, event_handler) \
	new wxButton((parent), obj_id, (text), wxPoint((x), (y)), wxSize((width), (height))); \
    Connect(obj_id, wxEVT_BUTTON, wxCommandEventHandler(event_handler)); \
    obj_id++;

#define CHECKBOX_E(parent, text, event_handler) \
    new wxCheckBox((parent), obj_id, text, wxPoint(10, y)); \
    Connect(obj_id, wxEVT_CHECKBOX, wxCommandEventHandler(event_handler)); \
    obj_id++; y+=20;

#define CHECKBOX(parent, text) \
    new wxCheckBox((parent), obj_id, text, wxPoint(10, y)); \
    obj_id++; y+=20;

#define COMBOBOX(parent, text) \
    new wxComboBox((parent), obj_id, text, wxPoint(10, y)); \
    y+=20;

#define CHECKLISTBOX(parent) \
    new wxCheckListBox((parent), obj_id, wxPoint(10, y), wxSize(180, 70)); \
    y+=70;

#define BITMAPBUTTON(parent, bitmap, x, event_handler)\
    new wxBitmapButton((parent), obj_id, bitmap, wxPoint((x), (y))); \
    Connect(obj_id, wxEVT_BUTTON, wxCommandEventHandler(event_handler));\
    obj_id++;
    

#define STATICLINE() \
    new wxStaticLine(this, -1, wxPoint(5, y+10), wxSize(180, 1)); y+=20;

#define STATICTEXT(text, indent) \
    new wxStaticText(this, -1, text, wxPoint(indent, y)); y+=20;

#ifdef MAC_OS   
    #define ABOUT_SIZE wxSize(800, 600)
    #define LICENSE_SIZE wxSize(550,700)
    #define MENU_SIZE wxSize(400, 350)
    #define SIM_SIZE wxSize(900, 800)
    #define STATS_SIZE wxSize(600, 570)   
#else    
    #define ABOUT_SIZE wxSize(1024, 800)
    #define LICENSE_SIZE wxSize(550,700)
    #define MENU_SIZE wxSize(550, 550)
    #define SIM_SIZE wxSize(900, 800)
    #define STATS_SIZE wxSize(600, 540) 
#endif