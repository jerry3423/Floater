/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Window
//  -----------------------------
//  
//  The render window class is really just a container
//  for tracking the visual hierarchy.  While it is
//  entirely possible to use Qt Creator, I try to avoid
//  over-commitment to it because I need to write code in
//  multiple environments, some of which are not well-suited
//  to IDEs in general, let alone Qt Creator
//
//  Also, by providing sample code, the didactic purpose of 
//  showing how things fit together is better served.
//
/////////////////////////////////////////////////////////////////

// include guard
#ifndef _RENDER_WINDOW_H
#define _RENDER_WINDOW_H

// include standard Qt widgets
#include <QtWidgets>

// include a custom arcball widget 
#include "ArcBallWidget.h"
// include the custom render widget
#include "RenderWidget.h"
// and the RGB Object class
#include "AttributedObject.h"

// a window that displays an geometric model with controls
class RenderWindow : public QWidget
    { // class RenderWindow
    private:
    // the geometric object being shown
    AttributedObject       		*attributedObject;
    
    // the values set in the interface
    RenderParameters            *renderParameters;

    // window layout    
    QGridLayout                 *windowLayout;

    // custom widgets
    ArcBallWidget               *modelRotator;
    RenderWidget                *renderWidget;

    // check boxes for modelling options
    QCheckBox					*renderWireframeBox;
    QCheckBox                   *useNormalBox;
    QCheckBox                   *useTexCoordsBox;
    QCheckBox                   *renderTextureBox;
    QCheckBox                   *renderNormalMapBox;

    // sliders for spatial manipulation
    QSlider                     *xTranslateSlider;
    QSlider                     *yTranslateSlider;
    QSlider                     *zoomSlider;

    // labels for sliders & arcballs
    QLabel                      *modelRotatorLabel;
    QLabel                      *yTranslateLabel;
    QLabel                      *zoomLabel;

    public:
    // constructor
    RenderWindow
        (
        // the object to be rendered
        AttributedObject         		*newAttributedObject, 
        // the model object storing render parameters
        RenderParameters        *newRenderParameters,
        // the title for the window (with default value)
        const char              *windowName = "Object Renderer"
        );  
    
    // routine to reset the interface
    // sets every visual control to match the model
    // gets called by the controller after each change in the model
    void ResetInterface();

    // declare the render controller class a friend so it can access the UI elements
    friend class RenderController;

    }; // class RenderWindow

// end of include guard
#endif
