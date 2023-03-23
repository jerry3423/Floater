/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Controller
//  -----------------------------
//  
//  We're using the Model-View-Controller pattern
//  so most of the control logic goes here
//  which means we need a slot for substantially
//  every possible UI manipulation
//
/////////////////////////////////////////////////////////////////

#include "RenderController.h"
#include <stdio.h>

// constructor
RenderController::RenderController
        (
        // the geometric object to show
        AttributedObject	*newAttributedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // the render window that it controls
        RenderWindow        *newRenderWindow
        )
    :
    attributedObject(newAttributedObject),
    renderParameters(newRenderParameters),
    renderWindow    (newRenderWindow),
    dragButton      (Qt::NoButton)
    { // RenderController::RenderController()
    
    // connect up signals to slots

    // signals for arcballs
    QObject::connect(   renderWindow->modelRotator,                 SIGNAL(RotationChanged()),
                        this,                                       SLOT(objectRotationChanged()));

    // signals for main widget to control arcball
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(BeginScaledDrag(int, float, float)),
                        this,                                       SLOT(BeginScaledDrag(int, float, float)));
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(ContinueScaledDrag(float, float)),
                        this,                                       SLOT(ContinueScaledDrag(float, float)));
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(EndScaledDrag(float, float)),
                        this,                                       SLOT(EndScaledDrag(float, float)));

    // signal for zoom slider
    QObject::connect(   renderWindow->zoomSlider,                   SIGNAL(valueChanged(int)),
                        this,                                       SLOT(zoomChanged(int)));

    // signal for x translate sliders
    QObject::connect(   renderWindow->xTranslateSlider,             SIGNAL(valueChanged(int)),
                        this,                                       SLOT(xTranslateChanged(int)));

    // signal for y translate slider
    QObject::connect(   renderWindow->yTranslateSlider,             SIGNAL(valueChanged(int)),
                        this,                                       SLOT(yTranslateChanged(int)));

    // signal for check box for normals
    QObject::connect(   renderWindow->renderWireframeBox,           SIGNAL(stateChanged(int)),
                        this,                                       SLOT(useWireframeCheckChanged(int)));

    // signal for check box for normals
    QObject::connect(   renderWindow->useNormalBox,                 SIGNAL(stateChanged(int)),
                        this,                                       SLOT(useNormalCheckChanged(int)));

    // signal for check box for texture coordinates
    QObject::connect(   renderWindow->useTexCoordsBox,              SIGNAL(stateChanged(int)),
                        this,                                       SLOT(useTexCoordsCheckChanged(int)));

    // signal for check box for texture render
    QObject::connect(   renderWindow->renderTextureBox,             SIGNAL(stateChanged(int)),
                        this,                                       SLOT(renderTextureCheckChanged(int)));

    // signal for check box for normal map render
    QObject::connect(   renderWindow->renderNormalMapBox,           SIGNAL(stateChanged(int)),
                        this,                                       SLOT(renderNormalMapCheckChanged(int)));

    // copy the rotation matrix from the widgets to the model
    renderParameters->rotationMatrix = renderWindow->modelRotator->RotationMatrix();
    } // RenderController::RenderController()

// slot for responding to arcball rotation for object
void RenderController::objectRotationChanged()
    { // RenderController::objectRotationChanged()
    // copy the rotation matrix from the widget to the model
    renderParameters->rotationMatrix = renderWindow->modelRotator->RotationMatrix();
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::objectRotationChanged()

// slot for responding to zoom slider
void RenderController::zoomChanged(int value)
    { // RenderController::zoomChanged()
    // compute the new scale
    float newZoomScale = pow(10.0, (float) value / 100.0);

    // clamp it
    if (newZoomScale < ZOOM_SCALE_MIN)
        newZoomScale = ZOOM_SCALE_MIN;
    else if (newZoomScale > ZOOM_SCALE_MAX)
        newZoomScale = ZOOM_SCALE_MAX;

    // and reset the value  
    renderParameters->zoomScale = newZoomScale;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::zoomChanged()

// slot for responding to x translate sliders
void RenderController::xTranslateChanged(int value)
    { // RenderController::xTranslateChanged()
    // reset the model's x translation (slider ticks are 1/100 each)
    renderParameters->xTranslate = (float) value / 100.0;

    // clamp it
    if (renderParameters->xTranslate < TRANSLATE_MIN)
        renderParameters->xTranslate = TRANSLATE_MIN;
    else if (renderParameters->xTranslate > TRANSLATE_MAX)
        renderParameters->xTranslate = TRANSLATE_MAX;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::xTranslateChanged()

// slot for responding to y translate slider
void RenderController::yTranslateChanged(int value)
    { // RenderController::tTranslateChanged()
    // reset the model's y translation (slider ticks are 1/100 each)
    renderParameters->yTranslate =  (float) value / 100.0;

    // clamp it
    if (renderParameters->yTranslate < TRANSLATE_MIN)
        renderParameters->yTranslate = TRANSLATE_MIN;
    else if (renderParameters->yTranslate > TRANSLATE_MAX)
        renderParameters->yTranslate = TRANSLATE_MAX;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::yTranslateChanged()

// slot for toggling wireframe
void RenderController::useWireframeCheckChanged(int state)
    { // RenderController::useWireframeCheckChanged()
    // reset the model's flag
    renderParameters->useWireframe = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::useWireframeCheckChanged()

// slot for toggling normal mapping
void RenderController::useNormalCheckChanged(int state)
    { // RenderController::useNormalCheckChanged()
    // reset the model's flag
    renderParameters->useNormal = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::useNormalCheckChanged()
    
// slot for toggling texture coords
void RenderController::useTexCoordsCheckChanged(int state)
    { // RenderController::useTexCoordsCheckChanged()
    // reset the model's flag
    renderParameters->useTexCoords = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::useTexCoordsCheckChanged()
    
// slot for toggling texture render
void RenderController::renderTextureCheckChanged(int state)
    { // RenderController::renderTextureCheckChanged()
    // reset the model's flag
    renderParameters->renderTexture = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::renderTextureCheckChanged()
    
// slot for toggling normal map render
void RenderController::renderNormalMapCheckChanged(int state)
    { // RenderController::renderNormalMapCheckChanged()
    // reset the model's flag
    renderParameters->renderNormalMap = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::renderNormalMapCheckChanged()
    
// slots for responding to arcball manipulations
// these are general purpose signals which pass the mouse moves to the controller
// after scaling to the notional unit sphere
void RenderController::BeginScaledDrag(int whichButton, float x, float y)
    { // RenderController::BeginScaledDrag()
    // depends on which button was depressed, so save that for the duration
    dragButton = whichButton;

    // now switch on it to determine behaviour
    switch (dragButton)
        { // switch on the drag button
        // left button drags the model
        case Qt::LeftButton:
            renderWindow->modelRotator->BeginDrag(x, y);
            break;
        } // switch on the drag button

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::BeginScaledDrag()
    
// note that Continue & End assume the button has already been set
void RenderController::ContinueScaledDrag(float x, float y)
    { // RenderController::ContinueScaledDrag()
    // switch on the drag button to determine behaviour
    switch (dragButton)
        { // switch on the drag button
        // left button drags the model
        case Qt::LeftButton:
            renderWindow->modelRotator->ContinueDrag(x, y);
            break;
        } // switch on the drag button

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::ContinueScaledDrag()

void RenderController::EndScaledDrag(float x, float y)
    { // RenderController::EndScaledDrag()
    // now switch on it to determine behaviour
    switch (dragButton)
        { // switch on the drag button
        // left button drags the model
        case Qt::LeftButton:
            renderWindow->modelRotator->EndDrag(x, y);
            break;
        } // switch on the drag button

    // and reset the drag button
    dragButton = Qt::NoButton;

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::EndScaledDrag()

