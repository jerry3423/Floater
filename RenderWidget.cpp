//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Widget
//  -----------------------------
//  
//  Since the render code is in the geometric object class
//  this widget primarily sets up the transformation matrices and 
//  lighting.
//
//  It implements the UI for an arcball controller, but in an
//  abstractable way that allows it to share an arcball with another
//  widget (the arcball controller) - thus you can manipulate the
//  object either directly in the widget or indirectly through the
//  arcball displayed visually.
//  
//  It also supports translation by visual dragging.  This results in
//  needing three mouse buttons, and not all systems have them, so it
//  will be up to the user to make sure that they can use it on their
//  own machine (especially on OSX)
//
//  Since the controls are (potentially) shared with other widgets, 
//  this widget is only responsible for scaling the x,y of mouse events
//  then passing them to the controller
//  
////////////////////////////////////////////////////////////////////////

#include <math.h>

// include the header file
#include "RenderWidget.h"

// constructor
RenderWidget::RenderWidget
        (   
        // the geometric object to show
        AttributedObject 	*newAttributedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // parent widget in visual hierarchy
        QWidget             *parent
        )
    // the : indicates variable instantiation rather than arbitrary code
    // it is considered good style to use it where possible
    : 
    // start by calling inherited constructor with parent widget's pointer
    QOpenGLWidget(parent),
    // then store the pointers that were passed in
    attributedObject(newAttributedObject),
    renderParameters(newRenderParameters)
    { // constructor
    // leaves nothing to put into the constructor body
    } // constructor    

// destructor
RenderWidget::~RenderWidget()
    { // destructor
    // empty (for now)
    // all of our pointers are to data owned by another class
    // so we have no responsibility for destruction
    // and OpenGL cleanup is taken care of by Qt
    } // destructor                                                                 

// called when OpenGL context is set up
void RenderWidget::initializeGL()
    { // RenderWidget::initializeGL()
    // set lighting parameters (may be reset later)
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    
    // background is yellowish-grey
    glClearColor(0.8, 0.8, 0.6, 1.0);

    // enable depth-buffering
	glEnable(GL_DEPTH_TEST);
    } // RenderWidget::initializeGL()

// called every time the widget is resized
void RenderWidget::resizeGL(int w, int h)
    { // RenderWidget::resizeGL()
    // reset the viewport
    glViewport(0, 0, w, h);
    
    // set projection matrix to be glOrtho based on zoom & window size
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // compute the aspect ratio of the widget
    float aspectRatio = (float) w / (float) h;
    
    // we want to capture a sphere of radius 1.0 without distortion
    // so we set the ortho projection based on whether the window is portrait (> 1.0) or landscape
    // portrait ratio is wider, so make bottom & top -1.0 & 1.0
    if (aspectRatio > 1.0)
        glOrtho(-aspectRatio, aspectRatio, -1.0, 1.0, -1.0, 1.0);
    // otherwise, make left & right -1.0 & 1.0
    else
        glOrtho(-1.0, 1.0, -1.0/aspectRatio, 1.0/aspectRatio, -1.0, 1.0);

    } // RenderWidget::resizeGL()
    
// called every time the widget needs painting
void RenderWidget::paintGL()
    { // RenderWidget::paintGL()
    // clear the buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set model view matrix based on stored translation, rotation &c.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
        
    // start with lighting turned off
    glDisable(GL_LIGHTING);

    // translate by the visual translation
    glTranslatef(renderParameters->xTranslate, renderParameters->yTranslate, 0.0f);

    // apply rotation matrix from arcball
    glMultMatrixf(renderParameters->rotationMatrix.columnMajor().coordinates);

    // tell the object to draw itself, 
    // passing in the render parameters for reference
    attributedObject->Render(renderParameters);
    } // RenderWidget::paintGL()
    
// mouse-handling
void RenderWidget::mousePressEvent(QMouseEvent *event)
    { // RenderWidget::mousePressEvent()
    // store the button for future reference
    int whichButton = event->button();
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;
    
    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    int modifiers = event->modifiers();
    
    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;
    
    // send signal to the controller for detailed processing
    emit BeginScaledDrag(whichButton, x,y);
    } // RenderWidget::mousePressEvent()
    
void RenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // RenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // RenderWidget::mouseMoveEvent()
    
void RenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // RenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // RenderWidget::mouseReleaseEvent()
