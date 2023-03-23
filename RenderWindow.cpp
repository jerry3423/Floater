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

#include "RenderWindow.h"
#include "RenderParameters.h"

// constructor
RenderWindow::RenderWindow
        (
        // the object to be rendered
        AttributedObject         *newAttributedObject, 
        // the model object storing render parameters
        RenderParameters        *newRenderParameters,
        // the title for the window (with default value)
        const char              *windowName
        )
    // call the inherited constructor
    // NULL indicates that this widget has no parent
    // i.e. that it is a top-level window
    :
    // member instantiation
    QWidget(NULL),
    attributedObject(newAttributedObject),
    renderParameters(newRenderParameters)
    { // RenderWindow::RenderWindow()
    // set the window's title
    setWindowTitle(QString(windowName));
    
    // initialise the grid layout
    windowLayout = new QGridLayout(this);
    
    // create all of the widgets, starting with the custom render widgets
    renderWidget                = new RenderWidget              (newAttributedObject,     		newRenderParameters,        this);

    // construct custom arcball Widgets
    modelRotator                = new ArcBallWidget             (                       this);

    // construct standard QT widgets
    
    // modelling options
    renderWireframeBox          = new QCheckBox                 ("Wireframe",           this);  
    useNormalBox                = new QCheckBox                 ("Nrm -> RGB",          this);  
    useTexCoordsBox             = new QCheckBox                 ("UVW -> RGB",          this);
    renderTextureBox            = new QCheckBox                 ("Texture",             this);
    renderNormalMapBox          = new QCheckBox                 ("Normal Map",          this);
    
    // spatial sliders
    xTranslateSlider            = new QSlider                   (Qt::Horizontal,        this);
    yTranslateSlider            = new QSlider                   (Qt::Vertical,          this);
    zoomSlider                  = new QSlider                   (Qt::Vertical,          this);
    
    // labels for sliders and arcballs
    modelRotatorLabel           = new QLabel                    ("Model",               this);
    yTranslateLabel             = new QLabel                    ("Y",                   this);
    zoomLabel                   = new QLabel                    ("Zm",                  this);
    
    // add all of the widgets to the grid               Row         Column      Row Span    Column Span
    
    int nStacked = 7;
    
    windowLayout->addWidget(renderWidget,               0,          1,          nStacked,   1           );
    windowLayout->addWidget(yTranslateSlider,           0,          2,          nStacked,   1           );
    windowLayout->addWidget(zoomSlider,                 0,          4,          nStacked,   1           );

    // the stack in the middle

    windowLayout->addWidget(modelRotator,               0,          3,          1,          1           );
    windowLayout->addWidget(modelRotatorLabel,          1,          3,          1,          1           );
    windowLayout->addWidget(renderWireframeBox,         2,          3,          1,          1           );
    windowLayout->addWidget(useNormalBox,               3,          3,          1,          1           );
    windowLayout->addWidget(useTexCoordsBox,            4,          3,          1,          1           );
    windowLayout->addWidget(renderTextureBox,           5,          3,          1,          1           );
    windowLayout->addWidget(renderNormalMapBox,         6,          3,          1,          1           );

    // Translate Slider Row
    windowLayout->addWidget(xTranslateSlider,           nStacked,   1,          1,          1           );
    windowLayout->addWidget(yTranslateLabel,            nStacked,   2,          1,          1           );
    // nothing in column 3
    windowLayout->addWidget(zoomLabel,                  nStacked,   4,          1,          1           );
    
    // now reset all of the control elements to match the render parameters passed in
    ResetInterface();
    } // RenderWindow::RenderWindow()

// routine to reset interface
// sets every visual control to match the model
// gets called by the controller after each change in the model
void RenderWindow::ResetInterface()
    { // RenderWindow::ResetInterface()
    // set check boxes
    renderWireframeBox      ->setChecked        (renderParameters   ->  useWireframe);
    useNormalBox            ->setChecked        (renderParameters   ->  useNormal);
    useTexCoordsBox         ->setChecked        (renderParameters   ->  useTexCoords);
    renderTextureBox        ->setChecked        (renderParameters   ->  renderTexture);
    renderNormalMapBox      ->setChecked        (renderParameters   ->  renderNormalMap);
    
    // set sliders
    // x & y translate are scaled to notional unit sphere in render widgets
    // but because the slider is defined as integer, we multiply by a 100 for all sliders
    xTranslateSlider        ->setMinimum        ((int) (TRANSLATE_MIN                               * PARAMETER_SCALING));
    xTranslateSlider        ->setMaximum        ((int) (TRANSLATE_MAX                               * PARAMETER_SCALING));
    xTranslateSlider        ->setValue          ((int) (renderParameters -> xTranslate              * PARAMETER_SCALING));
    
    yTranslateSlider        ->setMinimum        ((int) (TRANSLATE_MIN                               * PARAMETER_SCALING));
    yTranslateSlider        ->setMaximum        ((int) (TRANSLATE_MAX                               * PARAMETER_SCALING));
    yTranslateSlider        ->setValue          ((int) (renderParameters -> yTranslate              * PARAMETER_SCALING));

    // zoom slider is a logarithmic scale, so we want a narrow range
    zoomSlider              ->setMinimum        ((int) (ZOOM_SCALE_LOG_MIN                          * PARAMETER_SCALING));
    zoomSlider              ->setMaximum        ((int) (ZOOM_SCALE_LOG_MAX                          * PARAMETER_SCALING));
    zoomSlider              ->setValue          ((int) (log10(renderParameters -> zoomScale)        * PARAMETER_SCALING));

    // now flag them all for update 
    renderWidget            ->update();
    modelRotator            ->update();
    xTranslateSlider        ->update();
    yTranslateSlider        ->update();
    zoomSlider              ->update();
    renderWireframeBox      ->update();
    useNormalBox            ->update();
    useTexCoordsBox         ->update();
    renderTextureBox        ->update();
    renderNormalMapBox      ->update();
    } // RenderWindow::ResetInterface()
