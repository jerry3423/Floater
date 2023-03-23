/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Parameters
//  -----------------------------
//  
//  This is part of the "model" in the MVC paradigm
//  We separate out the render parameters from the object rendered
//
/////////////////////////////////////////////////////////////////

// include guard
#ifndef _RENDER_PARAMETERS_H
#define _RENDER_PARAMETERS_H

#include "Matrix4.h"

// class for the render parameters
class RenderParameters
    { // class RenderParameters
    public:
    
    // we have a widget with an arcball which stores the rotation
    // we'll be lazy and leave it there instead of moving it here
    
    // we store x & y translations
    float xTranslate, yTranslate;

    // and a zoom scale
    float zoomScale;
    
    Matrix4 rotationMatrix;
    
    // and the booleans
    bool useWireframe;
	bool useNormal;
	bool useTexCoords;
	bool renderTexture;
	bool renderNormalMap;

    // constructor
    RenderParameters()
        :
        xTranslate(0.0), 
        yTranslate(0.0),
        zoomScale(1.0),
        useWireframe(false),
        useNormal(false),
        useTexCoords(false),
        renderTexture(false),
        renderNormalMap(false)
        { // constructor

        // because we are paranoid, we will initialise the matrices to the identity
        rotationMatrix.SetIdentity();
        } // constructor

    // accessor for scaledXTranslate

    }; // class RenderParameters

// now define some macros for bounds on parameters
#define TRANSLATE_MIN -1.0
#define TRANSLATE_MAX 1.0

#define ZOOM_SCALE_LOG_MIN -2.0
#define ZOOM_SCALE_LOG_MAX 2.0
#define ZOOM_SCALE_MIN 0.01
#define ZOOM_SCALE_MAX 100.0

// this is to scale to/from integer values
#define PARAMETER_SCALING 100

// end of include guard
#endif