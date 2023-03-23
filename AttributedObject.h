///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.h
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//  
///////////////////////////////////////////////////

// include guard for AttributedObject
#ifndef _ATTRIBUTED_OBJECT_H
#define _ATTRIBUTED_OBJECT_H

// include the C++ standard libraries we need for the header
#include <vector>
#include<map>
#include<set>
#include<math.h>
#include<algorithm>
#include <iostream>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"

// define a macro for "not used" flag
#define NO_SUCH_ELEMENT -1
#define M_PI 3.1415926

// use macros for the "previous" and "next" IDs
#define PREVIOUS_EDGE(x) ((x) % 3) ? ((x) - 1) : ((x) + 2)
#define NEXT_EDGE(x) (((x) % 3) == 2) ? ((x) - 2) : ((x) + 1)

class AttributedObject
    { // class AttributedObject
    public:
    // vector of vertices
    std::vector<Cartesian3> vertices;

    // vector of colours stored as cartesian triples in float
    std::vector<Cartesian3> colours;

    // vector of normals
    std::vector<Cartesian3> normals;

    // vector of texture coordinates (stored as triple to simplify code)
    std::vector<Cartesian3> textureCoords;

    // vector for the "first" directed edge for each vertex
    std::vector<long> firstDirectedEdge;

    // vector of faces - doubles as the "to" array for edges
    std::vector<unsigned int> faceVertices;

    // vector of the "other half" of the directed edges
    std::vector<long> otherHalf;

    struct HalfEdge{
        long edge_index;
        long face_index;
        long vertex_from;
        long vertex_to;
        long next_edge;
        long opposite_edge;
        bool boundary_flag = false;
        HalfEdge():
            edge_index(-1),
            face_index(-1),
            vertex_from(-1),
            vertex_to(-1),
            next_edge(-1),
            opposite_edge(-1){}
    };

    //store half edge structure
    std::vector<AttributedObject::HalfEdge>halfedge;

    std::map<std::pair<long, long>, long> halfedge2faceindex;

    std::map<std::pair<long,long>,long> de2edgeindex;

    std::map<std::pair<unsigned int, unsigned int>, unsigned int> halfedgemap;

    std::map<unsigned long,bool>boundaryVertex;

    //store unordered 1-ring neighbor vertex
    std::vector<std::vector<unsigned long>>neighborVertex;

    //store ordered 1-ring neighbor vertex
    std::vector<std::vector<unsigned long>>orderedneighborVertex;

    //generate texture map
    Cartesian3 *uvMap;

    //generate normal map
    Cartesian3 *normalMap;



    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    // constructor will initialise to safe values
    AttributedObject();
   
    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    // routine to render
    void Render(RenderParameters *renderParameters);

    void edgecreate();

    void BuildhalfEdge();

    long checkboundaryEdge(unsigned long x,unsigned long y);

    long findoppositeEdge(unsigned long x,unsigned long y);

    void findboundaryVertex();

    std::vector<unsigned long> findneighborVertex(unsigned long index);

    double CalculateAngle(Cartesian3 v1, Cartesian3 v2);

    void TextureAssign();

    void TextureMap();

    void NormalMap();
    }; // class AttributedObject

// end of include guard for AttributedObject
#endif
