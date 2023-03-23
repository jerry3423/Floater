///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
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

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))

#define N_ITERATIONS 100000

unsigned long numofedge;

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
	otherHalf.resize(0);
    } // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();
        
//         std::cout << "Read: " << firstChar << std::endl;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;
                
            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();
                
                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        geometryStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;                  
                        } // tex coord
                    default:
                        break;
                    } // switch on second character 
                break;
                } // some sort of vertex data
                
            case 'f':       // face data
                { // face
				// make a hard assumption that we have a single triangle per line
                unsigned int vertexID;
                
                // read in three vertices
				for (unsigned int vertex = 0; vertex < 3; vertex++)
					{ // per vertex
					// read a vertex ID
					geometryStream >> vertexID;

					// subtract one and store them (OBJ uses 1-based numbering)
					faceVertices.push_back(vertexID-1);
					} // per vertex
				break;
                } // face
                
            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
                
            } // per vertex
        } // non-empty vertex set

// 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
// 	std::cout << "Object Size:       " << objectSize << std::endl;

    BuildhalfEdge();
    TextureAssign();
    // return a success code
    return true;
	} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        geometryStream << "f";
        
        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
			} // per vertex
		// end the line
        geometryStream << std::endl;
        } // per face
    
    } // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
    { // Render()
	// make sure that textures are disabled
    //glDisable(GL_TEXTURE_2D);

	float scale = renderParameters->zoomScale;
	scale /= objectSize;
	// Scale defaults to the zoom setting
	glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);

		
	if (renderParameters->useWireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // start rendering
    glBegin(GL_TRIANGLES);
	
    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        
		// now do a loop over three vertices
		for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
            if(renderParameters->useTexCoords){
                glColor3f
                    (
                    textureCoords[faceVertices[face+vertex]].x * 255,
                    textureCoords[faceVertices[face+vertex]].y * 255,
                    textureCoords[faceVertices[face+vertex]].z * 255
                    );
            }
            else if(renderParameters->useNormal || renderParameters->renderNormalMap){
                glColor3f
                    (
                    0.5 + normals[faceVertices[face+vertex]].x * 0.5,
                    0.5 + normals[faceVertices[face+vertex]].y * 0.5,
                    0.5 + normals[faceVertices[face+vertex]].z * 0.5
                    );
            }
            else{
                // set colour using vertex ID
                glColor3f
                    (
                    colours[faceVertices[face+vertex]].x,
                    colours[faceVertices[face+vertex]].y,
                    colours[faceVertices[face+vertex]].z
                    );
            }

            if(renderParameters->renderTexture){
                glVertex3f
                    (
                    scale * textureCoords[faceVertices[face+vertex]].x,
                    scale * textureCoords[faceVertices[face+vertex]].y,
                    scale * textureCoords[faceVertices[face+vertex]].z
                    );
            }
            else{
              // use scaled xyz for vertex position
                glVertex3f
                    (
                    scale * vertices[faceVertices[face+vertex]].x,
                    scale * vertices[faceVertices[face+vertex]].y,
                    scale * vertices[faceVertices[face+vertex]].z
                    );
            }
			} // per vertex
        } // per face

    // close off the triangles
    glEnd();

    // revert render mode  
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    //I tried, but not working
    /*if(renderParameters->renderNormalMap){
        glDrawPixels(1024,1024, GL_RGBA, GL_UNSIGNED_BYTE,normalMap);
    }
    if(renderParameters->renderTexture){
        glDrawPixels(1024,1024, GL_RGBA, GL_UNSIGNED_BYTE,uvMap);
    }*/

    } // Render()

long AttributedObject::checkboundaryEdge(unsigned long x,unsigned long y){
    std::map<std::pair<long, long>, long>::iterator it = halfedge2faceindex.find(std::make_pair(x, y));
        if (it == halfedge2faceindex.end()) {
            return -1;
        }
        return it->second;
}

long AttributedObject::findoppositeEdge(unsigned long x, unsigned long y){
    std::map<std::pair<long, long>, long>::iterator it = de2edgeindex.find(std::make_pair(x, y));
        if (it == de2edgeindex.end()) {
            return -1;
        }
        return it->second;
}

void AttributedObject::findboundaryVertex(){
    for(unsigned int i=0;i<halfedge.size();i++){
        if(halfedge[i].boundary_flag==true){
            boundaryVertex[halfedge[i].vertex_from] = true;
            boundaryVertex[halfedge[i].vertex_to] = true;
        }
    }
}

std::vector<unsigned long> AttributedObject::findneighborVertex(unsigned long index){
    std::vector<unsigned long>temp;
    for(unsigned long i=0;i<halfedge.size();i++){
        if(halfedge[i].vertex_to==(index)){
            temp.push_back(halfedge[i].vertex_from);
        }
        if(halfedge[i].vertex_from==(index)){
            temp.push_back(halfedge[i].vertex_to);
        }
    }

    sort(temp.begin(), temp.end());
    auto it = unique(temp.begin(), temp.end());
    temp.erase(it, temp.end());
    return temp;
}

//Creating a set of undirected edges
void AttributedObject::edgecreate() {
    typedef std::set<std::pair<long, long>> edgeset;
    edgeset edgess;
    for (long i = 0; i < faceVertices.size(); i+=3) {
        edgess.insert(std::make_pair(std::min(faceVertices[i], faceVertices[i+1]), std::max(faceVertices[i], faceVertices[i+1])));
        edgess.insert(std::make_pair(std::min(faceVertices[i+1], faceVertices[i+2]), std::max(faceVertices[i+1], faceVertices[i+2])));
        edgess.insert(std::make_pair(std::min(faceVertices[i+2], faceVertices[i]), std::max(faceVertices[i+2], faceVertices[i])));
    }
    numofedge = edgess.size();
}

void AttributedObject::BuildhalfEdge(){
    //create a map to store all the directed edge based on faceindex in triangles
        for (unsigned long i = 0; i < faceVertices.size(); i+=3) {
            halfedge2faceindex[std::make_pair(faceVertices[i], faceVertices[i+1])] = i%3;
            halfedge2faceindex[std::make_pair(faceVertices[i+1], faceVertices[i+2])] = i%3;
            halfedge2faceindex[std::make_pair(faceVertices[i+2], faceVertices[i])] = i%3;
        }

        //store the fromvertex and tovertex based on the halfedge index
        for (unsigned long i = 0; i < faceVertices.size(); i+=3) {
            HalfEdge he0;
            HalfEdge he1;
            HalfEdge he2;
            he0.face_index = checkboundaryEdge(faceVertices[i+2], faceVertices[i]);
            he0.vertex_to = faceVertices[i];
            he0.vertex_from = faceVertices[i+2];
            he0.edge_index = i;
            halfedge.push_back(he0);
            he1.face_index = checkboundaryEdge(faceVertices[i], faceVertices[i+1]);
            he1.vertex_to = faceVertices[i+1];
            he1.vertex_from = faceVertices[i];
            he1.edge_index = i+1;
            halfedge.push_back(he1);
            he2.face_index = checkboundaryEdge(faceVertices[i+1], faceVertices[i+2]);
            he2.vertex_to = faceVertices[i+2];
            he2.vertex_from = faceVertices[i+1];
            he2.edge_index = i+2;
            halfedge.push_back(he2);

        }

        //Determine if the opposite halfedge is on the triangle, if not, then it is the boundary halfedge
        //store the fromvertex and tovertex based on the halfedge index
        unsigned long count = 0;
        for (unsigned long i = 0; i < faceVertices.size(); i+=3) {
            if (checkboundaryEdge(faceVertices[i], faceVertices[i+2]) == -1) {
                HalfEdge temp;
                temp.face_index = -1;
                temp.vertex_from = faceVertices[i];
                temp.vertex_to = faceVertices[i+2];
                temp.edge_index = count + faceVertices.size();
                temp.boundary_flag = true;
                halfedge.push_back(temp);
                count++;
                if (count == 2 * numofedge - faceVertices.size()) {
                    break;
                }
            }
            if (checkboundaryEdge(faceVertices[i+1], faceVertices[i]) == -1) {
                HalfEdge temp;
                temp.face_index = -1;
                temp.vertex_from = faceVertices[i+1];
                temp.vertex_to = faceVertices[i];
                temp.edge_index = count + faceVertices.size();
                temp.boundary_flag = true;
                halfedge.push_back(temp);
                count++;
                if (count == 2 * numofedge - faceVertices.size()) {
                    break;
                }
            }
            if (checkboundaryEdge(faceVertices[i+2], faceVertices[i+1]) == -1) {
                HalfEdge temp;
                temp.face_index = -1;
                temp.vertex_from = faceVertices[i+2];
                temp.vertex_to = faceVertices[i+1];
                temp.edge_index = count +faceVertices.size();
                temp.boundary_flag = true;
                halfedge.push_back(temp);
                count++;
                if (count == 2 * numofedge - faceVertices.size()) {
                    break;
                }
            }
        }
        //make boundary half edge ordered
        std::vector<HalfEdge>unorderhalfedges;
        for(unsigned int i=faceVertices.size();i<halfedge.size();i++){
            unorderhalfedges.push_back(halfedge[i]);
        }
        unsigned int cnt = faceVertices.size();
        unsigned int a = unorderhalfedges[0].vertex_to;
        for(unsigned int i=1;i<unorderhalfedges.size();i++){
            cnt++;
            for(unsigned int j=0;j<unorderhalfedges.size();j++){
                unsigned int b = unorderhalfedges[j].vertex_from;
                if(a==b){
                    halfedge[cnt] = unorderhalfedges[j];
                    a = unorderhalfedges[j].vertex_to;
                    break;
                }
            }
        }

        //find the nexthalfedge in triangles based on the halfedge struction
        for (unsigned long i = 0; i < faceVertices.size(); i += 3) {
            halfedge.at(i).next_edge = halfedge.at(i + 1).edge_index;
            halfedge.at(i + 1).next_edge = halfedge.at(i + 2).edge_index;
            halfedge.at(i + 2).next_edge = halfedge.at(i).edge_index;
        }
        //find the nexthalfedge outside the triangles based on the halfedge struction
        for (unsigned long i = 0; i < halfedge.size() - faceVertices.size() - 1; i++) {
            halfedge.at(i + faceVertices.size()).next_edge = halfedge.at(i + faceVertices.size() + 1).edge_index;
        }
        halfedge.at(halfedge.size() - 1).next_edge = halfedge.at(faceVertices.size() - 1).edge_index;

        //create a map to store all the directed edge based on halfedge index in triangles
        for (unsigned long i = 0; i < faceVertices.size(); i+=3) {
            de2edgeindex[std::make_pair(faceVertices[i+2], faceVertices[i])] = i;
            de2edgeindex[std::make_pair(faceVertices[i], faceVertices[i+1])] = i + 1;
            de2edgeindex[std::make_pair(faceVertices[i+1], faceVertices[i+2])] = i + 2;
        }

        for (unsigned long i = faceVertices.size(); i < 2 * numofedge; i++) {
            de2edgeindex[std::make_pair(halfedge.at(i).vertex_from, halfedge.at(i).vertex_to)] = i;
        }

        //find the otherhalfedge index
        for (unsigned long i = 0; i < halfedge.size(); i++) {
            halfedge.at(i).opposite_edge = findoppositeEdge(halfedge.at(i).vertex_to, halfedge.at(i).vertex_from);
        }

        firstDirectedEdge.resize(vertices.size(),-1);
        for(unsigned int i=0;i<vertices.size();i++){
            for(unsigned int j=0;j<halfedge.size();j++){
                if(i==halfedge[j].vertex_to){
                    firstDirectedEdge[i] = halfedge[j].edge_index;
                    break;
                }
            }
        }
        findboundaryVertex();
        for(unsigned int i=0;i<vertices.size();i++){
            neighborVertex.push_back(findneighborVertex(i));
            orderedneighborVertex.push_back(findneighborVertex(i));
        }

        for(unsigned int i=0;i<vertices.size();i++){
            textureCoords.push_back(vertices[i]);
        }
}

double AttributedObject::CalculateAngle(Cartesian3 v1, Cartesian3 v2){
    Cartesian3 temp1 = v1.unit();
    Cartesian3 temp2 = v2.unit();
    double cosangle = temp1.dot(temp2) / (temp1.length() * temp2.length());
    double angle = acos(cosangle);
    return angle;
}

void AttributedObject::TextureAssign(){
    //assign new boundary vertex coordinations to square
    for(unsigned int i=faceVertices.size();i<halfedge.size();i++){
        double a = (i - faceVertices.size()) / double(halfedge.size() - faceVertices.size());
        double x,y,z = 0;
        if(a >=0 && a< 0.25){
            x = 0.5 - 4 * a;
            y = 1;
        }
        if(a >=0.25 && a< 0.5){
            x = -0.5;
            y = 2 - 4 * a;
        }
        if(a >=0.5 && a< 0.75){
            x = 4 * a - 2.5;
            y = 0;
        }
        if(a >=0.75 && a< 1){
            x = 0.5;
            y = 4 * a - 3;
        }
        textureCoords[halfedge[i].vertex_from] = Cartesian3(x,y,z);
    }

    //make the neighbor vertex ordered
    for(unsigned int i=0;i<vertices.size();i++){
        if(!boundaryVertex[i] && neighborVertex[i].size()!=0){
            orderedneighborVertex[i].at(0) = neighborVertex[i].at(0);
            for(unsigned int j = 0;j<neighborVertex[i].size() - 1;j++){
                unsigned int neighborj = orderedneighborVertex[i].at(j);
                unsigned int p = halfedge[halfedge[de2edgeindex[std::make_pair(i,neighborj)]].next_edge].vertex_to;
                orderedneighborVertex[i].at(j + 1) = p;
            }
        }
    }


    //average weight
    for(unsigned int k=0;k<500;k++){
        for(unsigned int i=0;i<vertices.size();i++){
            if(!boundaryVertex[i] && neighborVertex[i].size()!=0){
                Cartesian3 sum = Cartesian3(0,0,0);
                unsigned int num = neighborVertex[i].size();
                for(unsigned int j = 0;j<neighborVertex[i].size();j++){
                    sum = sum + textureCoords[neighborVertex[i].at(j)];
                }
                textureCoords[i] = sum / float(num);
                textureCoords[i].z = 0;
            }
        }
    }

    /*
    //floater algorithm
    for(unsigned int i=0;i<vertices.size();i++){
        if(!boundaryVertex[i] && neighborVertex[i].size()!=0){
            double sum_weight = 0;
            Cartesian3 sum = Cartesian3(0,0,0);
            for(unsigned int j = 0;j<neighborVertex[i].size();j++){
                unsigned int neighborj = neighborVertex[i].at(j);
                //two neighbor vertex of neighborj
                unsigned int p = halfedge[halfedge[de2edgeindex[std::make_pair(i,neighborj)]].next_edge].vertex_to;
                unsigned int q = halfedge[halfedge[de2edgeindex[std::make_pair(neighborj,i)]].next_edge].vertex_to;
                Cartesian3 i_p = vertices[p] - vertices[i];
                Cartesian3 i_q = vertices[q] - vertices[i];
                Cartesian3 i_neighborj = vertices[neighborj] - vertices[i];
                //calculate two angles
                double angle1 = CalculateAngle(i_p,i_neighborj);
                double angle2 = CalculateAngle(i_q,i_neighborj);
                Cartesian3 temp1 = Cartesian3(textureCoords[neighborj].x,textureCoords[neighborj].y,0);
                Cartesian3 temp2 = Cartesian3(textureCoords[i].x,textureCoords[i].y,0);
                Cartesian3 real_dis = textureCoords[neighborj] - textureCoords[i];
                double weight = (1 / real_dis.length()) * (tan(angle1 / 2.0f) + tan(angle2 / 2.0f));
                sum_weight+=weight;
                sum  = sum + weight * vertices[neighborj];
            }
            textureCoords[i] = -sum_weight * sum;
            textureCoords[i].z = 0;
        }
    }*/

    //calculate normal vector
    normals.resize(vertices.size());
    for(unsigned int i=0;i<faceVertices.size();i+=3){
        Cartesian3 ab = vertices[faceVertices[i+1]] - vertices[faceVertices[i]];
        Cartesian3 ac = vertices[faceVertices[i+2]] - vertices[faceVertices[i]];
        normals[faceVertices[i]] = ab.cross(ac).unit();
        Cartesian3 bc = vertices[faceVertices[i+2]] - vertices[faceVertices[i+1]];
        Cartesian3 ba = vertices[faceVertices[i]] - vertices[faceVertices[i+1]];
        normals[faceVertices[i+1]] = bc.cross(ba).unit();
        Cartesian3 ca = vertices[faceVertices[i]] - vertices[faceVertices[i+2]];
        Cartesian3 cb = vertices[faceVertices[i+1]] - vertices[faceVertices[i+2]];
        normals[faceVertices[i+2]] = ca.cross(cb).unit();

    }

    TextureMap();
    NormalMap();

}

void AttributedObject::TextureMap(){
    uvMap = (Cartesian3 *)calloc(1024 * 1024,sizeof (Cartesian3));
    Cartesian3 color(1,1,1);
    for(unsigned int i=0;i<vertices.size();i++){
        unsigned int v = (textureCoords[i].y) * 1024;
        unsigned int u = (textureCoords[i].x + 0.5) * 1024;

        uvMap[u * 1023 + v] = 255 * color;

    }

}

void AttributedObject::NormalMap(){
    normalMap = (Cartesian3 *)calloc(1024 * 1024,sizeof (Cartesian3));
    for(unsigned int i=0;i<vertices.size();i++){
        unsigned int v = (textureCoords[i].y) * 1024;
        unsigned int u = (textureCoords[i].x + 0.5) * 1024;

        normalMap[u * 1023 + v].x = 0.5 + normals[i].x * 0.5;
        normalMap[u * 1023 + v].y = 0.5 + normals[i].y * 0.5;
        normalMap[u * 1023 + v].z = 0.5 + normals[i].z * 0.5;

    }
}



