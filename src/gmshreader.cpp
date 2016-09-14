#include "gmshreader.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

/*    GMSH supported element types
 *
 *  ( source: http://gmsh.info/doc/texinfo/gmsh.pdf )
 *
 *  1 : 2-node line.
 *  2 : 3-node triangle.
 *  3 : 4-node quadrangle.
 *  4 : 4-node tetrahedron.
 *  5 : 8-node hexahedron.
 *  6 : 6-node prism.
 *  7 : 5-node pyramid.
 *  8 : 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
 *  9 : 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
 * 10 : 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
 * 11 : 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
 * 12 : 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
 * 13 : 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
 * 14 : 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
 * 15 : 1-node point.
 * 16 : 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
 * 17 : 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
 * 18 : 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
 * 19 : 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
 * 20 : 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
 * 21 : 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
 * 22 : 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
 * 23 : 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
 * 24 : 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
 * 25 : 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
 * 26 : 4-node  third  order  edge  (2  nodes  associated  with  the  vertices,  2 internal to the edge)
 * 27 : 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
 * 28 : 6-node  fifth  order  edge  (2  nodes  associated  with  the  vertices,  4 internal to the edge)
 * 29 : 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
 * 30 : 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
 * 31 : 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
 * 92 : 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)
 * 93 : 125-node  fourth  order  hexahedron  (8  nodes  associated  with  the vertices, 36 with the edges, 54 with the faces, 27 in the volume)
 */


using namespace std;
using namespace simit;

int read(string gmsh_file, Set& node, Set& line, Set& triangle, bool verbose) {

    int numb_nodes[100]; // number of nodes for a given element type, see description above
    numb_nodes[ 1] =  2;
    numb_nodes[ 2] =  3;
    numb_nodes[ 3] =  4;
    numb_nodes[ 4] =  4;
    numb_nodes[ 5] =  8;
    numb_nodes[ 6] =  6;
    numb_nodes[ 7] =  5;
    numb_nodes[ 8] =  3;
    numb_nodes[ 9] =  6;
    numb_nodes[10] =  9;
    numb_nodes[11] =  10;
    numb_nodes[12] =  27;
    numb_nodes[13] =  18;
    numb_nodes[14] =  14;
    numb_nodes[15] =  1;
    numb_nodes[16] =  8;
    numb_nodes[17] =  20;
    numb_nodes[18] =  15;
    numb_nodes[19] =  13;
    numb_nodes[20] =  9;
    numb_nodes[21] =  10;
    numb_nodes[22] =  12;
    numb_nodes[23] =  15;
    numb_nodes[24] =  15;
    numb_nodes[25] =  21;
    numb_nodes[26] =  4;
    numb_nodes[27] =  5;
    numb_nodes[28] =  6;
    numb_nodes[29] =  20;
    numb_nodes[30] =  35;
    numb_nodes[31] =  56;
    numb_nodes[92] =  64;
    numb_nodes[93] =  125;

    ifstream is(gmsh_file);
    if(!is.good()) {
        cerr << "Error: reading file: " << gmsh_file << endl;
        return 1;
    }

    FieldRef<int>      ref_node_i   = node.addField<int>("i");
    FieldRef<bool>     ref_node_edg = node.addField<bool>("edge");
    FieldRef<double,3> ref_node_x   = node.addField<double,3>("x");
    FieldRef<double>   ref_node_u   = node.addField<double>("u");

    FieldRef<int> ref_line_i = line.addField<int>("i");

    FieldRef<int> ref_tri_i  = triangle.addField<int>("i");

    if(verbose)
        cout << "Reading file " << gmsh_file << endl;
    vector<ElementRef> ref_node;
    while(is.good()) {
        int i;
        double x1,x2, x3;
        string tag;
        is >> tag;
        if(tag.compare("$MeshFormat") == 0) {
            double version;
            int filetype, datasize;
            is >> version >> filetype >> datasize;
            is >> tag;
            if(verbose) {
                cout << "GMSH version  : " << version << endl;
                cout << "     file type: " << ((filetype==0)?"ASCII":"BINARY") << endl;
                cout << "     data size: " << datasize << endl;
            }
                
            if( tag.compare("$EndMeshFormat") != 0) {
                cerr << "Error: malformed gmsh file" << endl;
                return 2;
            }
        } else if(tag.compare("$Nodes") == 0) {
            int n, i;
            double x,y,z;
            is >> n;
            for(int j=0; j<n; j++)
                ref_node.push_back(node.add());
            for(int j=0; j<n; j++) {
                // j is the order which they appear in the file (0-index)
                // i is the listed index of this node           (1-index)
                // Why would you ever create a file with i+1 != j? If this ever 
                // happens, then there is no guarantee that the ref_node array
                // will be fully populated, and may contain invalid refernces
                is >> i >> x >> y >> z;
                ref_node_x(  ref_node[j]) = {x,y,z};
                ref_node_u(  ref_node[j]) = 0.0;
                ref_node_edg(ref_node[j]) = false;
                ref_node_i(  ref_node[j]) = i;
                if(verbose)
                    cout << "Read node #" << i << ": (" << x << ", " << y << ", " << z << ")" << endl;
            }
            is >> tag;
            if( tag.compare("$EndNodes") != 0) {
                cerr << "Error: malformed gmsh file" << endl;
                return 2;
            }
        } else if(tag.compare("$Elements") == 0) {
            int n;
            int i, type, n_tags;
            int indices[125]; // largest supported element type is 125-node fourth order hexahedron
            is >> n;
            for(int j=0; j<n; j++) {
                is >> i >> type >> n_tags;
                int tags[n_tags];
                for(int k=0; k<n_tags; k++)
                    is >> tags[k];
                for(int k=0; k<numb_nodes[type]; k++)
                    is >> indices[k];
                if(verbose)
                    cout << "Read element #" << i << " of type = " << type << " (" << numb_nodes[type] << "-gon)" << endl;

                if(type == 1) {// line piece
                    ElementRef node1 = ref_node[indices[0]-1]; // gmsh is using 1-indexing
                    ElementRef node2 = ref_node[indices[1]-1];
                    ElementRef el    = line.add(node1, node2);
                } else if(type == 2) { // triangle
                    ElementRef node1 = ref_node[indices[0]-1]; // gmsh is using 1-indexing
                    ElementRef node2 = ref_node[indices[1]-1];
                    ElementRef node3 = ref_node[indices[2]-1];
                    ElementRef el    = triangle.add(node1, node2, node3);
                }
            }
            is >> tag;
            if( tag.compare("$EndElements") != 0) {
                cerr << "Error: malformed gmsh file" << endl;
                return 2;
            }
        }
    }
} 
