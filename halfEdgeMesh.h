//
// Created by bill on 2017-09-05.
//

#include <map>
#include <GL/glut.h>
#include "tdMathLibNew.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>

#ifndef GSCEPT_LAB_ENV_HALFEDGEMESH_H
#define GSCEPT_LAB_ENV_HALFEDGEMESH_H
struct HE_edge;
struct HE_vertex;
struct HE_face;

struct HE_vertex
{
    u_int idx;
    fdVector pos;

    HE_edge* edge;
};

struct HE_face
{
    HE_edge* edge;
    int idx;
};

struct HE_edge
{
    HE_vertex* vert;
    HE_edge* pair = nullptr;
    HE_face* face;
    HE_edge* nxt;
    bool split = false;
    fdVector uv;

    int idx = 9999999999999999;
    int splitidx;
};



class halfEdgeMesh
{


public:

    halfEdgeMesh(){};
    ~halfEdgeMesh()
    {
        for (int i = 0; i < nr_of_e; ++i) {
            delete edges[i];
        }
        for (int j = 0; j < nr_of_f; ++j) {
            delete faces[j];
        }
        for (int k = 0; k < nr_of_v; ++k) {
            delete verts[k];
        }
        delete[] verts;
        delete[] edges;
        delete[] faces;
    };

    int steps = 0;

    int nrofsub = 0;

    float* get_vertices();
    float* get_uvs();


    void convert_str2(char* str);
    void parseobj(char* filepath);

    void generate_HE();
    void generate_Faces(int ,int ,int* ,float* , std::vector<fdVector>,float*);

    void edgesplit(HE_edge** edges, HE_face** faces, HE_vertex** vertices, int e, int loops, HE_edge* temp);
    void cmc_subdivide();

    void draw_mesh();

    HE_edge** edges;
    HE_vertex** verts;
    HE_face** faces;

    int nr_of_e = 0;
    int nr_of_f = 0;
    int nr_of_v = 0;
};


#endif //GSCEPT_LAB_ENV_HALFEDGEMESH_H
