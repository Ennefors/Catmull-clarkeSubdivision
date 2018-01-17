//
// Created by bill on 2017-09-05.
//

#include <map>
#include <GL/glut.h>
#include "tdMathLibNew.h"
#include <sstream>
#include <fstream>
#include <string.h>
#include <chrono>

#ifndef GSCEPT_LAB_ENV_HALFEDGEMESH_H
#define GSCEPT_LAB_ENV_HALFEDGEMESH_H
#include <vector>

#define POW 12
#define BLOCKSIZE (1<<POW)
template <class T>
class mem_pool
{
private:
    struct block
    {
        T data[BLOCKSIZE];
    };
    std::vector<block*> blocks;
    int next_free;
public:
    mem_pool();
    void new_block();
    void clear();
    T* new_element();
    T* operator[](int);
};

template <class T>
inline mem_pool<T>::mem_pool()
{
    next_free = -1;
};

template<class T>
inline void mem_pool<T>::new_block()
{
    blocks.push_back(new block);
}

template<class T>
inline T* mem_pool<T>::new_element()
{
    next_free++;
    if((next_free&(BLOCKSIZE-1))==0)
    {
        new_block();
    }
    return &blocks[next_free>>POW]->data[(next_free&(BLOCKSIZE-1))];
};

template<class T>
inline void mem_pool<T>::clear()
{
    for (int i = 0; i < blocks.size(); ++i) {
        delete blocks[i];
    }
}

template<class T>
inline T* mem_pool<T>::operator[](int i)
{
    return &blocks[i>>POW]->data[(i&(BLOCKSIZE-1))];
};
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



class memhalfEdgeMesh
{


public:

    memhalfEdgeMesh(){};
    ~memhalfEdgeMesh()
    {
        edges.clear();
        faces.clear();
        verts.clear();
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

    mem_pool<HE_edge> edges;
    mem_pool<HE_vertex> verts;
    mem_pool<HE_face> faces;

    int nr_of_e = 0;
    int nr_of_f = 0;
    int nr_of_v = 0;
};


#endif //GSCEPT_LAB_ENV_HALFEDGEMESH_H
