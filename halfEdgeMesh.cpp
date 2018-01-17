// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
//
// Created by bill on 2017-09-05.
//

#include "halfEdgeMesh.h"
#include <chrono>

float* convert_str(char* str)
{
    char* strtoken;
    strtoken = strtok(str, "/");
    float vec[4];

    int i = 0;
    while (strtoken != NULL)
    {
        vec[i] = stof(strtoken);
        strtoken = strtok (NULL, "/");
        i++;
    }

    return vec;
}

void halfEdgeMesh::parseobj(char *filepath)
{
    FILE* obj = fopen(filepath, "r");

    bool hasNorm;
    bool hasUV;



    int match;

    if (obj)
    {
        std::vector<unsigned int> vertexIndices, uvIndices, normIndices;
        std::vector<fdVector> temp_ver, temp_uv, temp_norm;

        while (1)
        {
            char line[128];
            int res = fscanf(obj, "%s", line);
            if(res == EOF)
                break;
            if(strcmp(line, "v") == 0)
            {
                float f[3];
                fscanf(obj, "%f %f %f\n", &f[0], &f[1], &f[2]);
                fdVector p(f[0], f[1], f[2]);
                temp_ver.push_back(p);
            }
            else if(strcmp(line, "vt") == 0)
            {
                float f[2];
                fscanf(obj, "%f %f %f\n", &f[0], &f[1], &f[2]);
                fdVector p;
                p.set_Vector(f[0], f[1], 0);
                temp_uv.push_back(p);
            }
            else if(strcmp(line, "vn") == 0)
            {
                float f[3];
                fscanf(obj, "%f %f %f\n", &f[0], &f[1], &f[2]);
                fdVector p(f[0], f[1], f[2]);
                temp_norm.push_back(p);
            }
            else if(strcmp(line, "f") == 0)
            {
                nr_of_e++;

                if (temp_uv.size() == 0)
                    hasUV = false;
                else
                    hasUV = true;
                if (temp_norm.size() == 0)
                    hasNorm = false;
                else
                    hasNorm = true;

                unsigned int vertIndex[4], uvIndex[4], normIndex[4];



                if (hasNorm && !hasUV)
                {
                    match = fscanf(obj, " %d%*[/]%d %d%*[/]%d %d%*[/]%d %d%*[/]%d\n", &vertIndex[0], &normIndex[0],
                                   &vertIndex[1],  &normIndex[1],
                                   &vertIndex[2], &normIndex[2],
                                   &vertIndex[3], &normIndex[3]);
                }

                else if (!hasNorm && hasUV)
                {
                    match = fscanf(obj, " %d%*[/]%d%*[/] %d%*[/]%d%*[/] %d%*[/]%d%*[/] %d%*[/]%d%*[/]\n", &vertIndex[0], &uvIndex[0],
                                   &vertIndex[1], &uvIndex[1],
                                   &vertIndex[2], &uvIndex[2],
                                   &vertIndex[3], &uvIndex[3]);
                }
                else if (hasNorm && hasUV)
                {
                    match = fscanf(obj, "%d%*[/]%d%*[/]%d %d%*[/]%d%*[/]%d %d%*[/]%d%*[/]%d %d%*[/]%d%*[/]%d\n",
                                   &vertIndex[0], &uvIndex[0], &normIndex[0], &vertIndex[1], &uvIndex[1], &normIndex[1],
                                   &vertIndex[2], &uvIndex[2], &normIndex[2], &vertIndex[3], &uvIndex[3],
                                   &normIndex[3]);
                } else
                {
                    match = fscanf(obj, "%d %d %d %d\n", &vertIndex[0], &vertIndex[1], &vertIndex[2], &vertIndex[3]);
                }

                switch(match) {
                    case 12:
                    {
                        vertexIndices.push_back(vertIndex[0] - 1);
                        vertexIndices.push_back(vertIndex[1] - 1);
                        vertexIndices.push_back(vertIndex[2] - 1);
                        vertexIndices.push_back(vertIndex[3] - 1);

                        uvIndices.push_back(uvIndex[0] - 1);
                        uvIndices.push_back(uvIndex[1] - 1);
                        uvIndices.push_back(uvIndex[2] - 1);
                        uvIndices.push_back(uvIndex[3] - 1);

                        normIndices.push_back(normIndex[0] - 1);
                        normIndices.push_back(normIndex[1] - 1);
                        normIndices.push_back(normIndex[2] - 1);
                        normIndices.push_back(normIndex[3] - 1);
                        break;
                    }
                    case 9:
                    {
                        vertexIndices.push_back(vertIndex[0] - 1);
                        vertexIndices.push_back(vertIndex[1] - 1);
                        vertexIndices.push_back(vertIndex[2] - 1);

                        uvIndices.push_back(uvIndex[0] - 1);
                        uvIndices.push_back(uvIndex[1] - 1);
                        uvIndices.push_back(uvIndex[2] - 1);

                        normIndices.push_back(normIndex[0] - 1);
                        normIndices.push_back(normIndex[1] - 1);
                        normIndices.push_back(normIndex[2] - 1);
                        break;
                    }
                    case 6:
                    {
                        if (hasNorm)
                        {
                            vertexIndices.push_back(vertIndex[0] - 1);
                            vertexIndices.push_back(vertIndex[1] - 1);
                            vertexIndices.push_back(vertIndex[2] - 1);

                            normIndices.push_back(normIndex[0] - 1);
                            normIndices.push_back(normIndex[1] - 1);
                            normIndices.push_back(normIndex[2] - 1);
                            break;
                        }
                        if (hasUV)
                        {
                            vertexIndices.push_back(vertIndex[0] - 1);
                            vertexIndices.push_back(vertIndex[1] - 1);
                            vertexIndices.push_back(vertIndex[2] - 1);

                            uvIndices.push_back(uvIndex[0] - 1);
                            uvIndices.push_back(uvIndex[1] - 1);
                            uvIndices.push_back(uvIndex[2] - 1);
                            break;
                        }
                    }
                    case 4:
                    {
                        vertexIndices.push_back(vertIndex[0] - 1);
                        vertexIndices.push_back(vertIndex[1] - 1);
                        vertexIndices.push_back(vertIndex[2] - 1);
                        vertexIndices.push_back(vertIndex[3] - 1);
                        break;

                    }
                    default:
                    {
                        vertexIndices.push_back(vertIndex[0] - 1);
                        vertexIndices.push_back(vertIndex[1] - 1);
                        vertexIndices.push_back(vertIndex[2] - 1);
                        break;
                    }
                }
            }
        }

        if (match % 4 == 0) {

            steps = 4;

            int IBlength = vertexIndices.size();
            int UVLength = uvIndices.size();
            int nLength = normIndices.size();
            int vec_length = temp_ver.size();

            int *IB = new int[IBlength];
            float *UVB = new float[UVLength];
            float *NB = new float[nLength];
            verts = new HE_vertex*[vec_length];

            nr_of_f = nr_of_e;
            nr_of_e *= steps;

            faces = new HE_face*[nr_of_f];
            edges = new HE_edge*[nr_of_e];

            for (int i = 0; i < IBlength; i++) {
                IB[i] = vertexIndices[i];
            }
            for (int i = 0; i < UVLength; i++) {
                UVB[i] = uvIndices[i];
            }
            for (int i = 0; i < nLength; i++) {
                NB[i] = normIndices[i];
            }
            for (int i = 0; i < vec_length; ++i) {
                nr_of_v++;
                verts[i] = new HE_vertex;
                verts[i]->pos = temp_ver[i];
            }

            fclose(obj);

            generate_Faces(IBlength, steps, IB, UVB,temp_uv, NB);

            delete[] IB;
            delete[] UVB;
            delete[] NB;
        } else
        {
            steps = 3;

            int IBlength = vertexIndices.size();
            int UVLength = uvIndices.size();
            int nLength = normIndices.size();
            int vec_length = temp_ver.size();

            int *IB = new int[IBlength];
            float *UVB = new float[UVLength];
            float *NB = new float[nLength];
            verts = new HE_vertex*[vec_length];

            nr_of_f = nr_of_e;
            nr_of_e *= steps;

            faces = new HE_face*[nr_of_f];
            edges = new HE_edge*[nr_of_e];

            for (int i = 0; i < IBlength; i++) {
                IB[i] = vertexIndices[i];
            }
            for (int i = 0; i < UVLength; i++) {
                UVB[i] = uvIndices[i];
            }
            for (int i = 0; i < nLength; i++) {
                NB[i] = normIndices[i];
            }
            for (int i = 0; i < vec_length; i++) {
                nr_of_v++;
                verts[i]->pos = temp_ver[i];
            }

            fclose(obj);

            generate_Faces(IBlength, steps, IB, UVB, temp_uv, NB);

            delete[] IB;
            delete[] UVB;
            delete[] NB;
        }

    }
}




void halfEdgeMesh::generate_Faces(int IBlength, int steps, int* indexBuffer, float* uvBuffer,vector<fdVector> uv, float* normBuffer) {
    if (steps == 4){
        int j = 0;
        for (int i = 0; i < IBlength; i += 4)
        {
            faces[j] = new HE_face;

            edges[i] = new HE_edge;
            edges[i + 1] = new HE_edge;
            edges[i + 2] = new HE_edge;
            edges[i + 3] = new HE_edge;


            edges[i]->face = faces[j];
            edges[i]->nxt = edges[i+1];
            edges[i]->vert = verts[indexBuffer[i]];
            edges[i]->uv = uv[uvBuffer[i]];
            verts[indexBuffer[i]]->edge = edges[i];

            edges[i+1]->face = faces[j];
            edges[i+1]->nxt = edges[i+2];
            edges[i+1]->vert = verts[indexBuffer[i + 1]];
            edges[i+1]->uv = uv[uvBuffer[i + 1]];
            verts[indexBuffer[i + 1]]->edge = edges[i + 1];

            edges[i+2]->face = faces[j];
            edges[i+2]->nxt = edges[i+3];
            edges[i+2]->vert = verts[indexBuffer[i + 2]];
            edges[i+2]->uv = uv[uvBuffer[i + 2]];
            verts[indexBuffer[i + 2]]->edge = edges[i + 2];

            edges[i+3]->face = faces[j];
            edges[i+3]->nxt = edges[i];
            edges[i+3]->vert = verts[indexBuffer[i + 3]];
            edges[i+3]->uv = uv[uvBuffer[i + 3]];
            verts[indexBuffer[i + 3]]->edge = edges[i + 3];

            faces[j]->edge = edges[i];

            j++;
        }
    } else
    {
        int j = 0;
        for (int i = 0; i < IBlength; i += 3)
        {
            HE_face *f1 = new HE_face;

            HE_edge* e1 = new HE_edge;
            HE_edge* e2 = new HE_edge;
            HE_edge* e3 = new HE_edge;

            e1->face = f1;
            e1->nxt = e2;
            e1->vert = verts[indexBuffer[i]];

            e2->face = f1;
            e2->nxt = e3;
            e2->vert = verts[indexBuffer[i + 1]];

            e3->face = f1;
            e3->nxt = e1;
            e3->vert = verts[indexBuffer[i + 2]];

            f1->edge = e1;

            faces[j] = f1;

            edges[i] = e1;
            edges[i + 1] = e2;
            edges[i + 2] = e3;

            j++;
        }
    }
    int check = 0;
    for (int i = 0; i < nr_of_e ; i++) {
        for (int j = 0; j < nr_of_e; j++) {
            if (edges[i]->vert->pos == edges[j]->nxt->vert->pos && edges[i]->nxt->vert->pos == edges[j]->vert->pos)
            {

                edges[j]->pair = edges[i];
                edges[i]->pair = edges[j];
                check++;
            }
        }
    }
}



float* halfEdgeMesh::get_vertices()
{
    float *VB;
    if (steps == 3) {
        VB = new float[nr_of_e * 3];

        int j = 0;
        bool f = true;

        for (int i = 0; i < nr_of_f; i++) {
            HE_edge *cur_e = faces[i]->edge;
            fdVector current = faces[i]->edge->vert->pos;
            fdVector first;
            f = true;
            while (current != first) {
                if (f) {
                    first = current;
                    f = false;
                }
                VB[j] = current.get_X() * 0.2;
                VB[j + 1] = current.get_Y() * 0.2;
                VB[j + 2] = current.get_Z() * 0.2;
                cur_e = cur_e->nxt;
                current = cur_e->vert->pos;
                j += 3;
            }
        }
    } else {
        VB = new float[nr_of_e * 4];

        int j = 0;
        bool f = true;

        for (int i = 0; i < nr_of_f; i++) {
            HE_edge *cur_e = faces[i]->edge;
            fdVector current = faces[i]->edge->vert->pos;
            fdVector first;
            f = true;
            while (!((fabs(current.get_X() - first.get_X()) <= 0.0001) && (fabs(current.get_Y() - first.get_Y()) <= 0.0001) && (fabs(current.get_Z() - first.get_Z()) <= 0.0001))) {
                if (f) {
                    first = current;
                    f = false;
                }
                VB[j] = current.get_X()*0.6;
                VB[j + 1] = current.get_Y()*0.6;
                VB[j + 2] = current.get_Z()*0.6;
                cur_e = cur_e->nxt;
                current = cur_e->vert->pos;
                j += 3;
            }
        }
    }
    return VB;
};


float* halfEdgeMesh::get_uvs()
{
    float *VB;
    if (steps == 3) {
        VB = new float[nr_of_e * 3];

        int j = 0;
        bool f = true;

        for (int i = 0; i < nr_of_f; i++) {
            HE_edge *cur_e = faces[i]->edge;
            fdVector current = faces[i]->edge->vert->pos;
            fdVector first;
            f = true;
            while (current != first) {
                if (f) {
                    first = current;
                    f = false;
                }
                VB[j] = current.get_X()*0.6;
                VB[j + 1] = current.get_Y();
                VB[j + 2] = current.get_Z();
                cur_e = cur_e->nxt;
                current = cur_e->vert->pos;
                j += 3;
            }
        }
    } else {
        VB = new float[nr_of_e * 2];

        int j = 0;
        bool f = true;

        for (int i = 0; i < nr_of_f; i++) {
            HE_edge *cur_e = faces[i]->edge;
            fdVector current = faces[i]->edge->vert->pos;
            fdVector first(9999999,9999999,9999999,1);
            f = true;
            while (current != first) {
                if (f) {
                    first = current;
                    f = false;
                }
                VB[j] = cur_e->uv.get_X();
                VB[j + 1] = cur_e->uv.get_Y();
                cur_e = cur_e->nxt;
                current = cur_e->vert->pos;
                j += 2;
            }
        }
    }
    return VB;
};


void halfEdgeMesh::cmc_subdivide()
{
    nrofsub++;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::time_point<std::chrono::system_clock> startf, endf;
    std::chrono::time_point<std::chrono::system_clock> startv, endv;
    std::chrono::time_point<std::chrono::system_clock> starte, ende;
    std::chrono::time_point<std::chrono::system_clock> startd, endd;
    start = std::chrono::system_clock::now();

    HE_vertex** newvertices = new HE_vertex*[nr_of_f+nr_of_e/2+nr_of_v];

    int j = 0;

    startf = std::chrono::system_clock::now();

    for (int i = 0; i < nr_of_f; i++)
    {
        bool f = true;

        fdVector temp_f;

        HE_edge* cur_e = nullptr;
        int bajs = 0;
        while (cur_e != faces[i]->edge) {
            if (f)
            {
                f = false;
                cur_e  = faces[i]->edge;
            }
            temp_f = temp_f + cur_e->vert->pos;
            bajs++;
            cur_e = cur_e->nxt;
        }

        newvertices[j] = new HE_vertex;
        newvertices[j]->pos = temp_f * 0.25;
        newvertices[j]->idx = j;

        faces[i]->idx = j;

        j++;
    }
    endf = std::chrono::system_clock::now();
    double elapsed_seconds4 = (double)std::chrono::duration_cast<std::chrono::milliseconds>(endf-startf).count();

    std::cout << "face fixing " << nrofsub <<" took: " << elapsed_seconds4*0.001 << " seconds!"<< '\n';


    starte = std::chrono::system_clock::now();
    for (int i = 0; i < nr_of_e; i++)
    {
        fdVector ep = (edges[i]->vert->pos + edges[i]->nxt->vert->pos + newvertices[edges[i]->face->idx]->pos + newvertices[edges[i]->pair->face->idx]->pos)*0.25;

        if (edges[i]->pair->idx < i+nr_of_f)
        {
            edges[i]->idx = edges[i]->pair->idx;

        } else
        {
            newvertices[j] = new HE_vertex;
            newvertices[j]->pos = ep;
            newvertices[j]->idx = j;
            edges[i]->idx = j;
            j++;
        }

    }

    int n = 0;

    ende = std::chrono::system_clock::now();
    double elapsed_seconds3 = (double)std::chrono::duration_cast<std::chrono::milliseconds>(ende-starte).count();

    std::cout << "edge fixing " << nrofsub <<" took: " << elapsed_seconds3*0.001 << " seconds!"<< '\n';

    startv = std::chrono::system_clock::now();
    for (int i = 0; i < nr_of_v; i++)
    {
        n = 0;
        float rx = 0, ry = 0, rz = 0;
        float finX = 0, finY = 0, finZ = 0;
        float fx = 0, fy = 0, fz = 0;
        fdVector temp = verts[i]->pos;

        bool f = true;

        HE_edge* cur_e = nullptr;

        while (cur_e != verts[i]->edge) {
            if (f)
            {
                f = false;
                cur_e  = verts[i]->edge;
            }
            rx = rx + (cur_e->vert->pos.get_X()+cur_e->nxt->vert->pos.get_X())/2;
            ry = ry + (cur_e->vert->pos.get_Y()+cur_e->nxt->vert->pos.get_Y())/2;
            rz = rz + (cur_e->vert->pos.get_Z()+cur_e->nxt->vert->pos.get_Z())/2;

            fx = fx + newvertices[cur_e->face->idx]->pos.get_X();
            fy = fy + newvertices[cur_e->face->idx]->pos.get_Y();
            fz = fz + newvertices[cur_e->face->idx]->pos.get_Z();
            n++;
            cur_e = cur_e->pair->nxt;
        }

        rx /= n;
        ry /= n;
        rz /= n;
        fx /= n;
        fy /= n;
        fz /= n;

        finX = (fx + 2*rx + (n - 3)*temp.get_X())/n;
        finY = (fy + 2*ry + (n - 3)*temp.get_Y())/n;
        finZ = (fz + 2*rz + (n - 3)*temp.get_Z())/n;
        fdVector fin(finX, finY, finZ);

        newvertices[j] = verts[i];
        newvertices[j]->pos = fin;
        newvertices[j]->idx = j;

        j++;
    }

    endv = std::chrono::system_clock::now();
    double elapsed_seconds2 = (double)std::chrono::duration_cast<std::chrono::milliseconds>(endv-startv).count();

    std::cout << "vertice fixing " << nrofsub <<" took: " << elapsed_seconds2*0.001 << " seconds!"<< '\n';

    HE_edge** newedges = new HE_edge*[nr_of_f*4*4];
    HE_face** newfaces = new HE_face*[nr_of_f*4];

    int loops = 0;
    int e = 0;

    startd = std::chrono::system_clock::now();

    HE_edge* temp;
    HE_edge* temp2;
    HE_edge* temp3;
    HE_edge* temp4;

    for (int l = 0; l < nr_of_f; ++l) {
        temp = faces[l]->edge;
        temp2 = faces[l]->edge->nxt;
        temp3 = faces[l]->edge->nxt->nxt;
        temp4 = faces[l]->edge->nxt->nxt->nxt;

        newfaces[loops] = faces[l];
        newfaces[loops + 1] = new HE_face;
        newfaces[loops + 2] = new HE_face;
        newfaces[loops + 3] = new HE_face;

        newedges[e] = new HE_edge;
        newedges[e + 1] = new HE_edge;

        newedges[e + 2] = new HE_edge;
        newedges[e + 3] = new HE_edge;

        newedges[e + 4] = new HE_edge;
        newedges[e + 5] = new HE_edge;

        newedges[e + 6] = new HE_edge;
        newedges[e + 7] = new HE_edge;

        newedges[e + 8] = new HE_edge;
        newedges[e + 9] = new HE_edge;
        newedges[e + 10] = new HE_edge;
        newedges[e + 11] = new HE_edge;
        newedges[e + 12] = new HE_edge;
        newedges[e + 13] = new HE_edge;
        newedges[e + 14] = new HE_edge;
        newedges[e + 15] = new HE_edge;

////----First Edge----////
        newedges[e]->vert = newvertices[temp->vert->idx];
        newedges[e]->uv = temp->uv;
        newvertices[temp->vert->idx]->edge = newedges[e];
        if (temp->pair->split)
        {
            newedges[e]->pair = newedges[temp->pair->splitidx]->nxt->pair->nxt;
            newedges[temp->pair->splitidx]->nxt->pair->nxt->pair = newedges[e];
            newedges[e + 1]->pair = newedges[temp->pair->splitidx];
            newedges[temp->pair->splitidx]->pair = newedges[e + 1];
        }

        newedges[e + 1]->vert = newvertices[temp->idx];
        newedges[e + 1]->nxt = newedges[e + 2];
        newedges[e + 1]->uv = (temp->uv + temp->nxt->uv)*0.5;
        newvertices[temp->idx]->edge = newedges[e + 1];


        temp->split = true;
        temp->splitidx = e;
        newedges[e+1]->splitidx = e;

////----Second Edge----////
        newedges[e + 2]->vert = newvertices[temp2->vert->idx];
        newedges[e + 2]->uv = temp2->uv;
        newvertices[temp2->vert->idx]->edge = newedges[e + 2];
        if (temp2->pair->split)
        {
            newedges[e + 2]->pair = newedges[temp2->pair->splitidx]->nxt->pair->nxt;
            newedges[temp2->pair->splitidx]->nxt->pair->nxt->pair = newedges[e + 2];
            newedges[e + 3]->pair = newedges[temp2->pair->splitidx];
            newedges[temp2->pair->splitidx]->pair = newedges[e + 3];
        }

        newedges[e + 3]->vert = newvertices[temp2->idx];
        newedges[e + 3]->nxt = newedges[e + 4];
        newedges[e + 3]->uv = (temp2->uv + temp2->nxt->uv)*0.5;
        newvertices[temp2->idx]->edge = newedges[e + 3];


        temp2->split = true;
        temp2->splitidx = e+2;
        newedges[e+3]->splitidx = e+2;

        ////----Third Edge----////
        newedges[e + 4]->vert = newvertices[temp3->vert->idx];
        newedges[e + 4]->uv = temp3->uv;
        newvertices[temp3->vert->idx]->edge = newedges[e + 4];
        if (temp3->pair->split)
        {
            newedges[e + 4]->pair = newedges[temp3->pair->splitidx]->nxt->pair->nxt;
            newedges[temp3->pair->splitidx]->nxt->pair->nxt->pair = newedges[e + 4];
            newedges[e + 5]->pair = newedges[temp3->pair->splitidx];
            newedges[temp3->pair->splitidx]->pair = newedges[e + 5];
        }

        newedges[e + 5]->vert = newvertices[temp3->idx];
        newedges[e + 5]->nxt = newedges[e + 6];
        newedges[e + 5]->uv = (temp3->uv + temp3->nxt->uv)*0.5;
        newvertices[temp3->idx]->edge = newedges[e + 5];


        temp3->split = true;
        temp3->splitidx = e+4;
        newedges[e+5]->splitidx = e+4;

        ////----Fourth Edge----////
        newedges[e + 6]->vert = newvertices[temp4->vert->idx];
        newedges[e + 6]->nxt = newedges[e + 7];
        newedges[e + 6]->uv = temp4->uv;
        newvertices[temp4->vert->idx]->edge = newedges[e + 6];
        if (temp4->pair->split)
        {
            newedges[e + 6]->pair = newedges[temp4->pair->splitidx]->nxt->pair->nxt;
            newedges[temp4->pair->splitidx]->nxt->pair->nxt->pair = newedges[e + 6];
            newedges[e + 7]->pair = newedges[temp4->pair->splitidx];
            newedges[temp4->pair->splitidx]->pair = newedges[e + 7];
        }

        newedges[e + 7]->vert = newvertices[temp4->idx];
        newedges[e + 7]->nxt = newedges[e];
        newedges[e + 7]->uv = (temp4->uv + temp4->nxt->uv)*0.5;
        newvertices[temp4->idx]->edge = newedges[e + 7];


        temp4->split = true;
        temp4->splitidx = e+6;
        newedges[e+7]->splitidx = e+6;

        newedges[e]->nxt = newedges[e + 8];
        newedges[e]->face = newfaces[loops];
        newfaces[loops]->edge = newedges[e];

        newedges[e + 8]->vert = newvertices[temp->idx];
        newedges[e + 8]->face = newfaces[loops];
        newedges[e + 8]->nxt = newedges[e + 9];
        newedges[e + 8]->pair = newedges[e + 11];
        newedges[e + 8]->uv = (temp->uv + temp->nxt->uv)*0.5;

        newedges[e + 9]->vert = newvertices[temp->face->idx];
        newedges[e + 9]->face = newfaces[loops];
        newedges[e + 9]->nxt = newedges[e + 7];
        newedges[e + 9]->pair = newedges[e + 14];
        newedges[e + 9]->uv = (temp->uv + temp2->uv + temp3->uv + temp4->uv) * 0.25;
        newvertices[temp->face->idx]->edge = newedges[e + 9];

        newedges[e + 7]->face = newfaces[loops];

        newedges[e + 2]->nxt = newedges[e + 10];
        newedges[e + 2]->face = newfaces[loops + 1];
        newfaces[loops + 1]->edge = newedges[e + 2];

        newedges[e + 10]->vert = newvertices[temp2->idx];
        newedges[e + 10]->face = newfaces[loops + 1];
        newedges[e + 10]->nxt = newedges[e + 11];
        newedges[e + 10]->pair = newedges[e + 13];
        newedges[e + 10]->uv = (temp2->uv + temp3->uv)*0.5;

        newedges[e + 11]->vert = newvertices[temp->face->idx];
        newedges[e + 11]->face = newfaces[loops + 1];
        newedges[e + 11]->nxt = newedges[e + 1];
        newedges[e + 11]->pair = newedges[e + 8];
        newedges[e + 11]->uv = (temp->uv + temp2->uv + temp3->uv + temp4->uv) * 0.25;

        newedges[e + 1]->face = newfaces[loops + 1];

        newedges[e + 4]->nxt = newedges[e + 12];
        newedges[e + 4]->face = newfaces[loops + 2];
        newfaces[loops + 2]->edge = newedges[e + 4];

        newedges[e + 12]->vert = newvertices[temp3->idx];
        newedges[e + 12]->face = newfaces[loops + 2];
        newedges[e + 12]->nxt = newedges[e + 13];
        newedges[e + 12]->pair = newedges[e + 15];
        newedges[e + 12]->uv = (temp3->uv + temp4->uv)*0.5;

        newedges[e + 13]->vert = newvertices[temp->face->idx];
        newedges[e + 13]->face = newfaces[loops + 2];
        newedges[e + 13]->nxt = newedges[e + 3];
        newedges[e + 13]->pair = newedges[e + 10];
        newedges[e + 13]->uv = (temp->uv + temp2->uv + temp3->uv + temp4->uv) * 0.25;

        newedges[e + 3]->face = newfaces[loops + 2];

        newedges[e + 6]->nxt = newedges[e + 14];
        newedges[e + 6]->face = newfaces[loops + 3];
        newfaces[loops + 3]->edge = newedges[e + 6];

        newedges[e + 14]->vert = newvertices[temp4->idx];
        newedges[e + 14]->face = newfaces[loops + 3];
        newedges[e + 14]->nxt = newedges[e + 15];
        newedges[e + 14]->pair = newedges[e + 9];
        newedges[e + 14]->uv = (temp4->uv + temp->uv)*0.5;

        newedges[e + 15]->vert = newvertices[temp->face->idx];
        newedges[e + 15]->face = newfaces[loops + 3];
        newedges[e + 15]->nxt = newedges[e + 5];
        newedges[e + 15]->pair = newedges[e + 12];
        newedges[e + 15]->uv = (temp->uv + temp2->uv + temp3->uv + temp4->uv) * 0.25;

        newedges[e + 5]->face = newfaces[loops + 3];

        newfaces[loops]->idx = loops;
        newfaces[loops + 1]->idx = loops + 1;
        newfaces[loops + 2]->idx = loops + 2;
        newfaces[loops + 3]->idx = loops + 3;

        e += 16;
        loops += 4;
    }

    endd = std::chrono::system_clock::now();
    double elapsed_seconds5 = (double)std::chrono::duration_cast<std::chrono::milliseconds>(endd-startd).count();

    std::cout << "Assembly " << nrofsub <<" took: " << elapsed_seconds5*0.001 << " seconds!"<< '\n';

    for (int i = 0; i < nr_of_e; i++)
    {
        delete edges[i];
    }


    delete[] verts;
    delete[] edges;
    delete[] faces;

    nr_of_v = nr_of_f+nr_of_e/2+nr_of_v;
    nr_of_e = nr_of_f*4*4;
    nr_of_f = nr_of_f*4;

    edges = newedges;
    faces = newfaces;
    verts = newvertices;
    end = std::chrono::system_clock::now();
    double elapsed_seconds = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();


    std::cout << "subdivision " << nrofsub <<" took: " << elapsed_seconds*0.001 << " seconds!"<< '\n' << endl;
}

