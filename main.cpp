// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
#include <iostream>
#include "memHalfEdgeMesh.h"

int main(int argc, char* args[]) {
    memhalfEdgeMesh mesh;
    mesh.parseobj("/home/bill/wilenn-5/S0017D/Labb3/q.obj");
    int t;
    for (int j = 0; j < argc; ++j) {
        if (atoi(args[j]))
            t = atoi(args[j]);
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    for (int i = 0; i < 10; ++i)
    {
       mesh.cmc_subdivide();
    }


    end = std::chrono::system_clock::now();
    double elapsed_seconds = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();


    std::cout << "subdivisions took: " << elapsed_seconds*0.001 << " seconds!"<< '\n' << endl;

    return 0;
}