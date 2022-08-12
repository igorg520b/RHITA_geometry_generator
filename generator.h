#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include "mesh2d.h"


class Generator
{
public:
    Generator();

    std::string outputFileName = "default.msh";
    double elemSize = 0.05;
    double notchOffset = 1.1;
    double indenterRadius = 0.161925;
    double indentationDepth = 0.05;
    bool insertCZs = false;

    constexpr static double blockHeight = 1.;
    constexpr static double blockLength = 2.5;

    icy::Mesh2D mesh2d;

    void Generate();
};

#endif // GENERATOR_H
