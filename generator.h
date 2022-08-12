#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include "mesh2d.h"

#include "cohesivezone2d.h"

class Generator
{
public:
    Generator();

    std::string outputFileName = "default";
    double elemSize = 0.05;
    double notchOffset = 1.1;
    double indenterRadius = 0.161925;
    double indentationDepth = 0.05;
    bool insertCZs = false;

    constexpr static double blockHeight = 1.;
    constexpr static double blockLength = 2.5;
    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e5;
    constexpr static double czElasticity = 5e9;
    constexpr static double czEnergy = 50;
    constexpr static double indentationRate = 0.2;

    icy::Mesh2D mesh2d;

    void Generate();
};

#endif // GENERATOR_H
