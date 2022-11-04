#ifndef GENERATORALT1_H
#define GENERATORALT1_H

#include <string>
#include "mesh2d.h"
#include "cohesivezone2d.h"

class GeneratorAlt1
{
public:
    std::string outputFileName = "default";

    // block geometry
    constexpr static double blockHeight = 1.;
    constexpr static double blockLength = 2.5;

    // model settings
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 6e6;
    constexpr static double czElasticity = 1e11;
    constexpr static double czEnergy = 500;
    constexpr static double indentationRate = 0.2;

    // other parameters
    constexpr static double indenterRadius = 0.161925;
    constexpr static double indentationDepth = 0.05;
    constexpr static double elemSize = 0.03;
    constexpr static double interactionRadius = 0.01; // for exponential interaction property
    constexpr static double timeToRun = 10;
    constexpr static int nFrames = 2000;

    constexpr static int numberOfCores = 12;

    icy::Mesh2D mesh2d;
    void Generate();
    void CreateCDP(std::ofstream &s);
    void CreatePy();

};

#endif // GENERATORALT1_H
