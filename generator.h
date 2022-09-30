#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include "mesh2d.h"


class Generator
{
public:

    std::string outputFileName = "default";
    double elemSize = 0.025;
    double notchOffset = 1.1;
    double indenterRadius = 0.161925;
    double indentationDepth = 0.05;
    double indenterOffset = 0;

    double interactionRadius = 0.01; // for exponential interaction property

    bool insertCZs = false;
    bool plasticity = false;
    bool createCDP = false;

    double timeToRun = 10;
    int nFrames = 2000;

    bool loadWithIndenter = true;   // if false -> static load

    double blockHeight = 1.;
    double blockLength = 2.5;
    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e5;
    constexpr static double czElasticity = 5e9;
    constexpr static double czEnergy = 50;
    constexpr static double indentationRate = 0.2;

    icy::Mesh2D mesh2d;

    void Generate();
    void LoadFromFile(std::string MSHFileName);
    void CreatePyWithIndenter2D();
    void CreateCDP(std::ofstream &s);
};

#endif // GENERATOR_H
