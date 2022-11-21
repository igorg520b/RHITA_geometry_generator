#ifndef GENERATOR3D_H
#define GENERATOR3D_H


#include <string>
#include "mesh.h"

class Generator3D
{
public:
    std::string outputFileName = "default";
    double indenterRadius = 0.161925;
    double indentationDepth = 0.05;

    double interactionRadius = 0.01; // for exponential interaction property

    double timeToRun = 10;
    int nFrames = 2000;

    bool loadWithIndenter = true;   // if false -> static load

    double blockHeight = 1.;
    double blockLength = 2.5;
    double blockWidth = 1.5;
    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 6e6;
    constexpr static double czElasticity = 1e12;
    constexpr static double czEnergy = 30;
    constexpr static double indentationRate = 0.2;

    bool insertCZs = false;

    icy::Mesh mesh;

    void LoadFromFile(std::string MSHFileName);
    void CreatePy();
    void CreateCDP(std::ofstream &s);
};

#endif // GENERATOR3D_H
