#ifndef GENERATOR3D_H
#define GENERATOR3D_H


#include <string>
#include <unordered_set>
#include <algorithm>
#include "mesh.h"
#include "node.h"
#include "element.h"


class Generator3D
{
public:
    std::string outputFileName = "default";
    constexpr static double indenterRadius = 0.161925;
    constexpr static double indentationDepth = 0.05;

    double interactionRadius = 0.01; // for exponential interaction property

    double timeToRun = 10;
    int nFrames = 2000;

    bool loadWithIndenter = true;   // if false -> static load

    double blockHeight = 1.;
    double blockLength = 2.5;
    double blockWidth = 1.5;
    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e6;
    constexpr static double czElasticity = 1e11;
    constexpr static double czEnergy = 100;
    constexpr static double indentationRate = 0.2;

    bool insertCZs = false;
    bool isHalfSphere = false; // this is a different geometry - attach on top
    bool makeCutout = true; // remove elements to create an "indentation"
    constexpr static double cutoutX = 1;
    constexpr static double deltaX = 0.02;  // depth of layer for static force application

    icy::Mesh mesh;

    void LoadFromFile(std::string MSHFileName);
    void LoadFromFileWithCrop(std::string MSHFileName);
    void CreatePy();
    void CreateCDP(std::ofstream &s);
};

#endif // GENERATOR3D_H
