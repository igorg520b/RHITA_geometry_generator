#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>

class Generator
{
public:
    Generator();

    std::string outputFileName;
    double elemSize;
    double notchOffset;
    double indenterRadius;

    void Generate();
};

#endif // GENERATOR_H
