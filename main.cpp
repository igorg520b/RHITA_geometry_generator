#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>
#include <QFileInfo>

#include "gmsh.h"
#include "generator.h"


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QCoreApplication::setApplicationName("RHITA geometry generator");
    QCoreApplication::setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Batch generator of Abaqus .inp files");
    parser.addHelpOption();
    parser.addVersionOption();
    //    parser.addPositionalArgument("source", QCoreApplication::translate("main", "Configuration file."));

    QCommandLineOption outputFileOption(QStringList() << "o" << "out",
                                        QCoreApplication::translate("main", "Output file name"));
    parser.addOption(outputFileOption);


    QCommandLineOption indenterOption(QStringList() << "i" << "ind",
                                      QCoreApplication::translate("main", "Indenter radius"));
    parser.addOption(indenterOption);

    QCommandLineOption elemSizeOption(QStringList() << "e" << "elemsize",
                                      QCoreApplication::translate("main", "Element size"));
    parser.addOption(elemSizeOption);

    QCommandLineOption notchOffsetOption(QStringList() << "n" << "notch",
                                         QCoreApplication::translate("main", "Notch offset"));
    parser.addOption(notchOffsetOption);

    // -o test.msh -i 0.161925 -e 0.1 -n 1.1

    // Process the actual command line arguments given by the user
    parser.process(a);

    std::string outputFileName = "default.msh";
    double elemSize = 0.05;
    double notchOffset = 1.1;
    double indenterRadius = 0.161925;

    if(parser.isSet(outputFileOption)) outputFileName = parser.value(outputFileOption).toStdString();
    if(parser.isSet(indenterOption)) indenterRadius = parser.value(indenterOption).toDouble();
    if(parser.isSet(elemSizeOption)) elemSize = parser.value(elemSizeOption).toDouble();
    if(parser.isSet(notchOffsetOption)) notchOffset = parser.value(notchOffsetOption).toDouble();

    Generator g;
    g.elemSize = elemSize;
    g.indenterRadius = indenterRadius;
    g.notchOffset = notchOffset;
    g.outputFileName = outputFileName;

    gmsh::initialize();
    g.Generate();

}
