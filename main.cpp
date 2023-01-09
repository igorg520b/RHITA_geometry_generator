#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>
#include <QFileInfo>

#include "gmsh.h"
#include "generator.h"
#include "generator3d.h"


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
                                        QCoreApplication::translate("main", "Output file name"), "outputFile");

    QCommandLineOption loadFileOption(QStringList() << "m" << "meshfile",
                                         QCoreApplication::translate("main", "Load from MSH file"), "FileName");

    QCommandLineOption insertCZSOption(QStringList() << "czs", QCoreApplication::translate("main", "InsertCZS"));
    QCommandLineOption createCDPOption(QStringList() << "cdp", QCoreApplication::translate("main", "Create CDP material"));
    QCommandLineOption createSCollisionsOption(QStringList() << "sCollisions", QCoreApplication::translate("main", "Create self collisions"));
    QCommandLineOption use3DGeneratorOption(QStringList() << "3d", QCoreApplication::translate("main", "Use 3D Generator"));
    QCommandLineOption attachTopOption(QStringList() << "attachTop", QCoreApplication::translate("main", "Attach top side of the block"));

    parser.addOption(outputFileOption);
    parser.addOption(insertCZSOption);
    parser.addOption(loadFileOption);
    parser.addOption(createCDPOption);
    parser.addOption(createSCollisionsOption);
    parser.addOption(use3DGeneratorOption);
    parser.addOption(attachTopOption);

    // -o x500 -m x500.msh --cdp --czs --3d
    // -o t300 -m t300.msh --cdp --czs --sCollisions

    // Process the actual command line arguments given by the user
    parser.process(a);

    if(parser.isSet(use3DGeneratorOption))
    {
        // 3D geometry
        Generator3D g3d;
        if(parser.isSet(outputFileOption)) g3d.outputFileName = parser.value(outputFileOption).toStdString();
        if(parser.isSet(insertCZSOption)) g3d.insertCZs = true;
        g3d.LoadFromFile(parser.value(loadFileOption).toStdString());
    }
    else
    {
        // 2D
        Generator g;
        if(parser.isSet(outputFileOption)) g.outputFileName = parser.value(outputFileOption).toStdString();
        if(parser.isSet(insertCZSOption)) g.insertCZs = true;
        if(parser.isSet(createCDPOption)) g.createCDP = true;
        if(parser.isSet(createSCollisionsOption)) g.createSelfCollisions = true;
        if(parser.isSet(attachTopOption)) g.attachTop = true;
        if(parser.isSet(loadFileOption))
            g.LoadFromFile(parser.value(loadFileOption).toStdString());
        else
            g.Generate();
    }



}
