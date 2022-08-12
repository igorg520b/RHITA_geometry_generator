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
                                        QCoreApplication::translate("main", "Output file name"), "default.py");
    parser.addOption(outputFileOption);


    QCommandLineOption indenterOption(QStringList() << "i" << "ind",
                                      QCoreApplication::translate("main", "Indenter radius"), "0.161925");
    parser.addOption(indenterOption);

    QCommandLineOption elemSizeOption(QStringList() << "e" << "elemsize",
                                      QCoreApplication::translate("main", "Element size"), "0.05");
    parser.addOption(elemSizeOption);

    QCommandLineOption notchOffsetOption(QStringList() << "n" << "notch",
                                         QCoreApplication::translate("main", "Notch offset"), "1.1");
    parser.addOption(notchOffsetOption);

    QCommandLineOption indentationDepthOption(QStringList() << "d" << "depth",
                                         QCoreApplication::translate("main", "Indentation depth"), "0.05");
    parser.addOption(indentationDepthOption);

    QCommandLineOption insertCZSOption(QStringList() << "czs",
                                         QCoreApplication::translate("main", "InsertCZS"));
    parser.addOption(insertCZSOption);

    // -o test.msh -i 0.161925 -e 0.1 -n 1.1

    // Process the actual command line arguments given by the user
    parser.process(a);

    Generator g;


    if(parser.isSet(outputFileOption)) {
        g.outputFileName = parser.value(outputFileOption).toStdString();
    }
    if(parser.isSet(indenterOption)) {
        g.indenterRadius = parser.value(indenterOption).toDouble();
    }
    if(parser.isSet(elemSizeOption))
    {
        g.elemSize = parser.value(elemSizeOption).toDouble();
    }

    if(parser.isSet(notchOffsetOption)) g.notchOffset = parser.value(notchOffsetOption).toDouble();
    if(parser.isSet(indentationDepthOption)) g.indentationDepth = parser.value(indentationDepthOption).toDouble();
    if(parser.isSet(insertCZSOption)) g.insertCZs = true;


    g.Generate();

}
