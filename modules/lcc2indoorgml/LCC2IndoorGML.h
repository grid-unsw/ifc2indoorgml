#ifndef LCC2INDOORGML_H
#define LCC2INDOORGML_H

#include <QTime>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QStringList>
#include <QFileDialog>

//#include "MainWindow.h"
#include "IndoorGML_writer.h"
#include "LCC2IndoorGML_Ops.h"
#include "UI/addthematiclayer.h"

extern std::map<str, Dart_handle> ori3CellIDs2IndoorGMLIDs;
extern std::map<str, Dart_handle> cellspace_dart, cellboundary_dart;

class LCC2IndoorGML
{

public:

    // Constructors
    LCC2IndoorGML():
        isReady(false)
    {}

    bool myIndoorGMLisEmpty()
    { return (indoorGML1.getprimalSpaceFeatures()->getCellSpaceMember())->size() == 0; }
    void setIndoorFeatures( IndoorFeatures InFt )
    { indoorGML1 = InFt; }
    IndoorFeatures *getIndoorFeatures()
    { return &indoorGML1; }

    bool myIndoorGML2isEmpty()
    { return (indoorGML2.getThematicLayerMap()->size() == 0
              && indoorGML2.getInterLayerConnectionMap()->size() == 0 ); }
    void setIndoorFeatures2( IndoorGML2::IndoorFeatures InFt )
    { indoorGML2 = InFt; }
    IndoorGML2::IndoorFeatures *getIndoorFeatures2()
    { return &indoorGML2; }

    void create_New_ThematicLayer(LCC &, std::string, bool visibleCellsOnly = false);
    AddThematicLayer* newTLdialog;

    void add_ThematicLayer(LCC &);
    void read_IndoorGML(LCC &, MainWindow *);
    void write_IndoorGML1(LCC &, MainWindow *);
    void write_IndoorGML2(LCC &, MainWindow *);
    void load_LCC_to_IndoorGML1(LCC &);
    void load_LCC_to_IndoorGML2(LCC &);
    void load_IndoorGML(LCC &, str, int v=1);
    bool save_IndoorGML(LCC &, str, int v=1);

    bool isReady;
    void init(LCC&);

    IndoorGML::IndoorFeatures* getIndoorGML1(){
        return &indoorGML1;
    }

    IndoorGML2::IndoorFeatures* getIndoorGML2(){
        return &indoorGML2;
    }

private:
    IndoorGML::IndoorFeatures indoorGML1;
    IndoorGML2::IndoorFeatures indoorGML2;
};

extern LCC2IndoorGML lcc2igmlHandler;

std::size_t generateNetworkFromLCC(LCC&, MainWindow*, bool use3Links=true);

#endif // LCC2INDOORGML_H
