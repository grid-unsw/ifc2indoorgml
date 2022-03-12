#ifndef IFC2LCC_H
#define IFC2LCC_H

#include <QTime>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QStringList>
#include <QFileDialog>

#include "MainWindow.h"
#include "IFC_loader_ui.h"
//#include "IFC_Importer.h"

///
/// \brief load_IFC
/// \callgraph
void load_IFC(const QString &, LCC &);

#endif
