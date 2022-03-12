// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.3/Linear_cell_complex/demo/Linear_cell_complex/MainWindow.h $
// $Id: MainWindow.h cc99fd9 2021-02-19T16:02:12+01:00 Maxime Gimeno
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//                 Abdoulaye Diakite <diakite.abdoulaye@gmail.com>
//
#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"

#include <CGAL/Qt/DemosMainWindow.h>

#include <QDialog>
#include <QSlider>
#include <QLabel>
#include <QFileDialog>

#include <QDockWidget>
#include <QTableWidget>
#include <QCheckBox>

class QWidget;

template < class First, class Second, class Third > struct Triplet
{
  First first;
  Second second;
  Third third;

  Triplet(First first, Second second, Third third)
  {
    this->first = first;
    this->second = second;
    this->third = third;
  }

  Triplet()
  {}
};

class MainWindow : public CGAL::Qt::DemosMainWindow, public Ui::MainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* parent = nullptr);

public Q_SLOTS:
  // File menu
  void on_actionSave_triggered();
  void on_actionLoad_triggered();
  void on_actionClear_triggered();

  // IndoorGML menu
  void on_actionGenerate_IndoorGML_triggered();
  void on_actionExport_IndoorGML_v1_triggered();
  void on_actionExport_IndoorGML_v2_triggered();

  // Operations menu
  void on_actionSew3_same_facets_triggered();
  void on_actionUnsew3_all_triggered();
  void on_actionInsideOut_triggered();
  void on_actionMerge_coplanar_faces_triggered();
  void on_actionRemove_filled_volumes_triggered();
  void on_actionTriangulate_all_facets_triggered();

  // View menu
  void on_actionReset_scene_view_triggered();

  // Other slots
  void load_depend_on_extension(const QString& fileName, bool clear=false);
  void load(const QString& fileName, bool clear=false);
  void save(const QString& fileName);
  void load_off(const QString& fileName, bool clear=false);
  void load_ifc(const QString& fileName, bool clear=false);

  void onSceneChanged();

  void connectVolumeListHandlers();
  void onCellChanged(int, int);
  void onHeaderClicked(int);

  void open(QString);

Q_SIGNALS:
  void sceneChanged();

protected:
  void clear_all();
  void on_new_volume(Dart_handle adart);
  void on_delete_volume(Dart_handle adart);
  void init_all_new_volumes();
  void init_ifc_volumes();
  void mark_all_filled_and_visible_volumes(LCC::size_type amark);

  void connect_actions();
  void update_operations_entries(bool show);

  bool is_volume_in_list(LCC::Attribute_handle<3>::type ah);
  void recreate_whole_volume_list();
  void update_volume_list_all_ckeckstates();
  void update_volume_list_add(LCC::Attribute_handle<3>::type ah);
  void update_volume_list_remove(int);
  void update_volume_list_remove(LCC::Attribute_handle<3>::type ah);

  Scene scene;

  QLabel*      statusMessage;

  QDockWidget* volumeListDock;
  QTableWidget* volumeList;
};

#endif
