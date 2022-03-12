// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.3/Linear_cell_complex/demo/Linear_cell_complex/Viewer.h $
// $Id: Viewer.h fb6f703 2021-05-04T14:07:49+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//                 Abdoulaye Diakite <diakite.abdoulaye@gmail.com>
#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <CGAL/draw_linear_cell_complex.h>

// Functor used by SimpleLCCViewerQt to colorize of not elements.
struct MyDrawingFunctorLCC
{
  /// @return true iff the volume containing dh is drawn.
  template<typename LCC>
  bool draw_volume(const LCC& alcc,
                   typename LCC::Dart_const_handle dh) const
  { return alcc.template info<3>(dh).is_visible(); }
  /// @return true iff the face containing dh is drawn.
  template<typename LCC>
  bool draw_face(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the edge containing dh is drawn.
  template<typename LCC>
  bool draw_edge(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the vertex containing dh is drawn.
  template<typename LCC>
  bool draw_vertex(const LCC&,
                   typename LCC::Dart_const_handle) const
  { return true; }

  /// @return true iff the volume containing dh is drawn in wireframe.
  template<typename LCC>
  bool volume_wireframe(const LCC& alcc,
                        typename LCC::Dart_const_handle dh) const
  { return !(alcc.template info<3>(dh).is_filled()); }
  /// @return true iff the face containing dh is drawn in wireframe.
  template<typename LCC>
  bool face_wireframe(const LCC&,
                        typename LCC::Dart_const_handle) const
  { return false; }

  /// @return true iff the volume containing dh is colored.
  template<typename LCC>
  bool colored_volume(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the face containing dh is colored.
  ///  if we have also colored_volume(alcc, dh), the volume color is
  ///  ignored and only the face color is considered.
  template<typename LCC>
  bool colored_face(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the edge containing dh is colored.
  template<typename LCC>
  bool colored_edge(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the vertex containing dh is colored.
  template<typename LCC>
  bool colored_vertex(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return false; }

  /// @return the color of the volume containing dh
  ///  used only if colored_volume(alcc, dh) and !colored_face(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color volume_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  { return alcc.template info<3>(dh).color(); }
  /// @return the color of the face containing dh
  ///  used only if colored_face(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color face_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the edge containing dh
  ///  used only if colored_edge(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color edge_color(const LCC&,
                         typename LCC::Dart_const_handle) const
  { return CGAL::IO::Color(0, 0, 0); }
  /// @return the color of the vertex containing dh
  ///  used only if colored_vertex(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color vertex_color(const LCC&,
                           typename LCC::Dart_const_handle) const
  { return CGAL::IO::Color(0, 0, 0); }
};


class Viewer : public CGAL::SimpleLCCViewerQt<LCC, MyDrawingFunctorLCC>
{
  Q_OBJECT

  typedef CGAL::SimpleLCCViewerQt<LCC, MyDrawingFunctorLCC> Base;

public:
  Viewer(QWidget* parent);
  void setScene(Scene* scene_, bool doredraw=true);
  void keyPressEvent(QKeyEvent *e);
  virtual QString helpString() const;

  /// Abdou
  void refresh();
  // Sending to the viewer primitives to be drawn
  void add_vol3d( std::string &, Drawer::vol3d& );
  void add_vtx_layer( std::string &, Drawer::vtx_layer& );
  void add_edge_layer( std::string &, Drawer::edge_layer& );
  // Drawing primitives external to the LCC
  void draw_edge_layers();
  void draw_vertex_layers();
  void draw_volume_layers();

  // Handling the primitive containers
  void print_size_layer_maps()
  {
      std::cout << "Vertices to draw: " << vtx_to_draw.size() << std::endl;
      std::cout << "Edges to draw: " << edge_to_draw.size() << std::endl;
      std::cout << "Vol3d to draw: " << vol3d_to_draw.size() << std::endl;
  }
  void clear_vtx_layer_map(){vtx_to_draw.clear();}
  void clear_edge_layer_map(){edge_to_draw.clear();}
  void clear_vol3d_layer_map(){vol3d_to_draw.clear();}
  void clear_layer_maps()
  {
      vtx_to_draw.clear();
      edge_to_draw.clear();
      vol3d_to_draw.clear();
  }

public Q_SLOTS:
  void sceneChanged();

private:
  Scene* scene;
  bool m_previous_scene_empty;

  /// Abdou
  vtx_layer_map vtx_to_draw;
  edge_layer_map edge_to_draw;
  vol3d_layer_map vol3d_to_draw;
};

#endif
