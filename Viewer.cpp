// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//                 Abdoulaye Diakite <diakite.abdoulaye@gmail.com>

#include "Viewer.h"
#include <CGAL/Qt/vec.h>

Viewer::Viewer(QWidget* parent) :
    Base(parent, nullptr, ""),
    m_previous_scene_empty(true)
{}

void Viewer::setScene(Scene* scene_, bool doredraw)
{
  scene = scene_;
  set_lcc(scene->lcc, doredraw);
}

void Viewer::sceneChanged()
{
  Base::compute_elements();

  // Draw network if any
  draw_edge_layers();
  draw_vertex_layers();

  this->camera()->
      setSceneBoundingBox(CGAL::qglviewer::Vec(m_bounding_box.xmin(),
                                               m_bounding_box.ymin(),
                                               m_bounding_box.zmin()),
                          CGAL::qglviewer::Vec(m_bounding_box.xmax(),
                                               m_bounding_box.ymax(),
                                               m_bounding_box.zmax()));

  Base::redraw();
  if (m_previous_scene_empty)
  { this->showEntireScene(); }

  m_previous_scene_empty = scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
  // const Qt::KeyboardModifiers modifiers = e->modifiers();
  Base::keyPressEvent(e);
}

QString Viewer::helpString() const
{ return Base::helpString("LCC Demo"); }



/// Abdou
void Viewer::add_vtx_layer( std::string& k, Drawer::vtx_layer& vl )
{
//    vtx_to_draw[k] = vl;
    vtx_to_draw.insert( std::pair<std::string, Drawer::vtx_layer>(k, vl) );
//    std::cout << "\nJust added something: " << vtx_to_draw.size() << std::endl;
}

void Viewer::add_edge_layer( std::string& k, Drawer::edge_layer& el )
{
    edge_to_draw[k] = el;
}

void Viewer::add_vol3d( std::string& k, Drawer::vol3d& v3d )
{
    vol3d_to_draw[k] = &v3d;
}

/// Draw the layers of (visible) vertices
void Viewer::draw_vertex_layers()
{
    // For each layer of vertices
    for (auto const &it : vtx_to_draw )
    {
        // If the layer is visible
        if (it.second.is_visible)
        {
            // add all the vertices
            for( const auto &itx : it.second.vertices )
            {
                add_point( itx, it.second.layer_color );
                m_bounding_box += itx.bbox();
//                std::cout << "\t" << itx << std::endl;
            }
        }
    }
}


/// Draw the layers of (visible) edges
void Viewer::draw_edge_layers()
{
    // For each layer of vertices
    for (const auto &it : edge_to_draw )
    {
        // If the layer is visible
        if (it.second.is_visible)
        {
            // add all the edges
            for( const auto &ite : it.second.edges )
            {
                add_segment( ite.v1, ite.v2, it.second.layer_color );
                m_bounding_box += ite.v1.bbox();
                m_bounding_box += ite.v2.bbox();
//                std::cout << "\t" << ite.v1 << " - " << ite.v2 << std::endl;
            }
        }
    }
}

/// Draw the layers of (visible) volumes

///////////
