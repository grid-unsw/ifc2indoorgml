// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Abdoulaye Diakite <diakite.abdoulaye1@gmail.com>
//
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Combinatorial_map_save_load.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>


// Abdou's change
// Function to split strings based on selected charachter
#include <regex>
inline std::vector<std::string> strSplit(const std::string str, const std::string regex_str)
{
    std::regex regexz(regex_str);
    std::vector<std::string> list(std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
                                  std::sregex_token_iterator());
    return list;
}

// Global random
extern CGAL::Random myrandom;

// Use to define properties on volumes.
#define LCC_DEMO_VISIBLE 1 // if not visible => hidden
#define LCC_DEMO_FILLED  2 // if not filled, wireframe

class Volume_info
{
    friend void CGAL::read_cmap_attribute_node<Volume_info>
    (const boost::property_tree::ptree::value_type &v,Volume_info &val);

    friend void CGAL::write_cmap_attribute_node<Volume_info>(boost::property_tree::ptree & node,
                                                             const Volume_info& arg);
public:
    Volume_info() : m_color(CGAL::Color(myrandom.get_int(0,256),
                                        myrandom.get_int(0,256),
                                        myrandom.get_int(0,256))),
        m_status( LCC_DEMO_VISIBLE | LCC_DEMO_FILLED )
      , m_id("")
      , m_semClass("")
      , m_label("")
    {}

    CGAL::Color& color()
    { return m_color; }
    const CGAL::Color& color() const
    { return m_color; }

    std::string color_name() const
    {
        std::ostringstream ss;
        ss<<std::setfill('0');
        ss<<"#"<<std::hex<<std::setw(2)<<(int)m_color.red()
         <<std::setw(2)<<(int)m_color.green()<<std::setw(2)<<(int)m_color.blue();
        return ss.str();
    }

    bool is_visible() const
    { return (m_status & LCC_DEMO_VISIBLE)!=0; }
    bool is_filled() const
    { return (m_status & LCC_DEMO_FILLED)!=0; }
    bool is_filled_and_visible() const
    { return is_filled() && is_visible(); }

    void set_visible(bool val=true)
    {
        if ( is_visible()==val ) return;
        if ( val ) m_status = m_status | LCC_DEMO_VISIBLE;
        else       m_status = m_status ^ LCC_DEMO_VISIBLE;
    }
    void set_filled(bool val=true)
    {
        if ( is_filled()==val ) return;
        if ( val ) m_status = m_status | LCC_DEMO_FILLED;
        else       m_status = m_status ^ LCC_DEMO_FILLED;
    }

    void negate_visible()
    { set_visible(!is_visible()); }
    void negate_filled()
    { set_filled(!is_filled()); }


    /// /// /// /// /// ///
    /// Abdou's change  ///
    /// /// /// /// /// ///
    // Unique ID of a volume
    std::string id()
    { return m_id; }
    void set_id( std::string id )
    { m_id = id; }

    // semantic class
    std::string semClass()
    { return m_semClass; }
    void set_semClass(std::string sc)
    { m_semClass = sc; }

    // volume's label
    std::string label()
    { return m_label; }
    void set_label(std::string l)
    { m_label = l; }

    // list of IDs of volumes related to this one
    std::vector<std::string> relatedCells()
    { return m_relatedCellIDs; }
    void addRelatedCell( std::string id )
    { m_relatedCellIDs.push_back(id); }

    // Set specific color for a volume
    void set_color(int r, int g, int b)
    { m_color = CGAL::Color(r,g,b); }

private:
    CGAL::Color m_color;
    char        m_status;

    /// Abdou's change
    std::string m_id, m_label, m_semClass;
    std::vector<std::string> m_relatedCellIDs;
};

namespace CGAL
{

template<>
inline void read_cmap_attribute_node<Volume_info>
(const boost::property_tree::ptree::value_type &v,Volume_info &val)
{
    try
    {
        val.m_status = v.second.get<int>("status");
    }
    catch(const std::exception &  )
    {}

    try
    {
        char r = v.second.get<int>("color-r");
        char g = v.second.get<int>("color-g");
        char b = v.second.get<int>("color-b");
        val.m_color = CGAL::Color(r,g,b);
    }
    catch(const std::exception &  )
    {}

    /// Abdou's change
    try
    {
        val.m_id = v.second.get<std::string>("unique-id");
    }
    catch(const std::exception &  )
    {}

    try
    {
        val.m_semClass = v.second.get<std::string>("semantic_class");
    }
    catch(const std::exception &  )
    {}

    try
    {
        val.m_label = v.second.get<std::string>("label");
    }
    catch(const std::exception &  )
    {}

    try
    {
        val.m_relatedCellIDs = strSplit(v.second.get<std::string>("related_cells"),",");
    }
    catch(const std::exception &  )
    {}
}

// Definition of function allowing to save custon information.
template<>
inline void write_cmap_attribute_node<Volume_info>(boost::property_tree::ptree & node,
                                                   const Volume_info& arg)
{
    boost::property_tree::ptree & nValue = node.add("v","");
    nValue.add("status",(int)arg.m_status);
    nValue.add("color-r",(int)arg.m_color.r());
    nValue.add("color-g",(int)arg.m_color.g());
    nValue.add("color-b",(int)arg.m_color.b());

    /// Abdou's change
    nValue.add("unique-id",arg.m_id);
    nValue.add("semantic_class",arg.m_semClass);
    nValue.add("label",arg.m_label);
    std::string related = "";
    for(int i=0; i<arg.m_relatedCellIDs.size(); i++){
        if ( i != 0 ) related += ",";
        related += arg.m_relatedCellIDs[i];
    }
    nValue.add("related_cells",related);
}

}

class Myitems
{
public:
    template < class Refs >
    struct Dart_wrapper
    {
        typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attrib;
        typedef CGAL::Cell_attribute< Refs, Volume_info> Volume_attrib;

        typedef CGAL::cpp11::tuple<Vertex_attrib,void,void,
        Volume_attrib> Attributes;
    };
};

typedef CGAL::Linear_cell_complex_traits
<3,CGAL::Exact_predicates_inexact_constructions_kernel> Mytraits;

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,Mytraits,Myitems> LCC;
typedef LCC::Dart_handle      Dart_handle;
typedef LCC::Vertex_attribute Vertex;

typedef LCC::Point    Point_3;
typedef LCC::Vector   Vector_3;

typedef CGAL::Timer Timer;

struct Vertex_info
{
    Dart_handle dh;
    Vector_3 v;
};

struct Face_info {
    bool exist_edge[3];
    bool is_external;
    bool is_process;
};

typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

typedef CGAL::Triangulation_face_base_with_info_2<Face_info,P_traits> Fb1;

typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
// typedef CGAL::No_intersection_tag                                     Itag;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
Itag>              CDT;

struct Scene {
    LCC* lcc;
};


/// //////////////////////////////////////////////////////////////////////////////
/// Structures and containers to handle geometry drawing of non-LCC features
/// //////////////////////////////////////////////////////////////////////////////
namespace Drawer
{
struct vtx_layer
{
    bool is_visible;
    CGAL::Color layer_color;
    std::vector<Point_3> vertices;

    vtx_layer(): is_visible(true),
        layer_color(CGAL::Color(0,0,0))
    {}
};

struct edge
{
    Point_3 v1, v2;
};

struct edge_layer
{
    bool is_visible;
    CGAL::Color layer_color;
    std::vector<edge> edges;

    edge_layer(): is_visible(true),
        layer_color(CGAL::Color(0,0,0))
    {}
};

struct face3d
{
    std::vector<Point_3> vertices;
};

class vol3d
{
public:
    vol3d() :
        uid(-1),
        label(""),
        color( CGAL::Color(100,100,100) )
    {}
    ~vol3d(){}

    void set_label(std::string l)
    { label = l; }
    std::string get_label()
    { return label; }

    // Set specific color for a volume
    void set_color(int r, int g, int b)
    { color = CGAL::Color(r,g,b); }
    // Get color value of a volume
    CGAL::Color get_color()
    { return color; }

private:
    int uid;
    std::list<face3d> facelist;
    std::string label;
    CGAL::Color color;
};
}

///// Global variables to store entities to be drawn
typedef std::map<std::string, Drawer::vol3d*>       vol3d_layer_map;
typedef std::map<std::string, Drawer::vtx_layer>    vtx_layer_map;
typedef std::map<std::string, Drawer::edge_layer>   edge_layer_map;

//extern vtx_layer_map vtx_to_draw;
//extern edge_layer_map edges_to_draw;
//extern vol3d_layer_map vol3d_to_draw;

#endif
