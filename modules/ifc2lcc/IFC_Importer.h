// Copyright (c) 2016.
// All rights reserved.
//
// Author(s)     : Abdoulaye A. Diakite <diakite.abdoulaye1@gmail.com>


#ifndef IFC_IMPORTER_H
#define IFC_IMPORTER_H

//#include "IFC2LCC_typedefs.h"

#include "typedefs.h"
#include "LCC_SpecialOps.h"

/// Kernels
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 IK;
typedef CGAL::Exact_predicates_exact_constructions_kernel                   EK;
typedef typename LCC::Traits::FT                                            FT;
typedef EK::FT                                                              ET;

/// Polyhedron 3 & Nef
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/convex_decomposition_3.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>

/// Polygon_mesh_processing package
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

typedef CGAL::Polyhedron_3<EK>                                              Polyhedron_EK;
typedef Polyhedron_EK::HalfedgeDS                                           HDS_EK;
typedef CGAL::Polyhedron_3<IK>                                              Polyhedron_IK;
typedef Polyhedron_IK::HalfedgeDS                                           HDS_IK;
typedef CGAL::Nef_polyhedron_3<EK>                                          Nef_polyhedron;

/// Inexact Primitives
//typedef typename LCC::Base::Dart_handle                                     Dart_handle;
typedef typename LCC::Base::Dart_const_handle                               Dart_const_handle;
typedef typename LCC::Base::Dart_range                                      Dart_range;

typedef typename LCC::Traits::Point                                         Point;
typedef typename LCC::Traits::Vector                                        Vector;
typedef typename LCC::Traits::Plane_3                                       Plane;
typedef typename LCC::Traits::Segment_3                                     Segment;

typedef typename LCC::Traits::Point_2                                       Point2d;
typedef typename LCC::Traits::Vector_2                                      Vector2d;
typedef typename LCC::Traits::Segment_2                                     Segment2d;

/// Converters
typedef CGAL::Cartesian_converter<EK, IK>                                   Converter_1;
typedef CGAL::Cartesian_converter<IK ,EK>                                   Converter_2;

extern double size_scale;
extern uint volume_label;
extern uint vol_counter;

#ifdef IFCPP_ON

//        #include <ifcpp/model/BuildingModel.h>
//        #include <ifcpp/reader/IfcPPReaderSTEP.h>
#include <ifcpp/model/BuildingModel.h>
#include <ifcpp/model/BuildingException.h>
#include <ifcpp/reader/ReaderSTEP.h>
#include <ifcpp/IFC4/include/IfcProduct.h>
#include <ifcpp/IFC4/include/IfcGloballyUniqueId.h>
#include <ifcpp/IFC4/include/IfcLabel.h>
#include <ifcpp/IFC4/include/IfcText.h>

#include <ifcpp/IFC4/include/IfcDoor.h>
#include <ifcpp/IFC4/include/IfcSpace.h>
#include <ifcpp/IFC4/include/IfcWindow.h>
#include <ifcpp/IFC4/include/IfcElement.h>
#include <ifcpp/IFC4/include/IfcRelAggregates.h>
#include <ifcpp/IFC4/include/IfcBuildingStorey.h>
#include <ifcpp/IFC4/include/IfcOpeningElement.h>
#include <ifcpp/IFC4/include/IfcVirtualElement.h>
#include <ifcpp/IFC4/include/IfcRelFillsElement.h>
#include <ifcpp/IFC4/include/IfcRelSpaceBoundary.h>
#include "ifcpp/IFC4/include/IfcRelConnectsElements.h"
#include <ifcpp/IFC4/include/IfcSpatialStructureElement.h>
#include <ifcpp/IFC4/include/IfcRelContainedInSpatialStructure.h>


//#include <ifcpp/IFC4/include/IfcProject.h>
//        #include <ifcpp/geometry/GeometryConverter.h>
//        #include <ifcpp/geometry/GeometryInputData.h>


#include <ifcpp/geometry/Carve/GeometryConverter.h>
#include <ifcpp/geometry/Carve/ConverterOSG.h>
#include <ifcpp/geometry/Carve/GeomUtils.h>

//        #include <ifcpp/geometry/OCC/GeometryConverterOCC.h>
//        #include <ifcpp/geometry/OCC/SceneGraphConverterOCC.h>
//        #include <ifcpp/geometry/OCC/GeomUtilsOCC.h>

#include <ifcpp/geometry/SceneGraphUtils.h>
#include <ifcpp/geometry/GeometrySettings.h>

    extern shared_ptr<GeometryConverter> geomConverter;
    extern std::map<std::string, shared_ptr<ProductShapeData> > map_ifc_shapes;
    extern std::map<int, shared_ptr<BuildingEntity> > map_ifc_entities;
    extern std::vector< shared_ptr<IfcBuildingStorey> > ifc_storey_list;
    extern std::vector<std::vector< shared_ptr<IfcSpace> > > spaces_of_storeys;

#endif


namespace CGAL
{

    /// Simplif
    typedef std::vector<Point>                                              vec_pt3d;
    typedef std::vector <Dart_handle>                                       vec_dart;
    typedef std::pair<Point, Point>                                         pair_pt3d;
    typedef std::pair<Dart_handle, Dart_handle>                             pair_dart;
    typedef std::pair<std::string, std::string>                             string_pair;
    typedef std::pair<Dart_handle, Vector>                                  pair_dart_vec3d;
    typedef std::pair<pair_dart_vec3d, pair_dart_vec3d>                     double_pair_dart_vec3d;
    typedef std::map<double, std::pair<vec_dart, vec_dart> >                map_angles_pair_vec_darts;



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    ///                                             Convenient Functions
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */

    std::string get_unique_key_pt( Point& );
//    Vector poly_normal (const LCC&, Dart_handle, bool normalize = true);
    bool Bbox_is_degenerate(Bbox_3);
    Dart_handle Create_bbox_volume (LCC&, Bbox_3&);


    /* ////////////////////////////////////////////////////////////////////////////// */
    /// Return the angle between two 2D vectors v1 and v2, in the range [0Â°, 360Â°]
    /* ////////////////////////////////////////////////////////////////////////////// */
    template<class K> double compute_2d_angle(const Vector_2<K>& v1, const Vector_2<K>& v2)
    {
        double s1 = to_double(v1[0]), s2 = to_double(v1[1]),
               s3 = to_double(v2[0]), s4 = to_double(v2[1]);

        double a = atan2(s4, s3) - atan2(s2, s1);

        return ( ((a/M_PI)*180.0) + (a > 0.0 ? 0.0 : 360.0));
    }


    /* ////////////////////////////////////////////////////////////////////////////// */
    /// Creates a Polyhedron_3 from a 3-cell of an LCC
    /* ////////////////////////////////////////////////////////////////////////////// */
    class Build_exact_polyhedron_3_from_lcc : public CGAL::Modifier_base<HDS_EK>
    {
        public:
            LCC* local_lcc;
            Dart_handle local_d;
            uint v, f, nd;
            Converter_2 in2ex;

            Build_exact_polyhedron_3_from_lcc () {}
            Build_exact_polyhedron_3_from_lcc( LCC& alcc, Dart_handle d )
            {
                CGAL_assertion( d != LCC::null_handle );

                local_lcc = &alcc;
                local_d = d;

    //            std::cout << "The volume is a(n) " << local_lcc->template info<3>(d).label() << std::endl;

                // Number of vertices, faces and darts (edges)
                v = local_lcc->template one_dart_per_incident_cell<0,3>(d).size();
                std::cout << "Number of vertices: " << v << std::endl;

                f = local_lcc->template one_dart_per_incident_cell<2,3>(d).size();
                std::cout << "Number of faces: " << f << std::endl;

                nd = local_lcc->template darts_of_cell<3>(d).size();
            }
            void operator()( HDS_EK& hds)
            {
                // Postcondition: hds is a valid polyhedral surface.
                std::map <std::string, uint, LCCtools::cmp_string> pt_id;
                uint id = 0;

                Polyhedron_incremental_builder_3<HDS_EK> Pbuild(hds, true);
                Pbuild.begin_surface( v, f, nd );

                // Create a unique key map of the vertices (no vertex duplication)
                // and build vertices
                typename LCC::Base::template One_dart_per_incident_cell_range<0,3>::iterator
                        vitr = local_lcc->template one_dart_per_incident_cell<0,3>(local_d).begin(),
                        vitrend = local_lcc->template one_dart_per_incident_cell<0,3>(local_d).end();
                for (; vitr != vitrend; vitr++)
                {
                    pt_id[ get_unique_key_pt( local_lcc->point( vitr ) ) ] = id++;
                    //                    std::cout << /*id << " ---> " <<*/ local_lcc->point( vitr ) << std::endl;

                    Pbuild.add_vertex( in2ex(local_lcc->point( vitr )) );
                }

                // Build faces
                typename LCC::Base::template One_dart_per_incident_cell_range<2,3>::iterator
                        fitr = local_lcc->template one_dart_per_incident_cell<2,3>(local_d).begin(),
                        fitrend = local_lcc->template one_dart_per_incident_cell<2,3>(local_d).end();
                for (; fitr != fitrend; fitr++)
                {
                    Pbuild.begin_facet();

                    // go through all the points of the face
                    typename LCC::Base::template Dart_of_orbit_range<1>::iterator itr(*local_lcc, fitr);
                    for (; itr.cont(); itr++)
                        Pbuild.add_vertex_to_facet( pt_id[ get_unique_key_pt( local_lcc->point( itr ) ) ] );

                    Pbuild.end_facet();
                }

                Pbuild.end_surface();
            }
    };

    /* ////////////////////////////////////////////////////////////////////////////// */
    /// Creates a Polyhedron_3 from a 3-cell of an LCC
    /* ////////////////////////////////////////////////////////////////////////////// */
    class Build_inexact_polyhedron_3_from_lcc : public CGAL::Modifier_base<HDS_IK>
    {
        public:
            LCC* local_lcc;
            Dart_handle local_d;
            uint v, f, nd;

            Build_inexact_polyhedron_3_from_lcc () {}
            Build_inexact_polyhedron_3_from_lcc( LCC& alcc, Dart_handle d )
            {
                CGAL_assertion( d != LCC::null_handle );

                local_lcc = &alcc;
                local_d = d;

    //            std::cout << "The volume is a(n) " << alcc.template info<3>(d).label() << std::endl;

                // Number of vertices, faces and darts (edges)
                v = alcc.template one_dart_per_incident_cell<0,3>(d).size();
                f = alcc.template one_dart_per_incident_cell<2,3>(d).size();
                nd = alcc.template darts_of_cell<3>(d).size();
            }
            void operator()( HDS_IK& hds)
            {
                // Postcondition: hds is a valid polyhedral surface.
                std::map <std::string, uint, LCCtools::cmp_string> pt_id;
                uint id = 0;

                Polyhedron_incremental_builder_3<HDS_IK> Pbuild(hds, true);
                Pbuild.begin_surface( v, f, nd );

                // Create a unique key map of the vertices (no vertex duplication)
                // and build vertices
                typename LCC::Base::template One_dart_per_incident_cell_range<0,3>::iterator
                        vitr = local_lcc->template one_dart_per_incident_cell<0,3>(local_d).begin(),
                        vitrend = local_lcc->template one_dart_per_incident_cell<0,3>(local_d).end();
                for (; vitr != vitrend; vitr++)
                {
                    pt_id[ get_unique_key_pt( local_lcc->point( vitr ) ) ] = id++;
                    //                    std::cout << /*id << " ---> " <<*/ local_lcc->point( vitr ) << std::endl;

                    Pbuild.add_vertex( local_lcc->point( vitr ) );
                }

                // Build faces
                typename LCC::Base::template One_dart_per_incident_cell_range<2,3>::iterator
                        fitr = local_lcc->template one_dart_per_incident_cell<2,3>(local_d).begin(),
                        fitrend = local_lcc->template one_dart_per_incident_cell<2,3>(local_d).end();
                for (; fitr != fitrend; fitr++)
                {
                    Pbuild.begin_facet();

                    // go through all the points of the face
                    typename LCC::Base::template Dart_of_orbit_range<1>::iterator itr(*local_lcc, fitr);
                    for (; itr.cont(); itr++)
                        Pbuild.add_vertex_to_facet( pt_id[ get_unique_key_pt( local_lcc->point( itr ) ) ] );

                    Pbuild.end_facet();
                }

                Pbuild.end_surface();
            }
    };


    bool Convert_3cell_2_inexact_polyhedron3( LCC&, Dart_handle, CGAL::Polyhedron_3<IK>&);
    bool Convert_3cell_2_exact_polyhedron3( LCC&, Dart_handle, CGAL::Polyhedron_3<EK>&);
    bool Convert_3cell_2_Nef3( LCC&, Dart_handle, Nef_polyhedron&);
    Dart_handle Build_lcc_from_exact_polyhedron_3(LCC&, Polyhedron_EK&);
    Dart_handle Build_lcc_from_inexact_polyhedron_3(LCC&, Polyhedron_IK&);
    Dart_handle Build_lcc_from_nef_3(LCC&, Nef_polyhedron&);
    void Perform_corefinement_reconstruction (LCC &alcc, vec_dart &, vec_dart &, vec_dart &);


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    ///                                        End of Convenient Functions
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */





    /// Structures for Spaces and Sub-Spaces
    struct Sub_Space
    {
        std::string uid;
        uint space_id;
        //own dart and dart of the related O-Space
        Dart_handle dart, related_object;
        double buffer_size;

        Sub_Space() : uid(""),
                      space_id(100000),
                      dart(LCC::null_handle),
                      related_object(LCC::null_handle),
                      buffer_size(0.0)
        {}
    };
    typedef struct Sub_Space Sub_Space;

    struct O_Space
    {
        std::string uid;
        uint space_id;
        // own dart + attached b-space and f-space (if any)
        Dart_handle dart, b_space, f_space;
        //simplified objects aggregated in the O-Space
        vec_dart contained_obj;

        O_Space() : uid(""),
                    space_id(100000),
                    dart(LCC::null_handle),
                    b_space(LCC::null_handle),
                    f_space(LCC::null_handle)
        {}
    };
    typedef struct O_Space O_Space;

    struct R_Space
    {
        std::string uid;
        uint space_id;
        Dart_handle dart;
        Nef_polyhedron nef;

        R_Space() : uid(""),
                    space_id(100000),
                    dart(LCC::null_handle)
        {}
    };
    typedef struct R_Space R_Space;

    struct Space
    {
        std::string uid;
        Dart_handle dart;
        #ifdef IFCPP_ON
            shared_ptr<IfcSpace> to_IfcSpace;
        #endif
        //simplified objects contained in the space and surrounding openings
        vec_dart contained_elem,
                 surrounding_ops;

        std::vector<O_Space*> contained_O_Spaces;
        std::vector<Sub_Space*> contained_B_Spaces;
        std::vector<Sub_Space*> contained_F_Spaces;
        std::vector<R_Space*> contained_R_Spaces;

        Space() : uid(""),
                  dart(LCC::null_handle)
        {}
    };
    typedef struct Space Space;

    // Structure to keep/exchange the attributes of a volume
    struct Vol_Attribs
    {
        // Collect attributes from 3-cell(d1)
        std::string m_uID; // unique ID
        uint m_vol_label; // label common to all elements of a same volume
        std::string m_label; // label specific to a semantic class of 3-cell
        std::vector<std::string> m_IfcInfo; // Ifc Global Unique ID
        bool m_select; // Whether a 3-cell is selected or not
        bool m_init; // Whether a 3-cell is fully initialized or not
        CGAL::Color m_color; // Color of the 3-cell

        Vol_Attribs() : m_uID(""),
                        m_vol_label(0),
                        m_label(""),
                        m_IfcInfo({"","","","","",""}),
                        m_select(false),
                        m_init(false),
                        m_color( CGAL::Color(255, 255, 255) )
        {}
    };
    typedef Vol_Attribs Vol_Attribs;


    void Rescale_LCC_size (LCC& alcc, double scaling);
    void Copy_Attributes_from_3cell(LCC&, Dart_handle&, Vol_Attribs&);
    void Paste_Attributes_to_3cell(LCC&, Dart_handle&, Vol_Attribs&);
    void Transfer_Attributes_of_3cell(LCC&, Dart_handle&, Dart_handle&);
#ifdef IFCPP_ON
        Point Get_Osg_3dpt_from_Mesh_Vertex (carve::mesh::Vertex<3>&, Point shift = Point(ORIGIN));
        std::vector<Point> Get_Vertices_from_Mesh (carve::mesh::Mesh<3>*);
        std::vector<Point> Get_Vertices_from_Single_Meshset (shared_ptr<carve::mesh::MeshSet<3> >& );
        std::vector<Point> Get_Vertices_from_Meshsets (std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& );
#endif
        void Pointset_Simplification (std::vector<Point >&, std::pair<double, std::vector<Point> >&, int mode = 1);


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Applies a rotation angle to a given vector
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class K> Vector_2<K> Rotate_vector (Vector_2<K>& v, double ang)
    {
        // convert to radian
        ang = (ang*M_PI)/180.0;
        double x = to_double(v.x())*cos(ang) - to_double(v.y())*sin(ang),
               y = to_double(v.x())*sin(ang) + to_double(v.y())*cos(ang);
        return Vector_2<K>(x, y);
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Returns the area of a rectangular box
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class K> double Area_of_box (std::vector<Point_2<K> >& corner)
    {
        assert(corner.size() == 4);

        double width, length;

        // The area is width * length
        width = sqrt( pow((corner[1].x() - corner[0].x()), 2) + pow((corner[1].y() - corner[0].y()), 2) );
        length = sqrt( pow((corner[2].x() - corner[1].x()), 2) + pow((corner[2].y() - corner[1].y()), 2) );

        return width * length;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Find the intersection point between two lines
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class K> bool Get_lines_intersection_point( Line_2<K>& l1, Line_2<K>& l2, Point_2<K>& res )
    {
//                CGAL::cpp11::result_of<K::Intersect_2(Line_2, Line_2)>::type result = intersection(l1, l2);
        auto result = intersection(l1, l2);
        if (result)
        {
            if (const Line_2<K>* l = boost::get<Line_2<K> >(&*result))
            {
                std::cout << "Got a line xD... " << *l << std::endl;
                return false;
            }
            else
            {
                const Point_2<K>* p = boost::get<Point_2<K> >(&*result);
//                        std::cout << "Got an intersection point! " <<  *p << std::endl;
                res = *p;
                return true;
            }
        }

        return false;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Simplify 2D polygon's boundary by avoiding more than 2 aligned points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class K> void Remove_aligned_points_on_poly (std::vector<Point_2<K> >& poly)
    {
        std::vector<Point_2<K> > res;
        double ang, ang_sup = 181.0, ang_inf = 179.0;

        // for the first point
        ang = compute_2d_angle<K>( Vector_2<K>( poly[0], poly[1] ),
                                   Vector_2<K>( poly[0], poly[ poly.size()-1 ] ));
        if ( ang > ang_sup || ang < ang_inf )
            res.push_back( poly[0] );

        uint i;
        for(i=1; i<poly.size()-1; i++)
        {
            ang = compute_2d_angle<K>( Vector_2<K>( poly[i], poly[i+1] ),
                                       Vector_2<K>( poly[i], poly[i-1] ));
            if ( ang > ang_sup || ang < ang_inf )
                res.push_back( poly[i] );
        }

        std::cout << "size? " << poly.size() << "; i? " << i << std::endl;
        // for the last point
        ang = compute_2d_angle<K>( Vector_2<K>( poly[i], poly[0] ),
                                   Vector_2<K>( poly[i], poly[i-1] ));
        if ( ang > ang_sup || ang < ang_inf )
            res.push_back( poly[i] );

        std::cout << poly.size() - res.size() << " point(s) removed!" << std::endl;

        poly = res;

        std::cout << "\nNew poly after aligned points removed: " << std::endl;
        for(uint j=0; j<poly.size(); j++)
            std::cout << poly[j] << std::endl;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Find the smallest bounding rectangle of a (convex) 2D polygon using the rotating calipers approach
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template <class K> void Get_smallest_bounding_rectangle ( std::vector<Point_2<K> >& poly, std::vector<Point_2<K> >& res, double eps = 0.001 )
    {
        Remove_aligned_points_on_poly<K>(poly);
        Polygon_2<K> mypoly;

        // Determine the initial bbox of the polygon
        LCCtools::cmp_eps_double mycomp(eps);
        std::map<double, std::vector<uint>, LCCtools::cmp_eps_double > min_max_x (mycomp), min_max_y (mycomp);
        std::vector<bool> is_used(poly.size(), false);

        // Sort the point to get the extremities and store the "next_pt" of each point of the polygon
        std::map< uint, Point_2<K> > next_pt;
        for(uint i=0; i<poly.size(); i++)
        {
            min_max_x[ to_double(poly[i].x()) ].push_back( i );
            min_max_y[ to_double(poly[i].y()) ].push_back( i );

            if (i == poly.size()-1)
                next_pt[i] = poly[0];
            else
                next_pt[i] = poly[i+1];

            mypoly.push_back( poly[i] );
        }

        // Keep the indices of the corners of the bbox
//                std::cout << "min_max size: " << min_max_x.size() << " & "<< min_max_y.size() << std::endl;
        uint idx[4];
        idx[0] = min_max_x.rbegin()->second[0];
        idx[1] = min_max_y.rbegin()->second[0];
        idx[2] = min_max_x.begin()->second[0];
        idx[3] = min_max_y.begin()->second[0];
//                std::cout << "First pts: " << idx[0] << ", " << idx[1] << ", " << idx[2] << ", " << idx[3] << std::endl;

        // Set the first bbox throught the computed corners
        std::vector<Point_2<K> > corner;
        corner.push_back(Point_2<K>(min_max_x.rbegin()->first, min_max_y.rbegin()->first));
        corner.push_back(Point_2<K>(min_max_x.begin()->first, min_max_y.rbegin()->first));
        corner.push_back(Point_2<K>(min_max_x.begin()->first, min_max_y.begin()->first));
        corner.push_back(Point_2<K>(min_max_x.rbegin()->first, min_max_y.begin()->first));
        std::cout << "First bbox: " << std::endl;

        // Compute the area of the bbox
        std::map<double, std::vector<Point_2<K> > > box_area;
        double cur_area = Area_of_box(corner);
        box_area[ cur_area ] = corner;

        // We determine two vectors: v_p -> a vector from one point of the polygon to its next
        //                           v_l -> a vector from the point of the polygon to the follow one on the bbox
        // in order to compute the angle former by the polygon and the corners of the bbox.
        Point_2<K> pt[4];
        Vector_2<K> v_p[4];
        for(uint i=0; i<4; i++)
        {
            pt[i] = poly[ idx[i] ];
            is_used[ idx[i] ] = true;
            v_p[i] = Vector_2<K> ( pt[i], next_pt[ idx[i] ] );

//                    std::cout << "v_p" << i << ": " << v_p[i] << std::endl;
            std::cout << corner[i] << std::endl;
        }

        Vector_2<K> v_l[4];
        Point_2<K>  pt_1 = Point_2<K>( pt[0].x(), to_double(pt[0].y()+1.0) ),
                    pt_2 = Point_2<K>( to_double(pt[1].x()-1.0), pt[1].y() ),
                    pt_3 = Point_2<K>( pt[2].x(), to_double(pt[2].y()-1.0) ),
                    pt_4 = Point_2<K>( to_double(pt[3].x()+1.0), pt[3].y() );
        v_l[0] = Vector_2<K> ( pt[0], pt_1 );
        v_l[1] = Vector_2<K> ( pt[1], pt_2 ),
        v_l[2] = Vector_2<K> ( pt[2], pt_3 ),
        v_l[3] = Vector_2<K> ( pt[3], pt_4 );

//                std::cout << "v_p" << 0 << ": " << v_p[0] << std::endl;
//                std::cout << "v_l" << 0 << ": " << v_l[0] << std::endl;
//                std::cout << "\nv_p" << 1 << ": " << v_p[1] << std::endl;
//                std::cout << "v_l" << 1 << ": " << v_l[1] << std::endl;
//                std::cout << "\nv_p" << 2 << ": " << v_p[2] << std::endl;
//                std::cout << "v_l" << 2 << ": " << v_l[2] << std::endl;
//                std::cout << "\nv_p" << 3 << ": " << v_p[3] << std::endl;
//                std::cout << "v_l" << 3 << ": " << v_l[3] << std::endl;

        // Computation of the angles
        double ang[4];
        std::map<double, std::vector<uint> > sort_angles;
        std::map<double, std::vector<uint> > angle_2;
        for(uint i=0; i<4; i++)
        {
            ang[i] = compute_2d_angle<K>( v_l[i], v_p[i] );
            std::cout << "\tang" << i << ": " << ang[i] << std::endl;

            // set the angle to 0 if necessary, to avoid ambiguities.
            if (ang[i] < 5.0 || ang[i] > 355.0)
                ang[i] = 0.0;

            sort_angles[ ang[i] ].push_back( idx[i] );
            angle_2[ ang[i] ].push_back( i );
        }


        // The max number of iteration should not be more than the number of edges of the polygon
        uint nb_iter = poly.size();
        std::cout << "\nMAX ITERATION: " << nb_iter << std::endl;
        double rot_angle, new_angle;

        for( uint i=0; i<nb_iter; i++ )
        {
            std::map<double, std::vector<uint> >::iterator it_ang (sort_angles.begin()), it_ang_nxt( it_ang );
            for(it_ang_nxt++; it_ang!= sort_angles.end(); it_ang++)
            {
                // If the angle at the current point is null, the angle of the next rotation
                // needed to reach zero (meaning, to align an edge of the polygon with at least
                // one side of the caliper) has to be computed. It will correspond to the smallest one
                // in the direction of the rotation.
                if ( (it_ang->first == 0.0) )
                {
                    std::cout << sort_angles.begin()->second.size() << " points aligned: pt" << sort_angles.begin()->second[0] << std::endl;

                    // if the box fits more or less entirely the polygon, stop the process
                    if (nb_iter <= sort_angles.begin()->second.size())
                        break;

                    // Otherwise, remove the nb of edges aligned with the caliper's sides, to have less edges to process...
                    nb_iter -= sort_angles.begin()->second.size();
                    std::cout << "\nremaining iterations: " << nb_iter << std::endl;

                    std::map<double, std::vector<uint> > new_sort_angles;

                    for(uint j=0; j<sort_angles.begin()->second.size(); j++)
                    {
                        uint id = sort_angles.begin()->second[j];
                        uint id_local = angle_2.begin()->second[j];

                        // update the rotation points
                        pt[id_local] = next_pt[id];

//                                Vector_2<K> v_t (poly[id], next_pt[id]);
//                                Point_2<K> pt_t ( 2*( poly[id].x() + v_t.x() ),
//                                                  2*( poly[id].y() + v_t.y() ) );
//                                v_l[id_local] = Vector_2<K>( next_pt[id], pt_t );

                        v_l[id_local] = Vector_2<K> (poly[id], next_pt[id]);

                        if( id + 1 == poly.size() )
                        {
                            v_p[id_local] = Vector_2<K>( next_pt[id], next_pt[0] );
                            idx[id_local] = 0;
                        }
                        else
                        {
                            v_p[id_local] = Vector_2<K>( next_pt[id], next_pt[id+1] );
                            idx[id_local] = id+1;
                        }


//                                std::cout << "\nnew v_p" << id_local << ": " << v_p[id_local] << std::endl;
//                                std::cout << "new v_l" << id_local << ": " << v_l[id_local] << std::endl;

                        new_angle = compute_2d_angle<K>( v_l[id_local], v_p[id_local] );
//                                std::cout << "\tnew angle" << id_local << ": " << new_angle << std::endl;

                        new_sort_angles[ new_angle ].push_back( idx[id_local] );
                    }

                    // No zero angle should be found for the following rotation
                    // In other words, no 3 aligned points allowed on the convex input polygon.
                    assert ( new_sort_angles.begin()->first > 0.0 );

                    if ( it_ang_nxt != sort_angles.end()
                         && new_sort_angles.begin()->first > it_ang_nxt->first )
                    {
                        rot_angle = it_ang_nxt->first;
//                                std::cout << "Rotation angle 1: " << rot_angle;
//                                for( uint k=0; k<it_ang_nxt->second.size(); k++ )
//                                    std::cout << "\n\tfrom pt" << it_ang_nxt->second[k] << std::endl;
                    }
                    else
                    {
                        rot_angle = new_sort_angles.begin()->first;
//                                std::cout << "Rotation angle 2: " << rot_angle;
//                                for( uint k=0; k<new_sort_angles.begin()->second.size(); k++ )
//                                    std::cout << "\n\tfrom pt" << new_sort_angles.begin()->second[k] << std::endl;
                    }
                }

                else
                {
                    rot_angle = it_ang->first;
                    std::cout << "Not zero angle here: " << rot_angle << std::endl;
//                            std::cout << "Rotation angle 3: " << rot_angle;
//                            for( uint k=0; k<it_ang->second.size(); k++ )
//                                std::cout << "\n\tfrom pt" << it_ang->second[k] << std::endl;
                }

                // Perform the rotation of the caliper
                sort_angles.clear();
                angle_2.clear();
                Line_2<K> line[4];
                for(uint j=0; j<4; j++)
                {
                    v_l[j] = Rotate_vector<K>( v_l[j], rot_angle );
                    line[j] = Line_2<K>( pt[j], v_l[j] );

                    ang[j] = compute_2d_angle<K>( v_l[j], v_p[j] );

                    // set the angle to 0 if necessary, to avoid ambiguities.
                    if (ang[j] < 5.0 || ang[j] > 355.0)
                        ang[j] = 0.0;

                    sort_angles[ ang[j] ].push_back( idx[j] );
                    angle_2[ ang[j] ].push_back( j );
                }

                std::cout << i << " --- done over a number of iterations of: " << nb_iter << std::endl;
                Get_lines_intersection_point<K>(line[0], line[1], corner[0]);
                Get_lines_intersection_point<K>(line[1], line[2], corner[1]);
                Get_lines_intersection_point<K>(line[2], line[3], corner[2]);
                Get_lines_intersection_point<K>(line[3], line[0], corner[3]);

//                std::cout << "\tAngle corner 0: " << compute_2d_angle<K>( Vector_2<K>( corner[0], corner[1] ),
//                                                                          Vector_2<K>( corner[0], corner[3] ) ) << std::endl;
//                std::cout << "\tAngle corner 1: " << compute_2d_angle<K>( Vector_2<K>( corner[1], corner[2] ),
//                                                                          Vector_2<K>( corner[1], corner[0] ) ) << std::endl;
//                std::cout << "\tAngle corner 2: " << compute_2d_angle<K>( Vector_2<K>( corner[2], corner[3] ),
//                                                                          Vector_2<K>( corner[2], corner[1] ) ) << std::endl;
//                std::cout << "\tAngle corner 3: " << compute_2d_angle<K>( Vector_2<K>( corner[3], corner[0] ),
//                                                                          Vector_2<K>( corner[3], corner[2] ) ) << std::endl;

                cur_area = Area_of_box(corner);
                box_area[ cur_area ] = corner;

                break;
            }
        }

        double init_area = mypoly.area();
        for(auto it = box_area.begin(); it != box_area.end(); it++)
        {
            // Filter out the boxes smaller than the original shape
            if ( it->first >= init_area )
            {
                res.push_back( it->second[0] );
                res.push_back( it->second[1] );
                res.push_back( it->second[2] );
                res.push_back( it->second[3] );

                break;
            }
        }
    }

    void Z_extrude_2cell (LCC&, Dart_handle, double&);
    Dart_handle Z_extrude_rectangle (LCC&, std::pair<double, std::vector<Point> >&);
    Dart_handle Z_extrude_polygon (LCC&, std::vector<Point>&, double);
    void Z_extrude_upper_lower_box_faces (LCC&, Dart_handle, double, bool upper = true, bool lower = false, bool dir = true);
    Dart_handle OBB_Pointset_Simplification (LCC&, std::vector<Point>&);
    Dart_handle Convex_hull_Pointset_Simplification (LCC&, std::vector<Point>&);
    Dart_handle Concave_hull_Pointset_Simplification (LCC&, std::vector<Point>&);
    Bbox_3 Get_Bbox_Vector_Point (std::vector<Point>&);
//    Bbox_3 Get_Bbox_vol (LCC&, Dart_handle);
    void Rescale_Bbox (Bbox_3 &, double &, bool increase = true);
    bool Boxes_intersect (LCC&, Dart_handle, Dart_handle);
    void Remove_selected_3_cells(LCC&, vec_dart&);

#ifdef IFCPP_ON
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates a 3-cell in the LCC given a (closed) mesh of OSG (assumed to potentially contain half-edge links)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template < class LCC > vec_dart osg_mesh2lcc(LCC& alcc, carve::mesh::Mesh<3>*& mesh, bool sew2 = true)
    {
        vec_dart res, one_res, final_res;

        // each key is unique for one edge
        std::map< std::string, vec_dart, LCCtools::cmp_string > edge_map;
        std::map< std::string, vec_dart, LCCtools::cmp_string >::iterator it;

        Dart_handle d = LCC::null_handle,
                    prev = LCC::null_handle;
        Dart_handle firstFacet = LCC::null_handle;
        Point pt1, pt2;

        uint nb_pt = 0;

//                std::wcout << "Hypothetical nb of edges: " << mesh->closed_edges.size() + mesh->open_edges.size() << std::endl;
//                std::wcout << "Included open edges: " << mesh->open_edges.size() << std::endl;

        // First traversal to build the darts and link them.
        // Get the faces of the mesh
        std::vector<carve::mesh::Face<3>* >& vec_faces = mesh->faces;
        for( size_t i_face = 0; i_face < vec_faces.size(); ++i_face )
        {
            carve::mesh::Face<3>* face = vec_faces[i_face];
            prev = LCC::null_handle;

            // iterate through edges:
            carve::mesh::Edge<3>* edge = face->edge;
            do
            {
                pt1 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v1() );
//                        pt1 = Scale_3D_point(pt1, size_scale);
                pt2 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v2() );
//                        pt2 = Scale_3D_point(pt2, size_scale);

//                        if ( Points_are_eps_close( pt1, pt2 ) )
//                        {
//                            std::wcout << "Very close points..." << std::endl;
//                            std::cout << "\t" << pt1 << std::endl;
//                            std::cout << "\t" << pt2 << std::endl;
//                        }

                {
                    d = alcc.create_dart();

                    //                        std::string key = get_unique_key_edge( pt1, pt2 );
                    std::string key ( get_unique_key_pt(pt1) + get_unique_key_pt(pt2) );
                    edge_map[ key ].push_back(d);

                    //                        std::wcout << "\tkey: " << key.c_str() << key.c_str() << std::endl;

                    if ( prev != LCC::null_handle ){
                        if (alcc.template is_sewable<1>(prev, d)){
                            alcc.template sew<1>(prev, d);
                        }
                        else if (alcc.template is_sewable<1>(d, prev)){
                            alcc.template sew<1>(d, prev);
//                        alcc.template link_beta<1>(prev, d);
                        }
                    }
                    else
                        firstFacet = d;

                    // If the opposite edge is already in the map, sew2 them
                    if (sew2)
                    {
                        it = edge_map.find( std::string( get_unique_key_pt(pt2) +
                                                         get_unique_key_pt(pt1) ) );
                        if ( it != edge_map.end() && it->second.size() == 1 ){
                            if (alcc.template is_sewable<2>(d, it->second[0])){
                                alcc.template sew<2>(d, it->second[0]);
                            }
                            else if (alcc.template is_sewable<2>(it->second[0], d)){
                                alcc.template sew<2>(it->second[0], d);
    //                        alcc.template link_beta<2>(it->second[0], d);
                            }
                        }
                    }

                    prev = d;
                }

                // increment edge
                edge = edge->next;
                nb_pt++;
            }
            while( edge != face->edge );

            if( alcc.template is_sewable<1>(prev, firstFacet) )
                alcc.template sew<1>(prev, firstFacet);
            else if( alcc.template is_sewable<1>(firstFacet, prev) )
                alcc.template sew<1>(firstFacet, prev);
        }

//        std::wcout << "expected nb of points: " << vec_faces.size() * 3 << std::endl;
//        std::wcout << "nb points obtained: " << nb_pt << std::endl;
//        std::wcout << "nb darts created: " << edge_map.size() << std::endl;
//        if (nb_pt != edge_map.size())
//            std::cout << "Something went wrong with the previous object... "<< std::endl;

        // Second traversal to update the geometry.
        // We run once again through the facets of the HalfEdge data structure.
        for( size_t i_face = 0; i_face < vec_faces.size(); ++i_face )
        {
            carve::mesh::Face<3>* face = vec_faces[i_face];
            // iterate through edges:
            carve::mesh::Edge<3>* edge = face->edge;
            do
            {
                pt1 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v1() );
//                        pt1 = Scale_3D_point(pt1, size_scale);
                pt2 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v2() );
//                        pt2 = Scale_3D_point(pt2, size_scale);

//                        if ((pt1.x() < 1.e-10) &&
//                                (pt1.y() < 1.e-10) &&
//                                (pt1.z() < 1.e-10) )
//                            std::cout << "WTFFFFFF??? " << pt1 << "\t" << pt2 << std::endl;

                // Get the dart associated to the Halfedge
//                        if (!Points_are_eps_close( pt1, pt2 ) )
                {
                    vec_dart ds = edge_map[ std::string( get_unique_key_pt(pt1) +
                                                         get_unique_key_pt(pt2) ) ];

                    for (uint i=0; i<ds.size(); i++)
                    {
                        d = ds[i];

                        if (alcc.template attribute<3>(d) == LCC::null_handle )
                        {
                            alcc.set_vertex_attribute (d, alcc.create_vertex_attribute( pt1 ));
                            res.push_back(d);
                        }
                    }
                }

                // increment edge
                edge = edge->next;
            }
            while( edge != face->edge );
        }

        for (uint i=0; i<res.size(); i++)
        {
            if ( res[i] == LCC::null_handle || alcc.template attribute<0>(res[i]) == LCC::null_handle )
            {
                std::cout << "\tDAMMIT STILL ONE NULL DART HERE!!!" << std::endl;
                alcc.erase_dart( res[i]);
            }
            else
                final_res.push_back( res[i] );
        }


        if (sew2 && final_res.size()>0)
        {
            one_res.push_back( final_res[0] );
            return one_res;
        }
        else
            return final_res;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extracts and simplify the geometry of an IFC object (if it can be casted as an IfcProduct).
    /// method 1 = AABB, method 2 = OBB.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class AnyIfcObject> vec_dart Get_Simplified_Geometry_from_IfcObject(LCC& alcc, shared_ptr<AnyIfcObject>& ifcobject,
                                                                                 std::map<std::string, shared_ptr<ProductShapeData> >& map_of_ifc_shapes,
                                                                                 int method)
    {
        vec_dart res;

        std::string semClass = std::string("Spl_") + ifcobject->className();

        shared_ptr<IfcGloballyUniqueId> iGUID = ifcobject->m_GlobalId;
        std::string GUID (iGUID->m_value.begin(), iGUID->m_value.end());
//        std::cout << "\nThe product will be a: " << label << ", GUID: " << GUID << std::endl;

        std::string Name = "",
                    LongName = "";
        if (ifcobject->m_Name)
        {
            shared_ptr<IfcLabel> iName = ifcobject->m_Name;
            Name = std::string (iName->m_value.begin(), iName->m_value.end());
        }

        shared_ptr<IfcSpace> space = dynamic_pointer_cast<IfcSpace> (ifcobject);
        if (space && space->m_LongName)
        {
            shared_ptr<IfcLabel> iName = space->m_LongName;
            LongName = std::string (iName->m_value.begin(), iName->m_value.end());
        }

        // Check if the "AnyIfcObject" can be casted into an IfcProduct
        shared_ptr<IfcProduct> product = dynamic_pointer_cast<IfcProduct> (ifcobject);

//                idx 0 = EntityId,
//                idx 1 = ClassName,
//                idx 2 = GlobalUId,
//                idx 3 = Name,
//                idx 4 = LongName,
//                idx 5 = ContainingStorey
        std::vector<std::string> ifcinfo(6, "");
        ifcinfo[0] = std::to_string(ifcobject->m_entity_id);
        ifcinfo[1] = semClass;
        ifcinfo[2] = GUID;
        ifcinfo[3] = Name;
        ifcinfo[4] = LongName;
        ifcinfo[5] = Get_storey_containing_IfcObject(ifcobject);

        double r = 130,
               g = 130,
               b = 130;

        if (product)
        {
//            shared_ptr<BuildingModel> ifc_model( new BuildingModel() );
//            shared_ptr<GeometryConverter> geometry_converter(new GeometryConverter( ifc_model ));
//            shared_ptr<ProductShapeData> product_shape ( new ProductShapeData() );
//            shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(product);
//            product_shape->m_ifc_object_definition = obj_def;

//            geometry_converter->convertIfcProductShape( product_shape );

            std::vector<shared_ptr<carve::mesh::MeshSet<3> > > all_meshsets;
            // Dart of the resulting simplified 3-cell
            Dart_handle d = LCC::null_handle;

            std::string pGuid;
            if (product->m_GlobalId){
                std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converterX;
                pGuid = converterX.to_bytes(product->m_GlobalId->m_value);
            }

            auto it_shape = map_of_ifc_shapes.find( pGuid );
            if ( it_shape != map_of_ifc_shapes.end() )
            {
                shared_ptr<ProductShapeData> product_shape = it_shape->second;
//                std::wcout << "vec shapes: " << product_shape->m_vec_representations.size() << std::endl;


                // Collect all the geometries of the IfcProduct to make one simplified shape
                // for each IfcProduct, there can be mulitple geometric representation items:
                std::vector<shared_ptr<RepresentationData> >& vec_rep = product_shape->m_vec_representations;
                for( size_t i_rep = 0; i_rep < vec_rep.size(); ++i_rep )
                {
                    shared_ptr<RepresentationData>& rep_data = vec_rep[i_rep];
                    std::vector<shared_ptr<ItemShapeData> >& vec_item_data = rep_data->m_vec_item_data;
//                    std::wcout << "\t# items: " << vec_item_data.size() << std::endl;

                    for( size_t i_item = 0; i_item < vec_item_data.size(); ++i_item )
                    {
                        // every item can have several meshsets...
                        shared_ptr<ItemShapeData>& item_data = vec_item_data[i_item];
                        item_data->applyTransformToItem( product_shape->getTransform() );

                        // Try to get one color if any (the last color found will be the final one)
                        if (item_data->m_vec_item_appearances.size()>0)
                        {
                            shared_ptr<AppearanceData> appdata = item_data->m_vec_item_appearances[0];
                            r = appdata->m_color_ambient.m_r * 255;
                            g = appdata->m_color_ambient.m_g * 255;
                            b = appdata->m_color_ambient.m_b * 255;
                        }

                        // closed meshsets
                        std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets = item_data->m_meshsets;
                        //                    std::wcout << "\t# meshesets: " << vec_item_meshsets.size() << std::endl;
                        all_meshsets.insert( all_meshsets.end(),  vec_item_meshsets.begin(), vec_item_meshsets.end());

                        // open meshsets
                        std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets_op = item_data->m_meshsets_open;
                        //                    std::wcout << "\t# open meshesets: " << vec_item_meshsets_op.size() << std::endl;
                        all_meshsets.insert( all_meshsets.end(),  vec_item_meshsets_op.begin(), vec_item_meshsets_op.end());
                    }
                }
            }

//            std::wcout << "Total number of Meshsets obtained: " << all_meshsets.size() << std::endl;

            // Simplify all points of meshsets in a single BB
            if (all_meshsets.size() > 0)
            {
                std::vector<Point> all_pts = Get_Vertices_from_Meshsets( all_meshsets );
                assert(all_pts.size() > 0);

                /// AABB based simplification
                if (method == 1)
                {
                    Bbox_3 bb_all_pts = Get_Bbox_Vector_Point( all_pts );
                    if ( !Bbox_is_degenerate(bb_all_pts) )
                    {
//                        std::wcout << "Good Bbox " << std::endl;
                        d = Create_bbox_volume(alcc, bb_all_pts );
                    }
                    else
                        std::wcout << "Degenerate volume! Furniture " << product->m_Name->m_value << std::endl;
                }
                /// OBB based simplification
                else if (method == 2)
                    d = OBB_Pointset_Simplification(alcc, all_pts);
//                    d = Convex_hull_Pointset_Simplification(alcc, all_pts);


                /// Keep the dart of the resulting volume
                if ( d != LCC::null_handle )
                    res.push_back( d );
            }

            // Simplify each meshet separately
//            for (uint i_m = 0; i_m < all_meshsets.size(); i_m++)
//            {
//                std::vector<Point> pts = Get_Vertices_from_Single_Meshset( all_meshsets[i_m] );
//                assert(pts.size() > 0);

//                /// AABB based simplification
//                if (method == 1)
//                {
//                    Bbox_3 bb_all_pts = Get_Bbox_Vector_Point( pts );
//                    if ( !Bbox_is_degenerate(bb_all_pts) )
//                    {
//                        std::wcout << "Good Bbox " << std::endl;
//                        d = Create_bbox_volume(alcc, bb_all_pts );
//                    }
//                    else
//                        std::wcout << "Degenerate volume! Furniture " << product->m_Name->m_value << std::endl;
//                }
//                /// OBB based simplification
//                else if (method == 2)
//                    d = OBB_Pointset_Simplification(alcc, pts);


//                /// Keep the dart of the resulting volume
//                if ( d != LCC::null_handle )
//                    res.push_back( d );
//            }


        }
        else
            std::cout << "Not a valid IfcProduct! No geometry extracted!" << std::endl;

        for(uint i=0; i<res.size(); i++)
        {
            alcc.set_attribute<3>( res[i], alcc.create_attribute<3>() );
            alcc.info<3>(res[i]).set_semClass(semClass);
            alcc.info<3>(res[i]).set_label(LongName + "_" + Name);
//            alcc.info<3>(res[i]).set_vol_label(volume_label);
            alcc.info<3>(res[i]).set_id(std::to_string(vol_counter++));
            alcc.info<3>(res[i]).set_color( int(r), int(g), int(b) );

        }

        volume_label++;
        return res;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Identify in which storey a a given IfcProduct is contained
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class AnyIfcObject> std::string Get_storey_containing_IfcObject (shared_ptr<AnyIfcObject> &myIfcObject)
    {
        if (myIfcObject)
        {
            if (strcmp(myIfcObject->className(), "IfcSpace") == 0)
            {
                std::vector<weak_ptr<IfcRelAggregates> > vec_decomposes = myIfcObject->m_Decomposes_inverse;
//                std::wcout << "This " << myIfcObject->className() << " is contained in "
//                           << vec_decomposes.size() << " Structure(s)." << std::endl;
                for (uint i=0; i<vec_decomposes.size(); i++)
                {
                    shared_ptr<IfcRelAggregates> rel_agg (vec_decomposes[i]);
                    if(rel_agg)
                    {
                        shared_ptr<IfcObjectDefinition> objdef = rel_agg->m_RelatingObject;
                        if (strcmp(objdef->className(), "IfcSpace") == 0)
                        {
                            return Get_storey_containing_IfcObject(objdef);
                        }
                        else
                        {
//                            std::wcout << "\tOne structure is a(n) " << objdef->className();
                            if ( objdef->m_Name )
                            {
//                                std::wcout << " named " << objdef->m_Name->m_value;
//        //                        if ( objdef->m_LongName )
//        //                            std::wcout << " and (long)named " << objdef->m_LongName->m_value;
//        //                        if ( objdef->m_Description )
//        //                            std::wcout << " with the description: " << objdef->m_Description->m_value;
//                                std::wcout << std::endl;
                                return std::string(objdef->m_Name->m_value.begin(), objdef->m_Name->m_value.end());
                            }
                        }
                    }
                }
            }
            else
            {
                shared_ptr<IfcElement> myelem = dynamic_pointer_cast<IfcElement>( myIfcObject );
                if( myelem )
                {
                    std::vector<weak_ptr<IfcRelContainedInSpatialStructure> > vec_containedIn = myelem->m_ContainedInStructure_inverse;
//                    std::wcout << "This " << myIfcObject->className() << " is contained in "
//                               << vec_containedIn.size() << " Structure(s)." << std::endl;

                    for(uint i=0; i<vec_containedIn.size(); i++)
                    {
                        shared_ptr<IfcRelContainedInSpatialStructure> OneStruct(vec_containedIn[i]);
                        if (OneStruct && OneStruct->m_RelatingStructure)
                        {
                            shared_ptr<IfcSpatialElement> myContainer = OneStruct->m_RelatingStructure;
                            if (strcmp(myContainer->className(), "IfcSpace") == 0)
                            {
                                return Get_storey_containing_IfcObject(myContainer);
                            }
                            else
                            {
//                                std::wcout << "\tOne structure is a(n) " << myContainer->className();
                                if ( myContainer->m_Name )
                                {
//                                    std::wcout << " named " << myContainer->m_Name->m_value;
//        //                            if ( myContainer->m_LongName )
//        //                                std::wcout << " and (long)named " << myContainer->m_LongName->m_value;
//        //                            if ( myContainer->m_Description )
//        //                                std::wcout << " with the description: " << myContainer->m_Description->m_value;
//                                    std::wcout << std::endl;
                                    return std::string(myContainer->m_Name->m_value.begin(), myContainer->m_Name->m_value.end());
                                }
                            }
                        }
//                        else
//                            std::wcout << "\tCould not get one structure..." << std::endl;
                    }
                }
//                else
//                    std::wcout << "This " << myIfcObject->className() << " could not be casted into an IfcElement!" << std::endl;
            }
        }

        return "";
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extracts the geometry of an Ifc object as long as the latter can be casted into an IfcProduct.
    /// Returns false if the casting fails.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    template<class AnyIfcObject> vec_dart Get_Geometry_from_IfcObject( LCC& alcc, shared_ptr<AnyIfcObject>& anyifcobj,
                                                                       std::map<std::string, shared_ptr<ProductShapeData> >& map_of_ifc_shapes,
                                                                       bool simplify = true, bool no_sew2 = false, bool only_closed_meshsets = true,
                                                                       bool correct_wrong_spaces = false)
    {
        vec_dart res, cur_res;
//                int label = 4;

//        __try
//        {
//            try
//            {
                // Check if the "AnyIfcObject" can be casted into an IfcProduct
                std::string semClass = anyifcobj->className();
                shared_ptr<IfcGloballyUniqueId> iGUID = anyifcobj->m_GlobalId;
                std::string GUID (iGUID->m_value.begin(), iGUID->m_value.end());

                std::string Name = "",
                            LongName = "";
                if (anyifcobj->m_Name)
                {
                    shared_ptr<IfcLabel> iName = anyifcobj->m_Name;
                    Name = std::string (iName->m_value.begin(), iName->m_value.end());
                }

                shared_ptr<IfcSpace> space = dynamic_pointer_cast<IfcSpace>(anyifcobj);
                if (space && space->m_LongName)
                {
                    shared_ptr<IfcLabel> iName = space->m_LongName;
                    LongName = std::string (iName->m_value.begin(), iName->m_value.end());
                }

                shared_ptr<IfcProduct> product = dynamic_pointer_cast<IfcProduct> (anyifcobj);

//                idx 0 = EntityId,
//                idx 1 = ClassName,
//                idx 2 = GlobalUId,
//                idx 3 = Name,
//                idx 4 = LongName,
//                idx 5 = ContainingStorey
                std::vector<std::string> ifcinfo(6, "");
                ifcinfo[0] = std::to_string(anyifcobj->m_entity_id);
                ifcinfo[1] = semClass;
                ifcinfo[2] = GUID;
                ifcinfo[3] = Name;
                ifcinfo[4] = LongName;
                ifcinfo[5] = Get_storey_containing_IfcObject(anyifcobj);


//                std::cout << "\nThe product is a: " << label << ", GUID: " << GUID << std::endl;


                /// If the product at this stage is valid
                if (product)
                {
                    /// In case of simplification required
                    int isdoor = strcmp(semClass.c_str(), "IfcDoor"),
                        iswin = strcmp(semClass.c_str(), "IfcWindow"),
                        isVirtualElem = strcmp(semClass.c_str(), "IfcVirtualElement");

                    if ( simplify && (isdoor == 0 || iswin == 0 || isVirtualElem == 0) )
                    {
                        semClass = std::string("Spl_") + semClass;
                        ifcinfo[1] = semClass;
//                        return Get_Simplified_Geometry_from_IfcObject(alcc, anyifcobj, map_of_ifc_shapes, 2);
                        std::vector<weak_ptr<IfcRelFillsElement> > fillings;
                        if (isdoor == 0)
                        {
                            shared_ptr<IfcDoor> door = dynamic_pointer_cast<IfcDoor> (anyifcobj);
                            fillings = door->m_FillsVoids_inverse;
                            //                        label = 5;
                        }
                        else if (iswin == 0)
                        {
                            shared_ptr<IfcWindow> window = dynamic_pointer_cast<IfcWindow> (anyifcobj);
                            fillings = window->m_FillsVoids_inverse;
                            //                        label = 6;
                        }
                        else if (isVirtualElem == 0)
                        {
                            shared_ptr<IfcVirtualElement> virtualElem = dynamic_pointer_cast<IfcVirtualElement> (anyifcobj);

                            std::vector<weak_ptr<IfcRelVoidsElement>> voidElems = virtualElem->m_HasOpenings_inverse;
                            std::wcout << "\tNumber of related opening(s) found: " << voidElems.size() << std::endl;
                            if ( voidElems.size() > 0 ){
                                shared_ptr<IfcRelVoidsElement> vElem = voidElems[0].lock();
                                shared_ptr<IfcFeatureElementSubtraction> vElemOps = vElem->m_RelatedOpeningElement;

                                fillings = vElemOps->m_FillsVoids_inverse;
//                                fillings = virtualElem->m_FillsVoids_inverse;
                            }
                        }

                        // if there is a filling, product
                        if (fillings.size() == 1)
                        {
//                            std::wcout << "\tThe " << semClass.c_str() << " got a filling space..." << std::endl;
                            for (uint i=0; i<fillings.size(); i++)
                            {
                                shared_ptr<IfcRelFillsElement> fill_elem (fillings[i]);
                                shared_ptr<IfcOpeningElement> opening = fill_elem->m_RelatingOpeningElement;
                                product = dynamic_pointer_cast<IfcProduct> (opening);
                            }
                        }
                        else{
                            if ( fillings.size() > 0 )
                                std::wcout << "\tThe " << semClass.c_str() << " has MANY FILLINGS... How?? xD" << std::endl;

                            std::wcout << "\tAttempting simplification of " << semClass.c_str() << " " << LongName.c_str() << " " << Name.c_str() << std::endl;
                            if ( isVirtualElem != 0 )
                                return Get_Simplified_Geometry_from_IfcObject(alcc, anyifcobj, map_of_ifc_shapes, 2);
                        }
                    }


//                    std::wcout << "\tSo let's keep going!" << std::endl;

////                    shared_ptr<BuildingModel> ifc_model( new BuildingModel() );
////                    shared_ptr<GeometryConverter> geometry_converter(new GeometryConverter( ifc_model ));
//                    shared_ptr<ProductShapeData> product_shape ( new ProductShapeData() );
//                    shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(product);
//                    product_shape->m_ifc_object_definition = obj_def;

//                    geometry_converter->convertIfcProductShape( product_shape );
//                    //            std::wcout << "vec shapes: " << product_shape->m_vec_representations.size() << std::endl;

                    if (strcmp(semClass.c_str(), "IfcVirtualElement") == 0)
                        std::wcout << "IfcVirtualElement could be casted into IfcProduct" << std::endl;

                    std::string pGuid;
                    if (product->m_GlobalId){
                        std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converterX;
                        pGuid = converterX.to_bytes(product->m_GlobalId->m_value);
                    }


                    // If the shape data is found in the geometry converter used it
                    shared_ptr<ProductShapeData> product_shape (new ProductShapeData() );
                    auto it_shape = map_of_ifc_shapes.find( pGuid );
                    if ( it_shape != map_of_ifc_shapes.end() ){
                        product_shape = it_shape->second;
                    }
                    // otherwise try to get the shape data
                    else{
                        shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(product);
                        product_shape->m_ifc_object_definition = obj_def;
                        geomConverter->convertIfcProductShape( product_shape );
                    }


                    if ( product_shape ){
//                        unsigned int sz = product_shape->getAppearances().size();
//                        if ( sz == 0 )
//                            std::wcout << "No appearance data in ProductShape!" << std::endl;

                        // for each IfcProduct, there can be mulitple geometric representation items:
                        std::vector<shared_ptr<RepresentationData> >& vec_rep = product_shape->m_vec_representations;
                        for( size_t i_rep = 0; i_rep < vec_rep.size(); ++i_rep )
                        {
                            shared_ptr<RepresentationData>& rep_data = vec_rep[i_rep];

                            std::vector<shared_ptr<ItemShapeData> >& vec_item_data = rep_data->m_vec_item_data;
                            //                std::wcout << "\t# items: " << vec_item_data.size() << std::endl;

                            for( size_t i_item = 0; i_item < vec_item_data.size(); ++i_item )
                            {
                                // every item can have several meshsets...
                                shared_ptr<ItemShapeData>& item_data = vec_item_data[i_item];
                                item_data->applyTransformToItem( product_shape->getTransform() );

//                                std::wcout << "Appearance data?? " << item_data->m_vec_item_appearances.size() << std::endl;


                                // closed meshsets
                                std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets = item_data->m_meshsets;
                                //                    std::wcout << "\t# meshesets: " << vec_item_meshsets.size() << std::endl;
                                for( size_t i_meshset = 0; i_meshset < vec_item_meshsets.size(); ++i_meshset )
                                {
                                    shared_ptr<carve::mesh::MeshSet<3> >& meshset = vec_item_meshsets[i_meshset];
                                    std::vector<carve::mesh::Mesh<3>* >& vec_meshes = meshset->meshes;
//                                    std::wcout << "\t#closed meshes: " << vec_meshes.size() << std::endl;

                                    for( size_t i_mesh = 0; i_mesh < vec_meshes.size(); ++i_mesh )
                                    {
                                        carve::mesh::Mesh<3>* mesh = vec_meshes[i_mesh];
                                        cur_res = osg_mesh2lcc<LCC>( alcc, mesh, !no_sew2 );

                                        for(uint i=0; i<cur_res.size(); i++)
                                        {
                                            alcc.set_attribute<3>( cur_res[i], alcc.create_attribute<3>());
                                            alcc.info<3>(cur_res[i]).set_label(LongName + "_" + Name);
                                            alcc.info<3>(cur_res[i]).set_semClass(semClass);
                                            alcc.info<3>(cur_res[i]).set_id(std::to_string(vol_counter++));

                                            double r = 130,
                                                   g = 130,
                                                   b = 130;

                                            if (item_data->m_vec_item_appearances.size()>0)
                                            {
                                                shared_ptr<AppearanceData> appdata = item_data->m_vec_item_appearances[0];
                                                r = appdata->m_color_ambient.m_r * 255;
                                                g = appdata->m_color_ambient.m_g * 255;
                                                b = appdata->m_color_ambient.m_b * 255;

            //                                    std::wcout << "color found: " << r << ", " << g << ", " << b << std::endl;
                                            }
                                            alcc.info<3>(cur_res[i]).set_color( int(r), int(g), int(b) );
                                        }

                                        res.insert( res.end(), cur_res.begin(), cur_res.end() );
                                    }
                                }

                                // open meshsets
                                if (!only_closed_meshsets)
                                {
                                    std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets_op = item_data->m_meshsets_open;
                                    //                        std::wcout << "\t# meshesets open?: " << vec_item_meshsets_op.size() << std::endl;
                                    for( size_t i_meshset = 0; i_meshset < vec_item_meshsets_op.size(); ++i_meshset )
                                    {
                                        shared_ptr<carve::mesh::MeshSet<3> >& meshset = vec_item_meshsets_op[i_meshset];
                                        std::vector<carve::mesh::Mesh<3>* >& vec_meshes = meshset->meshes;
//                                        std::wcout << "\t#open meshes: " << vec_meshes.size() << std::endl;

                                        for( size_t i_mesh = 0; i_mesh < vec_meshes.size(); ++i_mesh )
                                        {
                                            carve::mesh::Mesh<3>* mesh = vec_meshes[i_mesh];
                                            cur_res = osg_mesh2lcc<LCC>( alcc, mesh, !no_sew2 );

                                            for(uint i=0; i<cur_res.size(); i++)
                                            {
                                                alcc.set_attribute<3>( cur_res[i], alcc.create_attribute<3>());
                                                alcc.info<3>(cur_res[i]).set_label(LongName + "_" + Name);
                                                alcc.info<3>(cur_res[i]).set_semClass(semClass);
                                                alcc.info<3>(cur_res[i]).set_id(std::to_string(vol_counter++));

                                                double r = 130,
                                                       g = 130,
                                                       b = 130;

                                                if (item_data->m_vec_item_appearances.size()>0)
                                                {
                                                    shared_ptr<AppearanceData> appdata = item_data->m_vec_item_appearances[0];
                                                    r = appdata->m_color_ambient.m_r * 255;
                                                    g = appdata->m_color_ambient.m_g * 255;
                                                    b = appdata->m_color_ambient.m_b * 255;

                //                                    std::wcout << "color found: " << r << ", " << g << ", " << b << std::endl;
                                                }
                                                alcc.info<3>(cur_res[i]).set_color( int(r), int(g), int(b) );

                                            }

                                            res.insert( res.end(), cur_res.begin(), cur_res.end() );
                                        }
                                    }
                                }
                            }
                        }

                        map_ifc_shapes[ pGuid ] = product_shape;
                    }
                }
                else
                    std::cout << "Not a valid IfcProduct! No geometry extracted!" << std::endl;

//                // Tries to repair multiple volumes and obtain only a single one instead
//                if (correct_wrong_spaces && (res.size() > 1) && strcmp(label.c_str().c_str(), "IfcSpace") == 0)
//                {
//                    std::cout << res.size() << " IfcSpace volumes detected as one single " << label << std::endl;
//                    std::cout << "**************** Beginning of volume repairing ****************" << label << std::endl;
//                    //            res = CGAL::BSP_volume_healing_poly_soup(alcc, res);
//                    res = BSP_volume_healing_poly_soup_ex(alcc, res);
//                    std::cout << res.size() << " volume(s) obtained after repairing!" << label << std::endl;

//                    for(uint i=0; i<res.size(); i++)
//                    {
//                        alcc.set_attribute<3>( res[i], alcc.create_attribute<3>());
////                        alcc.info<3>(cur_res[i]).set_label(LongName + "_" + Name);
////                        alcc.info<3>(cur_res[i]).set_semClass(semClass);
////                        alcc.info<3>(res[i]).set_id(vol_counter++);
////                        alcc.info<3>(res[i]).set_IfcInfo(ifcinfo);
//                    }
//                }
//                else
//                {
//                    std::cout << res.size() << " volumes detected as one single " << label << std::endl;
////                    std::cout << "No volume repairing..." << label << std::endl;
//                }


                // All the 3-cells resulting from this geometry should belong to the same volume
//                for(uint i=0; i<res.size(); i++)
//                {
//                    alcc.info<3>(res[i]).set_vol_label(volume_label);
//                    alcc.info<3>(res[i]).set_IfcInfo(ifcinfo);
//                }

                volume_label++;
//            }
//            catch(...)
//            {
//                std::cout << "******************** Something went wrong on this building part..." << std::endl;
//            }
//        }
//        __except(GetExceptionCode())
//        {
//            std::cout << "******************** Something went wrong on this building part..." << std::endl;
//        }

        return res;
    }


    bool Get_Building_Storeys (shared_ptr<BuildingModel>&, std::vector< shared_ptr<IfcBuildingStorey> >&);

    bool Get_IfcSpaces_of_Storey (shared_ptr<IfcBuildingStorey>&, std::vector< shared_ptr<IfcSpace> >&);

    void Get_IfcElements_around_IfcSpace (shared_ptr<IfcSpace>&, std::map< const char*,
                                          std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& );
    void Get_IfcElements_around_IfcDoor (shared_ptr<IfcDoor>&, std::map< const char*,
                                          std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& );

    void Get_all_IfcProducts_related_to_IfcBuildingStorey (shared_ptr<IfcBuildingStorey>&,
                                                           std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >& );

    void Get_all_elements_in_IfcSpace (shared_ptr<IfcSpace>&, std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >&,
                                       std::map< const char*, bool, LCCtools::cmp_const_char >&);
    void Get_all_elements_in_IfcSpace (shared_ptr<IfcSpace>&, std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >&);

    void Get_Openings_around_IfcSpace (shared_ptr<IfcSpace>&, std::vector< shared_ptr<IfcElement> >&,
                                          std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& );


#endif
    bool Extract_Real_Free_Space ( LCC&, Dart_handle&, vec_dart&, R_Space&/*, Nef_polyhedron&*/ );

    bool Extract_fspace_from_free_Space (LCC&, Dart_handle&,  Nef_polyhedron&, vec_dart&, vec_dart&, std::string );

    bool Extract_fspace_from_free_Space_without_furnitures (LCC&, Dart_handle&, vec_dart&, vec_dart&, std::string);

    void Get_Buffers_for_Openings(LCC&, R_Space&, vec_dart&, std::vector<Sub_Space>&, double buffer_size = 0.6);

    void Get_Buffers_or_Functional_for_O_Spaces (LCC&, R_Space&, std::vector<O_Space*> &, std::vector<double> &,
                                                 std::vector<Sub_Space>&, std::string&, double buffer_size = 0.6);

//    vec_dart Rebuild_faces_after_Arrangement_Correction(LCC&, arr_tuple&);

    Dart_handle Create_3cell_from_set_of_2cells(LCC&, vec_dart&);

    void Clean_Vol_Boundaries( LCC&, Dart_handle& );

    double_pair_dart_vec3d Extrude_opening_in_thickeness_direction(LCC&, Dart_handle&, double, bool increase = true);

    void Extrude_face_in_vector_direction( LCC&, Dart_handle&, Vector&, double);

    void Extrude_box_volume( LCC&, Dart_handle&, double, bool increase = true);

    void Make_contact_between_volA_and_volB( LCC&, Dart_handle&, Dart_handle&, double dist = 0.1, bool handle_floating = false);

    void Aggregate_Set_of_3cells (LCC&, vec_dart&, vec_dart&, bool buffer = false, double buffer_size = 0.6);

    void Aggregate_Set_of_O_Spaces (LCC&, std::vector<O_Space>&, std::vector<O_Space>&, int loop = 0,
                                    bool buffer = false, double buffer_size = 0.6, int mode = 0);

#ifdef IFCPP_ON
    void ifc2lcc(LCC&, std::wstring&);

    // UI functions
    bool pre_load_ifc(std::wstring&);

    void fully_load_ifc2lcc(LCC&, std::map<int, shared_ptr<BuildingEntity> >&,
                            bool simplify = false, bool no_sew2 = false, bool only_closed_meshes = false,
                            bool correct_wrong_spaces = false);

    void load_all_in_storey_ifc2lcc(LCC&, shared_ptr<IfcBuildingStorey>&,
                                    bool simplify = false, bool no_sew2 = false, bool only_closed_meshes = false,
                                    bool correct_wrong_spaces = false);

    void load_selected_in_storey_ifc2lcc(LCC&, shared_ptr<IfcBuildingStorey>&, std::vector<const char*>&,
                                         bool simplify = false, bool no_sew2 = false, bool only_closed_meshes = false,
                                         bool correct_wrong_spaces = false);
#endif


    void my_load_off(LCC&, std::istream&);
}

#endif
