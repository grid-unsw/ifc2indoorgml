#include "IFC_Importer.h"

// Global variables
Vector V_height (0.0, 0.0, 1.0);
//Point global_shift_pt (0.0, 0.0, 0.0) /*(638990.0, 240268.0, 10480.0)*/;
double scale_f = 1.0,
       size_scale = 1.0;
ET    eps1 = 0.001;
double      eps2_default = 0.001,
            coplanarity_default = 0.01;

//extern double eps2,
//       eps_angle;


uint volume_label;
uint vol_counter = 0;

Converter_1 exa2inexa;
Converter_2 inexa2exa;


#ifdef IFCPP_ON
    /// Abdou's change
    shared_ptr<GeometryConverter> geomConverter;
    std::map<std::string, shared_ptr<ProductShapeData> > map_ifc_shapes;
    std::map<int, shared_ptr<BuildingEntity> > map_ifc_entities;
    std::vector< shared_ptr<IfcBuildingStorey> > ifc_storey_list;
    std::vector<std::vector< shared_ptr<IfcSpace> > > spaces_of_storeys;
#endif

namespace CGAL
{

/* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
///                                             Convenient Functions
/* /////////////////////////////////////////////////////////////////////////////////////////////////////// */

Point Scale_3D_point (Point& pt, double ss)
{
    double x = pt.x() * ss,
           y = pt.y() * ss,
           z = pt.z() * ss;

    return Point(x,y,z);
}

std::string get_unique_key_pt( Point& pt )
{
    return  std::string ( std::to_string( int(pt.x()/eps2) )
                        + std::to_string( int(pt.y()/eps2) )
                        + std::to_string( int(pt.z()/eps2) ));
}


/* ////////////////////////////////////////////////////////////////////////////// */
/// Return the angle between two 2D vectors v1 and v2, in the range [0°, 180°]
/* ////////////////////////////////////////////////////////////////////////////// */
double compute_3d_angle(const Vector& v1, const Vector& v2)
{
    double a = to_double( (v1*v2) / ( sqrt(v1.squared_length()) * sqrt(v2.squared_length()) ) ) ;

    if (isnan(a) || isinf(a))
        a = 1.0;
//        std::cout << "\nRaw angle:\t" << a <<std::endl;

    if (a < -1.0)
        return acos(-1.0)/M_PI*180.0;
    else if (a > 1.0)
        return acos(1.0)/M_PI*180.0;
    else
        return acos(a)/M_PI*180.0;
}

/* //////////////////////////////////////////////////////////////////////////////// */
/// Returns true if v1 and v2 are (eps_angle) colinear
/* //////////////////////////////////////////////////////////////////////////////// */
bool vectors_are_eps_colinear (Vector v1, Vector v2)
{
    double angle = compute_3d_angle (v1, v2);
    if ( (angle < (0.0 + eps_angle) || angle > (360.0 - eps_angle) ) ||
         (angle < (180.0 + eps_angle) && angle > (180.0 - eps_angle) ) )
        return true;
    return false;
}

/* ////////////////////////////////////////////////////////// */
/// Returns true if the Bbox is degenerate (flat volume)
/* ////////////////////////////////////////////////////////// */
bool Bbox_is_degenerate(Bbox_3 bb)
{
    if ( fabs(bb.xmax() - bb.xmin()) <= eps2
         || fabs(bb.ymax() - bb.ymin()) <= eps2
            || fabs(bb.zmax() - bb.zmin()) <= eps2 )
        return true;

    return false;
}

/* ///////////////////////////////////////////////////////////////// */
/// Resize the Bbox by increasing/decreasing it in all directions
/* ///////////////////////////////////////////////////////////////// */
void Rescale_Bbox (Bbox_3 &bbox, double &scale, bool increase)
{
    if (!increase)
        scale = -1 * scale;

    double xmin, xmax,
           ymin, ymax,
           zmin, zmax;

    xmin = bbox.xmin() - scale;
    ymin = bbox.ymin() - scale;
    zmin = bbox.zmin() - scale;

    xmax = bbox.xmax() + scale;
    ymax = bbox.ymax() + scale;
    zmax = bbox.zmax() + scale;

    bbox = Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
}

/* ////////////////////////////////////////////////////////// */
/// Create a given Bounding Box as a 3-cell in the LCC
/* ////////////////////////////////////////////////////////// */
Dart_handle Create_bbox_volume (LCC& alcc, Bbox_3& bb)
{
    Point   p0 (bb.xmax(), bb.ymax(), bb.zmin()),
            p1 (bb.xmin(), bb.ymax(), bb.zmin()),
            p2 (bb.xmin(), bb.ymin(), bb.zmin()),
            p3 (bb.xmax(), bb.ymin(), bb.zmin()),
            p4 (bb.xmax(), bb.ymin(), bb.zmax()),
            p5 (bb.xmax(), bb.ymax(), bb.zmax()),
            p6 (bb.xmin(), bb.ymax(), bb.zmax()),
            p7 (bb.xmin(), bb.ymin(), bb.zmax());

    return alcc.make_hexahedron(p0,p1,p2,p3,p4,p5,p6,p7);
}


/* ////////////////////////////////////////////////////////// */
/// Collect one dart per face of a given volume in one vector
/* ////////////////////////////////////////////////////////// */
vec_dart one_dart_per_face_in_volume (LCC& alcc, Dart_handle vol)
{
    vec_dart res;

    // For each edge of the volume
    typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
            itr = alcc.one_dart_per_incident_cell<2,3>(vol).begin(),
            itrend = alcc.one_dart_per_incident_cell<2,3>(vol).end();
    for (; itr != itrend; itr++)
    {
        res.push_back( itr );
    }

    return res;
}


/* ////////////////////////////////////////////////////////// */
/// Compute the volume of the 3-cell containing d
/* ////////////////////////////////////////////////////////// */
double Get_volume_magnitude(LCC& alcc, Dart_handle d)
{
    double vol = 0.0;

    typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
            itr = alcc.one_dart_per_incident_cell<2,3>(d).begin(),
            itrend = alcc.one_dart_per_incident_cell<2,3>(d).end();
    for (; itr != itrend; itr++)
    {
        Point pt = alcc.point(itr);
        Vector N = LCCtools::poly_normal(alcc, itr, 1);

        // plane in which lies the polygon
        Plane plan(pt, N);

        // project the points on the plane of the polygon to get their 2D coordinates
        std::vector<Point2d> poly_2d;
        typename LCC::Base::Dart_of_orbit_range<1>::iterator it(alcc, itr);
        for ( ; it.cont(); ++it )
        {
            poly_2d.push_back( plan.to_2d( alcc.point(it) ) );
        }
        // the point n is the point 0
        poly_2d.push_back( plan.to_2d( alcc.point(itr) ) );

        double a, A = 0.0, Cx = 0.0, Cy = 0.0;
        for (uint j = 0; j<poly_2d.size()-1; j++)
        {
            a = to_double( poly_2d[j].x() * poly_2d[j+1].y() )
                    - to_double( poly_2d[j+1].x() * poly_2d[j].y() );

            A += a;
            Cx += a * ( to_double(poly_2d[j].x() + poly_2d[j+1].x() ) );
            Cy += a * ( to_double(poly_2d[j].y() + poly_2d[j+1].y() ) );
        }

        A = A/2.0;
        Cx = Cx / (6.0*A);
        Cy = Cy / (6.0*A);

        Point2d C (Cx, Cy);
        Point C_3d = plan.to_3d(C);

        vol += to_double(C_3d.x() * N.x() + C_3d.y() * N.y() + C_3d.z() * N.z()) * A;

        poly_2d.clear();
    }

    return fabs( vol/(3.0 /** pow(scale_f, 3)*/ ) );
}


/* ////////////////////////////////////////////////////////// */
/// Collect all the faces of a given volume having the same
/// normal direction than the face containing the dart d
/* ////////////////////////////////////////////////////////// */
std::vector <vec_dart > Collect_parallel_faces_in_volume (LCC& alcc, Dart_handle& vol, bool coplanar = false )
{
    vec_dart odpf_vol = one_dart_per_face_in_volume(alcc, vol);

    int mark = alcc.get_new_mark();
    CGAL_assertion( mark!=-1 );

    std::vector < vec_dart > res;
    vec_dart cur;

    Vector N_cur, N_ref;

    for (uint i=0; i<odpf_vol.size(); i++)
    {
        if ( !alcc.is_marked(odpf_vol[i], mark) )
        {
            N_ref = LCCtools::poly_normal(alcc, odpf_vol[i]);
            alcc.mark(odpf_vol[i], mark);
            cur.push_back(odpf_vol[i]);

            // plane containing the current ref face
            Plane p_ref(alcc.point( odpf_vol[i] ), N_ref);

            // Check for each face of the volume
            for (uint j=i+1; j<odpf_vol.size(); j++)
            {
                if ( !alcc.is_marked(odpf_vol[j], mark) )
                {
                    N_cur = LCCtools::poly_normal(alcc, odpf_vol[j]);
                    Vector V( alcc.point(odpf_vol[j]), p_ref.projection( alcc.point(odpf_vol[j]) ) );

                    // if the face is coplanar to the ref face
                    if ( coplanar && vectors_are_eps_colinear(N_ref, N_cur)
                            && V.squared_length() <= (coplanarity*coplanarity) )
                    {
                        cur.push_back(odpf_vol[j]);
                        alcc.mark(odpf_vol[j], mark);
                    }
                    // if just parallel
                    else if ( !coplanar && vectors_are_eps_colinear(N_ref, N_cur) )
                    {
                        cur.push_back(odpf_vol[j]);
                        alcc.mark(odpf_vol[j], mark);
                    }
                }
            }

            res.push_back(cur);
            cur.clear();
        }
    }


    alcc.unmark_all(mark);
    alcc.free_mark(mark);

    return res;
}

/* ////////////////////////////////////////////////////////// */
/// Create a new face in the LCC based on a vector of points
/* ////////////////////////////////////////////////////////// */
Dart_handle Insert_new_2cell (LCC& alcc, vec_pt3d& one_face, bool orientation = true)
{
//        std::cout << "\nface size: " << one_face.size() << std::endl;
    if (one_face.size() != 0)
    {
        Dart_handle d = alcc.make_combinatorial_polygon(one_face.size());

         Dart_handle cur_d = d;

         // normal orientation
         if (orientation)
         {
             for(uint i=0; i<one_face.size(); i++)
             {
                 alcc.set_vertex_attribute_of_dart(cur_d, alcc.create_vertex_attribute(one_face[i]) );
                 cur_d=alcc.beta(cur_d,1);
             }
         }

         // opposite orientation
         else
         {
             for(int i=one_face.size()-1; i >= 0; i--)
             {
//                     std::cout << "i = " << i << std::endl;
                 alcc.set_vertex_attribute_of_dart(cur_d, alcc.create_vertex_attribute(one_face[i]) );
                 cur_d=alcc.beta(cur_d,1);
             }
         }

        return d;
    }

    return LCC::null_handle;
}

/* ///////////////////////////////////////////////////////////////////// */
/// Create a new face in the LCC based on a vector of darts (overloaded)
/* ///////////////////////////////////////////////////////////////////// */
Dart_handle Insert_new_2cell (LCC& alcc, vec_dart& one_face, bool orientation = true)
{
//        std::cout << "\nface size: " << one_face.size() << std::endl;
    if (one_face.size() != 0)
    {
        Dart_handle d = alcc.make_combinatorial_polygon(one_face.size());
        Dart_handle cur_d = d;

        // normal orientation
        if (orientation)
        {
            for(uint i=0; i<one_face.size(); i++)
            {
                alcc.set_vertex_attribute_of_dart(cur_d, alcc.create_vertex_attribute(alcc.point(one_face[i])) );
                cur_d=alcc.beta(cur_d,1);
            }
        }

        // opposite orientation
        else
        {
            for(int i=one_face.size()-1; i >= 0; i--)
            {
                //                     std::cout << "i = " << i << std::endl;
                alcc.set_vertex_attribute_of_dart(cur_d, alcc.create_vertex_attribute(alcc.point(one_face[i])) );
                cur_d=alcc.beta(cur_d,1);
            }
        }

        return d;
    }

    return LCC::null_handle;
}

/* //////////////////////////////////////////////////////////////////////////////// */
/// Remove a selected set of 3-cells from the LCC
/* //////////////////////////////////////////////////////////////////////////////// */
void Remove_selected_3_cells (LCC& alcc, vec_dart& darts)
{
    alcc.set_update_attributes(false);
    for(uint i=0; i<darts.size(); i++)
    {
        if (alcc.is_dart_used(darts[i])
                && alcc.is_removable<3>(darts[i]))
        {
//                std::cout << "YES! DART IN USE! (3-removal)" << std::endl;
            alcc.remove_cell<3>( darts[i]);
        }
//            else
//                std::cout << "NO!! DART UNUSED! (3-removal)" << std::endl;
    }
    alcc.set_update_attributes(true);
}


/* ////////////////////////////////////////////////////////////////////////////// */
/// Converts an existing 3-cell into a Polyhedron_3
/* ////////////////////////////////////////////////////////////////////////////// */
bool Convert_3cell_2_exact_polyhedron3( LCC& alcc, Dart_handle d, Polyhedron_3<EK>& P )
{
    try
    {
        CGAL::Build_exact_polyhedron_3_from_lcc PolyLcc( alcc, d );
        P.delegate( PolyLcc );
        if ( P.is_valid() )
        {
            //            std::cout << "\nPolyhedron valid!" << std::endl;
            return true;
        }
        else
        {
            std::cout << "\nInvalid Polyhedron!" << std::endl;
            return false;
        }
    }
    catch (...)
    {
        std::cout << "\nInvalid Polyhedron!" << std::endl;
        return false;
    }
}

/* ////////////////////////////////////////////////////////////////////////////// */
/// Converts an existing 3-cell into a Polyhedron_3
/* ////////////////////////////////////////////////////////////////////////////// */
bool Convert_3cell_2_inexact_polyhedron3( LCC& alcc, Dart_handle d, Polyhedron_3<IK>& P )
{
    try
    {
        Build_inexact_polyhedron_3_from_lcc PolyLcc( alcc, d );
        P.delegate( PolyLcc );
        if ( P.is_valid() )
        {
            //            std::cout << "\nPolyhedron valid!" << std::endl;
            return true;
        }
        else
        {
            std::cout << "\nInvalid Polyhedron!" << std::endl;
            return false;
        }
    }
    catch (...)
    {
        std::cout << "\nInvalid Polyhedron!" << std::endl;
        return false;
    }
}

/* ////////////////////////////////////////////////////////////////////////////// */
/// Converts a 3-cell into a Nef_Polyhedron_3 object
/* ////////////////////////////////////////////////////////////////////////////// */
bool Convert_3cell_2_Nef3( LCC& alcc, Dart_handle d, Nef_polyhedron& N )
{
    try
    {
        Polyhedron_EK P;
        if ( Convert_3cell_2_exact_polyhedron3(alcc, d, P) )
        {
            N += Nef_polyhedron( P );
            return N.is_simple();
        }
        else
            return false;
    }
    catch (...)
    {
        std::cout << "\nCouldn't convert to Nef!" << std::endl;
        return false;
    }
}


/* ////////////////////////////////////////////////////////////////////////////// */
/// Creates a 3-cell out of an exact kernel Polyhedron_3 input
/* ////////////////////////////////////////////////////////////////////////////// */
Dart_handle Build_lcc_from_exact_polyhedron_3(LCC& alcc, Polyhedron_EK& apoly)
{
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==3 );

    typedef typename Polyhedron_EK::Halfedge_const_handle  Halfedge_handle;
    typedef typename Polyhedron_EK::Facet_const_iterator   Facet_iterator;
    typedef typename Polyhedron_EK::Halfedge_around_facet_const_circulator
            HF_circulator;

    typedef std::map < Halfedge_handle, typename Dart_handle>
            Halfedge_handle_map;
    typedef typename Halfedge_handle_map::iterator itmap_hds;
    Halfedge_handle_map TC;

    itmap_hds it;
    typename Dart_handle d = LCC::null_handle, prev = LCC::null_handle;
    typename Dart_handle firstFacet = LCC::null_handle, firstAll = LCC::null_handle;

    try
    {
        // First traversal to build the darts and link them.
        for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
        {
            HF_circulator j = i->facet_begin();
            prev = LCC::null_handle;
            do
            {
                d = alcc.create_dart();
                TC[j] = d;

                if (prev != LCC::null_handle) alcc.template link_beta<1>(prev, d);
                else firstFacet = d;
                it = TC.find(j->opposite());
                if (it != TC.end())
                    alcc.template link_beta<2>(d, it->second);
                prev = d;
            }
            while (++j != i->facet_begin());
            alcc.template link_beta<1>(prev, firstFacet);
            if (firstAll == LCC::null_handle) firstAll = firstFacet;
        }

        // Second traversal to update the geometry.
        // We run one again through the facets of the HDS.
        for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
        {
            HF_circulator j = i->facet_begin();
            do
            {
                d = TC[j]; // Get the dart associated to the Halfedge
                if (alcc.attribute<0>(d) == LCC::null_handle)
                {
                    alcc.set_vertex_attribute
                            (d, alcc.create_vertex_attribute( exa2inexa( j->opposite()->vertex()->point() )));
                }
            }
            while (++j != i->facet_begin());
        }
    }
    catch (...)
    {
        std::cout << "\nInvalid Polyhedron!" << std::endl;
        return firstAll;
    }

    return firstAll;
}


/* ////////////////////////////////////////////////////////////////////////////// */
/// Creates a 3-cell out of an inexact kernel Polyhedron_3 input
/* ////////////////////////////////////////////////////////////////////////////// */
Dart_handle Build_lcc_from_inexact_polyhedron_3(LCC& alcc, Polyhedron_IK& apoly)
{
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==3 );

    typedef typename Polyhedron_IK::Halfedge_const_handle  Halfedge_handle;
    typedef typename Polyhedron_IK::Facet_const_iterator   Facet_iterator;
    typedef typename Polyhedron_IK::Halfedge_around_facet_const_circulator
            HF_circulator;

    typedef std::map < Halfedge_handle, typename Dart_handle>
            Halfedge_handle_map;
    typedef typename Halfedge_handle_map::iterator itmap_hds;
    Halfedge_handle_map TC;

    itmap_hds it;
    typename Dart_handle d = LCC::null_handle, prev = LCC::null_handle;
    typename Dart_handle firstFacet = LCC::null_handle, firstAll = LCC::null_handle;

    // First traversal to build the darts and link them.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
        HF_circulator j = i->facet_begin();
        prev = LCC::null_handle;
        do
        {
            d = alcc.create_dart();
            TC[j] = d;

            if (prev != LCC::null_handle) alcc.template link_beta<1>(prev, d);
            else firstFacet = d;
            it = TC.find(j->opposite());
            if (it != TC.end())
                alcc.template link_beta<2>(d, it->second);
            prev = d;
        }
        while (++j != i->facet_begin());
        alcc.template link_beta<1>(prev, firstFacet);
        if (firstAll == LCC::null_handle) firstAll = firstFacet;
    }

    // Second traversal to update the geometry.
    // We run one again through the facets of the HDS.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
        HF_circulator j = i->facet_begin();
        do
        {
            d = TC[j]; // Get the dart associated to the Halfedge
            if (alcc.attribute<0>(d) == LCC::null_handle)
            {
                alcc.set_vertex_attribute
                        (d, alcc.create_vertex_attribute( j->opposite()->vertex()->point() ));
            }
        }
        while (++j != i->facet_begin());
    }

    return firstAll;
}

/* ////////////////////////////////////////////////////////////////////////////// */
/// Creates a (set of) 3-cell(s) out of a Nef_Polyhedron_3
/* ////////////////////////////////////////////////////////////////////////////// */
Dart_handle Build_lcc_from_nef_3(LCC& alcc, Nef_polyhedron& N)
{
    std::cout << "\nBuilding an LCC from a NEF: " << std::endl;
    Dart_handle d = LCC::null_handle;
    uint total = 0, invalid = 0;

    if(!N.is_simple())
        std::cout << "(The NEF is not simple though ...)" << std::endl;
    else
        std::cout << "(The NEF is simple ...)" << std::endl;

    // the first volume is the outer volume, which is
    // ignored in the decomposition
    Nef_polyhedron::Volume_const_iterator ci = ++N.volumes_begin();
    for( ; ci != N.volumes_end(); ++ci)
    {
        if(ci->mark())
        {
            Polyhedron_EK P;
            N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);

            // A valid polyhedron may not be a valid manifold...
            if (P.is_valid())
            {
//                    std::string output = "../../data/Outputs/From_Nef" + std::to_string(total) + ".off";
//                    std::ofstream out (output);
//                    write_off(out, P);

//                    std::cout << P << "\n\n" << std::endl;


//                    // The generated convex volume are often pointing at the opposite of their original volume
//                    Polygon_mesh_processing::reverse_face_orientations(P);
                d = Build_lcc_from_exact_polyhedron_3(alcc, P);
            }
            else
            {
                invalid++;
                std::cout << "\tInvalid volume detected from the Nef (Polyhedron) ... xD" << std::endl;
            }
        }
        total++;
    }

    std::cout << "\nTotal 3-Cell(s) detected: " << total << " (with " << invalid << " found invalid)." << std::endl;
    return d;
}


/* ////////////////////////////////////////////////////////////////////////////// */
/// Gets the convex decomposition of a Nef_Polyhedron_3 as a set of Polyhedron
/* ////////////////////////////////////////////////////////////////////////////// */
void Get_convex_decomposition(Nef_polyhedron& N, std::vector<Polyhedron_EK>& convex_P )
{
    CGAL::convex_decomposition_3(N);

    // the first volume is the outer volume, which is
    // ignored in the decomposition
    Nef_polyhedron::Volume_const_iterator ci = ++N.volumes_begin();
    for( ; ci != N.volumes_end(); ++ci)
    {
        if(ci->mark())
        {
            Polyhedron_EK P;
            N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);

            // A valid polyhedron may not be a valid manifold...
            if (P.is_valid())
            {
                // The generated convex volume are often pointing at the opposite of their original volume
                Polygon_mesh_processing::reverse_face_orientations(P);

                convex_P.push_back(P);
            }
            else
                std::cout << "\nInvalid polyhedron from a convex part... xD" << std::endl;
        }
    }
    std::cout << "\nDecomposition into " << convex_P.size() << " convex parts " << std::endl;
}

/* ////////////////////////////////////////////////////////////////////////////////////////// */
/// Performs the corefinement method to retrieve the sew3 relations between 3-cells
/* ////////////////////////////////////////////////////////////////////////////////////////// */
void Perform_corefinement_reconstruction (LCC &alcc, vec_dart &vol_subset, vec_dart &corefined_vol_subset, vec_dart &vol_to_remove)
{
    std::cout << "\nPerforming Corefinement ..." << std::endl;

    std::vector<Bbox_3> vol_bboxes;
    std::vector<Polyhedron_EK> vol_poly;
    // collect all the volumes
    // If there are selected volumes
    if (vol_subset.size() > 0)
    {
        for(uint i=0; i<vol_subset.size(); i++)
        {
            Polyhedron_EK poly;
            Dart_handle vol = vol_subset[i];
            // Get the bboxes to check intersecting volumes only
            Bbox_3 bbox_vol = LCCtools::Get_Bbox_vol(alcc, vol);
            if ( !Bbox_is_degenerate( bbox_vol )
                 && Convert_3cell_2_exact_polyhedron3(alcc, vol, poly) )
            {
                vol_poly.push_back(poly);

                // slightly increase the bboxes
                Rescale_Bbox( bbox_vol, coplanarity );
                vol_bboxes.push_back(bbox_vol);

                vol_to_remove.push_back( vol );
            }
            else
            {
                std::cout << "\n\tDegenerate volume or couldn't convert to Polyhedron..." << std::endl;
            }
        }
    }

    // If no subset is specified, process the entire LCC
    else
    {
        typename LCC::Base::One_dart_per_cell_range<3>::iterator vol = alcc.one_dart_per_cell<3>().begin(),
                last_vol = alcc.one_dart_per_cell<3>().end();
        for(; vol!=last_vol; vol++)
        {
            Polyhedron_EK poly;
            // Get the bboxes to check intersecting volumes only
            Bbox_3 bbox_vol = LCCtools::Get_Bbox_vol(alcc, vol);
            if ( !Bbox_is_degenerate( bbox_vol )
                 && Convert_3cell_2_exact_polyhedron3(alcc, vol, poly) )
            {
                vol_poly.push_back(poly);

                // slightly increase the bboxes
                Rescale_Bbox( bbox_vol, coplanarity );
                vol_bboxes.push_back(bbox_vol);

                vol_to_remove.push_back( vol );
            }
            else
            {
                std::cout << "\n\tDegenerate volume or couldn't convert to Polyhedron..." << std::endl;
            }
        }
    }

    assert(vol_poly.size() == vol_bboxes.size());
    assert(vol_bboxes.size() == vol_to_remove.size());

    for(uint i=0; i<vol_poly.size(); i++)
    {
        for(uint j=i+1; j<vol_poly.size(); j++)
        {
            CGAL::Polygon_mesh_processing::corefine(vol_poly[i], vol_poly[j]);
        }

        Dart_handle crf_vol = Build_lcc_from_exact_polyhedron_3(alcc, vol_poly[i]);
        corefined_vol_subset.push_back(crf_vol);
    }

}


/* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
///                                        End of Convenient Functions
/* /////////////////////////////////////////////////////////////////////////////////////////////////////// */






    void Show_face(LCC& alcc, Dart_handle d)
    {
        int c = 0;
        typename LCC::Base::Dart_of_orbit_range<1>::const_iterator itr( alcc, d );
        for (; itr.cont(); ++itr )
        {
            std::cout << alcc.point(itr) << std::endl;
            c++;
        }
        std::cout << alcc.point(d) << std::endl;
        std::cout << "number of points: " << c << std::endl;
    }

    void Show_volume(LCC& alcc, Dart_handle d)
    {
        int c = 0;
        typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
                itr = alcc.one_dart_per_incident_cell<2,3>(d).begin(),
                itrend = alcc.one_dart_per_incident_cell<2,3>(d).end();
        for (; itr != itrend; itr++)
            Show_face(alcc,itr);
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Applies a size scaling factor to all the points of the LCC
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Rescale_LCC_size (LCC& alcc, double scaling)
    {
        // Go trhough all the 0-cells
        typename Dart_range::iterator it(alcc.darts().begin()),
                                      itend(alcc.darts().end());
        for(; it!=itend; it++)
        {
//                Point pt = alcc.point( it );
            alcc.point( it ) = Scale_3D_point( alcc.point( it ), scaling );
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Copies the attributes from 3-cell(d) to an attribute structure
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Copy_Attributes_from_3cell(LCC& alcc, Dart_handle& d, Vol_Attribs& attrib)
//    {
//        assert( alcc.attribute<3>(d) != LCC::null_handle );

//        // Collect attributes from 3-cell(d)
//        uint m_uID = alcc.info<3>(d).id(); // unique ID
//        uint m_vol_label = alcc.info<3>(d).vol_label(); // label common to all elements of a same volume
//        std::string m_label = alcc.info<3>(d).label(); // label specific to a semantic class of 3-cell
//        std::vector<std::string> m_IfcInfo = alcc.info<3>(d).IfcInfo(); // Ifc Global Unique ID
//        bool m_select = alcc.info<3>(d).is_selected(); // Whether a 3-cell is selected or not
//        bool m_init = alcc.info<3>(d).is_init(); // Whether a 3-cell is fully initialized or not
//        CGAL::Color m_color = alcc.info<3>(d).color(); // Color of the 3-cell

//        // copy to the structure
//        attrib.m_uID = m_uID;
//        attrib.m_vol_label = m_vol_label;
//        attrib.m_label = m_label;
//        attrib.m_IfcInfo = m_IfcInfo;
//        attrib.m_select = m_select;
//        attrib.m_init = m_init;
//        attrib.m_color = m_color;
//    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Copies the attributes from an attribute structure to 3-cell(d)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Paste_Attributes_to_3cell(LCC& alcc, Dart_handle& d, Vol_Attribs& attrib)
//    {
//        if ( alcc.attribute<3>(d)==LCC::null_handle )
//            alcc.set_attribute<3>(d, alcc.create_attribute<3>());

//        // copy to 3-cell(d)
//        alcc.info<3>(d).set_id( attrib.m_uID );
//        alcc.info<3>(d).set_vol_label( attrib.m_vol_label );
//        alcc.info<3>(d).set_label( attrib.m_label );
//        alcc.info<3>(d).set_selection( attrib.m_select );
//        alcc.info<3>(d).set_init( attrib.m_init );
//        alcc.info<3>(d).set_color(attrib.m_color.r(),
//                                  attrib.m_color.g(),
//                                  attrib.m_color.b());
//    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Transfers the attributes of 3-cell(d1) to 3-cell(d2)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Transfer_Attributes_of_3cell(LCC& alcc, Dart_handle& d1, Dart_handle& d2)
//    {
//        Vol_Attribs attrib;
//        Copy_Attributes_from_3cell(alcc, d1, attrib);
//        Paste_Attributes_to_3cell(alcc, d2, attrib);
//    }

#ifdef IFCPP_ON

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Converts a 3D point from OSG to LCC, given a mesh vertex structure of OSG.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    LCC::Point Get_Osg_3dpt_from_Mesh_Vertex (carve::mesh::Vertex<3>& vertex, Point global_shift_pt)
    {
        // coordinates of vertex:
        carve::geom::vector<3>& vertex_coords = vertex.v;
        double x = (vertex_coords.x * scale_f) - global_shift_pt.x(),
               y = (vertex_coords.y * scale_f) - global_shift_pt.y(),
               z = (vertex_coords.z * scale_f) - global_shift_pt.z();

//            std::wcout << "\t# point: [" << x << ", " << y << ", " << z << "]" << std::endl;

        return LCC::Point( x, y, z );
    }



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the certices of a set of meshsets into a vector of 3D points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    std::vector<Point> Get_Vertices_from_Mesh (carve::mesh::Mesh<3>* mesh)
    {
        std::vector<Point> res;
        Point pt, pt_;
        int step = 10;

        std::vector<carve::mesh::Face<3>* >& vec_faces = mesh->faces;
        for( size_t i_face = 0; i_face < vec_faces.size(); ++i_face )
        {
            carve::mesh::Face<3>* face = vec_faces[i_face];
            // iterate through edges:
            carve::mesh::Edge<3>* edge = face->edge;
            do
            {
                pt = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v1() );
                res.push_back(pt);

                /// Adding additional points between every segment extremity to optimize Alpha-shape simplification
//                    pt_ = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v2() );
//                    if ( Vector(pt, pt_).squared_length() > (0.1*0.1) )
//                    {
//                        int cx=1, cy=1, cz=1;
//                        if ( pt.x() > pt_.x() )
//                            cx = -1;
//                        if ( pt.y() > pt_.y() )
//                            cy = -1;
//                        if ( pt.z() > pt_.z() )
//                            cz = -1;

//                        double dx = fabs(pt_.x() - pt.x()) / step,
//                                dy = fabs(pt_.y() - pt.y()) / step,
//                                dz = fabs(pt_.z() - pt.z()) / step;

//                        for(uint i=1; i<=step; i++)
//                        {
//                            res.push_back( Point( pt.x() + (cx*dx*i),
//                                                  pt.y() + (cy*dy*i),
//                                                  pt.z() + (cz*dz*i)) );
//                        }
//                    }

                // increment edge
                edge = edge->next;
            }
            while( edge != face->edge );
        }

        std::cout << "\n\tNumber of points got from the mesh: " << res.size() << std::endl;
        return res;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the vertices of a single meshset into a vector of 3D points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    std::vector<Point> Get_Vertices_from_Single_Meshset (shared_ptr<carve::mesh::MeshSet<3> >& meshset)
    {
        std::vector<Point> res;
        Point pt1, pt2;
        int step = 10;

        std::vector<carve::mesh::Mesh<3>* >& vec_meshes = meshset->meshes;
        for( size_t i_mesh = 0; i_mesh < vec_meshes.size(); ++i_mesh )
        {
            carve::mesh::Mesh<3>* mesh = vec_meshes[i_mesh];
            std::vector<carve::mesh::Face<3>* >& vec_faces = mesh->faces;
            for( size_t i_face = 0; i_face < vec_faces.size(); ++i_face )
            {
                carve::mesh::Face<3>* face = vec_faces[i_face];
                // iterate through edges:
                carve::mesh::Edge<3>* edge = face->edge;
                do
                {
                    pt1 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v1() );

                    res.push_back(pt1);
//                            res.push_back(pt2);


                    /// Adding additional points between every segment extremity to optimize Alpha-shape simplification
//                        pt2 = Get_Osg_3dpt_from_Mesh_Vertex( *edge->v2() );
//                        if ( Vector(pt1, pt2).squared_length() > (0.1*0.1) )
//                        {
//                            int cx=1, cy=1, cz=1;
//                            if ( pt1.x() > pt2.x() )
//                                cx = -1;
//                            if ( pt1.y() > pt2.y() )
//                                cy = -1;
//                            if ( pt1.z() > pt2.z() )
//                                cz = -1;

//                            double dx = fabs(pt2.x() - pt1.x()) / step,
//                                    dy = fabs(pt2.y() - pt1.y()) / step,
//                                    dz = fabs(pt2.z() - pt1.z()) / step;

//                            for(uint i=1; i<=step; i++)
//                            {
//                                res.push_back( Point( pt1.x() + (cx*dx*i),
//                                                      pt1.y() + (cy*dy*i),
//                                                      pt1.z() + (cz*dz*i)) );
//                            }
//                        }

                    // increment edge
                    edge = edge->next;
                }
                while( edge != face->edge );
            }
        }

        return res;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the vertices of a set of meshsets into a vector of 3D points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    std::vector<Point> Get_Vertices_from_Meshsets (std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& meshsets)
    {
        std::vector<Point> res, res_;
        Point pt1, pt2;

        for( size_t i_meshset = 0; i_meshset < meshsets.size(); ++i_meshset )
        {
            res_ = Get_Vertices_from_Single_Meshset(meshsets[i_meshset]);
            res.insert(res.end(), res_.begin(), res_.end());
            res_.clear();
        }

        return res;
    }
#endif

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Computes the 2D shape simplification of a set of points from a group of meshsets (representing one IfcProduct)
    /// Returns the 3D points of the simplified shape (res.first) and the highest Z coordinate (res.second).
    /// mode 0: convex hull
    /// mode 1: OBB
    /// mode 2: Alpha-shape
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Pointset_Simplification (std::vector<Point>& all_points, std::pair<double,
                                                    std::vector<Point> >& res, int mode)
    {
        Plane Pl(Point(0.0, 0.0, 0.0),
                 V_height);
        std::map<double, bool> highest_z;
        std::vector<Point2d> input_ch, output_ch, mini_box, unique_pts_2d;
        std::vector<Point> unique_pts_input, unique_pts;

//            std::wcout << "\nInput for Convex Hull:" << std::endl;


        //Filter out duplicated point
        std::map<std::string, Point> all_points_no_dup_map;
        for(uint i=0; i<all_points.size(); i++)
        {
            all_points_no_dup_map[ get_unique_key_pt( all_points[i] ) ] = all_points[i];
        }

        for(auto it = all_points_no_dup_map.begin(); it!=all_points_no_dup_map.end(); it++)
            unique_pts_input.push_back( it->second );

        std::cout << "\nNumber of points before cleaning: " << all_points.size() <<std::endl;
        std::cout << "Number of points after cleaning: " << unique_pts_input.size() <<std::endl;


//            std::cout << "\nPoints used to compute the alpha shape: " << all_points.size() <<std::endl;
        // Collect all the points to compute Convex Hull
        for(uint i=0; i<unique_pts_input.size(); i++)
        {
//                std::cout << all_points[i] << std::endl;
            input_ch.push_back( Pl.to_2d( unique_pts_input[i] ) );
            highest_z[ to_double(unique_pts_input[i].z()) ] = true;
//                std::cout<< unique_pts_input[i] << std::endl;
        }

        /// CONVEX HULL + OBB APPROACH
        if (mode == 0 || mode == 1)
        {
            //Computes the 2D convex hull
            convex_hull_2( input_ch.begin(), input_ch.end(), std::back_inserter(output_ch) );

            std::cout << "\nFound Convex Hull! (original pointset: " << unique_pts_input.size()
                      << " and convex_hull: " << output_ch.size()<< ")" << std::endl;

            // Get rid of redundant points of the CH
            std::map<std::string, bool> unique_pts_map;
            for( uint i=0; i<output_ch.size(); i++ )
            {
                // TODO: find a way to make this part generic, in case V_height is not Z...
                Point cur_pt = Point( to_double(output_ch[i].x()),
                                      to_double(output_ch[i].y()),
                                      highest_z.begin()->first   );
                //                Point cur_pt = Pl.to_3d( output_ch[i] );
                //                cur_pt = Point( cur_pt.x(), cur_pt.y(), highest_z.begin()->first );
                std::string unique_key = get_unique_key_pt( cur_pt );
                if ( !unique_pts_map[ unique_key ] )
                {
                    unique_pts_map[ unique_key ] = true;
                    //                    unique_pts.push_back( cur_pt );
                    unique_pts_2d.push_back( Point2d( cur_pt.x(), cur_pt.y() ) );
                    //                    std::cout << cur_pt.x()<< ", " << cur_pt.y() << std::endl;
                }
            }

            if (mode == 0)
                mini_box = unique_pts_2d;
            else
                Get_smallest_bounding_rectangle(unique_pts_2d, mini_box, eps2);

            std::cout << "\nFinal simplified (2D) shape: " << std::endl;
            for( uint i=0; i<mini_box.size(); i++ )
            {
                std::cout << mini_box[i] << std::endl;
                Point cur_pt = Point( to_double(mini_box[i].x()),
                                      to_double(mini_box[i].y()),
                                      highest_z.begin()->first   );
                unique_pts.push_back( cur_pt );
            }

            std::cout << "\nFinal simplified (3D) shape: " << std::endl;
            for( uint i=0; i<unique_pts.size(); i++ )
                std::cout << unique_pts[i] << std::endl;

        }
        /// ALPHA-SHAPES (CONCAVE-HULL) APPRAOCH
//        else if (mode == 2)
//        {
//            Alpha_shape2d A(input_ch.begin(), input_ch.end(),
//                            FT(0.4),
//                            Alpha_shape2d::GENERAL);
//            std::vector<Segment2d> segments;
//            alpha_edges(A, std::back_inserter(segments));
//            std::cout << "\nFound An ALPHA SHAPE (original pointset: " << all_points.size()
//                      << " and Alpha-Shape: " << segments.size()<< " segments)" << std::endl;
//            std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1)<<std::endl;

//            std::map<std::string, std::vector<Point2d> > sorted_seg;
//            for(uint i=0; i<segments.size(); i++)
//            {
//                std::string uk2d = get_unique_key_pt_2d(segments[i].source(), 0.0000001);
//                // More than one target stored at this key means as many different segments
//                sorted_seg[ uk2d ].push_back( segments[i].target() );

////                    // If that vertex is already stored
////                    if ( sorted_seg.find( uk2d ) != sorted_seg.end() )
////                    {
////                        // Try to store the other extremity instead
////                        std::string uk2d_ = get_unique_key_pt_2d(segments[i].target());
////                        if ( sorted_seg.find( uk2d_ ) != sorted_seg.end() )
////                            sorted_seg[ uk2d_ ] = segments[i].source();
////                        else
////                            std::cout << "Well... This should not happen actually!" <<std::endl;
////                    }
////                    else
////                        sorted_seg[ uk2d ] = segments[i].target();
//            }

//            std::cout << "Sorted segments: " << sorted_seg.size()  <<std::endl;

//            Point2d ref_pt = segments[0].source(), cur_pt = segments[0].target(), prev_pt = ref_pt;
//            std::vector<Point2d> boundary_seq;
//            boundary_seq.push_back( ref_pt );
//            uint count = 0;
//            while ( (cur_pt != ref_pt) && count <= segments.size() )
//            {
//                boundary_seq.push_back( cur_pt );

//                std::string uk2d = get_unique_key_pt_2d(cur_pt, 0.0000001);
//                std::cout << "Number of vertices at " << uk2d
//                          << " (" << cur_pt << "): "
//                          << sorted_seg[ uk2d ].size() <<std::endl;
//                assert (sorted_seg[ uk2d ].size() <= 2);

//                bool got_one = false;
//                for (uint i=0; i<sorted_seg[ uk2d ].size(); i++)
//                {
//                    if (sorted_seg[ uk2d ][i] != prev_pt )
//                    {
//                        prev_pt = cur_pt;
//                        cur_pt = sorted_seg[ uk2d ][i];
//                        got_one = true;
//                        break;
//                    }
//                }

//                if (!got_one)
//                    std::cout << "This should definitely not happen actually!" <<std::endl;

//                count++;
//                if (count == segments.size() )
//                    std::cout << "The loop did not end up properly..." <<std::endl;
//            }

//            std::cout << "Final Sorted segments: " << boundary_seq.size()  <<std::endl;

//            // Get rid of redundant points of the AS
//            std::map<std::string, bool> unique_pts_map;
//            for( uint i=0; i<boundary_seq.size(); i++ )
//            {
//                // TODO: find a way to make this part generic, in case V_height is not Z...
//                Point cur_pt = Point( to_double(boundary_seq[i].x()),
//                                      to_double(boundary_seq[i].y()),
//                                      highest_z.begin()->first   );
//                //                Point cur_pt = Pl.to_3d( output_ch[i] );
//                //                cur_pt = Point( cur_pt.x(), cur_pt.y(), highest_z.begin()->first );
//                std::string unique_key = get_unique_key_pt( cur_pt );
//                if ( !unique_pts_map[ unique_key ] )
//                {
//                    unique_pts_map[ unique_key ] = true;
//                    //                    unique_pts.push_back( cur_pt );
//                    unique_pts_2d.push_back( Point2d( cur_pt.x(), cur_pt.y() ) );
//                    //                    std::cout << cur_pt.x()<< ", " << cur_pt.y() << std::endl;
//                }
//            }


//            std::cout << "@@@@@@@@@ New boundary for the Alpha-Shape" <<std::endl;
//            Polygon_2<IK> my_poly( unique_pts_2d.begin(), unique_pts_2d.end() );

//            if (my_poly.is_simple() && my_poly.is_counterclockwise_oriented() )
//            {
//                std::cout << "The Poly is CCW!" <<std::endl;
//            }
//            else
//                std::cout << "The Poly is either not simple or CW!" <<std::endl;

//            for( uint i=0; i<unique_pts_2d.size(); i++ )
//            {
//                std::cout << unique_pts_2d[i] <<std::endl;

//                Point cur_pt = Point( to_double(unique_pts_2d[i].x()),
//                                      to_double(unique_pts_2d[i].y()),
//                                      highest_z.begin()->first   );
//                unique_pts.push_back( cur_pt );
//            }
//        }



        // Keep the highest z and create the 3D projection of the 2D convex hull
        // by affecting the lowest z coordinate to all the resulting 2D points.
        res.first = highest_z.rbegin()->first;
        res.second = unique_pts;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Z-extrusion of a 2-cell containing the dart d
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Z_extrude_2cell (LCC& alcc, Dart_handle d, double& Z)
    {
        std::map<std::string, Dart_handle, LCCtools::cmp_string> existing_darts;
        std::map<std::string, Dart_handle, LCCtools::cmp_string>::iterator it;
        std::vector<Point> upper_face;
        int total_sew2 = 0;

        typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, d);
        for (; itr.cont(); itr++)
        {
            Point pt1 = alcc.point(itr),
                  pt2 = alcc.point(alcc.beta<1>(itr));
            Point pt1_up ( pt1.x(), pt1.y(), Z ),
                  pt2_up ( pt2.x(), pt2.y(), Z );

            // Add the dart of the base face in the map
            existing_darts[ std::string( get_unique_key_pt(pt1) +
                                         get_unique_key_pt(pt2) ) ] = itr;

            // Create the latteral face corresponding to the current edge
            std::vector<Point> latteral_face;
            latteral_face.push_back( pt2 );
            latteral_face.push_back( pt1 );
            latteral_face.push_back( pt1_up );
            latteral_face.push_back( pt2_up );

            Dart_handle l_face = Insert_new_2cell(alcc, latteral_face);

            std::wcout << std::endl;
            // Insert the darts of the face in the map and check for possible beta2 links
            typename LCC::Base::Dart_of_orbit_range<1>::iterator itr_(alcc, l_face);
            for (; itr_.cont(); itr_++)
            {
                Point pt1_ = alcc.point(itr_),
                      pt2_ = alcc.point(alcc.beta<1>(itr_));

                existing_darts[ std::string( get_unique_key_pt(pt1_) +
                                             get_unique_key_pt(pt2_) ) ] = itr_;

                it = existing_darts.find( std::string( get_unique_key_pt(pt2_) +
                                                       get_unique_key_pt(pt1_) ) );
                if ( it != existing_darts.end()
                     && alcc.template is_sewable<2>(itr_, it->second))
                {
                    std::wcout << "sew2!" << std::endl;
                    alcc.template link_beta<2>(itr_, it->second);
                    std::cout << "\t" << pt1_ << "\t" << pt2_ << std::endl;
                    std::cout << "\t" << alcc.point(alcc.beta<1>(it->second)) << "\t" << alcc.point(it->second) << std::endl;
                    total_sew2++;
                }
            }

            // collect points for the upper face creation
            upper_face.push_back( pt1_up );
        }

        // Create the upper face (opposite direction of the base face) and link it to the existing ones
        Dart_handle u_face = Insert_new_2cell(alcc, upper_face, false);
        std::wcout << "size of upper face: " << upper_face.size() << std::endl;
        typename LCC::Base::Dart_of_orbit_range<1>::iterator itr_(alcc, u_face);
        for (; itr_.cont(); itr_++)
        {
            Point pt1_ = alcc.point(itr_),
                  pt2_ = alcc.point(alcc.beta<1>(itr_));

            // The opposite darts should exist already!
            assert( existing_darts.find( std::string( get_unique_key_pt(pt2_) + get_unique_key_pt(pt1_) ) )
                    != existing_darts.end() );

            alcc.template link_beta<2>(itr_, existing_darts[ std::string( get_unique_key_pt(pt2_) +
                                                                          get_unique_key_pt(pt1_) ) ]);
            total_sew2++;
        }

        if ( total_sew2 == (upper_face.size() * 3) )
            std::wcout << "GOOD sew2!" << total_sew2  << std::endl;
        else
            std::wcout << "BAD sew2!!!!" << total_sew2  << std::endl;

    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Z-extrusion of a the rectangle described by the vector of points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    Dart_handle Z_extrude_rectangle (LCC& alcc, std::pair<double, std::vector<Point> >& rect_and_h)
    {
        assert(rect_and_h.second.size() == 4);

        double height = rect_and_h.first;
        std::vector<Point> p_down = rect_and_h.second;
        Point p_up0 ( p_down[3].x(), p_down[3].y(), height ),
              p_up1 ( p_down[0].x(), p_down[0].y(), height ),
              p_up2 ( p_down[1].x(), p_down[1].y(), height ),
              p_up3 ( p_down[2].x(), p_down[2].y(), height );

        Polygon_2<IK> pgn;
        pgn.push_back( Point2d (p_down[3].x(), p_down[3].y()) );
        pgn.push_back( Point2d (p_down[0].x(), p_down[0].y()) );
        pgn.push_back( Point2d (p_down[1].x(), p_down[1].y()) );
        pgn.push_back( Point2d (p_down[2].x(), p_down[2].y()) );

        if (height != 0.0
                && pgn.area() > (eps2*eps2) )
            return alcc.make_hexahedron(p_down[0], p_down[1], p_down[2], p_down[3],
                                        p_up0, p_up1, p_up2, p_up3);
        else
        {
            std::wcout << "The Height is null (flat volume): no OBB generated!" << std::endl;
            return LCC::null_handle;
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Z-extrusion of a the rectangle described by the vector of points
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    Dart_handle Z_extrude_polygon (LCC& alcc, std::vector<Point>& poly_and_h, double h)
//    {
//        assert(poly_and_h.size() > 3);

//        vec_dart new_f;
//        vec_pt3d new_f_vec, upper;

//        // Creating the base polygon in the LCC
//        new_f.push_back( Insert_new_2cell(alcc, poly_and_h) );

//        // Creating the side polygons
//        poly_and_h.push_back( poly_and_h[0] );
//        for( uint i=0; i<poly_and_h.size()-1; i++ )
//        {
//            Point cur = poly_and_h[i],
//                  next = poly_and_h[i+1];

//            new_f_vec.push_back( next );
//            new_f_vec.push_back( cur );
//            new_f_vec.push_back( Point(cur.x(), cur.y(), h) );
//            new_f_vec.push_back( Point(next.x(), next.y(), h) );

//            new_f.push_back( Insert_new_2cell(alcc, new_f_vec) );
//            new_f_vec.clear();
//        }
//        poly_and_h.pop_back();

//        // Creating the upper polygon
//        for( uint i=0; i<poly_and_h.size(); i++ )
//        {
//            upper.push_back( Point(poly_and_h[i].x(),
//                                   poly_and_h[i].y(),
//                                   h) );
//        }
//        new_f.push_back( Insert_new_2cell(alcc, upper, false) );

//        if ( Perform_Simple_Volume_Reconstruction_from_Polygon_Soup(alcc, new_f) )
//            return new_f[0];
//        else
//            return LCC::null_handle;
//    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Z-extrusion of a 2-cell containing the dart d
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Z_extrude_upper_lower_box_faces (LCC& alcc, Dart_handle d, double Z, bool upper, bool lower, bool dir)
//    {
//        std::vector <vec_dart > parallel_face = Collect_parallel_faces_in_volume(alcc, d);
//        assert(parallel_face.size() == 3);
//        Vector V_height_opp = V_height * (-1);

//        for(uint i=0; i<parallel_face.size(); i++)
//        {
//            assert(parallel_face[i].size() == 2);
//            double angle = compute_3d_angle( LCCtools::poly_normal( alcc, parallel_face[i][0] ), V_height );

//            // if it is the lower face that is found first
//            if (  (angle < (0.0 + eps_angle) || angle > (360.0 - eps_angle)) )
//            {
//                // and it is the one intended to be extruded
//                if (lower)
//                {
//                    if (dir)
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][0], V_height, Z );
//                    else
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][0], V_height_opp, Z );
//                }

//                // if the upper face is, as well.
//                if ( upper )
//                {
//                    if (dir)
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][1], V_height, Z );
//                    else
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][1], V_height_opp, Z );
//                }
//            }

//            // if it is the upper face that is found first, same!
//            else if ( (angle < (180.0 + eps_angle) && angle > (180.0 - eps_angle) ) )
//            {
//                // and it is the one intended to be extruded
//                if (upper)
//                {
//                    if (dir)
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][0], V_height, Z );
//                    else
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][0], V_height_opp, Z );
//                }

//                // if the upper face is, as well.
//                if ( lower )
//                {
//                    if (dir)
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][1], V_height, Z );
//                    else
//                        Extrude_face_in_vector_direction( alcc, parallel_face[i][1], V_height_opp, Z );
//                }
//            }
//        }
//    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Oriented Bounding Box of a point set
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    Dart_handle OBB_Pointset_Simplification (LCC& alcc, std::vector<Point>& pointset)
    {
        std::pair<double, std::vector<Point> > CH3d_and_height;
        Pointset_Simplification( pointset, CH3d_and_height );
        return Z_extrude_rectangle(alcc, CH3d_and_height);
//            return Z_extrude_polygon(alcc, CH3d_and_height.second, CH3d_and_height.first);
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Oriented Bounding Box of a point set
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    Dart_handle Convex_hull_Pointset_Simplification (LCC& alcc, std::vector<Point>& pointset)
//    {
//        std::pair<double, std::vector<Point> > CH3d_and_height;
//        Pointset_Simplification( pointset, CH3d_and_height, 0 );
//        return Z_extrude_polygon(alcc, CH3d_and_height.second, CH3d_and_height.first);
//    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Creates the 3-cell corresponding to the Oriented Bounding Box of a point set
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    Dart_handle Concave_hull_Pointset_Simplification (LCC& alcc, std::vector<Point>& pointset)
//    {
//        std::pair<double, std::vector<Point> > CH3d_and_height;
//        Pointset_Simplification( pointset, CH3d_and_height, 2 );
//        return Z_extrude_polygon(alcc, CH3d_and_height.second, CH3d_and_height.first);
//    }

    /* ////////////////////////////////////////////////////////// */
    /// Compute the AABB of a 3D Point set
    /* ////////////////////////////////////////////////////////// */
    Bbox_3 Get_Bbox_Vector_Point (std::vector<Point>& points)
    {
        assert (points.size()>0);

        Bbox_3 bb = points[0].bbox();

        for (uint i=1; i<points.size(); i++)
            bb = bb + points[i].bbox();

        return bb;
    }


#ifdef IFCPP_ON
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the IfcBuildingStorey of the IFC model.
    /// Returns false if there is no IfcBuildingStorey collected.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Get_Building_Storeys (shared_ptr<BuildingModel>& ifc_model, std::vector< shared_ptr<IfcBuildingStorey> >& building_storeys)
    {
        const std::map<int, shared_ptr<BuildingEntity> >& map_ifc_entities = ifc_model->getMapIfcEntities();
        std::cout << "Number of entities found: " << map_ifc_entities.size() << std::endl;

        uint nb_spaces = 0;

        for( auto it = map_ifc_entities.begin(); it != map_ifc_entities.end(); ++it )
        {
            const shared_ptr<BuildingEntity>& ifcpp_entity = it->second;
            shared_ptr<IfcProduct> ifc_product = dynamic_pointer_cast<IfcProduct>( ifcpp_entity );

            if( ifc_product )
            {
                if ( strcmp(ifc_product->className(), "IfcBuildingStorey") == 0 )
                {
                    shared_ptr<IfcBuildingStorey> storey = dynamic_pointer_cast<IfcBuildingStorey>( ifc_product );
                    std::wcout << "\t\nGot one storey!" <<  std::endl;
                    building_storeys.push_back( storey );
                }
                else if ( strcmp(ifc_product->className(), "IfcSpace") == 0 )
                    nb_spaces++;
            }
        }

        if ( building_storeys.size() > 0 )
        {
            std::wcout << "\nTotal number of storeys found: " << building_storeys.size() <<  std::endl;
            std::wcout << "\nTotal number of spaces found: " << nb_spaces <<  std::endl;
            return true;
        }

        return false;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the IfcSpaces related to a given IfcBuildingStorey.
    /// Returns false if there is no IfcSpace collected.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Get_IfcSpaces_of_Storey (shared_ptr<IfcBuildingStorey>& storey, std::vector< shared_ptr<IfcSpace> >& storey_spaces)
    {
        std::vector<shared_ptr<IfcObjectDefinition> >::iterator it_object_def;

        if( storey->m_IsDecomposedBy_inverse.size() > 0 )
        {
            std::vector<weak_ptr<IfcRelAggregates> >& vec_IsDecomposedBy = storey->m_IsDecomposedBy_inverse;
            std::vector<weak_ptr<IfcRelAggregates> >::iterator it_;
            for( it_=vec_IsDecomposedBy.begin(); it_!=vec_IsDecomposedBy.end(); ++it_ )
            {
                shared_ptr<IfcRelAggregates> rel_agg( *it_ );
                std::vector<shared_ptr<IfcObjectDefinition> >& vec = rel_agg->m_RelatedObjects;

                for( it_object_def=vec.begin(); it_object_def!=vec.end(); ++it_object_def )
                {
                    shared_ptr<IfcObjectDefinition> child_obj_def = (*it_object_def);
//                    std::wcout << "What's in there? ---> " << child_obj_def->className() << std::endl;

                    if ( strcmp(child_obj_def->className(), "IfcSpace") == 0 )
                    {
                        shared_ptr<IfcSpace> myspace = dynamic_pointer_cast<IfcSpace>( child_obj_def );
                        if (myspace)
                            storey_spaces.push_back( myspace );
                    }
                }

                if (storey_spaces.size() > 0)
                    std::wcout << "\nTotal number of IfcSpaces found for storey " << storey->m_Name->m_value << ": " << storey_spaces.size() <<  std::endl;
                else
                    std::wcout << "No IfcSpace found in this storey." << std::endl;

                return true;
            }
        }
        else if (storey->m_Name)
        {
            std::wcout << "\nTotal number of IfcSpaces found for storey " << storey->m_Name->m_value << ": " << storey_spaces.size() <<  std::endl;
            std::wcout << "(The storey is decomposed by nothing...??)" << std::endl;
        }

        return false;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the building elements (as IfcElement) related to a given IfcSpace
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Get_IfcElements_around_IfcSpace (shared_ptr<IfcSpace>& space, std::map< const char*,
                                          std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& elem_around_space )
    {
        std::vector<weak_ptr<IfcRelSpaceBoundary> >& around_myspace = space->m_BoundedBy_inverse;
        std::wcout << "\tCurrent IfcSpace bounded by " << around_myspace.size() << " elements!" << std::endl;

        std::vector<weak_ptr<IfcRelSpaceBoundary> >::iterator it_1;
        for( it_1=around_myspace.begin(); it_1!=around_myspace.end(); ++it_1 )
        {
            shared_ptr<IfcRelSpaceBoundary> rel_space( *it_1 );
            shared_ptr<IfcElement> rel_elem = rel_space->m_RelatedBuildingElement;

            if (rel_elem)
            {
//                std::wcout << "\t\t" << rel_elem->className() << " ---> " << rel_elem->m_entity_id << std::endl;
                elem_around_space[ rel_elem->className() ].push_back( rel_elem );
            }
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the building elements (as IfcElement) related to a given IfcDoor
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Get_IfcElements_around_IfcDoor (shared_ptr<IfcDoor>& door, std::map< const char*,
                                          std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& elemAroundDoor )
    {
        std::vector<weak_ptr<IfcRelConnectsElements> >& aroundMyDoor = door->m_ConnectedTo_inverse;
        std::wcout << "\tCurrent IfcDoor is connected to " << aroundMyDoor.size() << " elements!" << std::endl;

        std::vector<weak_ptr<IfcRelConnectsElements> >::iterator it_1;
        for( it_1=aroundMyDoor.begin(); it_1!=aroundMyDoor.end(); ++it_1 )
        {
            shared_ptr<IfcRelConnectsElements> conDoor( *it_1 );
            shared_ptr<IfcElement> rel_elem = conDoor->m_RelatedElement;

            if (rel_elem)
            {
                std::wcout << "\t\t" << rel_elem->className() << " ---> " << rel_elem->m_entity_id << std::endl;
                elemAroundDoor[ rel_elem->className() ].push_back( rel_elem );
            }
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the openings (Doors & Windows) related to a given IfcSpace
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Get_Openings_around_IfcSpace (shared_ptr<IfcSpace>& space, std::vector< shared_ptr<IfcElement> >& res,
                                          std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >& elem_around_space )
    {
        // If the elements surrounding the space are not collected yet
        if (elem_around_space.size() == 0)
            Get_IfcElements_around_IfcSpace(space, elem_around_space);

        // Collect doors if any
        for(uint i=0; i<elem_around_space["IfcDoor"].size(); i++)
            res.push_back( elem_around_space["IfcDoor"][i] );

        // Collect windows if any
        for(uint i=0; i<elem_around_space["IfcWindow"].size(); i++)
            res.push_back( elem_around_space["IfcWindow"][i] );

        // Collect IfcVirtualElement if any
        for(uint i=0; i<elem_around_space["IfcVirtualElement"].size(); i++)
            res.push_back( elem_around_space["IfcVirtualElement"][i] );
    }



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the building elements (as IfcProduct) related to a given IfcBuildingStorey
    /// (excluding the IfcSpaces)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Get_all_IfcProducts_related_to_IfcBuildingStorey (shared_ptr<IfcBuildingStorey>& storey,
                                                           std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >& elem_of_storey )
    {
        std::vector<weak_ptr<IfcRelContainedInSpatialStructure> > elem_in_storey = storey->m_ContainsElements_inverse;
        if (storey->m_Name)
            std::wcout << "\nThe storey " << storey->m_Name->m_value  << " contains " << elem_in_storey.size() << " elements:" << std::endl;

        if (elem_in_storey.size() > 0)
        {
            for(uint i=0; i<elem_in_storey.size(); i++)
            {
                shared_ptr<IfcRelContainedInSpatialStructure> one_elem_in_storey ( elem_in_storey[i] );
                if (one_elem_in_storey)
                {
                    std::wcout << "Here I have one " << one_elem_in_storey->className() << " related to "
                               << one_elem_in_storey->m_RelatedElements.size() << " element(s):" << std::endl;

                    std::vector<shared_ptr<IfcProduct> > related_prod = one_elem_in_storey->m_RelatedElements;
                    for(uint j=0; j<related_prod.size(); j++)
                    {

                        shared_ptr<IfcProduct> one_related_prod = dynamic_pointer_cast<IfcProduct>( related_prod[j] );
//                                std::wcout << "\t" << one_related_prod->className() << std::endl;
                        elem_of_storey[ one_related_prod->className() ].push_back( one_related_prod );
                    }
                }

                std::map < const char*, std::vector< shared_ptr<IfcProduct> > >::iterator it_elem_map (elem_of_storey.begin());
                for(; it_elem_map != elem_of_storey.end(); it_elem_map++)
                    std::wcout << "\t" << it_elem_map->first << " (" << it_elem_map->second.size() <<  ")" << std::endl;

            }
            std::wcout << std::endl;
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Collects all the building elements (as IfcProduct) contained in a given IfcSpace
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Get_all_elements_in_IfcSpace (shared_ptr<IfcSpace> &space,
                                       std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > &elem_of_space,
                                       std::map<const char *, bool, LCCtools::cmp_const_char> &unallowed_classes)
    {
        std::vector<weak_ptr<IfcRelContainedInSpatialStructure> > elem_in_space = space->m_ContainsElements_inverse;
        std::wcout << "\nThe space " << space->m_Name->m_value  << " contains " << elem_in_space.size() << " element(s)";
        if (elem_in_space.size() > 0)
        {
            for(uint i=0; i<elem_in_space.size(); i++)
            {
                shared_ptr<IfcRelContainedInSpatialStructure> one_elem_in_space ( elem_in_space[i] );
                if (one_elem_in_space)
                {
                    std::wcout << " " << one_elem_in_space->className() << " related to "
                               << one_elem_in_space->m_RelatedElements.size() << " elements:" << std::endl;

                    std::vector<shared_ptr<IfcProduct> > related_prod = one_elem_in_space->m_RelatedElements;

                    // If there is a predefined list of undesired classes, exclude them
                    if (unallowed_classes.size() > 0)
                    {
                        for(uint j=0; j<related_prod.size(); j++)
                        {
                            shared_ptr<IfcProduct> one_related_prod = dynamic_pointer_cast<IfcProduct>( related_prod[j] );
                            if (!unallowed_classes[one_related_prod->className()])
                                elem_of_space[ one_related_prod->className() ].push_back( one_related_prod );
                        }
                    }
                    // Otherwise, load everything
                    else
                    {
                        for(uint j=0; j<related_prod.size(); j++)
                        {
                            shared_ptr<IfcProduct> one_related_prod = dynamic_pointer_cast<IfcProduct>( related_prod[j] );
                            elem_of_space[ one_related_prod->className() ].push_back( one_related_prod );
                        }
                    }
                }

                std::map < const char*, std::vector< shared_ptr<IfcProduct> > >::iterator it_elem_map (elem_of_space.begin());
                for(; it_elem_map != elem_of_space.end(); it_elem_map++)
                    std::wcout << "\t" << it_elem_map->first << " (" << it_elem_map->second.size() <<  ")" << std::endl;

            }
            std::wcout << std::endl;
        }
    }

    void Get_all_elements_in_IfcSpace (shared_ptr<IfcSpace> &space,
                                       std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > &elem_of_space)
    {
        std::map<const char *, bool, LCCtools::cmp_const_char> unallowed_classes;

        // Enter the banned classes here
        unallowed_classes[ "IfcBuildingElementProxy" ] = true;

        Get_all_elements_in_IfcSpace (space, elem_of_space, unallowed_classes);
    }


#endif

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extracts and creates in the LCC the real free space of an IfcSpace after applying boolean
    /// operations to substract the shapes of the furnitures out of the original space.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Extract_Real_Free_Space ( LCC& alcc, Dart_handle& space, vec_dart& o_spaces, R_Space& res/*, Nef_polyhedron& n_spaces*/ )
    {
        assert( space != LCC::null_handle );
        Dart_handle d = LCC::null_handle;

        // Convert the initial space to Nef
        Nef_polyhedron N_space;
        if (Convert_3cell_2_Nef3(alcc, space, N_space))
        {
            std::cout << "Number of points of the Nef: " << N_space.number_of_vertices() << std::endl;

            // Subtract from the space each inner object found
            for (uint i=0; i<o_spaces.size(); i++)
            {
                Nef_polyhedron N_obj;
                if ( Convert_3cell_2_Nef3(alcc, o_spaces[i], N_obj) )
                    N_space -= N_obj;
            }

            Polyhedron_EK P_space;
            N_space.convert_to_Polyhedron( P_space );
//                std::cout << "Number of point of the polyhedron: " << P_space.size_of_vertices() << std::endl;

            if (!P_space.is_empty() && P_space.is_valid())
                d = Build_lcc_from_exact_polyhedron_3(alcc, P_space);

            if (d != LCC::null_handle)
            {
                alcc.set_attribute<3>(d, alcc.create_attribute<3>());
                alcc.info<3>(d).set_label("R-Space");
                alcc.info<3>(d).set_id(std::to_string(vol_counter++));
//                alcc.info<3>(d).set_vol_label(volume_label++);
//                    Clean_Vol_Boundaries(alcc, d);
//                    vec_dart odpf_vol = one_dart_per_face_in_volume(alcc, d);
//                    Merge_coplanar_2sewed_2cells_set( alcc, odpf_vol );

                // collect the output
//                res.dart = d;
//                res.nef = N_space;
//                res.uid = alcc.info<3>(d).id();
//                res.space_id = alcc.info<3>(space).id();

                return true;
            }
        }

        return false;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extracts and creates in the LCC the real free space of an IfcSpace after applying boolean
    /// operations to substract the shapes of the furnitures out of the original space.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Extract_fspace_from_free_Space ( LCC& alcc, Dart_handle& space, Nef_polyhedron& N_space, vec_dart& fspace, vec_dart& res, std::string label)
    {
        if (fspace.size() > 0)
        {
            assert( space != LCC::null_handle );

//                int label = 0;
//                if ( alcc.attribute<3>(space[0])!=LCC::null_handle )
//                    label = alcc.info<3>(space[0]).label();

            Bbox_3 bb_space = LCCtools::Get_Bbox_vol(alcc, space);
//                Merge_coplanar_2sewed_2cells_in_volume(alcc, space);


            Polyhedron_EK P_space;
//                Convert_3cell_2_polyhedron3(alcc, space, P_space);
//                Nef_polyhedron N_space( P_space );

            for (uint i=0; i<fspace.size(); i++)
            {
                Bbox_3 bb_op = LCCtools::Get_Bbox_vol(alcc, fspace[i]);

                if ( do_overlap(bb_space, bb_op) )
                {
                    Nef_polyhedron op_fspace, N_op;
                    Convert_3cell_2_Nef3(alcc, fspace[i], N_op);

                    // try to get the intersection as separated volume
                    op_fspace = N_space.intersection( N_op );

                    Polyhedron_EK P_fspace;
                    if (op_fspace.is_simple())
                    {
                        N_space -= N_op;
                        op_fspace.convert_to_Polyhedron( P_fspace );
                        if ( !P_fspace.is_empty() && P_fspace.is_valid() )
                        {
//                                std::cout << "Is there something in the poly? " << P_fspace.size_of_facets() << " faces and "
//                                          << P_fspace.size_of_vertices() << " vertices." << std::endl;
                            Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_fspace);
                            alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                            alcc.info<3>(d).set_label(label);

                            res.push_back(d);
                        }
                    }
                }
            }

            P_space.clear();
            N_space.convert_to_Polyhedron( P_space );

            if (!P_space.is_empty() && P_space.is_valid())
            {
                alcc.remove_cell<3>( space);
                space = Build_lcc_from_exact_polyhedron_3(alcc, P_space);
                alcc.set_attribute<3>(space, alcc.create_attribute<3>());
//                alcc.info<3>(space).set_label("Brut_free_space");
//                    Clean_Vol_Boundaries(alcc, d);
//                    vec_dart odpf_vol = one_dart_per_face_in_volume(alcc, d);
//                    Merge_coplanar_2sewed_2cells_set( alcc, odpf_vol );
            }
            else
                return false;
        }
        return true;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extracts and creates in the LCC the real free space of an IfcSpace after applying boolean
    /// operations to substract the shapes of the furnitures and openings out of the original space.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Extract_fspace_from_free_Space_without_furnitures ( LCC& alcc, Dart_handle& space, vec_dart& openings, vec_dart& res, std::string label )
    {
        if (openings.size() > 0)
        {
            assert( space != LCC::null_handle );

//                int label = 0;
//                if ( alcc.attribute<3>(space[0])!=LCC::null_handle )
//                    label = alcc.info<3>(space[0]).label();

            Bbox_3 bb_space = LCCtools::Get_Bbox_vol(alcc, space);
//                Merge_coplanar_2sewed_2cells_in_volume(alcc, space);

            Nef_polyhedron N_space;
            Convert_3cell_2_Nef3(alcc, space, N_space);

            for (uint i=0; i<openings.size(); i++)
            {
                Bbox_3 bb_op = LCCtools::Get_Bbox_vol(alcc, openings[i]);

                if ( do_overlap(bb_space, bb_op) )
                {
                    Nef_polyhedron op_fspace, N_op;
                    Convert_3cell_2_Nef3(alcc, openings[i], N_op);

                    // try to get the intersection as separated volume
                    op_fspace = N_space.intersection( N_op );

                    Polyhedron_EK P_fspace;
                    if (op_fspace.is_simple())
                    {
                        N_space -= N_op;
                        op_fspace.convert_to_Polyhedron( P_fspace );
                        if ( !P_fspace.is_empty() && P_fspace.is_valid() )
                        {
                            Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_fspace);
                            alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                            alcc.info<3>(d).set_label(label);

                            res.push_back(d);
                        }
                    }
                }
            }

            Polyhedron_EK P_space;
            N_space.convert_to_Polyhedron( P_space );

            if (!P_space.is_empty() && P_space.is_valid())
            {
                alcc.remove_cell<3>( space);
                space = Build_lcc_from_exact_polyhedron_3(alcc, P_space);
                alcc.set_attribute<3>(space, alcc.create_attribute<3>());
//                alcc.info<3>(space).set_label("Brut_free_space");
//                    Clean_Vol_Boundaries(alcc, d);
//                    vec_dart odpf_vol = one_dart_per_face_in_volume(alcc, d);
//                    Merge_coplanar_2sewed_2cells_set( alcc, odpf_vol );
            }
            else
                return false;
        }
        return true;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Gets the buffers of the openings from the R-Space and create their 3-cells
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Get_Buffers_for_Openings(LCC& alcc, R_Space& r_space, vec_dart& storey_openings, std::vector<Sub_Space>& res, double buffer_size)
//    {
//        // Extrude the opening to make their buffer inside the free spaces of the rooms
//        std::cout << "\tEXTRUDING OPENINGS TO CREATE BUFFERS" << std::endl;
//        for(uint op=0; op<storey_openings.size(); op++)
//        {
//            double_pair_dart_vec3d cur_op_extru;
//            Merge_coplanar_2sewed_2cells_in_volume(alcc, storey_openings[op]);
//            // Keep the extrusion direction and faces to get back to proper size afterward
//            cur_op_extru = Extrude_opening_in_thickeness_direction( alcc, storey_openings[op], buffer_size, false );

//            //force contact to the space boundaries to avoid thin spaces
//            Make_contact_between_volA_and_volB(alcc, r_space.dart, storey_openings[op]);

//            // Nefs of the buffer and the openings (for the time of the boolean operations)
//            Nef_polyhedron N_buffer, N_op;
//            Convert_3cell_2_Nef3(alcc, storey_openings[op], N_op);

//            // try to get the intersection as separated volume
//            N_buffer = r_space.nef.intersection( N_op );

//            Polyhedron_EK P_buffer;
//            if (N_buffer.is_simple())
//            {
//                // Remove the buffer from the R-Space
//                r_space.nef -= N_op;

//                // Convert to polyhedron to check validity
//                N_buffer.convert_to_Polyhedron( P_buffer );
//                // If valid, import as a 3-cell
//                if ( !P_buffer.is_empty() && P_buffer.is_valid() )
//                {
//                    Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_buffer);
//                    alcc.set_attribute<3>(d, alcc.create_attribute<3>());
////                    alcc.info<3>(d).set_label("B-Space");
////                    alcc.info<3>(d).set_id(vol_counter++);
////                    alcc.info<3>(d).set_vol_label(volume_label++);

//                    // collect the output
//                    Sub_Space bs;
//                    bs.dart = d;
////                    bs.uid = alcc.info<3>(d).id();
//                    bs.space_id = r_space.space_id;
//                    bs.related_object = storey_openings[op];
//                    res.push_back( bs );
//                }
//            }

//            // resize the openings by changing extrusion direction to opposite
//            cur_op_extru.first.second = cur_op_extru.first.second * -1;
//            cur_op_extru.second.second = cur_op_extru.second.second * -1;

//            Extrude_face_in_vector_direction( alcc, cur_op_extru.first.first, cur_op_extru.first.second, buffer_size );
//            Extrude_face_in_vector_direction( alcc, cur_op_extru.second.first, cur_op_extru.second.second, buffer_size );
//        }

//        // Update the 3-cell of the R-Space
//        Polyhedron_EK P_space;
//        r_space.nef.convert_to_Polyhedron( P_space );

//        if (!P_space.is_empty() && P_space.is_valid())
//        {
//            // by copying the attributes of the old one to the new one
//            Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_space);
//            Transfer_Attributes_of_3cell(alcc, r_space.dart, d);

//            // and replacing the two of them
//            alcc.remove_cell<3>( r_space.dart);
//            r_space.dart =  d;
////                    Clean_Vol_Boundaries(alcc, d);
////                    vec_dart odpf_vol = one_dart_per_face_in_volume(alcc, d);
////                    Merge_coplanar_2sewed_2cells_set( alcc, odpf_vol );
//        }
//    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Gets the buffers/functional spaces of objects (O-Spaces) from the R-Space and create their 3-cells
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Get_Buffers_or_Functional_for_O_Spaces(LCC& alcc, R_Space& r_space, std::vector<O_Space*>& o_spaces,
//                                                std::vector<double>& add_to_buffer, std::vector<Sub_Space>& res,
//                                                std::string& label, double buffer_size)
//    {
//        assert ( o_spaces.size() == add_to_buffer.size() );

//        // Extrude the object spaces of the furnitures and extract their buffers out of the free space
//        // starting from the biggest to the smallest
//        std::cout << "\tEXTRUDING O-SPACES TO CREATE BUFFERS" << std::endl;
//        std::map<double, O_Space*> volume_sorted_o_spaces;
//        for(uint os=0; os<o_spaces.size(); os++)
//        {
//            // extruding the O-Spaces for considering buffers
////                Extrude_box_volume(alcc, o_spaces[os]->dart, (buffer_size + add_to_buffer[os]) );
//            Extrude_box_volume(alcc, o_spaces[os]->dart, buffer_size );

//            //force contact to the space boundaries to avoid thin spaces (this messes up the oiginal size recovering)
////                Make_contact_between_volA_and_volB(alcc, r_space.dart, o_spaces[os]->dart);

//            // sort by biggest volume
//            volume_sorted_o_spaces[ Get_volume_magnitude(alcc, o_spaces[os]->dart) ] = o_spaces[os];
//        }

//        if (o_spaces.size() > 1)
//        {
//            o_spaces.clear();
//            std::map<double, O_Space*>::reverse_iterator it_os(volume_sorted_o_spaces.rbegin());
//            for( ; it_os!=volume_sorted_o_spaces.rend(); it_os++ )
//                o_spaces.push_back( it_os->second );
//        }


//        // Process the O-Spaces, by starting from the biggest
//        for(uint os=0; os<o_spaces.size(); os++)
//        {
//            // only if the O-Space doesn't have B/F-Space yet
//            if ( (label == "B-Space" && o_spaces[os]->b_space == LCC::null_handle )
//                 || (label == "F-Space" && o_spaces[os]->f_space == LCC::null_handle) )
//            {
//                // Nefs of the buffer and the O-Space (for the time of the boolean operations)
//                Nef_polyhedron N_buffer, N_obj;
//                Convert_3cell_2_Nef3(alcc, o_spaces[os]->dart, N_obj);

//                // try to get the intersection as separated volume
//                N_buffer = r_space.nef.intersection( N_obj );

//                Polyhedron_EK P_buffer;
//                if (N_buffer.is_simple())
//                {
//                    // Remove the buffer from the R-Space
//                    r_space.nef -= N_obj;

//                    // Convert to polyhedron to check validity
//                    N_buffer.convert_to_Polyhedron( P_buffer );
//                    // If valid, import as a 3-cell
//                    if ( !P_buffer.is_empty() && P_buffer.is_valid() )
//                    {
//                        Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_buffer);
//                        alcc.set_attribute<3>(d, alcc.create_attribute<3>());
////                        alcc.info<3>(d).set_label(label);
////                        alcc.info<3>(d).set_id(vol_counter++);
////                        alcc.info<3>(d).set_vol_label(volume_label++);

//                        // collect the output
//                        Sub_Space bs;
//                        bs.dart = d;
////                        bs.uid = alcc.info<3>(d).id();
//                        bs.space_id = r_space.space_id;
//                        bs.related_object = o_spaces[os]->dart;
//                        res.push_back( bs );

//                        //Store the link to the O-Space to which it belongs
//                        o_spaces[os]->b_space = bs.dart;
//                    }
//                }

//                // Update the 3-cell of the R-Space
//                Polyhedron_EK P_space;
//                r_space.nef.convert_to_Polyhedron( P_space );

//                if (!P_space.is_empty() && P_space.is_valid())
//                {
//                    // by copying the attributes of the old one to the new one
//                    Dart_handle d = Build_lcc_from_exact_polyhedron_3(alcc, P_space);
//                    Transfer_Attributes_of_3cell(alcc, r_space.dart, d);

//                    // and replacing the two of them
//                    alcc.remove_cell<3>( r_space.dart);
//                    r_space.dart =  d;
//                }
//            }

//            // back to the original size of the O-Spaces
////                Extrude_box_volume(alcc, o_spaces[os]->dart, (buffer_size + add_to_buffer[os]), false);
//            Extrude_box_volume(alcc, o_spaces[os]->dart, buffer_size, false);
//        }
//    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Check among the faces resulting from 2D arrangement which ones
    /// to keep to rebuild a volume without overlaping surfaces
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    vec_dart Rebuild_faces_after_Arrangement_Correction(LCC& alcc, arr_tuple& arr)
//    {
//        vec_dart res;
//        Plane Pl = std::get<0>(arr);
//        std::vector <Triangle_2<EK> > final_triz = std::get<1>(arr);
//        std::vector <Polygon_2<EK> > poly_ori_faces = std::get<2>(arr);

//        Polygon_2<EK> final_poly;
//        std::vector<Point> new_face;

//        for ( uint tr=0; tr < final_triz.size(); tr++ )
//        {
//            Triangle_2_Poly( final_triz[tr], final_poly );

//            // Try to detect faces specifically created by the arrangement (openings usually)
//            // a point lying inside such faces is not lying inside any face used to build the arrangement
////                if ( ( fabs( to_double( final_poly.area() ) ) / pow(scale_f, 2) )  > (eps2*eps2) )
//            {
//                Point_2d inside_pt;
//                if ( Get_point_in_2DPoly<EK> ( final_poly, inside_pt ) )
//                {
//                    bool go = true;
//                    uint count = 0;
//                    for(uint i=0; i<poly_ori_faces.size(); i++)
//                    {
//                        if ( poly_ori_faces[i].has_on_bounded_side( inside_pt ) )
//                        {
//                            count++;

//                            // if 2 original faces overlap, don't create the arrangement face
//                            if (count == 2)
//                            {
//                                go = false;
//                                break;
//                            }
//                        }
//                    }

//                    if (count == 0)
//                        go = false;

//                    if ( go )
//                    {
//                        std::cout << "\tSurface area: " << final_poly.area() << std::endl;
//                        typename Polygon_2<EK>::Vertex_iterator pt(final_poly.vertices_begin());
//                        for (; pt!= final_poly.vertices_end(); pt++)
//                            new_face.push_back(  Pl.to_3d( exa2inexa( *pt ) ) );

//                        Dart_handle d = Insert_new_2cell( alcc, new_face );
////                                    Print_face( alcc, d );

//                        res.push_back(d);
//                        new_face.clear();
//                    }

//                    else
//                        std::cout << "\n\tDIDN'T GO!!!!" << std::endl;
//                }

//                else
//                    std::cout << "\n\tDIDN'T GET NOTHIIIIIIIIIIINNNG!!!!" << std::endl;
//            }
//        }

//        return res;
//    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Perfoms 3-cell reconstruction from a set of unconnected 2-cells
    /// returns true if the 3-cell is closed (no 2-free cells).
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    Dart_handle Create_3cell_from_set_of_2cells(LCC& alcc, vec_dart& faces)
    {
        assert (faces.size() > 0);
        Dart_handle res = LCC::null_handle;

        // each key is unique for one edge
        std::map< std::string, Dart_handle, LCCtools::cmp_string > edge_map;
        std::map< std::string, Dart_handle, LCCtools::cmp_string >::iterator it, it2;
        Point pt1, pt2;

        for(uint i=0; i<faces.size(); i++)
        {
            typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, faces[i]);
            for (; itr.cont(); itr++)
            {
                pt1 = alcc.point(itr);
                pt2 = alcc.point( alcc.beta<1>(itr) );

                std::string key1 = get_unique_key_pt(pt1),
                            key2 = get_unique_key_pt(pt2);

                // insert the edge in the map
                std::string key ( key1 + key2 );
//                    std::cout << "\n" <<  key.c_str() << std::endl;
                std::cout << key << std::endl;
                edge_map[ key ] = itr;

                // If the opposite edge is already in the map, sew2 them
                it = edge_map.find( std::string( key2 + key1 ) );
                if ( it != edge_map.end() )
                  alcc.template link_beta<2>(itr, it->second);
            }
        }

        // Check if all the edges are 2-sewn
        for( it=edge_map.begin(); it!=edge_map.end(); it++ )
            if (alcc.is_free<2>(it->second))
            {
                std::string k1 = get_unique_key_pt( alcc.point(it->second) ),
                            k2 = get_unique_key_pt( alcc.point(alcc.beta<1>(it->second)) );
                std::string k (k2 + k1);

                std::cout << "\nNot 2-sewn to: " << k2 << " " << k1 << std::endl;
                std::cout << "\t" << alcc.point(it->second) << std::endl;
                std::cout << "\t" << alcc.point(alcc.beta<1>(it->second)) << std::endl;

                it2 = edge_map.find( k );

                if (it2 != edge_map.end())
//                    if ( edge_map.count( k ) > 0 )
                    std::cout << "BUT SHOULD!!!" << std::endl;

                return res;
            }

        std::cout << "Nb of points: " << faces.size() * 3 << std::endl;
        std::cout << "Nb of darts: " << edge_map.size() << std::endl;
        std::cout << "Created a closed 3-cell!" << std::endl;

        res = faces[0];
        return res;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Removes the flat extensions on volume boundaries resulting from boolean operations
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
//    void Clean_Vol_Boundaries( LCC& alcc, Dart_handle& vol)
//    {
//        std::cout << "Values of the global variables: " << std::endl;
//        std::cout << "Height: " << V_height << std::endl;
//        std::cout << "eps_angle: " << eps_angle << std::endl;

//        vec_dart odpf_cur, odpf_final;
//        std::vector <vec_dart > cop_f_in_vol = Collect_parallel_faces_in_volume(alcc, vol, true);

//        for (uint i=0; i<cop_f_in_vol.size(); i++)
//        {
//            arr_tuple arr = Arrangement_partitioning( alcc, cop_f_in_vol[i] );
//            odpf_cur = Rebuild_faces_after_Arrangement_Correction(alcc, arr);

//            odpf_final.insert( odpf_final.end(), odpf_cur.begin(), odpf_cur.end() );
//        }

//        Merge_close_points(alcc, odpf_final);
//        vol = Create_3cell_from_set_of_2cells(alcc, odpf_final);
//    }



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extrudes of length "l" a given face in the given direction
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Extrude_face_in_vector_direction( LCC& alcc, Dart_handle& face, Vector& dir, double l)
    {
//            std::cout << "\n\tI AM EXTRUDIIIIIING!!" << std::endl;
//            double l = sqrt(to_double( v.squared_length() ));
//            LCCtools::unit_normal(v);
        typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, face);
        for (; itr.cont(); itr++)
        {
            Point pt = alcc.point( itr );
            alcc.point( itr ) = Point( pt.x() + dir.x()*l,
                                       pt.y() + dir.y()*l,
                                       pt.z() + dir.z()*l);
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extrudes all faces of a box volume (all faces are perpendicular to each other) in a direction
    /// increasing (increase = true) the volume, or decreasing it (increase = false).
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Extrude_box_volume( LCC& alcc, Dart_handle& vol, double l, bool increase)
    {
        Vector dir;
        double orientation;

        if (increase)
            orientation = -1.0;
        else
            orientation= 1.0;

        // for each 2-cell of the box volume
        typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
                it = alcc.one_dart_per_incident_cell<2,3>(vol).begin(),
                it_end = alcc.one_dart_per_incident_cell<2,3>(vol).end();
        for (; it != it_end; it++)
        {
            dir = LCCtools::poly_normal( alcc, it ) * orientation;
            Extrude_face_in_vector_direction( alcc, it, dir, l );
        }
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Extrudes the part of furnitures that very close to room borders to force the contact.
    /// This allows to remove flat parts of volumes after boolean operations.
    /// Furnitures are assumed to be boxes. Works only if the faces of the 2 volumes are parallel
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Make_contact_between_volA_and_volB(LCC& alcc, Dart_handle &volA, Dart_handle &volB, double dist, bool handle_floating)
    {
        // For floating furnitures. If the distance is the same for two faces,
        // only one will be picked anyway, so no need for a vector of the pair.
        std::map<double, std::pair<Dart_handle, Vector> > closest_contact_face;
        bool floating = true;

        // Collect normals of the VolA faces once for all
        std::vector< Vector > room_faces_normals;
        vec_dart volA_vec;
        typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
                itA = alcc.one_dart_per_incident_cell<2,3>(volA).begin(),
                itA_end = alcc.one_dart_per_incident_cell<2,3>(volA).end();
        for (; itA != itA_end; itA++)
        {
            room_faces_normals.push_back( LCCtools::poly_normal( alcc, itA ) );
            // Keep one dart per face to avoid iterator mess...
            volA_vec.push_back( itA );
        }

        // for each 2-cell of volB
        typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
                it_furni = alcc.one_dart_per_incident_cell<2,3>(volB).begin(),
                itf_end = alcc.one_dart_per_incident_cell<2,3>(volB).end();
        for (; it_furni != itf_end; it_furni++)
        {
            Vector Nf = LCCtools::poly_normal( alcc, it_furni ),
                   Nf_opp = -Nf;
            Plane Pf (alcc.point(it_furni), Nf);

            // Compare its coplanarity and proximity to the room faces
            for (uint i=0; i<volA_vec.size(); i++)
            {
                if ( vectors_are_eps_colinear( Nf, room_faces_normals[i] ) )
                {
                    Vector Vr ( alcc.point(volA_vec[i]),
                                Pf.projection( alcc.point(volA_vec[i]) ) );
                    double l = sqrt(to_double( Vr.squared_length() ));
                    LCCtools::unit_normal(Vr);
                    double ang = compute_3d_angle( Nf, Vr );

                    closest_contact_face[ l ] = std::pair<Dart_handle, Vector>( it_furni, Nf_opp );

                    // If they are close
                    if ( l <= dist*scale_f )
                    {
                        floating = false;

                        // If Nf and Vr has the same direction, the furniture is completely inside the space
                        // Otherwise it the current face is slightly out of the space, so the contact is already sure.
                        if (ang < (0.0 + eps_angle) || ang > (360.0 - eps_angle) )
                            Extrude_face_in_vector_direction( alcc, it_furni, Nf_opp, l+coplanarity);

                        // Break the process for the given furniture face
                        break;
                    }
                }
            }
        }

        if (floating && handle_floating )
        {
            std::cout << "Floating object detected!!!" << std::endl;
            double l = closest_contact_face.begin()->first;
            Dart_handle d_furni = closest_contact_face.begin()->second.first;
            Vector Nf_opp = closest_contact_face.begin()->second.second;

            Extrude_face_in_vector_direction( alcc, d_furni, Nf_opp, l+coplanarity);
        }
    }



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Checks wich simplified furnitures are intersecting and aggregate them in a bigger OBB.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Aggregate_Set_of_3cells (LCC& alcc, vec_dart& all_furni, vec_dart& res, bool buffer, double buffer_size)
    {
        std::cout << "\n\n$$$$$$$$$$$$ New aggregation to go: " << all_furni.size() << std::endl;

        std::map<double, vec_dart> furni_by_volume;
        double vol;

        // Sort the furnitures (boxes) according to their volume
        for(uint i=0; i<all_furni.size(); i++)
        {
            vol = Get_volume_magnitude(alcc, all_furni[i]);
            furni_by_volume[vol].push_back( all_furni[i] );
        }

        int checked = alcc.get_new_mark();
        CGAL_assertion( checked!=-1 );

        vec_dart bigger_to_smaller_vol;
        std::vector< vec_dart > vol_to_aggregate;
        std::map<double, vec_dart>::reverse_iterator it( furni_by_volume.rbegin() );
        for(; it != furni_by_volume.rend(); it++)
            bigger_to_smaller_vol.insert( bigger_to_smaller_vol.end(), it->second.begin(), it->second.end() );

        std::cout << "bigger to smaller: " << bigger_to_smaller_vol.size() << std::endl;

        vec_dart one_aggregation;
        uint isolated = 0;
        for(uint i=0; i<bigger_to_smaller_vol.size();i++)
        {
            Dart_handle d = bigger_to_smaller_vol[i];
            if ( !alcc.is_marked(d, checked) )
            {
                std::cout << "Got one ref! " << std::endl;
                alcc.mark(d, checked);

                // extrude
                if (buffer)
                {
//                        std::cout << "Volume: " << Get_volume_magnitude(alcc, d) << std::endl;
                    Extrude_box_volume(alcc, d, buffer_size);
                }

                for(uint j=i+1; j<bigger_to_smaller_vol.size(); j++)
                {
                    Dart_handle d_cur = bigger_to_smaller_vol[j];
                    if ( !alcc.is_marked(d_cur, checked) )
                    {
//                            Polyhedron cur_poly;
//                            Convert_3cell_2_polyhedron3(alcc, d, cur_poly);
//                            Nef_polyhedron cur_nef (cur_poly);

//                            they_intersect = ref_nef.intersection( cur_nef );
//                            std::cout << "# Number of volumes of intersection: " << they_intersect.interior().number_of_volumes() << std::endl;
//                            std::cout << "# Number of volumes of intersection: " << they_intersect.interior().number_of_facets() << std::endl;
//                            std::cout << "# Number of volumes of intersection: " << they_intersect.interior().number_of_edges() << std::endl;
                        if (Boxes_intersect(alcc, d, d_cur))
                        {
                            std::cout << "################################## Intersection in here!!! " << std::endl;
                            alcc.mark(d_cur, checked);
                            one_aggregation.push_back( d_cur );
                        }
                    }
                }

                // back to the right size
                if (buffer)
                    Extrude_box_volume(alcc, d, buffer_size, false);

                if (one_aggregation.size() == 0)
                    isolated++;

                // add the ref volume
                one_aggregation.push_back( d );
                vol_to_aggregate.push_back(one_aggregation);
                one_aggregation.clear();
            }
        }

        alcc.unmark_all(checked);
        alcc.free_mark(checked);

        std::cout << "Number of aggregations: " << vol_to_aggregate.size() << std::endl;
        std::cout << "Number of isolated: " << isolated << std::endl;

        if ( vol_to_aggregate.size() < all_furni.size() )
        {
            /// In case there are already resulting functional spaces in res!
            res.clear();

            /// Perform the OBB simplification of the set of volumes to aggregate
            Dart_handle d_res;
            std::vector<Point> all_points;
            for(uint i=0; i<vol_to_aggregate.size(); i++)
            {
                std::cout << "Number of volumes to aggregate: " << vol_to_aggregate[i].size() << std::endl;
                // collect all the points of the current volume set to aggregate
                for( uint j=0; j<vol_to_aggregate[i].size(); j++ )
                {
                    // for each 0-cell of the volume
                    typename LCC::Base::One_dart_per_incident_cell_range<0,3>::iterator
                            itr = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j]).begin(),
                            itr_end = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j]).end();
                    for (; itr != itr_end; itr++)
                        all_points.push_back( alcc.point( itr ) );

                    // if the dart is the result of a previous aggregation, remove it.
                    if ( alcc.attribute<3>( vol_to_aggregate[i][j] ) != LCC::null_handle
//                         && alcc.info<3>( vol_to_aggregate[i][j] ).label() == "Functional_space"
                         )
                        alcc.remove_cell<3>( vol_to_aggregate[i][j]);
                }

                std::cout << "Number of collected points: " << all_points.size() << std::endl;
                d_res = OBB_Pointset_Simplification(alcc, all_points);
                res.push_back(d_res);
                alcc.set_attribute<3>(d_res, alcc.create_attribute<3>());
//                alcc.info<3>(d_res).set_label("Functional_space");
                all_points.clear();
            }


            /// Perform Minkowski Sum of volumes to aggregate
//                Dart_handle d_res;
//                std::cout << "Number of aggregations: " << vol_to_aggregate.size() << std::endl;
//                for(uint i=0; i<vol_to_aggregate.size(); i++)
//                {
//                    std::cout << "Number of volumes to aggregate: " << vol_to_aggregate[i].size() << std::endl;
//                    uint sz = vol_to_aggregate[i].size();

//                    // Create the nef of the last volume that corresponds to the biggest one
//                    Nef_polyhedron Mk_agg;
//                    Convert_3cell_2_Nef3(alcc, vol_to_aggregate[i][sz-1], Mk_agg);

//                    // extrude back to the right size
//                    Extrude_box_volume(alcc, vol_to_aggregate[i][sz-1], extru, false);

//                    // collect all the points of the current volume set to aggregate
//                    for( uint j=0; j<sz-1; j++ )
//                    {
//                        Nef_polyhedron Mk_agg_cur;
//                        Convert_3cell_2_Nef3(alcc, vol_to_aggregate[i][sz-1], Mk_agg_cur);

//                        // Minkowski sum of the big vol and the smaller
//                        Mk_agg = minkowski_sum_3(Mk_agg, Mk_agg_cur);
//                    }

//                    Polyhedron P_agg;
//                    Mk_agg.convert_to_Polyhedron( P_agg );

//                    if (!P_agg.is_empty() && P_agg.is_valid())
//                    {
//                        d_res = Build_lcc_from_exact_polyhedron_3(alcc, P_agg);
//                        res.push_back(d_res);
//                        alcc.set_attribute<3>(d_res, alcc.create_attribute<3>());
//                        alcc.info<3>(d_res).set_label(30);
//                    }
//                }

            /// If volumes were aggregated, then check intersection between aggregations until nothing intersects anymore
            Aggregate_Set_of_3cells(alcc, res, res);
        }
        else
            std::cout << "No more aggregations! remains " << res.size() << " including "
                      << isolated << " isolated." << std::endl;

//            goto
//                get_out;
//            get_out:
//                std::cout << "I'm out" << std::endl;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Checks wich simplified furnitures are intersecting and aggregate them in a bigger OBB.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void Aggregate_Set_of_O_Spaces (LCC& alcc, std::vector<O_Space>& all_obj, std::vector<O_Space>& res,
                                    int loop, bool buffer, double buffer_size, int mode)
    {
        std::cout << "\n\n$$$$$$$$$$$$ New aggregation to go: " << all_obj.size() << std::endl;
        loop++;
        std::cout << "\tloop " << loop << std::endl;

        std::map<double, std::vector<O_Space> > obj_by_volume_magnitude;
        double vol;

        // Sort the furnitures (boxes) according to their volume
        for(uint i=0; i<all_obj.size(); i++)
        {
            vol = Get_volume_magnitude(alcc, all_obj[i].dart);
            obj_by_volume_magnitude[vol].push_back( all_obj[i] );
        }

        int checked = alcc.get_new_mark();
        CGAL_assertion( checked!=-1 );

        std::vector<O_Space> bigger_to_smaller_vol;
        std::vector< std::vector<O_Space> > vol_to_aggregate;
        std::map<double, std::vector<O_Space> >::reverse_iterator it( obj_by_volume_magnitude.rbegin() );
        for(; it != obj_by_volume_magnitude.rend(); it++)
            bigger_to_smaller_vol.insert( bigger_to_smaller_vol.end(), it->second.begin(), it->second.end() );

        std::cout << "bigger to smaller: " << bigger_to_smaller_vol.size() << std::endl;

        std::vector<O_Space> one_aggregation;
        uint isolated = 0;
        for(uint i=0; i<bigger_to_smaller_vol.size();i++)
        {
            Dart_handle d = bigger_to_smaller_vol[i].dart;
            if ( !alcc.is_marked(d, checked) )
            {
                std::cout << "Got one ref! " << std::endl;
                alcc.mark(d, checked);

                // extrude
                if (buffer)
                {
//                        std::cout << "Volume: " << Get_volume_magnitude(alcc, d) << std::endl;
                    Extrude_box_volume(alcc, d, buffer_size);
                }

                for(uint j=i+1; j<bigger_to_smaller_vol.size(); j++)
                {
                    O_Space s_cur = bigger_to_smaller_vol[j];
                    if ( !alcc.is_marked(s_cur.dart, checked) )
                    {
                        if (Boxes_intersect(alcc, d, s_cur.dart))
                        {
                            std::cout << "################################## Intersection in here!!! " << std::endl;
                            alcc.mark(s_cur.dart, checked);
                            one_aggregation.push_back( s_cur );
                        }
                    }
                }

                // back to the right size
                if (buffer)
                    Extrude_box_volume(alcc, d, buffer_size, false);

                if (one_aggregation.size() == 0)
                    isolated++;

                // add the ref volume
                one_aggregation.push_back( bigger_to_smaller_vol[i] );
                vol_to_aggregate.push_back(one_aggregation);
                one_aggregation.clear();
            }
        }

        alcc.unmark_all(checked);
        alcc.free_mark(checked);

        std::cout << "Number of aggregations: " << vol_to_aggregate.size() << std::endl;
        std::cout << "Number of isolated: " << isolated << std::endl;

        /// OBB based aggregation
        if (mode == 0)
        {
            if ( loop == 1 || vol_to_aggregate.size() < all_obj.size() )
            {
                /// In case there are already resulting object spaces in res!
                res.clear();

                /// Perform the OBB simplification of the set of volumes to aggregate
                std::vector<Point> all_points;
                for(uint i=0; i<vol_to_aggregate.size(); i++)
                {
                    O_Space s_res;
                    std::cout << "Number of volumes to aggregate: " << vol_to_aggregate[i].size() << std::endl;
                    // collect all the points of the current volume set to aggregate
                    for( uint j=0; j<vol_to_aggregate[i].size(); j++ )
                    {
                        // for each 0-cell of the volume
                        typename LCC::Base::One_dart_per_incident_cell_range<0,3>::iterator
                                itr = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j].dart).begin(),
                                itr_end = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j].dart).end();
                        for (; itr != itr_end; itr++)
                            all_points.push_back( alcc.point( itr ) );

                        // if the dart is the result of a previous aggregation, remove it.
                        if ( alcc.attribute<3>( vol_to_aggregate[i][j].dart ) != LCC::null_handle
//                             && alcc.info<3>( vol_to_aggregate[i][j].dart ).label() == "O-Space"
                             )
                        {
                            alcc.remove_cell<3>( vol_to_aggregate[i][j].dart);
                            // Keep the list of contained objects of the O-Space
                            s_res.contained_obj.insert( s_res.contained_obj.end(),
                                                        vol_to_aggregate[i][j].contained_obj.begin(),
                                                        vol_to_aggregate[i][j].contained_obj.end() );
                        }

                        // Keep the contained elements
                        s_res.contained_obj.push_back( vol_to_aggregate[i][j].dart );
                    }

                    std::cout << "Number of collected points: " << all_points.size() << std::endl;
                    s_res.dart = OBB_Pointset_Simplification(alcc, all_points);
                    if (s_res.dart != LCC::null_handle)
                    {
                        s_res.space_id = vol_to_aggregate[0][0].space_id;
                        s_res.uid = vol_counter++;
                        alcc.set_attribute<3>(s_res.dart, alcc.create_attribute<3>());
                        alcc.info<3>(s_res.dart).set_label("O-Space");
                        alcc.info<3>(s_res.dart).set_id(s_res.uid);
//                        alcc.info<3>(s_res.dart).set_vol_label(volume_label++);

                        res.push_back(s_res);
                    }
                    all_points.clear();
                }

                /// If volumes were aggregated, then check intersection between aggregations until nothing intersects anymore
                Aggregate_Set_of_O_Spaces(alcc, res, res, loop, buffer, buffer_size);
            }
            else
                std::cout << "No more aggregations! remains " << res.size() << " including "
                          << isolated << " isolated." << std::endl;
        }
        /// Alpha-shape based aggregation
//        else if (mode == 1)
//        {
//            /// In case there are already resulting object spaces in res!
//            res.clear();
//            uint step = 10;

//            /// Perform the OBB simplification of the set of volumes to aggregate
//            std::vector<Point> all_points;
//            for(uint i=0; i<vol_to_aggregate.size(); i++)
//            {
//                O_Space s_res;
//                std::cout << "Number of volumes to aggregate: " << vol_to_aggregate[i].size() << std::endl;
//                // collect all the points of the current volume set to aggregate
//                for( uint j=0; j<vol_to_aggregate[i].size(); j++ )
//                {
//                    // for each 1-cell of the volume
//                    typename LCC::Base::One_dart_per_incident_cell_range<0,3>::iterator
//                            itr = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j].dart).begin(),
//                            itr_end = alcc.one_dart_per_incident_cell<0,3>(vol_to_aggregate[i][j].dart).end();
//                    for (; itr != itr_end; itr++)
//                    {
//                        /// Adding additional points between every segment extremity to optimize Alpha-shape simplification
//                        Point pt1 = alcc.point( itr ),
//                              pt2 = alcc.point( alcc.beta<1>(itr) );

//                        all_points.push_back(pt1);
//                        if ( Vector(pt1, pt2).squared_length() > (0.4*0.4) )
//                        {
//                            int cx=1, cy=1, cz=1;
//                            if ( pt1.x() > pt2.x() )
//                                cx = -1;
//                            if ( pt1.y() > pt2.y() )
//                                cy = -1;
//                            if ( pt1.z() > pt2.z() )
//                                cz = -1;

//                            double dx = fabs(pt2.x() - pt1.x()) / step,
//                                   dy = fabs(pt2.y() - pt1.y()) / step,
//                                   dz = fabs(pt2.z() - pt1.z()) / step;

//                            for(uint i=1; i<=step; i++)
//                            {
//                                all_points.push_back( Point( pt1.x() + (cx*dx*i),
//                                                             pt1.y() + (cy*dy*i),
//                                                             pt1.z() + (cz*dz*i)) );
//                            }
//                        }
//                    }

//                    // if the dart is the result of a previous aggregation, remove it.
//                    if ( alcc.attribute<3>( vol_to_aggregate[i][j].dart ) != LCC::null_handle
////                         && alcc.info<3>( vol_to_aggregate[i][j].dart ).label() == "O-Space"
//                         )
//                    {
//                        alcc.remove_cell<3>( vol_to_aggregate[i][j].dart);
//                        // Keep the list of contained objects of the O-Space
//                        s_res.contained_obj.insert( s_res.contained_obj.end(),
//                                                    vol_to_aggregate[i][j].contained_obj.begin(),
//                                                    vol_to_aggregate[i][j].contained_obj.end() );
//                    }

//                    // Keep the contained elements
//                    s_res.contained_obj.push_back( vol_to_aggregate[i][j].dart );
//                }

//                std::cout << "Number of collected points: " << all_points.size() << std::endl;

////                    for(uint i=0; i<all_points.size(); i++)
////                        std::cout << all_points[i] << std::endl;

//                s_res.dart = Concave_hull_Pointset_Simplification(alcc, all_points);
//                if (s_res.dart != LCC::null_handle)
//                {
//                    s_res.space_id = vol_to_aggregate[0][0].space_id;
//                    s_res.uid = vol_counter++;
//                    alcc.set_attribute<3>(s_res.dart, alcc.create_attribute<3>());
////                    alcc.info<3>(s_res.dart).set_label("O-Space");
////                    alcc.info<3>(s_res.dart).set_id(s_res.uid);
////                    alcc.info<3>(s_res.dart).set_vol_label(volume_label++);

//                    res.push_back(s_res);
//                }
//                all_points.clear();
//            }
//        }
    }



    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// As it names says... Thickness is considered to be the smallest distance among
    /// the farthest parallel faces of the given volume
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    double_pair_dart_vec3d Extrude_opening_in_thickeness_direction(LCC& alcc, Dart_handle& vol, double l, bool increase)
    {
        std::vector <vec_dart> pll_face_set;
        pll_face_set = Collect_parallel_faces_in_volume(alcc, vol);

        // sort by distance the indices of the pair of parallel faces
        std::map< double, uint > faces_idx_by_dist;

        //For now, we assume that openings are just simple rectangular boxes, so only 3 directions.
        assert( pll_face_set.size() == 3 );
        for (uint i=0; i<pll_face_set.size(); i++)
        {
            // Thus there are for each direction one pair of parallel faces.
            assert( pll_face_set[i].size() == 2 );

//                std::cout << "\nFace pair " << i+1 << std::endl;
//                Show_face(alcc, pll_face_set[i][0]);
//                Show_face(alcc, pll_face_set[i][1]);

            // Create the plane containing the first face
            Plane Pl ( alcc.point( pll_face_set[i][0] ),
                       LCCtools::poly_normal(alcc, pll_face_set[i][0]) );

            // And project a point of the second face on the plane to get their distance
            Vector dist ( alcc.point( pll_face_set[i][1] ),
                          Pl.projection( alcc.point( pll_face_set[i][1] ) ) );

            faces_idx_by_dist[ dist.squared_length() ] = i;

            Vector show = dist;
            LCCtools::unit_normal( show );
//                std::cout << "Distance " << dist.squared_length() << " for dir " <<  show << ", idx ---> " << i << std::endl;
        }

        // The first element of the map indicates the pair of faces to extrude
        uint idx = faces_idx_by_dist.begin()->second;
        Vector dir1, dir2;
        if (increase)
        {
            dir1 = -1 * LCCtools::poly_normal( alcc, pll_face_set[idx][0] ),
            dir2 = -1 * LCCtools::poly_normal( alcc, pll_face_set[idx][1] );
        }
        else
        {
            dir1 = LCCtools::poly_normal( alcc, pll_face_set[idx][0] ),
            dir2 = LCCtools::poly_normal( alcc, pll_face_set[idx][1] );
        }

//            std::cout << "Extrusion dir: " << dir1 << " and " << dir2 << " of " << l << std::endl;
        Extrude_face_in_vector_direction(alcc, pll_face_set[idx][0], dir1, l);
        Extrude_face_in_vector_direction(alcc, pll_face_set[idx][1], dir2, l);

        double_pair_dart_vec3d res;
        res.first.first = pll_face_set[idx][0]; res.first.second = dir1;
        res.second.first = pll_face_set[idx][1]; res.second.second = dir2;

        return res;
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Returns true if a corner of box2 lies inside each face of box1 when projected on it.
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Boxes_intersect (LCC& alcc, Dart_handle box1, Dart_handle box2)
    {
        std::cout << "\n\nNew box to box intersection check" << std::endl;
        // for each 2_cell of the box1
        typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
                it1 = alcc.one_dart_per_incident_cell<2,3>(box1).begin(),
                it1_end = alcc.one_dart_per_incident_cell<2,3>(box1).end();
        for (; it1 != it1_end; it1++)
        {
            uint count = 0;
            Plane Pl(alcc.point(it1), LCCtools::poly_normal(alcc, it1));
            Polygon_2<IK> cur_2cell_poly;

            // for each edge of the 2-cell
            typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, it1);
            for (; itr.cont(); itr++)
                cur_2cell_poly.push_back(Pl.to_2d( alcc.point( itr ) ) );

            // for each 0_cell of the box2 check if its projection lies inside the 2D poly
            typename LCC::Base::One_dart_per_incident_cell_range<0,3>::iterator
                    it2 = alcc.one_dart_per_incident_cell<0,3>(box2).begin(),
                    it2_end = alcc.one_dart_per_incident_cell<0,3>(box2).end();
            for (; it2 != it2_end; it2++)
            {
                // Add the mid point to higher the accuracy.... (still not enough though...^_^')
                Point mid_pt ( ((alcc.point(it2).x() + alcc.point(alcc.beta<1>(it2)).x())/2.0),
                               ((alcc.point(it2).y() + alcc.point(alcc.beta<1>(it2)).y())/2.0),
                               ((alcc.point(it2).z() + alcc.point(alcc.beta<1>(it2)).z())/2.0) );

                if (cur_2cell_poly.has_on_bounded_side( Pl.to_2d( alcc.point( it2 ) ) )
                        || cur_2cell_poly.has_on_bounded_side( Pl.to_2d( mid_pt )) )
                    count++;
            }

            // if no intersection found for that face, there is no volume intersection
            if (count == 0)
                return false;
        }

        return true;
    }


#ifdef IFCPP_ON
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Returns .
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    std::vector<Nef_polyhedron> Get_Nef_of_furniture (LCC& alcc, shared_ptr<IfcProduct>& furni, shared_ptr<GeometryConverter>& geom_conv)
    {
        double extru = 0.1;
        shared_ptr<ProductShapeData> product_shape ( new ProductShapeData() );
        shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(furni);
        product_shape->m_ifc_object_definition = obj_def;

        geom_conv->convertIfcProductShape( product_shape );
        std::wcout << "vec shapes: " << product_shape->m_vec_representations.size() << std::endl;

        // Convert each mesh of the furiture into a box and make a union of their Nefs
        std::vector<Nef_polyhedron> Nef_box_set_of_furni;
        vec_dart box_set_of_furni;

        // Dart of the resulting simplified 3-cell
        Dart_handle d = LCC::null_handle;
        Nef_polyhedron N;

        // for each IfcProduct, there can be mulitple geometric representation items:
        std::vector<shared_ptr<RepresentationData> >& vec_rep = product_shape->m_vec_representations;
        for( size_t i_rep = 0; i_rep < vec_rep.size(); ++i_rep )
        {
            shared_ptr<RepresentationData>& rep_data = vec_rep[i_rep];
            std::vector<shared_ptr<ItemShapeData> >& vec_item_data = rep_data->m_vec_item_data;
            std::wcout << "\t# items: " << vec_item_data.size() << std::endl;

            for( size_t i_item = 0; i_item < vec_item_data.size(); ++i_item )
            {
                // every item can have several meshsets...
                shared_ptr<ItemShapeData>& item_data = vec_item_data[i_item];
                item_data->applyTransformToItem( product_shape->getTransform() );

                // closed meshsets
                std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets = item_data->m_meshsets;
//                    std::wcout << "\t# meshesets: " << vec_item_meshsets.size() << std::endl;
                for( size_t i_meshset = 0; i_meshset < vec_item_meshsets.size(); ++i_meshset )
                {
                    shared_ptr<carve::mesh::MeshSet<3> >& meshset = vec_item_meshsets[i_meshset];
                    std::vector<carve::mesh::Mesh<3>* >& vec_meshes = meshset->meshes;
//                        std::wcout << "\t# meshes: " << vec_meshes.size() << std::endl;

                    for( size_t i_mesh = 0; i_mesh < vec_meshes.size(); ++i_mesh )
                    {
                        carve::mesh::Mesh<3>* mesh = vec_meshes[i_mesh];
                        std::vector<Point> mesh_pts = Get_Vertices_from_Mesh(mesh);
                        d = OBB_Pointset_Simplification(alcc, mesh_pts);
                        Show_volume(alcc,d);

                        if (d != LCC::null_handle)
                        {
//                                Extrude_box_volume(alcc, d, extru);

                            alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                            alcc.info<3>(d).set_label("Furniture_Nef");
                            box_set_of_furni.push_back( d );

                            Convert_3cell_2_Nef3(alcc, d, N);
                            Nef_box_set_of_furni.push_back(N);
                        }

                    }
                }

                // open meshsets
                std::vector<shared_ptr<carve::mesh::MeshSet<3> > >& vec_item_meshsets_op = item_data->m_meshsets_open;
//                    std::wcout << "\t# open meshesets: " << vec_item_meshsets_op.size() << std::endl;
                for( size_t i_meshset = 0; i_meshset < vec_item_meshsets_op.size(); ++i_meshset )
                {
                    shared_ptr<carve::mesh::MeshSet<3> >& meshset = vec_item_meshsets_op[i_meshset];
                    std::vector<carve::mesh::Mesh<3>* >& vec_meshes = meshset->meshes;
//                        std::wcout << "\t# meshes: " << vec_meshes.size() << std::endl;

                    for( size_t i_mesh = 0; i_mesh < vec_meshes.size(); ++i_mesh )
                    {
                        carve::mesh::Mesh<3>* mesh = vec_meshes[i_mesh];
                        std::vector<Point> mesh_pts = Get_Vertices_from_Mesh(mesh);
                        d = OBB_Pointset_Simplification(alcc, mesh_pts);
                        Show_volume(alcc,d);

                        if (d != LCC::null_handle)
                        {
//                                Extrude_box_volume(alcc, d, extru);

                            alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                            alcc.info<3>(d).set_label("Furniture_Nef");
                            box_set_of_furni.push_back( d );

                            Convert_3cell_2_Nef3(alcc, d, N);
                            Nef_box_set_of_furni.push_back(N);
                        }

                    }
                }
            }
        }

        return Nef_box_set_of_furni;


//            int checked = alcc.get_new_mark();
//            CGAL_assertion( checked!=-1 );

//            std::map<double, vec_dart> biggest_vol;
//            for(uint i=0; i<box_set_of_furni.size(); i++)
//            {
//                double vol = Get_volume_magnitude(alcc, box_set_of_furni[i]);
//                biggest_vol[ vol ].push_back( box_set_of_furni[i] );
//            }

//            vec_dart sorted_vols;
//            std::map<double, vec_dart>::iterator it_vol(biggest_vol.begin());
//            std::vector<Nef_polyhedron> sorted_nefs;
//            for(; it_vol!=biggest_vol.end(); it_vol++)
//            {
//                sorted_vols.insert( sorted_vols.end(), it_vol->second.begin(), it_vol->second.end() );

//                for(uint i=0; i<it_vol->second.size(); i++)
//                {
//                    Nef_polyhedron N;
//                    Convert_3cell_2_Nef3(alcc, it_vol->second[i], N);
//                    sorted_nefs.push_back( N );
//                }
//            }


//            assert(sorted_vols.size() > 0);
//            assert(sorted_vols.size() == sorted_nefs.size());

//            Nef_polyhedron N_global = sorted_nefs[0];
//            uint count = 0, total_loop = 0, cur_vol = N_global.number_of_edges();


//            std::cout << "\nTotal number of volumes for union: " << sorted_vols.size() << std::endl;
//            while ( (count < sorted_vols.size() - 1) || (total_loop < sorted_vols.size() - count) )
//            {
//                total_loop++;
//                uint nb_u = 0;
//                for(uint i=1; i<sorted_vols.size(); i++)
//                {
//                    if (!alcc.is_marked(sorted_vols[i], checked))
//                    {
//                        Nef_polyhedron N_cur = sorted_nefs[i];

//                        N_global += N_cur;
//                        uint new_vol = N_global.number_of_edges();

//                        // if there is union
//                        if ( cur_vol < new_vol )
//                        {
//                            std::cout << "\tYes UNION! " << new_vol - cur_vol << std::endl;
//                            alcc.mark(sorted_vols[i], checked);
//                            cur_vol = new_vol;
//                            count++;

//                            nb_u++;
//                        }
//                    }
//                }

//                std::cout << "\nloop? " << total_loop << std::endl;
//                std::cout << "union(s)? " << count << std::endl;

//                if (nb_u == 0)
//                    break;
//            }

//            if ( total_loop > sorted_vols.size() - 1 )
//                std::cout << "More than enouhg loop done... so WHAT???!" << total_loop << std::endl;
//            else if ( count == sorted_vols.size() - 1 )
//                std::cout << "seem WELL DONE!" << count << std::endl;


//            Dart_handle d_ref, d_cur;
//            for(uint i=0; i<box_set_of_furni.size(); i++)
//            {
//                d_ref = box_set_of_furni[i];
//                if (!alcc.is_marked( d_ref, checked ))
//                {
//                    Nef_polyhedron N_ref;
//                    Convert_3cell_2_Nef3(alcc, d_ref, N_ref);
//                    // mark the ref dart as processed
//                    alcc.mark( d_ref, checked );

//                    Bbox_3 bb_ref = LCCtools::Get_Bbox_vol(alcc, d_ref);
//                    d_cur = box_set_of_furni[j];
//                    for( uint j=i+1; j<box_set_of_furni.size(); j++ )
//                    {
//                        if (!alcc.is_marked( d_cur, checked ))
//                        {
//                            Bbox_3 bb_cur = LCCtools::Get_Bbox_vol(alcc, d_cur);
//                            // If the Bboxes intersect, check their union
//                            if (do_overlap(bb_ref, bb_cur))
//                            {
//                                Nef_polyhedron N_cur;
//                                Convert_3cell_2_Nef3(alcc, d_cur, N_cur);

//                                N_ref += N_cur;


//                            }
//                        }
//                    }
//                }
//            }


//            alcc.unmark_all(checked);
//            alcc.free_mark(checked);

//            Polyhedron_3<EK> P_global;
//            N_global.convert_to_Polyhedron( P_global );
//            d = Build_lcc_from_exact_polyhedron_3(alcc, P_global);

//            if (d != LCC::null_handle)
//            {
//                alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                alcc.info<3>(d).set_label(101);
//            }

//            res = N_global;
    }
#endif

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Returns .
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool Get_free_space_in_functional_space (LCC& alcc, Dart_handle& fspace, std::vector< Nef_polyhedron >& nef_furni)
    {
        // Get the Nef of the functional space
        Nef_polyhedron N_fspace;
        Dart_handle d = LCC::null_handle;
        Convert_3cell_2_Nef3(alcc, fspace, N_fspace);

        for(uint i=0; i<nef_furni.size(); i++)
        {
            N_fspace -= nef_furni[i];
        }

        if( N_fspace.is_simple() )
        {
            Polyhedron_EK P_fspace;
            N_fspace.convert_to_Polyhedron( P_fspace );
            d = Build_lcc_from_exact_polyhedron_3(alcc, P_fspace);

            if (d != LCC::null_handle)
            {
                alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                alcc.info<3>(d).set_label("Free_Space_in_Fspace");
//                    return true;
            }
        }
        else
        {
            std::cout << "\nNot simple volume... " << std::endl;
            std::vector<Polyhedron_EK> convex_poly;
            Get_convex_decomposition(N_fspace, convex_poly);

            vec_dart res;
            for(uint i=0; i<convex_poly.size(); i++)
            {
                d = Build_lcc_from_exact_polyhedron_3(alcc, convex_poly[i]);
                if (d != LCC::null_handle)
                {
                    alcc.set_attribute<3>(d, alcc.create_attribute<3>());
//                    alcc.info<3>(d).set_label("Free_Space_in_Fspace");
                }
                res.push_back(d);
            }

            return (res.size() > 0);
        }
        return false;
    }


#ifdef IFCPP_ON
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Import IFC file using the IFC++ library (http://www.ifcplusplus.com, visited on May 2016)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void ifc2lcc(LCC& alcc, std::wstring& in)
    {
        shared_ptr<BuildingModel> ifc_model( new BuildingModel() );
        shared_ptr<ReaderSTEP> reader( new ReaderSTEP() );
//            reader->loadModelFromFile( in, ifc_model );

        shared_ptr<GeometryConverter> geometry_converter(new GeometryConverter(ifc_model));
        std::wcout << "Number of shape before geometry conversion: "
                   << geometry_converter->getShapeInputData().size() << std::endl;

        geometry_converter->clearInputCache();
        reader->loadModelFromFile( in, geometry_converter->getBuildingModel() );
        // convert IFC geometric representations into Carve geometry
        geometry_converter->convertGeometry();
        ifc_model = geometry_converter->getBuildingModel();

        std::wcout << "Number of shape after geometry conversion: "
                   << geometry_converter->getShapeInputData().size() << std::endl;

//            const std::map<int, shared_ptr<BuildingEntity> >& map_ifc_entities = ifc_model->getMapIfcEntities();
//            std::wcout << "Number of entities found: " << map_ifc_entities.size() << std::endl;

        double buffer = 0.75;

        std::vector< shared_ptr<IfcBuildingStorey> > all_storeys;
        if ( Get_Building_Storeys( ifc_model, all_storeys ) )
        {
            for (uint i = all_storeys.size(); i>0; i--)
            {
                uint x = i-1;
                std::wcout << "i = " << x << std::endl;

                /// Load all elements (openings) of the storey
                vec_dart storey_openings, op_res, spaces_res, spaces_no_furni, aggreg_furni;
                std::vector<Nef_polyhedron> n_spaces;
                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > elem_of_storey;
                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >::iterator it_all_elem;

                Get_all_IfcProducts_related_to_IfcBuildingStorey(all_storeys[ x ], elem_of_storey);

                std::vector<const char*> key = {"IfcDoor", "IfcWindow", "IfcWall", "IfcWallStandardCase", "IfcSlab", "IfcCovering" };
                for(uint c=0; c<key.size(); c++)
                {
                    it_all_elem = elem_of_storey.find( key[c] );
                    if (it_all_elem != elem_of_storey.end())
                    {
                        std::cout << "\nGot " << it_all_elem->second.size() << " " << it_all_elem->first << " to go!" << std::endl;
                        for(uint k=0; k<it_all_elem->second.size(); k++)
                        {
                            std::wcout << "\t" << it_all_elem->first << " " << k << std::endl;
                            op_res = Get_Geometry_from_IfcObject(alcc, it_all_elem->second[k], map_ifc_shapes);

                            // store resulting darts in array
                            storey_openings.insert( storey_openings.end(), op_res.begin(), op_res.end() );
                            op_res.clear();
                        }
                    }
                }

                /// Load IfcSpaces of storeys
//                    std::vector< shared_ptr<IfcSpace> > spaces_of_storey;
//                    if ( Get_IfcSpaces_of_Storey( all_storeys[ x ], spaces_of_storey) )
//                    {
////                            std::wcout << "Number of IfcSpaces found in storey " << all_storeys.size() - i << ": " << spaces_of_storey.size() << std::endl;
////                        std::vector<uint> sp{7,8};
////                        for(uint s=0; s<sp.size(); s++)
////                        for (uint j=0; j<4; j++)
//                        for (uint j=0; j<spaces_of_storey.size(); j++)
//                        {
////                            uint j = sp[s];
////                            uint j = 0;

//                            vec_dart space, furnitures_cur, furnitures_all, aggregations;
//                            std::vector<Nef_polyhedron> nef_furni;

//                            std::wcout << "\nSpace " << j << "(" << spaces_of_storey[j]->m_Name->m_value << ")"
//                                       << " that contains " << spaces_of_storey[j]->m_ContainsElements_inverse.size() << std::endl;
//                            space = Get_Geometry_from_IfcObject(alcc, spaces_of_storey[j], geometry_converter, true);

//                            // more than space = ill IfcSpace
//                            if ( space.size() == 1 )
//                            {
//                                /// Get the furnitures (not limited to IfcFurnishingElements) contained in the current IfcSpace
////                                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > furnitures;
////                                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >::iterator it_f;

////                                Get_all_Furnitures_in_IfcSpace (spaces_of_storey[j], furnitures);
////                                for (it_f = furnitures.begin(); it_f!=furnitures.end(); it_f++)
////                                {
////                                    for (uint k=0; k<it_f->second.size(); k++)
////                                    {
////                                        std::wcout << "\n" << it_f->first << " " << k+1 << " ("
////                                                   << it_f->second[k]->m_Name->m_value << ")" << std::endl;
//////                                        Get_Geometry_from_IfcObject(alcc, it_f->second[k], geometry_converter);
////                                        Get_Geometry_from_IfcObject_as_PolySoup(alcc, it_f->second[k], geometry_converter);

////                                        furnitures_cur = Get_Geometry_from_IfcObject(alcc, it_f->second[k], geometry_converter,false,true);
////                                        furnitures_all.insert( furnitures_all.end(), furnitures_cur.begin(), furnitures_cur.end() );

////                                        std::vector<Nef_polyhedron> nef_furni_cur;
////                                        nef_furni_cur = Get_Nef_of_furniture(alcc, it_f->second[k], geometry_converter);

////                                        nef_furni.insert( nef_furni.end(), nef_furni_cur.begin(), nef_furni_cur.end() );
////                                    }
////                                }

//                                /// Create aggregation of furnitures if any and if needed
////                                std::cout << "\tCREATING O-SPACES" << std::endl;
////                                if (furnitures_all.size() > 0)
////                                {
////                                    Aggregate_Set_of_3cells(alcc, furnitures_all, aggregations);
////                                    for(uint k=0; k<aggregations.size(); k++)
////                                        Make_contact_between_volA_and_volB(alcc, space[0], aggregations[k], true);

////                                    if ( Extract_Real_Free_IfcSpace( alcc, space[0], aggregations, spaces_res, n_spaces ) )
////                                    {
////                                        std::cout << "\nFree space successfully extracted!" << std::endl;
//////                                        alcc.remove_cell<3>(space[0]);
////                                    }
////                                    else
////                                        std::cout << "\nThere are more than one space... probaly not a closed IfcSpace!" << std::endl;

////                                    aggreg_furni.insert(aggreg_furni.end(), aggregations.begin(), aggregations.end());
////                                }
////                                else
////                                    spaces_no_furni.push_back( space[0] );


//                                /// Try to get the difference between the fspaces and the furnitures.
////                                for(uint ag=0; ag<aggreg_furni.size(); ag++)
////                                {
////                                    if (Get_free_space_in_functional_space( alcc, aggreg_furni[ag], nef_furni ) )
////                                        std::cout << "\nSimple volume obtain for the fspace subdivision!" << std::endl;
////                                }



//                                /// for the extrusion of the furnitures to force contact with the room borders
////                                    for(uint k=0; k<furnitures_all.size(); k++)
////                                    {
////                                        Make_contact_between_volA_and_volB(alcc, space[0], furnitures_all[k], true);

////                                        // for contact with neighboring furnitures
////                                        for(uint l=k+1; l<furnitures_all.size(); l++)
////                                            Make_contact_between_volA_and_volB(alcc, furnitures_all[k], furnitures_all[l]);
////                                    }

////                                    if ( Extract_Real_Free_IfcSpace( alcc, space[0], furnitures_all, spaces_res ) )
////                                    {
////                                        std::cout << "\nFree space successfully extracted!" << std::endl;
////                                        alcc.remove_cell<3>(space[0]);
////                                    }
////                                    else
////                                        std::cout << "\nThere are more than one space... probaly not a closed IfcSpace!" << std::endl;


//                                /// Load components around IfcSpaces
//                                // test for storey 0
////                                    std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char > elem_around_space;
////                                    std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >::iterator it_elem;

////                                    Get_IfcElements_around_IfcSpace(spaces_of_storey[j], elem_around_space);
////    //                                for(it_elem = elem_around_space.begin(); it_elem != elem_around_space.end(); it_elem++)

////                                    const char* key[] = {"IfcDoor", "IfcWindow"};
////                                    for( int count = 0; count<2; count++ )
////                                    {
////                                        it_elem = elem_around_space.find( key[count] );
////                                        if (it_elem != elem_around_space.end())
////                                        {
////                                            std::cout << "\nGot " << it_elem->second.size() << " " << it_elem->first << " to go!" << std::endl;
////                                            for(uint k=0; k<it_elem->second.size(); k++)
////                                            {
////                                                std::wcout << "\t" << it_elem->first << " " << k << std::endl;
////                                                Get_Geometry_from_IfcObject(alcc, it_elem->second[k], geometry_converter);
////                                            }
////                                        }
////                                    }

////                                    if (furnitures_all.size() > 0)
////                                        break;
//                            }
//                        }


//                        /// Extrude the opening to make their buffer inside the free spaces of the rooms
////                        std::cout << "\tEXTRUDING OPENINGS TO CREATE BUFFERS" << std::endl;
////                        std::cout << "number of spaces to deal with: " << spaces_res.size() << " and " << spaces_no_furni.size() << std::endl;
////                        std::vector<double_pair_dart_vec3d> dir_to_extrude;
////                        for(uint op=0; op<storey_openings.size(); op++)
////                        {
////                            double_pair_dart_vec3d cur_op_extru;
////                            Merge_coplanar_2sewed_2cells_in_volume(alcc, storey_openings[op]);
////                            cur_op_extru = Extrude_opening_in_thickeness_direction( alcc, storey_openings[op], buffer );

////                            // Keep the extrusion direction and faces to get back to proper size afterward
////                            dir_to_extrude.push_back( cur_op_extru );

////                            //force contact to the space boundaries to avoid thin spaces
////                            for(uint sp=0; sp<spaces_res.size(); sp++)
////                                Make_contact_between_volA_and_volB(alcc, spaces_res[sp], storey_openings[op]);
////                            for(uint spn=0; spn<spaces_no_furni.size(); spn++)
////                                Make_contact_between_volA_and_volB(alcc, spaces_no_furni[spn], storey_openings[op]);
////                        }

////                        vec_dart op_fspaces;
////                        for(uint sp=0; sp<spaces_res.size(); sp++)
////                        {
////                            if (Extract_fspace_from_free_Space(alcc, spaces_res[sp], n_spaces[sp], storey_openings, op_fspaces, 10) )
////                                std::cout << "\tI'm gooood bruh!" << std::endl;
////                            else
////                                std::cout << "\tI'm NAH gooood bruh!" << std::endl;
////                        }
////                        for(uint spn=0; spn<spaces_no_furni.size(); spn++)
////                        {
////                            if (Extract_fspace_from_free_Space_without_furnitures(alcc, spaces_no_furni[spn], storey_openings, op_fspaces, 10) )
////                                std::cout << "\tI'm gooood (again) bruh!" << std::endl;
////                            else
////                                std::cout << "\tI'm NAH gooood (again) bruh!" << std::endl;
////                        }

////                        // resize the openings
////                        for(uint op=0; op<dir_to_extrude.size(); op++)
////                        {
////                            // change extrusion direction to opposite
////                            dir_to_extrude[op].first.second = dir_to_extrude[op].first.second * -1;
////                            dir_to_extrude[op].second.second = dir_to_extrude[op].second.second * -1;

////                            Extrude_face_in_vector_direction( alcc, dir_to_extrude[op].first.first, dir_to_extrude[op].first.second, buffer );
////                            Extrude_face_in_vector_direction( alcc, dir_to_extrude[op].second.first, dir_to_extrude[op].second.second, buffer );
////                        }



//                        /// Extrude the object spaces of the furnitures and extract their buffers out of the free space
//                        /// starting from the biggest to the smallest
////                        std::cout << "\tEXTRUDING O-SPACES TO CREATE BUFFERS" << std::endl;
////                        std::map<double, Dart_handle> sorted_furniture_fspaces;
////                        for(uint ag=0; ag<aggreg_furni.size(); ag++)
////                        {
////                            // extruding the object spaces for considering buffers
////                            Extrude_box_volume(alcc, aggreg_furni[ag], 0.5);

////                            //force contact to the space boundaries to avoid thin spaces
////                            for(uint sp=0; sp<spaces_res.size(); sp++)
////                                Make_contact_between_volA_and_volB(alcc, spaces_res[sp], aggreg_furni[ag]);

////                            // sort by biggest volume
////                            sorted_furniture_fspaces[ Get_volume_magnitude(alcc, aggreg_furni[ag]) ] = aggreg_furni[ag];
////                        }

////                        aggreg_furni.clear();
////                        std::map<double, Dart_handle>::reverse_iterator it_ag(sorted_furniture_fspaces.rbegin());
////                        for( ; it_ag!=sorted_furniture_fspaces.rend(); it_ag++ )
////                            aggreg_furni.push_back( it_ag->second );

////                        vec_dart furni_fspaces;
////                        for(uint sp=0; sp<spaces_res.size(); sp++)
////                        {
////                            if (Extract_fspace_from_free_Space(alcc, spaces_res[sp], n_spaces[sp], aggreg_furni, furni_fspaces, 20) )
////                                std::cout << "\tE'rthang alright!" << std::endl;
////                        }

////                        // back to the original size of the functional spaces
////                        for(uint ag=0; ag<aggreg_furni.size(); ag++)
////                            Extrude_box_volume(alcc, aggreg_furni[ag], 0.5, false);








//                        /// GET THE REAL F-SPACE ACTUALLY..... xD
////                        std::cout << "\tEXTRUDING F-SPACES TO CREATE BUFFERS" << std::endl;
////                        sorted_furniture_fspaces.clear();
////                        for(uint ag=0; ag<aggreg_furni.size(); ag++)
////                        {
////                            // extruding the functional spaces for considering buffers
////                            Extrude_box_volume(alcc, aggreg_furni[ag], 1.0);
////                            Z_extrude_upper_lower_box_faces(alcc, aggreg_furni[ag], 1.0);

////                            //force contact to the space boundaries to avoid thin spaces
////                            for(uint sp=0; sp<spaces_res.size(); sp++)
////                                Make_contact_between_volA_and_volB(alcc, spaces_res[sp], aggreg_furni[ag]);

////                            // sort by biggest volume
////                            sorted_furniture_fspaces[ Get_volume_magnitude(alcc, aggreg_furni[ag]) ] = aggreg_furni[ag];
////                        }

////                        aggreg_furni.clear();
////                        std::map<double, Dart_handle>::reverse_iterator it_ag2(sorted_furniture_fspaces.rbegin());
////                        for( ; it_ag2!=sorted_furniture_fspaces.rend(); it_ag2++ )
////                            aggreg_furni.push_back( it_ag2->second );

////                        furni_fspaces.clear();
////                        for(uint sp=0; sp<spaces_res.size(); sp++)
////                        {
////                            if (Extract_fspace_from_free_Space(alcc, spaces_res[sp], n_spaces[sp], aggreg_furni, furni_fspaces, 21) )
////                                std::cout << "\tE'rthang alright!" << std::endl;
////                        }

////                        // back to the original size of the functional spaces
////                        for(uint ag=0; ag<aggreg_furni.size(); ag++)
////                        {
////                            Extrude_box_volume(alcc, aggreg_furni[ag], 1.0, false);
////                            Z_extrude_upper_lower_box_faces(alcc, aggreg_furni[ag], 1.0,true, false, false);
////                        }







//                        // to load just one storey...
//                        break;
//                    }
            }
        }


        /// Try to get rid of the wrong darts...
        typename LCC::Base::One_dart_per_cell_range<3>::iterator it = alcc.one_dart_per_cell<3>().begin(),
                                                                 itend = alcc.one_dart_per_cell<3>().end();
        for(; it!=itend; it++)
        {
            if ( it == NULL || alcc.template attribute<0>(it)==NULL )
            {
                std::cout << "ONE NULL DART FOUND!!!" << std::endl;
//                    remove_cell<LCC,0>(alcc, it);
                alcc.erase_dart(it);
            }
        }

        /// count the 2-free darts
//                int c = 0;
//                vec_dart wrong_dart;
//                for (typename Dart_range::iterator it1(alcc.darts().begin()),
//                     itend(alcc.darts().end()); it1!=itend; ++it1 )
//                {
//                    if(!it1->is_free(2))
//                    {
////                        wrong_dart.push_back(it1);
////                        c++;
//                        alcc.template unsew<2>(it1);
//                    }
//                }
//                std::cout << "\n\t ------------ Number of 2-free dartsc: ------------ " << c
//                          << "(" << (c * 100.0) / alcc.darts().size() << "%)" << std::endl;

////                alcc.remove_cell<3>( wrong_dart[0]);
    }

    ///
    ///
    ///  UI FUNCTIONS
    ///
    ///


    void load_shapedata_of_openings()
    {
        for(auto it=map_ifc_entities.begin(); it!=map_ifc_entities.end(); it++)
        {
            if ( strcmp(it->second->className(), "IfcOpeningElement") == 0 )
            {
                shared_ptr<IfcProduct> product = dynamic_pointer_cast<IfcProduct>(it->second);
                if (product)
                {
                    shared_ptr<ProductShapeData> product_shape ( new ProductShapeData() );
                    shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(product);
                    product_shape->m_ifc_object_definition = obj_def;
                    geomConverter->convertIfcProductShape( product_shape );

                    std::string pGuid;
                    if (product->m_GlobalId){

//                        std::wcout << "product global ID: " << product->m_GlobalId->m_value << std::endl;
//                        std::wcout << "object def global ID: " << obj_def->m_GlobalId->m_value << std::endl;

                        std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converterX;
                        pGuid = converterX.to_bytes(product->m_GlobalId->m_value);
                        map_ifc_shapes[ pGuid ] = product_shape;
                    }
                    else
                        std::wcout << "the product has NO GUID!!!!" << std::endl;
                }
            }
            else if ( strcmp(it->second->className(), "IfcDoor") == 0 ){
                shared_ptr<IfcDoor> door = dynamic_pointer_cast<IfcDoor>(it->second);
                if(door){
                    std::map< const char*,std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char > elemAroundDoor;
                    Get_IfcElements_around_IfcDoor(door, elemAroundDoor);
                }
            }
        }
    }


    void loadShapedataOfSelectedEntities()
    {
        for(auto it=map_ifc_entities.begin(); it!=map_ifc_entities.end(); it++)
        {
            if ( strcmp(it->second->className(), "IfcOpeningElement") == 0 )
            {
                shared_ptr<IfcProduct> product = dynamic_pointer_cast<IfcProduct>(it->second);
                if (product)
                {
                    shared_ptr<ProductShapeData> product_shape ( new ProductShapeData() );
                    shared_ptr<IfcObjectDefinition> obj_def = dynamic_pointer_cast<IfcObjectDefinition>(product);
                    product_shape->m_ifc_object_definition = obj_def;
                    geomConverter->convertIfcProductShape( product_shape );

                    std::string pGuid;
                    if (product->m_GlobalId){

//                        std::wcout << "product global ID: " << product->m_GlobalId->m_value << std::endl;
//                        std::wcout << "object def global ID: " << obj_def->m_GlobalId->m_value << std::endl;

                        std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converterX;
                        pGuid = converterX.to_bytes(product->m_GlobalId->m_value);
                        map_ifc_shapes[ pGuid ] = product_shape;
                    }
                    else
                        std::wcout << "the product has NO GUID!!!!" << std::endl;
                }
            }
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Pre-load the entities of an IFC file using the IFC++ library
    /// This function will fill the global containers on which the rest of the code will rely
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    bool pre_load_ifc(std::wstring& in)
    {
        shared_ptr<BuildingModel> ifc_model( new BuildingModel() );
        shared_ptr<ReaderSTEP> reader( new ReaderSTEP() );

//            reader->loadModelFromFile( in, ifc_model );

        shared_ptr<GeometryConverter> geometry_converter(new GeometryConverter(ifc_model));
//            reader->loadModelFromFile( in, ifc_model );

        reader->loadModelFromFile( in, geometry_converter->getBuildingModel() );

        std::wcout << "Number of shapes before geometry conversion: "
                   << geometry_converter->getShapeInputData().size() << std::endl;

        // convert IFC geometric representations into Carve geometry
//        geometry_converter->convertGeometry();
        ifc_model = geometry_converter->getBuildingModel();

        std::wcout << "Number of shapes after geometry conversion: "
                   << geometry_converter->getShapeInputData().size() << std::endl;

        map_ifc_shapes = geometry_converter->getShapeInputData();
        map_ifc_entities = ifc_model->getMapIfcEntities();
        std::wcout << "Number of entities found (preloader): " << map_ifc_entities.size() << std::endl;

        geomConverter = geometry_converter;

        // Load shapedata of the opening elements
        load_shapedata_of_openings();

        // Get all the storeys
        std::vector< shared_ptr<IfcBuildingStorey> > all_storeys;
        if ( Get_Building_Storeys( ifc_model, all_storeys ) )
        {
            ifc_storey_list = all_storeys;

            std::vector< shared_ptr<IfcSpace> > cur_storey_spaces;
            for(uint i=0; i<ifc_storey_list.size(); i++)
            {
                Get_IfcSpaces_of_Storey(ifc_storey_list[i], cur_storey_spaces);
                spaces_of_storeys.push_back( cur_storey_spaces );
                cur_storey_spaces.clear();
            }
        }
        else
        {
            geometry_converter->clearInputCache();
            return false;
        }

        geometry_converter->clearInputCache();
        return true;
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Import IFC file using the IFC++ library (http://www.ifcplusplus.com, visited on May 2016)
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void fully_load_ifc2lcc(LCC& alcc, std::map<int, shared_ptr<BuildingEntity> >& map_ifc_entities,
                            bool simplify, bool no_sew2, bool only_closed_meshes, bool correct_wrong_spaces)
    {
        std::wcout << "Number of entities found: " << map_ifc_entities.size() << std::endl;

        for( auto it = map_ifc_entities.begin(); it != map_ifc_entities.end(); ++it )
        {
            const shared_ptr<BuildingEntity>& ifcpp_entity = it->second;
            shared_ptr<IfcProduct> ifc_product = dynamic_pointer_cast<IfcProduct>( ifcpp_entity );

            if (ifc_product /*&& strcmp(ifc_product->className(), "IfcBuilding") == 0*/ )
            {
//                    __try
                {
                    Get_Geometry_from_IfcObject(alcc, ifc_product, map_ifc_shapes, simplify, no_sew2, only_closed_meshes, correct_wrong_spaces);
                }
//                    __except(GetExceptionCode())
//                    {
//                        std::cout << "******************** Something went wrong on this building part..." << std::endl;
//                    }
            }
        }
    }

    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Load all element of a given storey of an IFC file
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void load_all_in_storey_ifc2lcc(LCC& alcc, shared_ptr<IfcBuildingStorey>& strorey,
                                    bool simplify, bool no_sew2, bool only_closed_meshes, bool correct_wrong_spaces)
    {
        std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > elem_of_storey;
        std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >::iterator it_all_elem;

        Get_all_IfcProducts_related_to_IfcBuildingStorey(strorey, elem_of_storey);

        for( it_all_elem = elem_of_storey.begin(); it_all_elem != elem_of_storey.end(); it_all_elem++)
        {
            for(uint i=0; i<it_all_elem->second.size(); i++)
            {
                Get_Geometry_from_IfcObject(alcc, it_all_elem->second[i], map_ifc_shapes, simplify, no_sew2, only_closed_meshes, correct_wrong_spaces);
            }
        }
    }


    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    /// Load selected elements of a given storey of an IFC file
    /* /////////////////////////////////////////////////////////////////////////////////////////////////////// */
    void load_selected_in_storey_ifc2lcc(LCC& alcc, shared_ptr<IfcBuildingStorey>& strorey,
                                         std::vector<const char*>& selected_elem, bool simplify,
                                         bool no_sew2, bool only_closed_meshes, bool correct_wrong_spaces)
    {
        std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > elem_of_storey;
        std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >::iterator it_all_elem;

        Get_all_IfcProducts_related_to_IfcBuildingStorey(strorey, elem_of_storey);

        for(uint i=0; i<selected_elem.size(); i++)
        {
            it_all_elem = elem_of_storey.find( selected_elem[i] );
            if (it_all_elem != elem_of_storey.end())
            {
                std::cout << "\nGot " << it_all_elem->second.size() << " " << it_all_elem->first << " to go!" << std::endl;
                for(uint j=0; j<it_all_elem->second.size(); j++)
                {
                    Get_Geometry_from_IfcObject(alcc, it_all_elem->second[j], map_ifc_shapes, simplify, no_sew2, only_closed_meshes, correct_wrong_spaces);
                }
            }
        }
    }
#endif
































//    void my_load_off(LCC& alcc, std::istream& in)
//    {
//        std::cout << "Loadind the OFF file: " << &in << std::endl;

//      File_header_OFF  m_file_header;
//      File_scanner_OFF scanner( in, m_file_header.verbose());
//      if (!in) return;
//      m_file_header = scanner;  // Remember file header after return.

//      std::cout << "Number of vertices: " << scanner.size_of_vertices() << std::endl;
//      std::cout << "Number of faces: " << scanner.size_of_facets() << std::endl;
//      std::cout << "Number of halfedges: " << scanner.size_of_halfedges() << std::endl;

//      Linear_cell_complex_incremental_builder_3<LCC> B(alcc);
//      B.begin_surface(scanner.size_of_vertices(),
//                      scanner.size_of_facets(),
//                      scanner.size_of_halfedges());

//      typedef typename LCC::Point Point;

//      // read in all vertices
//      std::size_t  i;
//      for (i = 0; i < scanner.size_of_vertices(); i++)
//      {
//        Point p;
//        file_scan_vertex(scanner, p);
//        B.add_vertex(p);
//        scanner.skip_to_next_vertex(i);
//      }
//      /* TODO rollback
//         if ( ! in  || B.error()) {
//         B.rollback();
//         in.clear( std::ios::badbit);
//         return;
//         }
//      */

//      // read in all facets
//      for (i=0; i<scanner.size_of_facets(); i++)
//      {
//        B.begin_facet();
//        std::size_t no;
//        scanner.scan_facet(no, i);
//        /* TODO manage errors
//           if( ! in || B.error() || no < 3) {
//           if ( scanner.verbose()) {
//           std::cerr << " " << std::endl;
//           std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
//           std::cerr << "operator()(): input error: facet " << i
//           << " has less than 3 vertices." << std::endl;
//           }
//           B.rollback();
//           in.clear( std::ios::badbit);
//           return;
//           } */
//        for (std::size_t j=0; j<no; j++)
//        {
//          std::size_t index;
//          scanner.scan_facet_vertex_index(index, i);
//          B.add_vertex_to_facet(index);
//        }
//        B.end_facet();
//        scanner.skip_to_next_facet(i);
//      }
//      /* TODO manage errors
//         if ( ! in  || B.error()) {
//         B.rollback();
//         in.clear( std::ios::badbit);
//         return;
//         }
//         if ( B.check_unconnected_vertices()) {
//         if ( ! B.remove_unconnected_vertices()) {
//         if ( scanner.verbose()) {
//         std::cerr << " " << std::endl;
//         std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
//         std::cerr << "operator()(): input error: cannot "
//         "succesfully remove isolated vertices."
//         << std::endl;
//         }
//         B.rollback();
//         in.clear( std::ios::badbit);
//         return;
//         }
//         }*/
//      B.end_surface();
//    }



}
