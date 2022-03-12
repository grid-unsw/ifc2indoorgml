#include "LCC_SpecialOps.h"


double eps2 = 0.001,
       eps_angle = 5.0,
       coplanarity = 0.01;
Point global_shift_pt (0.0, 0.0, 0.0);


/* ///////////////////////////////////////////// */
/// normalize a 3D vector such that ||v|| = 1
/* ///////////////////////////////////////////// */
void LCCtools::unit_normal(Vector& v)
{
    double s_length = v.squared_length();
    if(s_length != 0.0)
    {
        Vector n = v*(1/sqrt(s_length));

        double x = n.x(),
               y = n.y(),
               z = n.z();

        if (abs(x) < 0.0001)
            x = 0.0;
        if (abs(y) < 0.0001)
            y = 0.0;
        if (abs(z) < 0.0001)
            z = 0.0;

        v = Vector(x,y,z);
    }
}

/* /////////////////////////////////////////////////////////////////////////////////////// */
/// Return the (unit) normal of a 3D polygon using Newell's method (works for concave poly)
/* /////////////////////////////////////////////////////////////////////////////////////// */
Vector LCCtools::poly_normal (const LCC& alcc, Dart_handle d, bool normalize)
{
    double Nx = 0.0, Ny = 0.0, Nz = 0.0;
    LCC::Point pt1, pt2;
    unsigned int nb_seg = 0;

    typename LCC::Base::Dart_of_orbit_range<1>::const_iterator it1 (alcc, d), it2 = it1;
    for (it2++; it1.cont(); it1++, it2++)
    {
        if (it2 == alcc.darts_of_orbit<1>(d).end())
            it2 = alcc.darts_of_orbit<1>(d).begin();

        pt1 = alcc.point(it1);
        pt2 = alcc.point(it2);

        if (pt1 != pt2)
        {
            Nx += ( pt1.y() - pt2.y() ) * ( pt1.z() + pt2.z() );
            Ny += ( pt1.z() - pt2.z() ) * ( pt1.x() + pt2.x() );
            Nz += ( pt1.x() - pt2.x() ) * ( pt1.y() + pt2.y() );

            nb_seg++;
        }
    }
    Vector N;

    if ( nb_seg < 3 ){
//            std::cout << "\nBad face! nb of segments: " << nb_seg <<std::endl;
//            Print_face(alcc, d);
        N = Vector(0.0, 0.0, 0.0);
    }
    else{
        N = Vector(Nx/2.0, Ny/2.0, Nz/2.0);

        if (normalize)
            unit_normal(N);
    }

    return N;
}

/* ////////////////////////////////////////////////////////// */
/// Compute the AABB of a 3-cell
/* ////////////////////////////////////////////////////////// */
Bbox_3 LCCtools::Get_Bbox_vol (LCC& alcc, Dart_handle d)
{
    Bbox_3 bb = alcc.point(d).bbox();

    typename LCC::Base::One_dart_per_incident_cell_range<0,3>::iterator
            itr = alcc.one_dart_per_incident_cell<0,3>(d).begin(),
            itrend = alcc.one_dart_per_incident_cell<0,3>(d).end();
    for (itr++; itr != itrend; itr++)
    {
        bb = bb + alcc.point(itr).bbox() ;
    }

    return bb;
}

/* ///////////////////////////////////////////// */
/// Generate a unique ID based on current time
/* ///////////////////////////////////////////// */
std::string LCCtools::generate_unique_ID()
{
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);

    std::string res;
    res += std::to_string(now->tm_mday);
    res += std::to_string(now->tm_mon);
    res += std::to_string(now->tm_year + 1900);
    res += std::to_string(now->tm_hour);
    res += std::to_string(now->tm_min);
    res += std::to_string(now->tm_sec);

    CGAL::Random myrand;
    int final = std::stod( res ) / myrand.get_double(1000000000, 3000000000);
    return std::to_string(final);
}



/* ////////////////////////////////////////////////////////////////////////////// */
/// Return the angle between two 2D vectors v1 and v2, in the range [0°, 180°]
/* ////////////////////////////////////////////////////////////////////////////// */
double LCCtools::compute_3d_angle(const Vector& v1, const Vector& v2)
{
    double a = CGAL::to_double( (v1*v2) / ( sqrt(v1.squared_length()) * sqrt(v2.squared_length()) ) ) ;

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
bool LCCtools::vectors_are_eps_colinear (Vector v1, Vector v2)
{
    double angle = compute_3d_angle (v1, v2);
    if ( (angle < (0.0 + eps_angle) || angle > (360.0 - eps_angle) ) ||
         (angle < (180.0 + eps_angle) && angle > (180.0 - eps_angle) ) )
        return true;
    return false;
}


/* /////////////////////////////////////////////////////////////////////////////////////// */
/// Return the (unit) normal of a 3D polygon using Newell's method (works for concave poly)
/* /////////////////////////////////////////////////////////////////////////////////////// */
Vector LCCtools::poly_normal_vec3d (vec_pt3d& pts, bool normalize)
{
    double Nx = 0.0, Ny = 0.0, Nz = 0.0;
    Point pt1, pt2;

    for (uint i=0, j=1; i<pts.size(); i++, j++)
    {
        if (i == pts.size() - 1)
            j = 0;

        pt1 = pts[i];
        pt2 = pts[j];

        Nx += ( pt1.y() - pt2.y() ) * ( pt1.z() + pt2.z() );
        Ny += ( pt1.z() - pt2.z() ) * ( pt1.x() + pt2.x() );
        Nz += ( pt1.x() - pt2.x() ) * ( pt1.y() + pt2.y() );
    }

    Vector N(Nx/2.0, Ny/2.0, Nz/2.0);

    if (normalize)
        unit_normal(N);

    return N;
}

/* /////////////////////////////////////////////////////////////// */
/// Compute the surface area of a LCC 3D polygon using Newell's normal
/* /////////////////////////////////////////////////////////////// */
double LCCtools::poly_surface_area(LCC& alcc, Dart_handle d)
{
    // the area correspond to the half of the length of the un-normalized Newell's normal
    return ( sqrt( CGAL::to_double( poly_normal(alcc, d, false).squared_length() ) ) / (2.0 /** pow(scale_f, 2)*/ ) );
}

/* /////////////////////////////////////////////////////////////// */
/// Compute the surface area of a LCC 3D polygon using Newell's normal
/* /////////////////////////////////////////////////////////////// */
double LCCtools::poly_surface_area_vec3d(vec_pt3d& pts)
{
    // the area correspond to the half of the length of the un-normalized Newell's normal
    return ( sqrt( CGAL::to_double( poly_normal_vec3d(pts, false).squared_length() ) ) / (2.0 /** pow(scale_f, 2)*/ ) );
}


/* /////////////////////////////////////////////////////// */
/// Print all the vertices of the given face
/* /////////////////////////////////////////////////////// */
void LCCtools::Print_face(const LCC& alcc, Dart_handle d)
{
//        typename LCC::Base::Dart_of_cell_range<2>::iterator itr( alcc, d );
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

/* /////////////////////////////////////////////////////// */
/// Print all the vertices of the given face (vec_pt3d)
/* /////////////////////////////////////////////////////// */
void LCCtools::Print_face_vec(vec_pt3d &poly)
{
    for ( uint i=0; i<poly.size(); i++ )
        std::cout << poly[i] << std::endl;
    std::cout << poly[0] << std::endl;

    std::cout << "number of points: " << poly.size() << std::endl;
}

/* ////////////////////////////////////////////////////////////////////////////// */
/// Returns true if two 3D points are epsilon equal
/* ////////////////////////////////////////////////////////////////////////////// */
bool LCCtools::Points_are_eps_close (const Point& p1, const Point& p2)
{
    double x = fabs ( CGAL::to_double( (p1.x() - p2.x() ) ) ),
           y = fabs ( CGAL::to_double( (p1.y() - p2.y() ) ) ),
           z = fabs ( CGAL::to_double( (p1.z() - p2.z() ) ) );

//        std::cout << "dx = " << x << " \tdy = " << y << " \tdz = " << z << std::endl;

    if (x<=eps2 && y<=eps2 && z<=eps2)
//        if (x == 0.0 && y == 0.0 && z == 0.0)
        return true;
    return false;
}

/* ////////////////////////////////////////////////////////////////////////////////////////// */
/// Collect for each 1-cell of the LCC all the faces that are around and incident to it
/* ////////////////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Collect_faces_around_edges_dumb (LCC& alcc, std::vector< vec_dart >& edges, vec_dart& face_subset)
{
    size_t checked = alcc.get_new_mark();
    CGAL_assertion( checked !=-1 );
    vec_dart edge;

    // If a subset of faces is defined, limit the collection to them
    if (face_subset.size() > 0)
    {
//            std::cout << "\nNUMBER OF FACES TO PROCESS: " << face_subset.size() <<std::endl;
        vec_dart edges_of_faces;
        for(uint i=0; i<face_subset.size(); i++)
        {
//                Print_face(alcc, face_subset[i]);
            typename LCC::Base::One_dart_per_incident_cell_range<0,2>::iterator
                    itr = alcc.one_dart_per_incident_cell<0,2>(face_subset[i]).begin(),
                    itrend = alcc.one_dart_per_incident_cell<0,2>(face_subset[i]).end();
            for (; itr != itrend; itr++)
            {
                edges_of_faces.push_back( itr );
            }
        }

        for(uint i=0; i<edges_of_faces.size(); i++)
        {
            Dart_handle it = edges_of_faces[i];
            if ( !alcc.is_marked(it, checked) )
            {
                alcc.mark(it, checked);
                edge.push_back( it );
                Point pt1 = alcc.point( it ),
                      pt2 = alcc.point( alcc.beta<1>(it) );

                for(uint j=i+1; j<edges_of_faces.size(); j++)
                {
                    Dart_handle itp = edges_of_faces[j];
                    if ( !alcc.is_marked(itp, checked) )
                    {
                        Point ptp1 = alcc.point( itp ),
                              ptp2 = alcc.point( alcc.beta<1>(itp) );

                        if ( (Points_are_eps_close( pt1, ptp1 ) && Points_are_eps_close( pt2, ptp2 ))
                             || (Points_are_eps_close( pt1, ptp2 ) && Points_are_eps_close( pt2, ptp1 )) )
                        {
                            edge.push_back( itp );
                            alcc.mark(itp, checked);
                        }
                    }
                }
//                CGAL_assertion( edge.size() % 2 == 0 );
                edges.push_back( edge );
                edge.clear();
            }
        }
    }
    // Otherwise run through the enire LCC
    else
    {
        uint count = 0;
//            std::cout << "\n#Darts to be processed: " << alcc.darts().size() <<std::endl;
        typename LCC::Dart_range::iterator it(alcc.darts().begin()),
                itend(alcc.darts().end());
        for (; it!=itend; ++it )
        {
            if ( !alcc.is_marked(it, checked) )
            {
                count++;
                alcc.mark(it, checked);
                edge.push_back( it );
                Point pt1 = alcc.point( it ),
                      pt2 = alcc.point( alcc.beta<1>(it) );

                typename LCC::Dart_range::iterator itp(it),
                        itendp(alcc.darts().end());
                for (itp++; itp!=itendp; ++itp )
                {
                    if ( !alcc.is_marked(itp, checked) )
                    {
                        Point ptp1 = alcc.point( itp ),
                              ptp2 = alcc.point( alcc.beta<1>(itp) );

                        if ( (Points_are_eps_close( pt1, ptp1 ) && Points_are_eps_close( pt2, ptp2 ))
                             || (Points_are_eps_close( pt1, ptp2 ) && Points_are_eps_close( pt2, ptp1 )) )
                        {
                            edge.push_back( itp );
                            alcc.mark(itp, checked);
                        }
                    }
                }
                //                CGAL_assertion( edge.size() % 2 == 0 );
//                    if (edge.size() < 2)
//                        std::cout << "\tDANGLING EDGE ??? --> " << edge.size() <<std::endl;
//                    else if (edge.size() > 2)
//                        std::cout << "\tNON-MANIFOLD ZONE --> " << edge.size() << "[" << pt1 << "][" << pt2 << "]" <<std::endl;


                edges.push_back( edge );
                edge.clear();
            }
        }

//            std::cout << "#Darts processed: " << count <<std::endl;
    }

    alcc.unmark_all(checked);
    alcc.free_mark(checked);
}

/* ////////////////////////////////////////////////////////////////////////////////////////// */
/// Overloaded version of the above function
/* ////////////////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Collect_faces_around_edges_dumb (LCC& alcc, std::vector< vec_dart >& edges)
{
    vec_dart face_subset;
    Collect_faces_around_edges_dumb(alcc, edges, face_subset);
}


/* ////////////////////////////////////////// */
/// Counts the number of 2-free darts
/* ////////////////////////////////////////// */
bool LCCtools::There_are_2free_darts(LCC& alcc, vec_dart& set_of_2cells_to_process)
{
    int checked = alcc.get_new_mark();
    CGAL_assertion( checked !=-1 );

    uint count = 0, all = 0;

    // Process just the set of faces if defined
    if (set_of_2cells_to_process.size() > 0)
    {
        for (uint i=0; i<set_of_2cells_to_process.size(); i++)
        {
            Dart_handle d = set_of_2cells_to_process[i];
            // For each dart of the 2-cell
            typename LCC::Base::One_dart_per_incident_cell_range<0,2>::iterator
                    it1 = alcc.one_dart_per_incident_cell<0,2>(d).begin(),
                    it1end = alcc.one_dart_per_incident_cell<0,2>(d).end();
            for (; it1 != it1end; it1++)
            {
                if(alcc.is_free(it1, 2))
                {
//                        weird_edges.push_back(it1);
                    count++;
                }
            }
            all += alcc.one_dart_per_incident_cell<0,2>(d).size();
        }
    }
    // Otherwise process the entire LCC
    else
    {
        for (typename LCC::Dart_range::iterator it1(alcc.darts().begin()),
             itend(alcc.darts().end()); it1!=itend; ++it1 )
        {
            if(alcc.is_free(it1, 2))
            {
//                    weird_edges.push_back(it1);
                count++;
            }
        }

        all = alcc.darts().size();
    }

    std::cout << "\n\t ------------ Number of 2-free darts: ------------ " << count << "/"
              << all << "(" << (count * 100.0) / all << "%)" << std::endl;

    alcc.unmark_all(checked);
    alcc.free_mark(checked);

    return (count != 0);
}


/* ////////////////////////////////////////// */
/// Overloaded version of the above function
/* ////////////////////////////////////////// */
bool LCCtools::There_are_2free_darts(LCC& alcc)
{
    vec_dart empty_set_of_2cells_to_process;
    return There_are_2free_darts(alcc, empty_set_of_2cells_to_process);
}



/* ////////////////////////////////////////////////////////// */
/// Create a new face in the LCC based on a vector of points
/* ////////////////////////////////////////////////////////// */
Dart_handle LCCtools::Insert_new_2cell (LCC& alcc, vec_pt3d& one_face, bool orientation)
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

/* ////////////////////////////////////////////////////////// */
/// Generates a unique string key from a Point
/* ////////////////////////////////////////////////////////// */
std::string LCCtools::get_unique_key_pt( Point& pt )
{
    return  std::string ( std::to_string( int(pt.x()/eps2) )
                        + std::to_string( int(pt.y()/eps2) )
                        + std::to_string( int(pt.z()/eps2) ));
}

/* ////////////////////////////////////////////////////////// */
/// Computes the AABB of a list of faces
/* ////////////////////////////////////////////////////////// */
Bbox_3 LCCtools::Get_Bbox_face_list (LCC& alcc, vec_dart &face_list)
{
    Bbox_3 bb = alcc.point(face_list[0]).bbox();
    for(uint i=0; i<face_list.size(); i++)
    {
        typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, face_list[i]);
        for (; itr.cont(); itr++)
            bb += alcc.point(itr).bbox();
    }

    return bb;
}

/* ////////////////////////////////////////////////////////// */
/// Computes the AABB of a 3D face
/* ////////////////////////////////////////////////////////// */
Bbox_3 LCCtools::Get_Bbox_face (LCC& alcc, Dart_handle d)
{
    Bbox_3 bb = alcc.point(d).bbox();

    typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, d);
    for (itr++; itr.cont(); itr++)
    {
        bb = bb + alcc.point(itr).bbox();
    }

    return bb;
}

/* ////////////////////////////////////////////////////////// */
/// Computes the AABB of a 3D LCC
/* ////////////////////////////////////////////////////////// */
Bbox_3 LCCtools::Get_Bbox_LCC (LCC& alcc)
{
    Bbox_3 bb;
    typename LCC::Base::One_dart_per_cell_range<3>::iterator it = alcc.one_dart_per_cell<3>().begin(),
                                                             itend = alcc.one_dart_per_cell<3>().end();
    if (it != LCC::null_handle)
        bb = alcc.point(it).bbox();

    for(; it!=itend; it++)
    {
        bb = bb +  Get_Bbox_vol(alcc, it);
    }

    return bb;
}

/* ////////////////////////////////////////////////////////// */
/// Returns true if the Bbox is degenerate (flat volume)
/* ////////////////////////////////////////////////////////// */
bool LCCtools::Bbox_is_degenerate(Bbox_3 bb)
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
void LCCtools::Rescale_Bbox (Bbox_3 &bbox, double &scale, bool increase)
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





/* ////////////////////////////////////////////////////////////////////////////////////////// */
/// Reconstructs a bunch of unstrucured polygons into closed volume(s)
/// This method assumes that all the faces are well oriented
/* ////////////////////////////////////////////////////////////////////////////////////////// */
bool LCCtools::Perform_Simple_Volume_Reconstruction_from_Polygon_Soup(LCC& alcc, vec_dart& face_subset)
{
    if ( face_subset.size() > 3 )
    {
//            std::cout << "\nHere we go for a simple volume reconstruction: #faces = " << face_subset.size() <<std::endl;
        //        for(uint i = 0; i<face_subset.size(); i++)
        //            Print_face(alcc, face_subset[i]);

        std::vector <vec_dart> edges;
        // If face_subset is empty, all the edges of the LCC will be considered in the following function
        Collect_faces_around_edges_dumb(alcc, edges, face_subset);

//            std::cout << "\t\tNUMBER OF EDGES: " << edges.size() <<std::endl;

        uint pairs = 0, groups = 0, sewn1 = 0, sewn2 = 0;
        for ( uint i=0; i<edges.size(); i++ )
        {
            //            std::cout << "\tEdge: " << i << " (" << edges[i].size() << " faces)" <<std::endl;
            // If there are at most two 2-cells around the edge
            if ( edges[i].size() == 2 )
            {
                pairs+=2;
                if (alcc.is_free<2>(edges[i][0]) && alcc.is_free<2>(edges[i][1])
                        && alcc.is_sewable<2>(edges[i][0], edges[i][1]) )
                {
                    alcc.template sew<2>(edges[i][0], edges[i][1]);
                    sewn1+=2;
                }
            }

            // If there are more than two (but an even number of) 2-cells around the edge
//                else if ( edges[i].size() > 2 && edges[i].size() % 2 == 0 )
//                {
////                    std::cout << "Non-manifold zone here ======> " << edges[i].size() <<std::endl;
//                    groups += edges[i].size();
//                    for(uint j=0; j<edges[i].size(); j++)
//                        Print_face(alcc, edges[i][j]);
//                    //                sewn2 += angular_sorting_sew2_around_one_edge(alcc, edges[i]);
//                }

//                else if ( edges[i].size() < 2 || edges[i].size() > 2)
//                {
//                    std::cout << "Wrong #faces around edge ======> " << edges[i].size() <<std::endl;
//                    //                weird_edges.push_back(edges[i][0]);
//                    std::cout << "[" <<  alcc.point(edges[i][0]) << "] [" << alcc.point( alcc.beta<1>(edges[i][0]) ) << "]" << std::endl;
//                }

            //            if ( edges[i].size() > 8 )
            //                for (uint j=0; j<edges[i].size(); j++)
            //                {
            //                    std::cout << "\nFace " << j+1 <<std::endl;
            //                    Print_face(alcc, edges[i][j] );
            //                }
            //            c++;
        }

//            std::cout << "\nNUMBER OF ESTIMATED 2-SEWN DARTS PAIR: " << sewn1 << " - expected: " << pairs << std::endl;
//            std::cout << "                          DARTS GROUP: " << sewn2 << " - expected: " << groups << std::endl;
//            std::cout << "FINAL ESTIMATION OF 2-SEWN DARTS: " << sewn1+sewn2 << std::endl;

        return !There_are_2free_darts(alcc, face_subset);
    }
    else
        return false;
}

/* ////////////////////////////////////////////////////////////////////////////////////////// */
/// Overloaded version of the function above
/* ////////////////////////////////////////////////////////////////////////////////////////// */
bool LCCtools::Perform_Simple_Volume_Reconstruction_from_Polygon_Soup(LCC& alcc)
{
    vec_dart face_subset;
    return Perform_Simple_Volume_Reconstruction_from_Polygon_Soup(alcc, face_subset);
}

/* ///////////////////////////////////////////////////////////////////////// */
/// Tries to remove aligned points on edges by keeping only corner points
/* ///////////////////////////////////////////////////////////////////////// */
void LCCtools::Keep_corners_points_only (LCC& alcc, vec_dart& to_remove)
{
    typename LCC::Base::One_dart_per_cell_range<0>::iterator
            it = alcc.one_dart_per_cell<0>().begin(),
            itend = alcc.one_dart_per_cell<0>().end();
    for(; it!=itend; it++)
    {
        if (alcc.one_dart_per_incident_cell<2,0>(it).size() <= 2)
        {
            to_remove.push_back(it);
        }
    }
}

/* //////////////////////////////////////////////////////////////////////////////// */
/// Remove a selected set of 0-cells from the LCC
/* //////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Remove_selected_0_cells (LCC& alcc, vec_dart& darts)
{
    alcc.set_update_attributes(false);
    for(uint i=0; i<darts.size(); i++)
    {
        if (alcc.is_dart_used(darts[i])
                && alcc.is_removable<0>(darts[i]))
        {
//                std::cout << "YES! DART IN USE! (0-removal) " << i+1 << std::endl;
            alcc.remove_cell<0>( darts[i]);
        }
//            else
//                std::cout << "NO!! DART UNUSED! (0-removal)" << i+1 << std::endl;
    }
    alcc.set_update_attributes(true);
}

/* //////////////////////////////////////////////////////////////////////////////// */
/// Remove a selected set of 2-cells from the LCC
/* //////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Remove_selected_2_cells (LCC& alcc, vec_dart& darts)
{
    alcc.set_update_attributes(false);
    for(uint i=0; i<darts.size(); i++)
    {
        if (alcc.is_dart_used(darts[i])
                && alcc.is_removable<2>(darts[i]))
        {
//                std::cout << "YES! DART IN USE! (2-removal)" << std::endl;
            alcc.remove_cell<2>( darts[i]);
        }
//            else
//                std::cout << "NO!! DART UNUSED! (2-removal)" << std::endl;
    }
    alcc.set_update_attributes(true);
}

/* //////////////////////////////////////////////////////////////////////////////// */
/// Remove a selected set of 3-cells from the LCC
/* //////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Remove_selected_3_cells (LCC& alcc, vec_dart& darts)
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


/* ///////////////////////////////////////////////////////////////////////////////////////////// */
/// Shift the min corner of the bbox of the model to the origin if the coordinates are too big
/* ///////////////////////////////////////////////////////////////////////////////////////////// */
void LCCtools::Apply_global_shifting (LCC& alcc)
{
    if ( global_shift_pt == Point(0.0, 0.0, 0.0) )
    {
        double x_max, y_max, z_max,
                x_shift, y_shift, z_shift;

        Bbox_3 bbox = Get_Bbox_LCC(alcc);
        bool shift = false;
        std::cout << "\nMinimum point: ("<< bbox.xmin()
                  << ", " << bbox.ymin()
                  << ", " << bbox.zmin() << std::endl;
        std::cout << "\nMaximum point: ("<< bbox.xmax()
                  << ", " << bbox.ymax()
                  << ", " << bbox.zmax() << std::endl;

        // Find the real extremes (sign independent)
        x_max = std::max ( fabs(bbox.xmin()), fabs(bbox.xmax()) );
        y_max = std::max ( fabs(bbox.ymin()), fabs(bbox.ymax()) );
        z_max = std::max ( fabs(bbox.zmin()), fabs(bbox.zmax()) );

        // Shifting the centroid of the model
        //        x_shift = (bbox.xmax() + bbox.xmin())/2;
        //        y_shift = (bbox.ymax() + bbox.ymin())/2;
        //        z_shift = (bbox.zmax() + bbox.zmin())/2;

        // Shifting the min corner of the model
        x_shift = bbox.xmin();
        y_shift = bbox.ymin();
        z_shift = bbox.zmin();

        if( x_max > 1000.0
                || y_max > 1000.0
                || z_max > 1000.0 )
            shift = true;

        if (shift)
        {
            std::cout << "\nGlobal shifting required! ---> "
                      << x_shift << ", "
                      << y_shift << ", "
                      << z_shift << std::endl;

            typename LCC::Base::One_dart_per_cell_range<0>::iterator
                    it = alcc.one_dart_per_cell<0>().begin(),
                    itend = alcc.one_dart_per_cell<0>().end();
            for(; it!=itend; it++)
            {
                Point sh = alcc.point(it);

                alcc.point(it) = Point( (sh.x() - x_shift)/* / size_scale*/,
                                        (sh.y() - y_shift)/* / size_scale*/,
                                        (sh.z() - z_shift)/* / size_scale*/ );

//                    std::cout << "\nOld point: " << sh << std::endl;
//                    std::cout << "New point: " << alcc.point(it) << std::endl;
                global_shift_pt = Point(x_shift, y_shift, z_shift);
            }
        }
        else
            std::cout << "\nNo global shifting required!" << std::endl;
    }
}
