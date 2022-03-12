#ifndef INDOORGML_READER_H
#define INDOORGML_READER_H

#include "rapidxml/rapidxml.hpp"
#include "IndoorGML.h"

using namespace IndoorGML;

inline void remove_aligned_points(vec_pt3d &pts)
{
    std::cout<< "#Points Initial: " << pts.size() << std::endl;
    if ( pts.size() >= 3 )
    {
        vec_pt3d tmp;

        // First point
        Vector v1( pts[0], pts[pts.size()-1] ),
                v2( pts[0], pts[1] );
        if (!LCCtools::vectors_are_eps_colinear(v1,v2))
            tmp.push_back(pts[0]);

        // Middle points
        for(uint i=1; i<pts.size()-1; i++)
        {
            v1 = Vector( pts[i], pts[i-1] );
            v2 = Vector( pts[i], pts[i+1] );
            if (!LCCtools::vectors_are_eps_colinear(v1,v2))
                tmp.push_back(pts[i]);
        }

        // Last point
        v1 = Vector( pts[pts.size()-1], pts[pts.size()-2] );
        v2 = Vector( pts[pts.size()-1], pts[0] );
        if (!LCCtools::vectors_are_eps_colinear(v1,v2))
            tmp.push_back(pts[pts.size()-1]);

        pts = tmp;
        std::cout<< "#Points Removed: " << pts.size() << std::endl;
    }
}

// Can I do an uglier function?
inline str remove_irrelevant_around_string(str mystr, str toclean)
{
    str res, tmp_, tmp = mystr;
    // Remove all occurence of the "toclean" string until something else is reached
    for(str::iterator it = tmp.begin(); it!=tmp.end(); it++)
    {
        if ( *it != toclean[0] )
        {
            tmp.erase(tmp.begin(),it);
            break;
        }
    }
    // Flip the resulting string
    for(str::reverse_iterator it = tmp.rbegin(); it!=tmp.rend(); it++)
        tmp_.push_back(*it);

    // Reapeath the first operation
    for(str::iterator it = tmp_.begin(); it!=tmp_.end(); it++)
    {
        if ( *it != toclean[0] )
        {
            tmp_.erase(tmp_.begin(),it);
            break;
        }
    }
    // Flip again the resulting string (which is the final output)
    for(str::reverse_iterator it = tmp_.rbegin(); it!=tmp_.rend(); it++)
        res.push_back(*it);

    return res;
}


inline str clean_string(str mystr)
{
    str res = mystr;
    res = remove_irrelevant_around_string(res,"\n");
    res = remove_irrelevant_around_string(res," ");
    res = remove_irrelevant_around_string(res,"\n");

    return res;
}


// Return a string free from the "toclean" string
inline str remove_from_string(str mystr, str toclean)
{
    str res;
    for(str::iterator it = mystr.begin(); it!=mystr.end(); it++)
        if ( *it != toclean[0] )
            res += *it;

    return res;
}

// Split a string on the basis of the provided delimiter and return an array
inline std::vector<double> split_string_delimiter_double(str mystr, str delim)
{
    std::vector<double> res;
    std::size_t pos = 0;
    str tmp;
    while ((pos = mystr.find(delim)) != str::npos)
    {
        tmp = mystr.substr(0, pos);
        res.push_back(std::stod(tmp));
        mystr.erase(0, pos + delim.length());
    }

    res.push_back(std::stod(mystr)); //the last component
    return res;
}

// Return the seeked attribute name or the one of the first one
inline str_pair getNodeAttribute( rapidxml::xml_node<> *node, str name )
{
    str_pair res;
    for(const rapidxml::xml_attribute<>* a = node->first_attribute();
        a; a = a->next_attribute() )
    {
        res.first = a->name();
        res.second = a->value();
        if ( name == "" || std::strcmp (a->name(), name.c_str() ) == 0 )
            return res;
    }

    return res;
}


// Returns the Dart_handle of a face (SurfaceMember) of a CellSpace
// Can only handle SurfaceMember contents with gml:LinearRing
inline Dart_handle getSurfaceGeometry(LCC& alcc, rapidxml::xml_node<> *&sm)
{
    Dart_handle res = LCC::null_handle;

    rapidxml::xml_node<> *sm_chld = sm->first_node();
    while( std::strcmp (sm_chld->name(), "gml:LinearRing") != 0 )
        sm_chld = sm_chld->first_node();

    if ( std::strcmp (sm_chld->name(), "gml:LinearRing") == 0 )
    {
        vec_pt3d pts;
        for( rapidxml::xml_node<> *it_vertices = sm_chld->first_node();
             it_vertices; it_vertices = it_vertices->next_sibling() )
        {
            str pt_s(it_vertices->value());
            //std::cout << "Name --> " << it_vertices->name()
            //          << " Value --> " << it_vertices->value() << std::endl;

            std::vector<double> pt_vec = split_string_delimiter_double(clean_string(pt_s), " ");
//            std::cerr << "Nb of point in the poly: " << pt_vec.size() << std::endl;
            assert ( pt_vec.size() % 3 == 0 );
            for (uint i=0; i<pt_vec.size()-2; i+=3)
                pts.push_back( LCC::Point(pt_vec[i], pt_vec[i+1], pt_vec[i+2]) );
        }

        if ( LCCtools::Points_are_eps_close( pts[0], pts[ pts.size()-1 ] ) )
            pts.pop_back(); //Remove the last point as it is a copy of the first

        if ( LCCtools::poly_surface_area_vec3d(pts) < (0.01 * 0.01) )
        {
            std::cerr << "This polygon is degenerate!" << std::endl;
    //        remove_aligned_points(pts); //Remove aligned points in the list
        }
        else
            res = LCCtools::Insert_new_2cell(alcc, pts);
    }

    return res;
}


inline Dart_handle getCellSpaceGeometry(LCC& alcc, rapidxml::xml_node<> *&CS_chld)
{
    vec_dart faces_of_cell;

    std::cout << "\nGetting one cellspace:" << std::endl;
    rapidxml::xml_node<> *csg_chld = CS_chld->first_node();
    while( std::strcmp (csg_chld->name(), "gml:Shell") != 0 )
        csg_chld = csg_chld->first_node();

    if ( std::strcmp (csg_chld->name(), "gml:Shell") == 0 )
    {
//        std::cout<< "Found a Shell" << std::endl;
        // Looking for the first surfaceMember at the bottom of the branch
        rapidxml::xml_node<> *sm;
        while( csg_chld != 0 )
        {
            if (std::strcmp (csg_chld->name(), "gml:surfaceMember") == 0)
                sm = csg_chld;
            csg_chld = csg_chld->first_node();
        }

        // Go through the surfaceMember entities sibling by sibling
        if ( std::strcmp (sm->name(), "gml:surfaceMember") == 0 )
        {
//            std::cout<< "\tGot a surfaceMember!" << std::endl;
            for( ; sm; sm = sm->next_sibling() )
            {
                Dart_handle d = getSurfaceGeometry(alcc, sm);
                if ( d != LCC::null_handle )
                {
//                    LCCtools::Print_face(alcc, d);
                    faces_of_cell.push_back( d );
                }
            }

			std::cout << "\tNumber of faces for this cellSpace: " << faces_of_cell.size() << std::endl;
            if ( faces_of_cell.size() > 0
                 && LCCtools::Perform_Simple_Volume_Reconstruction_from_Polygon_Soup(alcc, faces_of_cell)
                 )
            {
                // Initialize the id attribute of the 3-cell
                alcc.set_attribute<3>( faces_of_cell[0], alcc.create_attribute<3>() );
                alcc.info<3>(faces_of_cell[0]).set_id( getNodeAttribute( CS_chld->parent(), "gml:id" ).second );

//                LCCtools::Print_face(alcc, faces_of_cell[0]);

                return faces_of_cell[0];
            }
        }
    }

    return LCC::null_handle;
}


inline Dart_handle getCellSpaceBoundaryGeometry(LCC& alcc, rapidxml::xml_node<> *&CS_chld)
{
    Dart_handle res = LCC::null_handle;

//    std::cout << "\nGetting one cellspaceBoundary:" << std::endl;
    rapidxml::xml_node<> *csg_chld = CS_chld->first_node();
    res = getSurfaceGeometry(alcc, csg_chld);

    // Initialize the id attribute of the 3-cell
    if ( res != LCC::null_handle )
    {
        alcc.set_attribute<3>( res, alcc.create_attribute<3>() );
        alcc.info<3>(res).set_id( getNodeAttribute( CS_chld->parent(), "gml:id" ).second );
    }

    return res;
}


inline CellSpace getCellSpace(LCC& alcc, rapidxml::xml_node<> *&CS)
{
    CellSpace res;
    res.className = CS->name();
    res.setId( getNodeAttribute(CS, "gml:id").second );

    for( rapidxml::xml_node<> *CS_chld = CS->first_node();
         CS_chld; CS_chld = CS_chld->next_sibling() )
    {
//        std::cout << "\tName: " << CS_chld->name() << " Value: " << CS_chld->value();
        if (std::strcmp (CS_chld->name(), "gml:description") == 0)
            res.description = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "gml:name") == 0)
            res.name = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "core:cellSpaceGeometry") == 0)
            res.setCellSpaceGeometry( getCellSpaceGeometry(alcc, CS_chld) );

        else if (std::strcmp (CS_chld->name(), "core:duality") == 0)
        {
            str xlink( getNodeAttribute(CS_chld, "xlink:href").second );
            if ( *(xlink.begin()) == '#'  )
                xlink.erase(0,1); // Get rid of the "#"
            res.setDuality( xlink );
        }

        else if (std::strcmp (CS_chld->name(), "core:partialboundedBy") == 0)
        {
            str xlink( getNodeAttribute(CS_chld, "xlink:href").second );
            if ( *(xlink.begin()) == '#'  )
                xlink.erase(0,1); // Get rid of the "#"
            std::map<std::string, bool, LCCtools::cmp_string> bounded;
            bounded[ xlink ] = true;
            res.setPartialboundedBy( bounded );
        }

        else if (std::strcmp (CS_chld->name(), "navi:class") == 0)
            res.naviclass = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "navi:function") == 0)
            res.function = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "navi:usage") == 0)
            res.usage = CS_chld->value();

    }

    return res;
}


inline CellSpaceBoundary getCellSpaceBoundary(LCC& alcc, rapidxml::xml_node<> *&CS)
{
    CellSpaceBoundary res;
    res.className = CS->name();
    res.setId( getNodeAttribute(CS, "gml:id").second );

    for( rapidxml::xml_node<> *CS_chld = CS->first_node();
         CS_chld; CS_chld = CS_chld->next_sibling() )
    {
//        std::cout << "\tName: " << CS_chld->name() << " Value: " << CS_chld->value();
        if (std::strcmp (CS_chld->name(), "gml:description") == 0)
            res.description = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "gml:name") == 0)
            res.name = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "core:cellSpaceBoundaryGeometry") == 0)
            res.setCellSpaceBoundaryGeometry( getCellSpaceBoundaryGeometry(alcc, CS_chld) );

        else if (std::strcmp (CS_chld->name(), "core:duality") == 0)
        {
            str xlink( getNodeAttribute(CS_chld, "xlink:href").second );
            if ( *(xlink.begin()) == '#'  )
                xlink.erase(0,1); // Get rid of the "#"
            res.setDuality( xlink );
        }

        else if (std::strcmp (CS_chld->name(), "navi:class") == 0)
            res.naviclass = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "navi:function") == 0)
            res.function = CS_chld->value();

        else if (std::strcmp (CS_chld->name(), "navi:usage") == 0)
            res.usage = CS_chld->value();

    }

    return res;
}


inline void getPrimalSpaceFeatures_node( LCC& alcc, IndoorFeatures& InFt, rapidxml::xml_node<> *&pSF )
{
    // Loop through the PrimalSpaceFeatures nodes (just one iteration normally)
//    for(rapidxml::xml_node<> *node_PSF = pSF->first_node();
//        node_PSF; node_PSF = node_PSF->next_sibling())
    PrimalSpaceFeatures psf;

    rapidxml::xml_node<> *node_PSF = pSF->first_node();
    if ( std::strcmp (node_PSF->name(), "core:PrimalSpaceFeatures") == 0 )
    {
//        std::cout << node_PSF->name() << std::endl;
        psf.setId( getNodeAttribute( node_PSF, "gml:id" ).second );

        // Loop through the children of the PrimalSpaceFeatures node
        // Correspond to core:cellSpaceMember and core:cellSpaceBoundaryMember nodes
        for( rapidxml::xml_node<> *node_cSM = node_PSF->first_node();
             node_cSM; node_cSM = node_cSM->next_sibling() )
        {
            if( std::strcmp( node_cSM->name(), "core:cellSpaceMember" ) == 0 )
            {
                for( rapidxml::xml_node<> *node_CS = node_cSM->first_node();
                     node_CS; node_CS = node_CS->next_sibling() )
                {
                    psf.addCellSpaceMember( getCellSpace(alcc, node_CS) );
                }
            }
            else if( std::strcmp( node_cSM->name(), "core:cellSpaceBoundaryMember" ) == 0 )
            {
                for( rapidxml::xml_node<> *node_CS = node_cSM->first_node();
                     node_CS; node_CS = node_CS->next_sibling() )
                {
                    psf.addCellSpaceBoundaryMember( getCellSpaceBoundary(alcc, node_CS) );
                }
            }
        }
    }
    else
        std::cerr << "Could not get the 'core:PrimalSpaceFeatures' node!" << std::endl;

    InFt.setprimalSpaceFeatures(psf);
}


inline State getStateMember(rapidxml::xml_node<> *&node_n)
{
    State st;

    rapidxml::xml_node<> *cur_st = node_n->first_node();
    if ( std::strcmp( cur_st->name(), "core:State" ) == 0 )
    {
        st.setId( getNodeAttribute(cur_st, "gml:id").second );

        for( rapidxml::xml_node<> *st_chld = cur_st->first_node();
             st_chld; st_chld = st_chld->next_sibling() )
        {
            if( std::strcmp( st_chld->name(), "gml:description" ) == 0 )
                st.description = clean_string(st_chld->value());

            else if ( std::strcmp( st_chld->name(), "gml:name" ) == 0 )
                st.name = clean_string(st_chld->value());

            else if ( std::strcmp( st_chld->name(), "core:duality" ) == 0 )
            {
                str xlink( getNodeAttribute(st_chld, "xlink:href").second );
                if ( *(xlink.begin()) == '#'  )
                    xlink.erase(0,1); // Get rid of the "#"
                st.setDuality( xlink );
            }

            else if (std::strcmp( st_chld->name(), "core:connects") == 0)
            {
                str xlink( getNodeAttribute(st_chld, "xlink:href").second );
                if ( *(xlink.begin()) == '#'  )
                    xlink.erase(0,1); // Get rid of the "#"
                st.addConnects( xlink );
            }

            else if (std::strcmp( st_chld->name(), "core:geometry") == 0)
            {
                rapidxml::xml_node<> *geom_chld = st_chld->first_node();
                while ( geom_chld != 0 )
                {
                    if ( std::strcmp( geom_chld->name(), "gml:pos") == 0 )
                    {
                        str pt_s(geom_chld->value());
                        std::vector<double> pt_vec = split_string_delimiter_double(clean_string(pt_s), " ");
                        st.setGeometry( LCC::Point(pt_vec[0], pt_vec[1], pt_vec[2]) );
                    }

                    geom_chld = geom_chld->first_node();
                }
            }
        }
    }
    else
        std::cerr << "Could not get a 'core:State' node in this stateMember!" << std::endl;

    return st;
}

inline Transition getTransitionMember(rapidxml::xml_node<> *&node_e)
{
    Transition trs;

    rapidxml::xml_node<> *cur_trs = node_e->first_node();
    if ( std::strcmp( cur_trs->name(), "core:Transition" ) == 0 )
    {
        trs.setId( getNodeAttribute(cur_trs, "gml:id").second );

        for( rapidxml::xml_node<> *trs_chld = cur_trs->first_node();
             trs_chld; trs_chld = trs_chld->next_sibling() )
        {
            if( std::strcmp( trs_chld->name(), "gml:description" ) == 0 )
                trs.description = clean_string(trs_chld->value());

            else if ( std::strcmp( trs_chld->name(), "gml:name" ) == 0 )
                trs.name = clean_string(trs_chld->value());

            else if ( std::strcmp( trs_chld->name(), "gml:weight" ) == 0 )
                trs.setWeight( trs_chld->value() );

            else if (std::strcmp( trs_chld->name(), "core:connects") == 0)
            {
                str xlink( getNodeAttribute(trs_chld, "xlink:href").second );
                if ( *(xlink.begin()) == '#'  )
                    xlink.erase(0,1); // Get rid of the "#"
                trs.addConnects( xlink );
            }

            else if ( std::strcmp( trs_chld->name(), "core:duality" ) == 0 )
            {
                str xlink( getNodeAttribute(trs_chld, "xlink:href").second );
                if ( *(xlink.begin()) == '#'  )
                    xlink.erase(0,1); // Get rid of the "#"
                trs.setDuality( xlink );
            }

            else if (std::strcmp( trs_chld->name(), "core:geometry") == 0)
            {
                rapidxml::xml_node<> *geom_chld = trs_chld->first_node();
                while ( geom_chld != 0 )
                {
                    if (std::strcmp( geom_chld->name(), "gml:pos") == 0
                            || std::strcmp( geom_chld->name(), "gml:posList") == 0)
                    {
                        rapidxml::xml_node<> *geom_chld_pos = geom_chld;
                        while ( geom_chld_pos != 0 )
                        {
                            str pt_s(geom_chld_pos->value());
                            std::vector<double> pt_vec = split_string_delimiter_double(clean_string(pt_s), " ");
                            for (uint i=0; i<pt_vec.size()-2; i+=3)
                                trs.addGeometry( LCC::Point(pt_vec[i], pt_vec[i+1], pt_vec[i+2]) );

                            geom_chld_pos = geom_chld_pos->next_sibling();
                        }
                    }
                    geom_chld = geom_chld->first_node();
                }
            }
        }
    }
    else
        std::cerr << "Could not get a 'core:Transition' node in this transitionMember!" << std::endl;

    return trs;
}

inline SpaceLayer getSpaceLayer( rapidxml::xml_node<> *SL )
{
    SpaceLayer res;
    res.setId( getNodeAttribute(SL, "gml:id").second );

    for( rapidxml::xml_node<> *node_sl = SL->first_node();
         node_sl; node_sl = node_sl->next_sibling() )
    {
        if( std::strcmp( node_sl->name(), "core:nodes" ) == 0 )
        {
            for( rapidxml::xml_node<> *node_n = node_sl->first_node();
                 node_n; node_n = node_n->next_sibling() )
            {
                if( std::strcmp( node_n->name(), "core:stateMember" ) == 0 )
                    res.addNode( getStateMember(node_n) );
            }
        }

        else if( std::strcmp( node_sl->name(), "core:edges" ) == 0 )
        {
            for( rapidxml::xml_node<> *node_e = node_sl->first_node();
                 node_e; node_e = node_e->next_sibling() )
            {
                if( std::strcmp( node_e->name(), "core:transitionMember" ) == 0 )
                    res.addTransition( getTransitionMember(node_e) );
            }
        }
    }

    return res;
}

inline void getMultiLayeredGraph_node( IndoorFeatures& InFt, rapidxml::xml_node<> *&mLG )
{
    MultiLayeredGraph mlg;

    rapidxml::xml_node<> *node_MLG = mLG->first_node();
    if ( std::strcmp (node_MLG->name(), "core:MultiLayeredGraph") == 0 )
    {
        mlg.setId( getNodeAttribute( node_MLG, "gml:id" ).second );

        // Getting to the core:SpaceLayers nodes
        for( rapidxml::xml_node<> *node_SLs = node_MLG->first_node();
             node_SLs; node_SLs = node_SLs->next_sibling() )
        {
            if( std::strcmp( node_SLs->name(), "core:spaceLayers" ) == 0 )
            {
                SpaceLayers SLs;
                SLs.setId( getNodeAttribute( node_SLs, "gml:id" ).second );

                // Getting to the core:SpaceLayer nodes (first node of the core:spaceLayerMember)
                for( rapidxml::xml_node<> *node_SL = node_SLs->first_node();
                     node_SL; node_SL = node_SL->next_sibling() )
                {
                    if( std::strcmp( node_SL->name(), "core:spaceLayerMember" ) == 0 )
                        SLs.addSpaceLayer( getSpaceLayer( node_SL->first_node() ) );
                }

                mlg.addSpaceLayers( SLs );
            }
        }
    }
    else
        std::cerr << "Could not get a 'core:MultiLayeredGraph' node!" << std::endl;

    InFt.setmultiLayeredGraph( mlg );
}


inline void getIndoorFeatures_node( LCC& alcc, rapidxml::xml_document<> &root, IndoorFeatures& InFt )
{
    rapidxml::xml_node<> *node = root.first_node(); // corresponds to the root node
    str myname = node->name();

    if ( std::strcmp( myname.c_str(), "core:IndoorFeatures" ) == 0 )
    {
        InFt.setId( getNodeAttribute( node, "gml:id" ).second );

        // get all the attributes in the header
        for(const rapidxml::xml_attribute<>* a = node->first_attribute();
            a; a = a->next_attribute() )
            InFt.addheader( str_pair( a->name(), a->value() ) );

        for(node = node->first_node(); node; node = node->next_sibling())
        {
            myname = node->name();
            if ( std::strcmp(myname.c_str(), "core:primalSpaceFeatures" ) == 0 )
                getPrimalSpaceFeatures_node(alcc, InFt, node);
            else if ( std::strcmp(myname.c_str(), "core:multiLayeredGraph" ) == 0 )
                getMultiLayeredGraph_node(InFt, node);
        }
    }


}

inline void readIndoorGML(LCC& alcc, IndoorFeatures& InFt, str myfile)
{
    std::ifstream aFile (myfile);
    std::vector<char> buffer((std::istreambuf_iterator<char>(aFile)), std::istreambuf_iterator<char>());
//    std::vector<char> aFile (myfile.c_str(), myfile.c_str() + myfile.size() + 1);
    buffer.push_back('\0');

    rapidxml::xml_document<> mydoc;
    mydoc.parse<0>(&buffer[0]);

//    std::cout << "Name of my first node is: " << mydoc.first_node()->name() << std::endl;
//    std::cout << "Node has value: " << mydoc.first_node()->first_attribute()->name() << " = "
//                                    << mydoc.first_node()->first_attribute()->value() << std::endl;
//    std::cout << "Node has value: " << mydoc.first_node()->first_node()->first_node()->name() << std::endl;

//    IndoorFeatures InFt;
    getIndoorFeatures_node(alcc, mydoc, InFt);




//    rapidxml::xml_node<> *PS = nullptr;
//    if (PS != nullptr)
//    {
//        std::cout << "Node PS has value: " << PS->first_attribute()->name() << " = "
//                  << PS->first_attribute()->value() << std::endl;
//    }



//    for(node = mydoc.first_node(); node!=NULL; node = mydoc.next_sibling())
//    {
//        str myname = node->name();
//        if (myname == "PrimalSpaceFeatures")
//            std::cout << "Value of PrimalSpaceFeatures: " << node->value() << std::endl;
//    }

}



#endif
