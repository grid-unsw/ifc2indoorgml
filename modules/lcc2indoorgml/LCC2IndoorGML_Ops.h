#ifndef LCC2INDOORGML_OPS_H
#define LCC2INDOORGML_OPS_H

#include "MainWindow.h"
#include "IndoorGML.h"

extern int Layers_num, Layer_num;
extern std::map<str, str> ori3CellIDs2IndoorGMLIds;

// Containers to ensure export with RapidXML
extern std::map<unsigned int, str> tmp_str_ops;
extern std::map<unsigned int, str_pair > tmp_ordered_pairs;
extern unsigned int i_tmp_str_ops, i_tmp_ordered_pairs;

inline str_pair ordered_str( str s1, str s2 )
{
    if ( s1.compare(s2) < 0 )
        return str_pair( s1, s2 );
    else
        return str_pair( s2, s1 );
}

inline Point getCellCentroid(LCC& alcc, Dart_handle cell)
{
    Bbox_3 bb = LCCtools::Get_Bbox_vol(alcc, cell);
    Point pt( (bb.xmax() + bb.xmin())/2.0,
              (bb.ymax() + bb.ymin())/2.0,
              (bb.zmax() + bb.zmin())/2.0);

    return pt;
}

inline Point getCellCentroid_face(LCC& alcc, Dart_handle face)
{
    double x = 0.0, y = 0.0, z = 0.0;
    int count = 0;
    typename LCC::Base::Dart_of_orbit_range<1>::iterator itr(alcc, face);
    for (; itr.cont(); itr++)
    {
        x += alcc.point( itr ).x();
        y += alcc.point( itr ).y();
        z += alcc.point( itr ).z();
        count++;
    }

    return Point( x/count, y/count, z/count );
}

inline IndoorGML::State createState( LCC& alcc, Dart_handle cell, str SL_id, unsigned int i, bool use3links = false )
{
    IndoorGML::State res;

    if ( cell != LCC::null_handle )
    {
        str csID = ori3CellIDs2IndoorGMLIds[alcc.info<3>(cell).id()];
        res.setId( str(SL_id + "_ST" + csID ) );
        res.setDuality( csID );

        // Run through all the 2-cells (faces) of cell to fill the 'connects' relationships
        if( use3links ){ // If there are some sew3 links, use them
            typename LCC::Base::template One_dart_per_incident_cell_range<2,3>::iterator
                    it = alcc.one_dart_per_incident_cell<2,3>(cell).begin(),
                    itend = alcc.one_dart_per_incident_cell<2,3>(cell).end();
            for (; it != itend; it++)
            {
                // If a 2-cell is 3-sewn to another one, then both cells are linked
                if ( !alcc.is_free(it, 3) )
                {
                    tmp_str_ops[ i_tmp_str_ops++ ] = alcc.info<3>(alcc.beta<3>(it)).id();
                    tmp_ordered_pairs[ i_tmp_ordered_pairs++ ] = ordered_str( ori3CellIDs2IndoorGMLIds[alcc.info<3>(cell).id()], tmp_str_ops[ i_tmp_str_ops-1 ] );
                    res.addConnects( str ( SL_id + "_T"
                                           + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first + "-"
                                           + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second ) );
                }
            }
        }
        else{ // Otherwise, use relatedCell attributes if available
            std::vector<std::string> relatedCells = alcc.info<3>(cell).relatedCells();
            for(uint i=0; i<relatedCells.size(); i++){
                // Looping through all cells recorded as related to the current cell
                std::string oriCellID = relatedCells[i];
                // The new (IndoorGML) IDs of the related cells can be retrieved with the ori3CellIDs2IndoorGMLIds map
                std::string newCellID = ori3CellIDs2IndoorGMLIds[ oriCellID ];

                tmp_str_ops[ i_tmp_str_ops++ ] = newCellID;
                tmp_ordered_pairs[ i_tmp_ordered_pairs++ ] = ordered_str( ori3CellIDs2IndoorGMLIds[alcc.info<3>(cell).id()], tmp_str_ops[ i_tmp_str_ops-1 ] );
                res.addConnects( str ( SL_id + "_T"
                                       + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first + "-"
                                       + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second ) );
            }
        }

        res.setGeometry( getCellCentroid(alcc, cell) );
    }

    return res;
}


inline std::map<str, IndoorGML::Transition, LCCtools::cmp_string> createTransitions( LCC& alcc, str SL_id, bool use3links = false )
{
    std::map<str, IndoorGML::Transition, LCCtools::cmp_string> res;

    if( use3links ){
        // Run through all the 2-cells (faces) of the LCC
        typename LCC::Base::template One_dart_per_cell_range<2>::iterator it = alcc.one_dart_per_cell<2>().begin(),
                                                                 itend = alcc.one_dart_per_cell<2>().end();
        for(; it!=itend; it++)
        {
            IndoorGML::Transition trs;
            if ( alcc.info<3>(it).is_visible() && !alcc.is_free(it, 3) )
            {
                tmp_str_ops[ i_tmp_str_ops++ ] = alcc.info<3>(alcc.beta<3>(it)).id();
                tmp_ordered_pairs[ i_tmp_ordered_pairs++ ] = ordered_str( ori3CellIDs2IndoorGMLIds[alcc.info<3>(it).id()], tmp_str_ops[ i_tmp_str_ops-1 ] );
                str id = str(SL_id + "_T"
                             + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first + "-"
                             + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second );
                trs.setId( id );

                trs.addConnects( ori3CellIDs2IndoorGMLIds[alcc.info<3>(it).id()] );
                trs.addConnects( tmp_str_ops[ i_tmp_str_ops-1 ] );

    //            Dart_handle d1 = cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first ],
    //                        d2 = cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second ];

    //            if ( d1 == LCC::null_handle )
    //                std::cout << "\td1 is a NULL DART! -- " << i_tmp_ordered_pairs-1
    //                          << " -- " << tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first <<std::endl;
    //            else
    //                std::cout << "\tid of d1: " << alcc.info<3>(d1).id() <<std::endl;

    //            if ( d2 == LCC::null_handle )
    //                std::cout << "\td2 is a NULL DART! -- " << i_tmp_ordered_pairs-1
    //                          << " -- " << tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second <<std::endl;
    //            else
    //                std::cout << "\tid of d2: " << alcc.info<3>(d2).id() <<std::endl;


                std::vector<Point_3> geom = { getCellCentroid(alcc, cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first ] ),
                                              getCellCentroid_face(alcc, it),
                                              getCellCentroid(alcc, cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second ] )
                                            };
                trs.setGeometry( geom );
                res[id] = trs;
            }
        }
    }
    // Otherwise, use relatedCell attributes if available
    else{
        // Run through all the 3-cells (volumes) of the LCC
        typename LCC::Base::template One_dart_per_cell_range<3>::iterator it = alcc.one_dart_per_cell<3>().begin(),
                                                                          itend = alcc.one_dart_per_cell<3>().end();
        for(; it!=itend; it++)
        {
            if ( alcc.info<3>(it).is_visible() )
            {
                // Loop through the cells related to it
                std::vector<std::string> relatedCells = alcc.info<3>(it).relatedCells();
                for(uint i=0; i<relatedCells.size(); i++){
                    IndoorGML::Transition trs;
                    std::string oriCellID = relatedCells[i];
                    // The new (IndoorGML) IDs of the related cells can be retrieved with the ori3CellIDs2IndoorGMLIds map
                    std::string newCellID = ori3CellIDs2IndoorGMLIds[ oriCellID ];

                    tmp_str_ops[ i_tmp_str_ops++ ] = newCellID;
                    tmp_ordered_pairs[ i_tmp_ordered_pairs++ ] = ordered_str( ori3CellIDs2IndoorGMLIds[alcc.info<3>(it).id()], tmp_str_ops[ i_tmp_str_ops-1 ] );
                    str id = str(SL_id + "_T"
                                 + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first + "-"
                                 + tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second );
                    trs.setId( id );

                    trs.addConnects( ori3CellIDs2IndoorGMLIds[alcc.info<3>(it).id()] );
                    trs.addConnects( tmp_str_ops[ i_tmp_str_ops-1 ] );


                    std::vector<Point_3> geom = { getCellCentroid(alcc, cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].first ] ),
//                                                  getCellCentroid_face(alcc, it),
                                                  getCellCentroid(alcc, cellspace_dart[ tmp_ordered_pairs[ i_tmp_ordered_pairs-1 ].second ] )
                                                };
                    trs.setGeometry( geom );
                    res[id] = trs;
                }
            }
        }
    }


    return res;
}


inline void createSpaceLayer(LCC& alcc, IndoorGML::SpaceLayer& res, str id)
{
    if ( id != "" )
        res.setId( id );
    else
        id = *(res.getId());

    /*Usage,
    TerminationDate,
    Function,
    CreationDate,
    Class;*/

    unsigned int count = 0;
    typename LCC::Base::One_dart_per_cell_range<3>::iterator it = alcc.one_dart_per_cell<3>().begin(),
                                                             itend = alcc.one_dart_per_cell<3>().end();
    for(; it!=itend; it++)
    {
        // Apply only to filled and visible volumes
        if ( alcc.info<3>(it).is_visible() ){
            res.addNode( createState( alcc, it, id, count++ ) );
            res.setEdges( createTransitions( alcc, id ) );
        }
    }
}

inline void createSpaceLayer(LCC& alcc, IndoorGML::SpaceLayer& res)
{
    createSpaceLayer(alcc, res, "");
}

/// Updating the Cell-Node correspondance in the IndoorFeatures class (IndoorGMLv1)
/// (necessary because nodes are generated based on the LCC, independently of the CellSpace classes)
inline void updateCellSpaceDuality( IndoorGML::IndoorFeatures &InFt, IndoorGML::SpaceLayer &SL )
{
    for( auto& it : *SL.getNodes() )
    {
        tmp_str_ops[ i_tmp_str_ops++ ] = *it.second.getDuality();
        // The line below is reeeeally not cool... (clearer alternative commented just below it)
        InFt.getprimalSpaceFeatures()->getCellSpaceMember()->operator[](tmp_str_ops[ i_tmp_str_ops-1 ]).setDuality( *it.second.getId() );

//        std::map<std::string, IndoorGML::CellSpace, LCCtools::cmp_string> *it_csm = InFt.getprimalSpaceFeatures()->getCellSpaceMember();
//        IndoorGML::CellSpace *cs = &((*it_csm)[ tmp_str_ops.back() ]);
//        cs->setDuality( *it.second.getId() );
    }
}


/// Updating the Cell-Node correspondance in the IndoorFeatures class (IndoorGMLv2)
/// (necessary because nodes are generated based on the LCC, independently of the CellSpace classes
inline void updateCellSpaceDuality2( IndoorGML2::ThematicLayer &TL )
{
    IndoorGML2::DualSpaceLayer *DSL = TL.getdualSpaceLayer();
    if ( TL.getprimalSpaceLayer()->getCellSpaceMember()->size() > 0 )
    {
        for( auto& it : *DSL->getNodes() )
        {
            tmp_str_ops[ i_tmp_str_ops++ ] = *it.second.getDuality();
            TL.getprimalSpaceLayer()->getCellSpaceMember()->operator[](tmp_str_ops[ i_tmp_str_ops-1 ]).setDuality( *it.second.getId() );
//            std::cout << "CellSpace " << tmp_str_ops[ i_tmp_str_ops-1 ] << " gets Node " << *it.second.getId() << std::endl;
        }
    }
}


inline void createDualLayer1( LCC& alcc, IndoorGML::IndoorFeatures &InFt )
{
    // cleaning the containers to ensure export with RapidXML
    i_tmp_str_ops = 0;
    tmp_str_ops.clear();
    i_tmp_ordered_pairs = 0;
    tmp_ordered_pairs.clear();

    IndoorGML::SpaceLayers SLs;
    SLs.setId( str("NRG" + std::to_string(Layers_num++)) );

    IndoorGML::SpaceLayer SL;
    createSpaceLayer( alcc, SL );
    updateCellSpaceDuality( InFt, SL );
    SLs.addSpaceLayer( SL );

    InFt.getmultiLayeredGraph()->addSpaceLayers( SLs );
}

inline void createDualLayer2( LCC& alcc, IndoorGML2::IndoorFeatures &InFt, IndoorGML2::ThematicLayer myTL, bool semantic = false )
{
    // cleaning the containers to ensure export with RapidXML
    i_tmp_str_ops = 0;
    tmp_str_ops.clear();
    i_tmp_ordered_pairs = 0;
    tmp_ordered_pairs.clear();

    IndoorGML2::DualSpaceLayer DSL;
    createSpaceLayer( alcc, DSL );
    myTL.setdualSpaceLayer(DSL);
    updateCellSpaceDuality2( myTL );

    myTL.semanticExtension = semantic;
    InFt.addThematicLayer(myTL);
}


//void view_Single_NRG_Layer(MainWindow* mw, IndoorGML::MultiLayeredGraph* MLG)
//{
//    // For each SpaceLayers (SLS) get the map of SpaceLayer (SL)
//    IndoorGML::SpaceLayers SLayers = it_SLS.second;
//    std::map<std::string, IndoorGML::SpaceLayer, LCCtools::cmp_string> *sl = SLayers.getSpaceLayer();
//    for ( const auto &it_sl : *sl )
//    {
//        // Get the id of the Layer
//        std::string layer_id = it_sl.first;

//        // Get the nodes
//        Drawer::vtx_layer vl;
//        IndoorGML::SpaceLayer SLayer = it_sl.second;
//        std::map<std::string, State, LCCtools::cmp_string> *nodes = SLayer.getNodes();
//        for ( const auto &it_nd : *nodes )
//        {
//            IndoorGML::State node = it_nd.second;
//            vl.vertices.push_back( *node.getGeometry() );
//        }
//        // color of the vertices of the layer
//        vl.layer_color = CGAL::Color( 0, 255, 0 );
//        // Add layer to the drawer containers
//        mw->viewer->add_vtx_layer( layer_id, vl );

//        // Get the edges
//        Drawer::edge_layer el;
//        std::map<std::string, IndoorGML::Transition, LCCtools::cmp_string> *edges = SLayer.getEdges();
//        for ( const auto &it_edg : *edges )
//        {
//            IndoorGML::Transition edge = it_edg.second;
//            Drawer::edge e;
//            std::vector<Point_3> *vtx = edge.getGeometry();
//            if ( vtx->size() > 1 )
//            {
//                for ( int i = 1; i<vtx->size(); i++ )
//                {
//                    e.v1 = vtx->operator[](i-1);
//                    e.v2 = vtx->operator[](i);
//                    el.edges.push_back( e );
//                }
//            }
//        }
//        // color of the edges of the layer
//        el.layer_color = CGAL::Color( 0, 0, 255 );
//        // Add layer to the drawer containers
//        mw->viewer->add_edge_layer( layer_id, el );
//    }
//    std::cout << "\tLayer in MLG " << it_SLS.first << " collected for drawing!" << std::endl;
//}

inline void view_All_NRG_Layers(MainWindow* mw, IndoorGML::MultiLayeredGraph* MLG)
{
    /// Trying to visualize imported network
    std::map<std::string, IndoorGML::SpaceLayers, LCCtools::cmp_string> *SLS = MLG->getSpaceLayers();
    for( const auto &it_SLS : *SLS )
    {
        // For each SpaceLayers (SLS) get the map of SpaceLayer (SL)
        IndoorGML::SpaceLayers SLayers = it_SLS.second;
        std::map<std::string, IndoorGML::SpaceLayer, LCCtools::cmp_string> *sl = SLayers.getSpaceLayer();
        for ( const auto &it_sl : *sl )
        {
            // Get the id of the Layer
            std::string layer_id = it_sl.first;

            // Get the nodes
            Drawer::vtx_layer vl;
            IndoorGML::SpaceLayer SLayer = it_sl.second;
            std::map<std::string, State, LCCtools::cmp_string> *nodes = SLayer.getNodes();
            for ( const auto &it_nd : *nodes )
            {
                IndoorGML::State node = it_nd.second;
                vl.vertices.push_back( *node.getGeometry() );
            }
            // color of the vertices of the layer
            vl.layer_color = CGAL::Color( 0, 255, 0 );
            // Add layer to the drawer containers
            mw->viewer->add_vtx_layer( layer_id, vl );

            // Get the edges
            Drawer::edge_layer el;
            std::map<std::string, IndoorGML::Transition, LCCtools::cmp_string> *edges = SLayer.getEdges();
            for ( const auto &it_edg : *edges )
            {
                IndoorGML::Transition edge = it_edg.second;
                Drawer::edge e;
                std::vector<Point_3> *vtx = edge.getGeometry();
                if ( vtx->size() > 1 )
                {
                    for ( int i = 1; i<vtx->size(); i++ )
                    {
                        e.v1 = vtx->operator[](i-1);
                        e.v2 = vtx->operator[](i);
                        el.edges.push_back( e );
                    }
                }
            }
            // color of the edges of the layer
            el.layer_color = CGAL::Color( 0, 0, 255 );
            // Add layer to the drawer containers
            mw->viewer->add_edge_layer( layer_id, el );
        }
        std::cout << "\tLayer in MLG " << it_SLS.first << " collected for drawing!" << std::endl;
    }
}

inline void view_All_NRG_Layers2(MainWindow* mw, IndoorGML2::IndoorFeatures& InFt)
{
    /// Trying to visualize imported network
    std::map<std::string, IndoorGML2::ThematicLayer> *TL = InFt.getThematicLayerMap();
    for( auto &it_TL : *TL )
    {
        // For each ThematicLayer (TL) get the DualSpaceLayer (DSL)
        IndoorGML2::DualSpaceLayer* DSL = it_TL.second.getdualSpaceLayer();

        // Get the id of the Layer
        std::string layer_id = *(DSL->getId());

        // Get the nodes
        Drawer::vtx_layer vl;
        std::map<std::string, State, LCCtools::cmp_string> *nodes = DSL->getNodes();
        for ( const auto &it_nd : *nodes )
        {
            IndoorGML::State node = it_nd.second;
            vl.vertices.push_back( *node.getGeometry() );
        }
        // color of the vertices of the layer
        vl.layer_color = CGAL::Color( 0, 255, 0 );
        // Add layer to the drawer containers
        mw->viewer->add_vtx_layer( layer_id, vl );

        // Get the edges
        Drawer::edge_layer el;
        std::map<std::string, IndoorGML::Transition, LCCtools::cmp_string> *edges = DSL->getEdges();
        for ( const auto &it_edg : *edges )
        {
            IndoorGML::Transition edge = it_edg.second;
            Drawer::edge e;
            std::vector<Point_3> *vtx = edge.getGeometry();
            if ( vtx->size() > 1 )
            {
                for ( int i = 1; i<vtx->size(); i++ )
                {
                    e.v1 = vtx->operator[](i-1);
                    e.v2 = vtx->operator[](i);
                    el.edges.push_back( e );
                }
            }
        }
        // color of the edges of the layer
        el.layer_color = CGAL::Color( 0, 0, 255 );
        // Add layer to the drawer containers
        mw->viewer->add_edge_layer( layer_id, el );

        std::cout << "\tLayer " << layer_id << " in TL " << it_TL.first << " collected for drawing!" << std::endl;
    }
}


#endif
