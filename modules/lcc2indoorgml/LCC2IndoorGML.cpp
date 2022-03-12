#include "LCC2IndoorGML.h"

int Layers_num = 0,Layer_num = 0;

std::map<str, str> ori3CellIDs2IndoorGMLIds;
std::map<str, Dart_handle> cellspace_dart, cellboundary_dart;

// Containers to ensure export with RapidXML
// To store temporary strings until they are written in the xml
std::map<unsigned int, str> tmp_str;
std::map<unsigned int, str> tmp_vec_str;
std::map<unsigned int, xml_node<>*> tmp_xml_node;
unsigned int i_tmp_str = 0, i_tmp_vec_str = 0, i_tmp_xml_node = 0;

std::map<unsigned int, str> tmp_str_ops;
std::map<unsigned int, str_pair > tmp_ordered_pairs;
unsigned int i_tmp_str_ops = 0, i_tmp_ordered_pairs = 0;

// A map of all the cells in the scene
std::unordered_map<std::string, Dart_handle> allCellsMap;
std::unordered_map<std::string, Dart_handle> allCellNodesDartMap;

// An LCC2IndoorGML instance
LCC2IndoorGML lcc2igmlHandler;

std::size_t generateNetworkFromLCC(LCC& alcc, MainWindow* mw, bool use3Links){
    alcc.set_update_attributes(false);

    if (!lcc2igmlHandler.isReady){
        lcc2igmlHandler.init(alcc);
        std::cout << "Lcc2IndoorGML Handler initialised." << std::endl;
    }

    // Collect all the 3-cells (CellSpaces) of the LCC in a container
    // because new 3-cells will be created
    vec_dart all3Cells;
    typename LCC::Base::template One_dart_per_cell_range<3>::iterator
            itCell = alcc.one_dart_per_cell<3>().begin(),
            itCellEnd = alcc.one_dart_per_cell<3>().end();
    for(; itCell!=itCellEnd; itCell++){
        if( itCell != LCC::null_handle ){
            all3Cells.push_back( itCell );
            allCellsMap[ alcc.info<3>(itCell).id() ] = itCell;
        }
    }

    // Vertex layer
    Drawer::vtx_layer vl;
    vl.layer_color = CGAL::Color( 0, 255, 0 );
    // Edges layer
    Drawer::edge_layer el;
    el.layer_color = CGAL::Color( 0, 0, 255 );

    for(uint c=0; c < all3Cells.size(); c++){
        // Add the centroid of the 3-cell to the vertex list
        Point ct = getCellCentroid(alcc, all3Cells[c]);
        vl.vertices.push_back( ct );
        std::vector<std::string> cellNeighbors = alcc.info<3>(all3Cells[c]).relatedCells();
        for( uint i=0; i<cellNeighbors.size(); i++ ){
            Drawer::edge e;
            e.v1 = ct;
            e.v2 = getCellCentroid(alcc, allCellsMap[ cellNeighbors[i] ]);
            el.edges.push_back( e );
        }
    }

    std::string layer_name = "Network_default";
    mw->viewer->add_vtx_layer( layer_name, vl );
    mw->viewer->add_edge_layer( layer_name, el );

    alcc.set_update_attributes(true);

    return el.edges.size();
}




void LCC2IndoorGML::load_LCC_to_IndoorGML1(LCC& alcc)
{
    IndoorFeatures InFt;
    int c = 0;

    typename LCC::Base::One_dart_per_cell_range<3>::iterator
            it = alcc.one_dart_per_cell<3>().begin(),
            itend = alcc.one_dart_per_cell<3>().end();

    for(; it!=itend; it++)
    {
        if ( it != LCC::null_handle )
        {
//            if( alcc.info<3>(it).id() != "" )
//                ori3CellIDs2IndoorGMLIds[ alcc.info<3>(it).id() ] = it;

            std::string st = std::to_string(c);
//            if ( alcc.info<3>(it).id() == "" )
//                alcc.info<3>(it).set_id( "vol" + st );

            CellSpace cs;
            std::string id = *(cs.getId()) + "_" + st;
            cs.setId( id );

            cs.setCellSpaceGeometry( it );
            // Give the cell id to the 3-cells
//            alcc.info<3>(it).set_id( id );
            if( ori3CellIDs2IndoorGMLIds.find( alcc.info<3>(it).id() ) == ori3CellIDs2IndoorGMLIds.end() )
                ori3CellIDs2IndoorGMLIds[ alcc.info<3>(it).id() ] = id;

            InFt.getprimalSpaceFeatures()->addCellSpaceMember( cs );
            c++;
        }
        else
            std::cout << "\nNULL DART (load LCC to IndoorGML)!" <<std::endl;
    }

    // Add the dual graph
    createDualLayer1(alcc, InFt);

    setIndoorFeatures( InFt );
    std::cout << "Model compatible with IndoorGML(v1.1): " <<  InFt.getprimalSpaceFeatures()->getCellSpaceMember()->size()
              << " CellSpace layer(s) and " << InFt.getmultiLayeredGraph()->getSpaceLayers()->size()
              << " NRG layer(s)." << std::endl;
}



///
/// \brief LCC_demo_IndoorGML2LCC_plugin::create_New_ThematicLayer_From_Visible_Scene
/// Adds
void LCC2IndoorGML::create_New_ThematicLayer(LCC &alcc, std::string id, bool visibleCellsOnly)
{
    IndoorGML2::ThematicLayer myTL;
    if ( strcmp(id.c_str(), "") != 0 )
        myTL.setId(id);
    IndoorGML2::IndoorFeatures *InFt = getIndoorFeatures2();
    int c = 0, sem = 0;

    vec_dart cells;
    typename LCC::Base::One_dart_per_cell_range<3>::iterator
            it = alcc.one_dart_per_cell<3>().begin(),
            itend = alcc.one_dart_per_cell<3>().end();
    for(; it!=itend; it++){
        if( it != LCC::null_handle ){
            if( visibleCellsOnly && alcc.info<3>(it).is_visible() ){
                cells.push_back(it);
            }
            else if (!visibleCellsOnly){
                cells.push_back(it);
            }
        }
        else
            std::cout << "\nNULL DART (load LCC to IndoorGML)!" <<std::endl;
    }

//    myTL.getprimalSpaceLayer()->setId( "PSL" + std::to_string(Layer_num++));
    // Get the PrimalSpace
    for( uint i=0; i<cells.size(); i++){
//            if (alcc.info<3>(it).id() != "" )
//                        ori3CellIDs2IndoorGMLIds[ alcc.info<3>(it).id() ] = it;

        // Handle CellSpaces
        std::string st = std::to_string(c);
//            if ( alcc.info<3>(it).id() == "" )
//                alcc.info<3>(it).set_id( "vol" + st );

        IndoorGML2::CellSpace cs;
        std::string id = *(cs.getId()) + "_" + st;
        cs.setId( id );

        // Get semantic info (if LCC is from IFC)
        std::string className = alcc.info<3>(cells[i]).semClass();
        if( (className.length() >= 8 && className.substr( className.length() - 8 ) == "IfcSpace") )
        {
            cs.className = "navi:GeneralSpace";
            cs.naviclass = "undefined";
            cs.function = "undefined";
            sem++;
        }
        else if ( (className.length() >= 4 && className.substr( className.length() - 4 ) == "Door")
                  || (className.length() >= 6 && className.substr( className.length() - 6 ) == "Window")
                  || (className.length() >= 14 && className.substr( className.length() - 14 ) == "OpeningElement"))
        {
            cs.className = "navi:TransitionSpace";
            cs.naviclass = "undefined";
            cs.function = "undefined";
            sem++;
        }

        else if( (className.length() >= 7 && className.substr( className.length() - 7 ) == "F-Space") )
        {
            cs.className = "navi:GeneralSpace";
            cs.naviclass = "undefined";
            cs.function = "undefined";
            sem++;
        }
        else if( (className.length() >= 7 && className.substr( className.length() - 7 ) == "O-Space") )
        {
            cs.isNavigable = false;
            cs.className = "navi:NonNavigableSpace";
            sem++;
        }

        // Get the geometry
        cs.setCellSpaceGeometry( cells[i] );
        // Give the cell id to the 3-cells
//                alcc.info<3>(it).set_id( id );
        if( ori3CellIDs2IndoorGMLIds.find( alcc.info<3>(cells[i]).id() ) == ori3CellIDs2IndoorGMLIds.end() )
            ori3CellIDs2IndoorGMLIds[ alcc.info<3>(cells[i]).id() ] = id;

        myTL.getprimalSpaceLayer()->addCellSpaceMember( cs );



        // Handle CellBoundaries

        // Run through all the 2-cells (faces) of the current volume
//            typename LCC::Base::One_dart_per_incident_cell_range<2,3>::iterator
//                    itf = alcc.one_dart_per_incident_cell<2,3>(it).begin(),
//                    itfend = alcc.one_dart_per_incident_cell<2,3>(it).end();
//            for(; itf!=itfend; itf++)
//            {
//                // If the face is 3-sewn it means that is will be an edge in the dual space
//                if( !alcc.is_free(itf, 3) ){
//                    IndoorGML2::CellBoundary cb;
//                    std::string id_cb = *(cb.getId()) + "_" + st;
//                    cb.setId( id_cb );


//                }
//            }


        c++;
    }

    // Get the DualSpace
    createDualLayer2(alcc, *InFt, myTL, (sem>0));

    setIndoorFeatures2( *InFt );
    std::cout << "Model compatible with IndoorGML(v2): " <<  InFt->getThematicLayerMap()->size()
              << " Thematic layer(s) and " << InFt->getInterLayerConnectionMap()->size()
              << " Interlayer connection(s)." << std::endl;

//    msgBox.setText("New thematic layer successfully created.");
//    msgBox.exec();
}


///
/// \brief LCC_demo_IndoorGML2LCC_plugin::on_actionLoad_LCC_to_IndoorGML2_triggered
/// Loads the current LCC as a fresh IndoorGML2 class (only one thematic layer so far)
///
void LCC2IndoorGML::load_LCC_to_IndoorGML2(LCC &alcc){
    create_New_ThematicLayer(alcc, "");
//        view_All_NRG_Layers2(mw, *getIndoorFeatures2());
}



void LCC2IndoorGML::read_IndoorGML(LCC &alcc, MainWindow* mw)
{
    QString fileName = QFileDialog::getOpenFileName (mw,
                                                     "Load IndoorGML files",
                                                     "./gml",
                                                     "IndoorGML files (*.gml)");

    if (!fileName.isEmpty ())
    {
        load_IndoorGML(alcc, fileName.toStdString());
        Q_EMIT(mw->sceneChanged());
    }
}

void LCC2IndoorGML::load_IndoorGML(LCC &alcc, str myfile, int v)
{
    if (v == 1)
    {
        IndoorFeatures InFt;
        readIndoorGML(alcc, InFt, myfile);
        std::cout << "The loaded LCC counts " << alcc.number_of_vertex_attributes() << " vertices!" << std::endl;
        setIndoorFeatures( InFt );

        // Shifts only when necessary
        global_shift_pt = Point( 0.0, 0.0, 0.0 );
        LCCtools::Apply_global_shifting(alcc);

        //    vec_dart to_clean;
        //    LCCtools::Keep_corners_points_only(alcc, to_clean);

        //    std::cout << "Number of 0-cell to be removed: " << to_clean.size() << std::endl;
        //    LCCtools::Remove_selected_0_cells(alcc, to_clean);
        //    std::cout << "\tDONE removing them all!" << std::endl;

        /// Trying to visualize imported network
//        view_All_NRG_Layers( mw, InFt.getmultiLayeredGraph() );

        std::cout << "Detected " << cellboundary_dart.size() << " cellSpaceBoundaries!" << std::endl;
    }
}

void LCC2IndoorGML::write_IndoorGML1(LCC &alcc, MainWindow *mw)
{

    if (!myIndoorGMLisEmpty())
    {
        CGAL::Timer t;
        QString fileName = QFileDialog::getSaveFileName (mw,
                                                         "Export as an IndoorGML file",
                                                         "myIndoorGML1",
                                                         "IndoorGML files (*.gml)");
        if (!fileName.isEmpty())
        {
            t.start();
            if ( !save_IndoorGML(alcc, fileName.toStdString()) )
                std::cout << "\nCould not generate an IndoorGMLv1.0.3 file!" << std::endl;
            t.stop();
            std::cout << "Time to export IndoorGMLv1.1: " << t.time() << std::endl;
        }
    }
    else
    {
        std::cout << "No IndoorGML entity to export!" << std::endl;
//        msgBox.setText("No IndoorGML entity to export!");
//        msgBox.exec();
    }
}

void LCC2IndoorGML::write_IndoorGML2(LCC &alcc, MainWindow *mw)
{
    if (!myIndoorGML2isEmpty())
    {
        QString fileName = QFileDialog::getSaveFileName (mw,
                                                         "Export as an IndoorGML file",
                                                         "myIndoorGML2",
                                                         "IndoorGML files (*.gml)");
        if (!fileName.isEmpty())
        {
            CGAL::Timer t;
            t.start();
            if ( !save_IndoorGML(alcc, fileName.toStdString(), 2) )
                std::cout << "\nCould not generate an IndoorGMLv2.0 file!" << std::endl;
            t.stop();
            std::cout << "Time to export IndoorGMLv2 (beta): " << t.time() << std::endl;
        }
    }
    else
    {
        std::cout << "No IndoorGML entity to export!" << std::endl;
//        msgBox.setText("No IndoorGML entity to export!");
//        msgBox.exec();
    }
}

bool LCC2IndoorGML::save_IndoorGML(LCC &alcc, str out, int v)
{
    std::ofstream outfile (out.c_str());
    if  (outfile.is_open())
    {
        if ( v == 1 )
        {
            writeIndoorGML1(alcc, *(getIndoorFeatures()), outfile );
            std::cout << "\nSuccessfully generated IndoorGMLv1.1 file: " << out << std::endl;

//            msgBox.setText("Successfully saved the IndoorGML(v1.1) file");
//            msgBox.exec();

            return true;
        }
        else if ( v == 2 )
        {
            writeIndoorGML2(alcc, *(getIndoorFeatures2()), outfile );
            std::cout << "\nSuccessfully generated IndoorGMLv2.0 file: " << out << std::endl;

//            msgBox.setText("Successfully saved the IndoorGML(v2.0 Beta) file");
//            msgBox.exec();

            return true;
        }
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
        return false;
    }
}


void LCC2IndoorGML::add_ThematicLayer(LCC& alcc)
{
//#ifdef CGAL_PROFILE_LCC_DEMO
//    CGAL::Timer timer;
//    timer.start();
//#endif

    std::vector<std::string> existingLayers;
    for( auto &it_TL : *(getIndoorGML2()->getThematicLayerMap()) )
        existingLayers.push_back( it_TL.first );

    newTLdialog = new AddThematicLayer(0, &("TL_" + LCCtools::generate_unique_ID()), &existingLayers);
//    std::cout << "Num of EXISITING Layers: " << existingLayers.size() << std::endl;

    newTLdialog->exec();


    if (newTLdialog->done)
    {
//        std::cout<<"Layer ID to be created: " << newTLdialog->layerID << std::endl;
        if ( alcc.number_of_vertex_attributes() > 0 )
        {
            // Creating a new ThematicLayer with only visible & filled entities
            create_New_ThematicLayer(alcc, newTLdialog->layerID, true);

            if ( newTLdialog->selectedLayers.size() > 0 )
            {
                // TL1 is the pre-existing ThematicLayer, while TL2 is the newly created one
                IndoorGML2::ThematicLayer* TL1, *TL2;
                std::map<std::string, IndoorGML2::ThematicLayer>* thematicLayerMap = getIndoorFeatures2()->getThematicLayerMap();

//                std::cout << "\nThematicLayers before? " << thematicLayerMap->size() << " -> Layer ID: " << newTLdialog->layerID << std::endl;

                //    std::map<std::string, IndoorGML2::InterLayerConnection>* interLayerConnectionMap = getIndoorFeatures2()->getInterLayerConnectionMap();
                TL2 = &(thematicLayerMap->operator[](newTLdialog->layerID));

                // For every selected layer (in the UI dialog) for ILC
                for ( auto &it_Layer : newTLdialog->selectedLayers )
                {
//                    std::cout << "Is this happening?" << std::endl;

                    TL1 = &(thematicLayerMap->operator[](it_Layer.first));
                    // Creating the corresponding ILC
                    IndoorGML2::InterLayerConnection newILC;
                    newILC.ConnectedLayers = std::make_pair(TL1, TL2);
                    newILC.typeOfTopoExpression = it_Layer.second;

                    // Add a new ILC linking the two TL
                    getIndoorFeatures2()->addInterLayerConnection(newILC);

                    // Get the ID of the newly created layer as key and a pointer to their ILC as the value
                    // to update the IL_connections map of the pre-existing ThematicLayer
                    TL1->IL_connections[ *(TL2->getId()) ] = &(getIndoorFeatures2()->getInterLayerConnectionMap()->operator[](*(newILC.getId())));
                }

//                std::cout << "\nThematicLayers after? " << thematicLayerMap->size() << std::endl;
//                std::cout << "Num of SELECTED layers after: " << newTLdialog->selectedLayers.size() << std::endl;
            }

            //    view_All_NRG_Layers( mw, getIndoorFeatures()->getmultiLayeredGraph() );

            QApplication::restoreOverrideCursor ();
        }
        else{
            std::cout << "No LCC loaded! Make sure the scene contains geometries." << std::endl;
//            msgBox.setText("No LCC loaded! Make sure the scene contains geometries.");
//            msgBox.exec();
        }
    }

//#ifdef CGAL_PROFILE_LCC_DEMO
//    timer.stop();
//    std::cout<<"Time to generate the NRG: "
//           <<timer.time()<<" seconds."<<std::endl;
//#endif
}

void LCC2IndoorGML::init(LCC& alcc){
    if ( alcc.number_of_vertex_attributes() > 0 ){
        load_LCC_to_IndoorGML1(alcc);
        load_LCC_to_IndoorGML2(alcc);
        isReady = true;
    }else{
        std::cout << "No feature detected! Make sure the scene contains geometries." << std::endl;
//        msgBox.setText("No feature detected! Make sure the scene contains geometries.");
//        msgBox.exec();
    }
}



#include "LCC2IndoorGML.moc"
