#include "IFC2LCC.h"

///
/// \brief load_IFC
/// \param fileName
/// \param alcc
///
void load_IFC(const QString & fileName, LCC& alcc)
{
#ifdef IFCPP_ON
    QApplication::setOverrideCursor (Qt::WaitCursor);

    std::cout << "Loading IFC file (using IFC++ library)..." << std::endl;
    //            CGAL::ifc2lcc(alcc, fileName.toStdWString());

    //            CGAL::Merge_all_close_points(alcc);
    //            CGAL::remove_duplicated_vertices(alcc);
    //            CGAL::Remove_degenerate_poly(alcc);

    CGAL::Timer t;
    t.start();

    // pre-load the data in tables
    map_ifc_shapes.clear();
    map_ifc_entities.clear();
    ifc_storey_list.clear();
    spaces_of_storeys.clear();

    std::wstring ifcfile = fileName.toStdWString();
    if ( CGAL::pre_load_ifc(ifcfile) )
    {
        std::vector<std::wstring> storey_names;
        for(uint i=0; i<ifc_storey_list.size(); i++)
            if (ifc_storey_list[i]->m_Name)
                storey_names.push_back( ifc_storey_list[i]->m_Name->m_value );

        IFC_Loader *ifc_loader = new IFC_Loader(0, &storey_names);
//        if(QDialog::Rejected == ifc_loader->exec()){
//            return;
//        }
//        else
//        {
//            if (ifc_loader->done_set)
//            {
//                if (ifc_loader->load_full)
//                {
//                    CGAL::fully_load_ifc2lcc(alcc, map_ifc_entities, ifc_loader->do_simplif, ifc_loader->no_sew2,
//                                             ifc_loader->only_closed_meshes, ifc_loader->correct_spaces);
//                }
//                else
//                {
                    // If all the storeys are selected
                    if (ifc_loader->load_all_storeys)
                    {
                        // If all the objects of the storey are requested
                        if (ifc_loader->on_all_objects)
                        {
                            for(uint i=0; i<ifc_storey_list.size(); i++)
                                CGAL::load_all_in_storey_ifc2lcc(alcc, ifc_storey_list[i],
                                                                 ifc_loader->do_simplif, ifc_loader->no_sew2,
                                                                 ifc_loader->only_closed_meshes, ifc_loader->correct_spaces);
                        }
                        // Otherwise load the picked ones
                        else
                        {
                            for(uint i=0; i<ifc_storey_list.size(); i++)
                                CGAL::load_selected_in_storey_ifc2lcc(alcc, ifc_storey_list[i], ifc_loader->selected_elem, ifc_loader->do_simplif,
                                                                      ifc_loader->no_sew2, ifc_loader->only_closed_meshes, ifc_loader->correct_spaces);
                        }
                    }
                    // Otherwise just get the selected storeys
                    else
                    {
                        assert(ifc_storey_list.size() == ifc_loader->checked_storeys.size());
                        for(uint i=0; i<ifc_storey_list.size(); i++)
                        {
                            if ( ifc_loader->checked_storeys[i] )
                            {
                                if (ifc_loader->on_all_objects)
                                    CGAL::load_all_in_storey_ifc2lcc(alcc, ifc_storey_list[i],
                                                                     ifc_loader->do_simplif, ifc_loader->no_sew2, ifc_loader->correct_spaces);
                                else
                                    CGAL::load_selected_in_storey_ifc2lcc(alcc, ifc_storey_list[i], ifc_loader->selected_elem, ifc_loader->do_simplif,
                                                                          ifc_loader->no_sew2, ifc_loader->only_closed_meshes, ifc_loader->correct_spaces);
                            }
                        }
                    }

                    /// UI form for the spaces
                    if (ifc_loader->load_spaces)
                    {
                        DialogLoadSpaces *dialogLoadSpaces = new DialogLoadSpaces(0);
                        dialogLoadSpaces->init();
//                        std::cout << "size of storey: " << spaces_of_storeys.size() << " and " << ifc_loader->checked_storeys.size() << std::endl;

                        dialogLoadSpaces->Set_space_list( spaces_of_storeys, ifc_loader->checked_storeys );

                        // If no space were found
                        if (dialogLoadSpaces->IFCSpaceLoader_ui->list_of_spaces->count() == 0)
                        {
                            QMessageBox msgBox;
                            msgBox.setText("No IfcSpace entities found in the model.");
                            if(t.is_running())
                                t.stop();
                            msgBox.exec();
                        }

                        else
                        {
                            if(t.is_running())
                                t.stop();
                            if(QDialog::Rejected == dialogLoadSpaces->exec()){
                                return;
                            }
                            t.start();

                            // Collect the openings to load, to avoid duplication (openings shared by several spaces)
                            std::map<int, shared_ptr<IfcElement> > map_op_elem_id;
                            // To know which opening belongs to which space(s)
                            std::map<int, std::vector<std::string> > map_op_to_space;

                            int totalOps = 0;
                            // The space(s)
                            for(uint i=0; i<dialogLoadSpaces->IFCSpaceLoader_ui->list_of_spaces->count(); i++)
                            {
                                if ( dialogLoadSpaces->IFCSpaceLoader_ui->list_of_spaces->item(i)->checkState() == Qt::Checked )
                                {
                                    CGAL::vec_dart space_res;
                                    space_res = CGAL::Get_Geometry_from_IfcObject(alcc, dialogLoadSpaces->space_set[i], map_ifc_shapes, false, ifc_loader->no_sew2,
                                                                                  ifc_loader->only_closed_meshes, ifc_loader->correct_spaces);

                                    /// Only deal with the single solid room volumes for now
                                    if (space_res.size() == 1)
                                    {
                                        if (alcc.attribute<3>(space_res[0]) != LCC::null_handle)
                                        {
                                            std::string vol_id = alcc.info<3>( space_res[0] ).id();
                                            dialogLoadSpaces->Set_space_uid_map( vol_id, space_res[0], dialogLoadSpaces->space_set[i] );

                                            // Building elements inside the space(s)
                                            if ( dialogLoadSpaces->IFCSpaceLoader_ui->inner_space_elem->isChecked() )
                                            {
                                                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char > elem_of_space;
                                                std::map< const char*, std::vector< shared_ptr<IfcProduct> >, LCCtools::cmp_const_char >::iterator itr;
                                                CGAL::Get_all_elements_in_IfcSpace(dialogLoadSpaces->space_set[i], elem_of_space);

                                                if (dialogLoadSpaces->simplif_elem)
                                                {
                                                    CGAL::vec_dart simple_inner_elems_darts;
                                                    for(itr=elem_of_space.begin(); itr!=elem_of_space.end(); itr++)
                                                    {
                                                        for(uint j=0; j<itr->second.size(); j++)
                                                        {
                                                            simple_inner_elems_darts = CGAL::Get_Simplified_Geometry_from_IfcObject(alcc, itr->second[j],
                                                                                                                                    map_ifc_shapes, dialogLoadSpaces->simplif_method);
                                                            // Insert all the inner elements (one dart per elem)
                                                            dialogLoadSpaces->space_set_uid[vol_id].contained_elem.insert
                                                                    ( dialogLoadSpaces->space_set_uid[vol_id].contained_elem.end(),
                                                                      simple_inner_elems_darts.begin(), simple_inner_elems_darts.end());

                                                            // Force contact between the objects and the space containing them
                                                            for(uint k=0; k<simple_inner_elems_darts.size(); k++)
                                                                CGAL::Make_contact_between_volA_and_volB(alcc, space_res[0], simple_inner_elems_darts[k] );
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    for(itr=elem_of_space.begin(); itr!=elem_of_space.end(); itr++)
                                                    {
                                                        for(uint j=0; j<itr->second.size(); j++)
                                                            CGAL::Get_Geometry_from_IfcObject(alcc, itr->second[j], map_ifc_shapes, true,
                                                                                              ifc_loader->no_sew2, ifc_loader->only_closed_meshes,
                                                                                              ifc_loader->correct_spaces);
                                                    }
                                                }
                                            }


                                            // Building elements/openings surrounding the space(s)
                                            if ( dialogLoadSpaces->IFCSpaceLoader_ui->surrounding_space_elem->isChecked()
                                                 || dialogLoadSpaces->IFCSpaceLoader_ui->void_of_openings->isChecked() )
                                            {
                                                std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char > elem_around_space;
                                                std::map< const char*, std::vector< shared_ptr<IfcElement> >, LCCtools::cmp_const_char >::iterator itr;
                                                CGAL::Get_IfcElements_around_IfcSpace(dialogLoadSpaces->space_set[i], elem_around_space);

                                                if (dialogLoadSpaces->IFCSpaceLoader_ui->surrounding_space_elem->isChecked())
                                                {
                                                    for(itr=elem_around_space.begin(); itr!=elem_around_space.end(); itr++)
                                                    {
                                                        for(uint j=0; j<itr->second.size(); j++)
                                                            CGAL::Get_Geometry_from_IfcObject(alcc, itr->second[j], map_ifc_shapes, ifc_loader->do_simplif,
                                                                                              ifc_loader->no_sew2, ifc_loader->only_closed_meshes,
                                                                                              ifc_loader->correct_spaces);
                                                    }
                                                }

                                                if (dialogLoadSpaces->IFCSpaceLoader_ui->void_of_openings->isChecked())
                                                {
                                                    std::vector< shared_ptr<IfcElement> > op_to_load;
                                                    CGAL::Get_Openings_around_IfcSpace( dialogLoadSpaces->space_set[i], op_to_load, elem_around_space );

                                                    totalOps += op_to_load.size();

                                                    for(uint j=0; j<op_to_load.size(); j++){
                                                        map_op_elem_id[ op_to_load[j]->m_entity_id ] = op_to_load[j];
                                                        map_op_to_space[ op_to_load[j]->m_entity_id ].push_back(vol_id);
                                                    }
                                                }
                                            }
                                        }
                                    }
//                                    else
//                                        std::cout << " @@@@@@@@@@@@@@ THERE ARE MORE THAN ONE (OR NOT AT ALL) VOLUME FOR ONE SPACE OBJECT!!! " << std::endl;
                                }
                            }

                            if( totalOps == 0 ){
                                QMessageBox msgBox;
                                msgBox.setText("No opening (door/window) detected. IfcRelSpaceBoundary relationships are probably missing in your IFC file.");
                                if (t.is_running())
                                    t.stop();
                                msgBox.exec();
                            }
                            else{
                                //Load voids of openings if selected
                                if (dialogLoadSpaces->IFCSpaceLoader_ui->void_of_openings->isChecked()){
                                    CGAL::vec_dart d_ops;
    //                                std::map<int,std::vector<std::string> > linkedByVirtualElem;

                                    std::map<int, shared_ptr<IfcElement> >::iterator it_ops( map_op_elem_id.begin() );
                                    for(; it_ops != map_op_elem_id.end(); it_ops++)
                                    {
    //                                    std::cout << "Trying to extract the geometry of an " << it_ops->second->className() << std::endl;
                                        d_ops = CGAL::Get_Geometry_from_IfcObject(alcc, it_ops->second, map_ifc_shapes);
                                        if(d_ops.size() > 0 && d_ops[0] != LCC::null_handle)
                                        {
    //                                        std::cout << "--- Got some geometry of opening..." << std::endl;
                                            // store the darts of the openings to their corresponding space structure
                                            for(uint i=0; i< map_op_to_space[ it_ops->first ].size(); i++){
                                                dialogLoadSpaces->space_set_uid[ map_op_to_space[ it_ops->first ][i] ].surrounding_ops.push_back( d_ops[0] );

                                                // Add the relation to the attributes of the 3-cells
                                                Dart_handle d_space = dialogLoadSpaces->space_set_uid[ map_op_to_space[ it_ops->first ][i] ].dart;
                                                alcc.info<3>(d_space).addRelatedCell( alcc.info<3>(d_ops[0]).id() );
                                                alcc.info<3>(d_ops[0]).addRelatedCell( alcc.info<3>(d_space).id() );
                                            }
                                        }
                                        else if( d_ops.size() > 1 )
                                            std::cout << " @@@@@@@@@@@@@@ THERE ARE MORE THAN ONE VOLUME FOR ONE OPENING OBJECT!!! "<< std::endl;
                                        else{
    //                                        std::cout << "--- Could not get a THANG!" << std::endl;
                                            if ( strcmp( it_ops->second->className(), "IfcVirtualElement" ) == 0 ){
    //                                            Dart_handle d_space = dialogLoadSpaces->space_set_uid[ map_op_to_space[ it_ops->first ][i] ].dart;
    //                                            linkedByVirtualElem[ it_ops->first ].push_back( alcc.info<3>(d_space).id() );

    //                                            for(uint i=0; i< map_op_to_space[ it_ops->first ].size(); i++){}

                                                if( map_op_to_space[ it_ops->first ].size() == 2 ){
                                                    Dart_handle d1 = dialogLoadSpaces->space_set_uid[ map_op_to_space[ it_ops->first ][0] ].dart,
                                                                d2 = dialogLoadSpaces->space_set_uid[ map_op_to_space[ it_ops->first ][1] ].dart;
                                                    alcc.info<3>(d1).addRelatedCell( alcc.info<3>(d2).id() );
                                                    alcc.info<3>(d2).addRelatedCell( alcc.info<3>(d1).id() );
                                                }
                                                std::cout << "--- one virtualElement for " << map_op_to_space[ it_ops->first ].size() << " 3cells." << std::endl;
                                            }
                                        }
                                    }

    //                                for( auto &it : linkedByVirtualElem ){
    //                                    if( it.second.size() == 2 ){
    //                                        Dart_handle d1 = dialogLoadSpaces->space_set_uid[ it.second[0] ].dart,
    //                                                    d2 = dialogLoadSpaces->space_set_uid[ it.second[1] ].dart;
    //                                        alcc.info<3>(d1).addRelatedCell( alcc.info<3>(d2).id() );
    //                                        alcc.info<3>(d2).addRelatedCell( alcc.info<3>(d1).id() );
    //                                    }
    //                                    std::cout << "--- one virtualElement for " << it.second.size() << " 3cells." << std::endl;
    //                                }
                                }
                            }
                        }
                    }
//                }
//            }
//        }

        map_ifc_shapes.clear();
        map_ifc_entities.clear();

//        CGAL::Bbox_3 bb;
//        for (LCC::One_dart_per_cell_range<0>::iterator
//             it(alcc.one_dart_per_cell<0>().begin());
//             it.cont(); ++it){
//            bb += alcc.point(it).bbox();
//        }

        for (LCC::Dart_range::iterator itD=alcc.darts().begin();
             itD!=alcc.darts().end(); ++itD){
            if( alcc.attribute<0>(itD) == LCC::null_handle ){
                std::cout << "Null dart detected!" << std::endl;
                alcc.erase_dart(itD);
            }
        }

//        std::cout << "Bbox extent of loaded model: \n"
//                  << "xmin= " << bb.xmin() << ", xmax= " << bb.xmax() << "\n"
//                  << "ymin= " << bb.ymin() << ", ymax= " << bb.ymax() << "\n"
//                  << "zmin= " << bb.zmin() << ", zmax= " << bb.zmax() << std::endl;
        if(t.is_running())
            t.stop();
        std::cout << "(real) Time to load IFC model: " << t.time() << std::endl;
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setText("Something seems wrong with the loaded IFC file... Sorry try another!");
        if(t.is_running())
            t.stop();
        msgBox.exec();
    }

//    //            CGAL::Correct_non_simple_poly_global(alcc);
//    Q_EMIT(mw->sceneChanged());

#else
    {
        QMessageBox msgBox;
        msgBox.setText("You don't seem to have IFC++ installed... So can't import your file!");
        msgBox.exec();
//        Q_EMIT(mw->sceneChanged());
    }
#endif
}


//#include "IFC2LCC.moc"
