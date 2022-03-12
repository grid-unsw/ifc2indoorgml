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
//                 Sylvain Brandel <sylvain.brandel@liris.cnrs.fr>
//                 Abdoulaye Diakite <diakite.abdoulaye@gmail.com>
//
#include "MainWindow.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include <QSettings>
#include <QHeaderView>
#include <CGAL/Timer.h>
#include <CGAL/ipower.h>

#include "modules/ifc2lcc/IFC2LCC.h"
#include "modules/lcc2indoorgml/LCC2IndoorGML.h"

#define DELAY_STATUSMSG 1500

MainWindow::MainWindow (QWidget * parent) : CGAL::Qt::DemosMainWindow (parent)
{
    setupUi (this);
    scene.lcc = new LCC;

    volumeListDock = new QDockWidget(QString(tr("Volume List")),this);
    volumeListDock->setAllowedAreas(Qt::RightDockWidgetArea |
                                    Qt::LeftDockWidgetArea);
    volumeList = new QTableWidget(0,4,volumeListDock);
    volumeList->verticalHeader()->hide();
    volumeList->setColumnHidden(3,true);
    QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                     this, SLOT(onCellChanged(int,int)));

    QStringList labels(QString(tr("CellSpace")));
    labels.append(QString(tr("Filled")));
    labels.append(QString(tr("Hidden")));
    volumeList->setHorizontalHeaderLabels(labels);
//    volumeList->resizeColumnsToContents();
//    volumeList->setFixedWidth(220);
//    volumeList->setDragEnabled(true);
//    volumeList->setMaximumWidth(300);
    /*  volumeList->setColumnWidth(0,85);
  volumeList->setColumnWidth(1,35);
  volumeList->setColumnWidth(2,35);*/

    volumeList->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    volumeList->setSelectionMode(QAbstractItemView::NoSelection);
    //volumeList->setSelectionBehavior(QAbstractItemView::SelectRows);
    volumeListDock->setWidget(volumeList);
    addDockWidget(Qt::RightDockWidgetArea,volumeListDock);
    menuView->addAction(volumeListDock->toggleViewAction());

    this->viewer->setScene(&scene, false);

    connect_actions ();
    this->addAboutDemo (":/cgal/help/about_Linear_cell_complex_3.html");
    this->addAboutCGAL ();

    this->addRecentFiles (this->menuFile, this->actionQuit);
    connect (this, SIGNAL (openRecentFile (QString)),
             this, SLOT (load_depend_on_extension(QString)));

    statusMessage = new QLabel
            ("Darts: 0,  Vertices: 0  (Points: 0),  Edges: 0, Facets: 0,"
             " Volume: 0 (Vol color: 0),  Connected components: 0");
    statusBar ()->addWidget (statusMessage);

}

void MainWindow::connect_actions ()
{
    QObject::connect (this->actionQuit, SIGNAL (triggered ()),
                      qApp, SLOT (quit ()));

    QObject::connect (this, SIGNAL (sceneChanged ()),
                      this, SLOT (onSceneChanged ()));

    QObject::connect(this->volumeList->horizontalHeader(),
                     SIGNAL(sectionClicked(int)),
                     this, SLOT(onHeaderClicked(int)));
}

void MainWindow::connectVolumeListHandlers()
{
    QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                     this, SLOT(onCellChanged(int,int)));
}

void MainWindow::update_operations_entries(bool show)
{
    actionClear->setEnabled(show);
    menuOperations->setEnabled(show);
}

void MainWindow::onSceneChanged ()
{
    QApplication::setOverrideCursor( Qt::WaitCursor );

    LCC::size_type mark = scene.lcc->get_new_mark ();
    scene.lcc->negate_mark (mark);

    std::vector<unsigned int> cells;
    cells.push_back(0);
    cells.push_back(1);
    cells.push_back(2);
    cells.push_back(3);
    cells.push_back(4);

    std::vector<unsigned int> res = scene.lcc->count_cells (cells);

    std::ostringstream os;
    os << "Darts: " << scene.lcc->number_of_darts ()
       << ",  Vertices:" << res[0]
       <<",  (Points:"<<scene.lcc->number_of_attributes<0>()<<")"
      << ",  Edges:" << res[1]
      << ",  Facets:" << res[2]
      << ",  Volumes:" << res[3]
      <<",  (Vol color:"<<scene.lcc->number_of_attributes<3>()<<")"
     << ",  Connected components:" << res[4]
     <<",  Valid:"<<(scene.lcc->is_valid()?"true":"FALSE");

    scene.lcc->negate_mark (mark);
    scene.lcc->free_mark (mark);

    // statusBar()->showMessage (QString ("Update OpenGL lists"), DELAY_STATUSMSG);

    viewer->sceneChanged ();

    statusMessage->setText (os.str().c_str ());
    QApplication::restoreOverrideCursor();
}

void MainWindow::clear_all()
{
    scene.lcc->clear();
    viewer->clear_layer_maps();

    volumeList->clearContents();
    volumeList->setRowCount(0);
}

void MainWindow::on_new_volume(Dart_handle adart)
{
    CGAL_assertion( scene.lcc->attribute<3>(adart)==LCC::null_handle);
    scene.lcc->set_attribute<3>(adart, scene.lcc->create_attribute<3>());
    update_volume_list_add(scene.lcc->attribute<3>(adart));
}

void MainWindow::init_all_new_volumes()
{
    for (LCC::One_dart_per_cell_range<3>::iterator
         it(scene.lcc->one_dart_per_cell<3>().begin());
         it.cont(); ++it)
        if ( scene.lcc->attribute<3>(it)==LCC::null_handle )
        { on_new_volume(it); }
}

void MainWindow::init_ifc_volumes()
{
    for (LCC::One_dart_per_cell_range<3>::iterator
         it(scene.lcc->one_dart_per_cell<3>().begin());
         it.cont(); ++it){

        if ( scene.lcc->attribute<3>(it)==LCC::null_handle )
            on_new_volume(it);
        else
            update_volume_list_add(scene.lcc->attribute<3>(it));
    }
}

void MainWindow::on_actionSave_triggered ()
{
    QString fileName = QFileDialog::getSaveFileName (this,
                                                     tr ("Save 3D Indoor File"),
                                                     "myIndoorModel.3map",
                                                     tr ("3-map files (*.3map)"));

    if (!fileName.isEmpty ())
    {
        save(fileName);
    }
}

void MainWindow::on_actionLoad_triggered ()
{
    QString fileName = QFileDialog::getOpenFileName (this,
                                                     tr ("Open a 3D Indoor file"),
                                                     "",
                                                     tr ("3D model (*.3map *.ifc)"));

    if (!fileName.isEmpty ())
    {
        load_depend_on_extension(fileName);
//        load(fileName, false);
    }
}


void MainWindow::load_depend_on_extension(const QString & fileName, bool clear)
{
    QString ext = QFileInfo(fileName).suffix();
    if ( ext=="3map")
    {
        load(fileName, clear);
    }
    else if (ext=="off")
    {
        load_off(fileName, clear);
    }
    else if (ext=="ifc"){
        load_ifc (fileName, clear);
    }
    else
    {
        std::cout<<"Extension not considered."<<std::endl;
    }
}

void MainWindow::open(QString fileName){
    load_depend_on_extension(fileName, true);
}

void MainWindow::load(const QString & fileName, bool clear)
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

    if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    bool res = load_combinatorial_map(fileName.toStdString().c_str(), *(scene.lcc));

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to load 3-map "<<qPrintable(fileName)<<": "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    init_all_new_volumes();
    recreate_whole_volume_list();

    this->addToRecentFiles(fileName);
    QApplication::restoreOverrideCursor ();

    if (res)
        statusBar ()->showMessage (QString ("3-map loaded ") + fileName,
                                   DELAY_STATUSMSG);
    else
        statusBar ()->showMessage (QString ("Problem: 3-map not loaded ") + fileName,
                                   DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
}

void MainWindow::save(const QString & fileName)
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    if ( save_combinatorial_map(*(scene.lcc), fileName.toStdString().c_str()) )
        statusBar ()->showMessage (QString ("3-map saved ") + fileName,
                                   DELAY_STATUSMSG);
    else
        statusBar ()->showMessage (QString ("Problem: 3-map not saved ") + fileName,
                                   DELAY_STATUSMSG);
    QApplication::restoreOverrideCursor ();

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to save 3-map "<<qPrintable(fileName)<<": "
            <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::load_off (const QString & fileName, bool clear)
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

    if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    std::ifstream ifs (qPrintable (fileName));

    CGAL::import_from_polyhedron_3_flux < LCC > (*scene.lcc, ifs);

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to load off "<<qPrintable(fileName)<<": "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    init_all_new_volumes();
    recreate_whole_volume_list();

    this->addToRecentFiles (fileName);
    QApplication::restoreOverrideCursor ();

    if (clear)
        statusBar ()->showMessage (QString ("Loaded off file") + fileName,
                                   DELAY_STATUSMSG);
    else
        statusBar ()->showMessage (QString ("Added off file") + fileName,
                                   DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
}


void MainWindow::load_ifc (const QString & fileName, bool clear)
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

    if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    load_IFC (fileName, *scene.lcc);

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to load ifc "<<qPrintable(fileName)<<": "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    init_ifc_volumes();

    this->addToRecentFiles (fileName);
    QApplication::restoreOverrideCursor ();

    if (clear)
        statusBar ()->showMessage (QString ("Loaded ifc file") + fileName,
                                   DELAY_STATUSMSG);
    else
        statusBar ()->showMessage (QString ("Added ifc file") + fileName,
                                   DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionClear_triggered()
{
    clear_all();
    statusBar ()->showMessage (QString ("Scene cleared"), DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
}


void MainWindow::on_actionSew3_same_facets_triggered()
{
    LCC::size_type mymark = scene.lcc->get_new_mark();
    mark_all_filled_and_visible_volumes(mymark);

    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    if ( scene.lcc->sew3_same_facets(mymark) > 0 )
    {
        statusBar()->showMessage
                (QString ("Same facets of visible and filled volume(s) are 3-sewn"),
                 DELAY_STATUSMSG);
        Q_EMIT (sceneChanged ());
    }
    else
        statusBar()->showMessage (QString ("No facets 3-sewn"), DELAY_STATUSMSG);

    scene.lcc->free_mark(mymark);

    QApplication::restoreOverrideCursor ();

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to sew3 all same facets: "
            <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionUnsew3_all_triggered()
{
    unsigned int nb=0;
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    for (LCC::Dart_range::iterator it=scene.lcc->darts().begin();
         it!=scene.lcc->darts().end(); ++it)
    {
        if ( !scene.lcc->is_free(it,3) &&
             scene.lcc->info<3>(it).is_filled_and_visible() &&
             scene.lcc->info<3>(scene.lcc->beta(it,3))
             .is_filled_and_visible())
        { scene.lcc->unsew<3>(it); ++nb; }
    }

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to unsew3 all filled volumes: "
            <<timer.time()<<" seconds."<<std::endl;
#endif
    QApplication::restoreOverrideCursor ();

    if ( nb > 0 )
    {
        statusBar()->showMessage
                (QString ("Darts between visible and filled volume(s) are 3-unsewn"),
                 DELAY_STATUSMSG);
        Q_EMIT (sceneChanged ());
    }
    else
        statusBar()->showMessage (QString ("No dart 3-unsewn"), DELAY_STATUSMSG);
}

void MainWindow::on_actionInsideOut_triggered()
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    LCC::size_type mymark=scene.lcc->get_new_mark();

    for (LCC::Attribute_range<3>::type::iterator
         it=scene.lcc->attributes<3>().begin(),
         itend=scene.lcc->attributes<3>().end(); it!=itend; )
    {
        LCC::Attribute_handle<3>::type cur = it++;
        if( !scene.lcc->is_marked(scene.lcc->get_attribute<3>(cur).dart(), mymark) &&
                scene.lcc->get_attribute<3>(cur).info().is_filled_and_visible() )
        {
            scene.lcc->reverse_orientation_connected_component
                    (scene.lcc->get_attribute<3>(cur).dart(), mymark);
        }
    }

    // unmark all the darts by iterating on all the darts
    // but we cannot do really better
    scene.lcc->free_mark(mymark);

    QApplication::restoreOverrideCursor ();
    Q_EMIT( sceneChanged());

    statusBar()->showMessage
            (QString("Orientation of visible and filled volume(s) reversed"),
             DELAY_STATUSMSG);

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to reverse the orientation of all filled volumes: "
            <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionRemove_filled_volumes_triggered()
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    unsigned int count = 0;
    for (LCC::Attribute_range<3>::type::iterator
         it=scene.lcc->attributes<3>().begin(),
         itend=scene.lcc->attributes<3>().end(); it!=itend; )
    {
        LCC::Attribute_handle<3>::type cur = it++;
        if( scene.lcc->get_attribute<3>(cur).info().is_filled_and_visible() )
        {
            scene.lcc->remove_cell<3>(scene.lcc->get_attribute<3>(cur).dart());
            ++count;
        }
    }

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to remove all filled volumes: "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    recreate_whole_volume_list();
    QApplication::restoreOverrideCursor ();
    Q_EMIT( sceneChanged());

    statusBar()->showMessage
            (QString::number(count)+QString("Visible and filled volume(s) removed"),
             DELAY_STATUSMSG);
}


double compute_angle3d(const Vector_3& v1, const Vector_3& v2)
{
    double a = CGAL::to_double( (v1*v2) /
                                ( sqrt(v1.squared_length()) * sqrt(v2.squared_length()) ) ) ;

    if (a < -1.0) return acos(-1.0)/CGAL_PI*180.0;
    else if (a > 1.0) return acos(1.0)/CGAL_PI*180.0;
    else return acos(a)/CGAL_PI*180.0;
}

void MainWindow::on_actionMerge_coplanar_faces_triggered()
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    scene.lcc->set_update_attributes(false);

    std::vector<Dart_handle> edges;
    LCC::size_type treated  = scene.lcc->get_new_mark();
    LCC::size_type treated2 = scene.lcc->get_new_mark();

    for ( LCC::Dart_range::iterator it= scene.lcc->darts().begin(),
          itend = scene.lcc->darts().end(); it!=itend; ++it )
    {
        if (!scene.lcc->is_marked(it, treated) )
        {
            if ( scene.lcc->is_removable<1>(it) )
            {
                LCC::Vector normal1 = CGAL::compute_normal_of_cell_2(*scene.lcc,it);
                LCC::Vector normal2 = CGAL::compute_normal_of_cell_2(*scene.lcc, scene.lcc->beta<2>(it) );
                double angle = compute_angle3d(normal1, normal2);

                if ( ((angle<5.0 || angle>355.0) || (angle<185.0 && angle>175.0)) )
                {
                    edges.push_back(it);
                }
            }
            CGAL::mark_cell<LCC, 1>(*scene.lcc, it, treated);
        }
    }


    for (std::vector<Dart_handle>::iterator it=edges.begin(),
         itend=edges.end(); it!=itend; ++it)
    {
        CGAL::mark_cell<LCC, 1>(*scene.lcc, *it, treated2);

        if ( scene.lcc->beta<0, 2>(*it)==*it || scene.lcc->beta<1, 2>(*it)==*it)
        { // To process dangling edges

            Dart_handle actu = *it, prev=nullptr;
            do
            {
                if ( scene.lcc->beta<0, 2>(actu)==actu ) prev = scene.lcc->beta<1>(actu);
                else prev = scene.lcc->beta<0>(actu);

                if (scene.lcc->is_marked(actu, treated2) &&
                        (scene.lcc->beta<0, 2>(actu)!=actu || scene.lcc->beta<1, 2>(actu)!=actu) )
                {
                    scene.lcc->remove_cell<1>(actu);
                    actu = prev;
                }
                else
                    actu = nullptr;
            }
            while (actu!=nullptr && (scene.lcc->beta<0, 2>(actu)==actu || scene.lcc->beta<1, 2>(actu)==actu));
        }
        else if ( !CGAL::belong_to_same_cell<LCC, 2>(*scene.lcc, *it,
                                                     scene.lcc->beta<2>(*it)) )
            scene.lcc->remove_cell<1>(*it);
    }

    assert(scene.lcc->is_whole_map_marked(treated));
    scene.lcc->free_mark(treated);
    scene.lcc->free_mark(treated2);

    scene.lcc->set_update_attributes(true);

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to merge all coplanar faces: "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    recreate_whole_volume_list();

    QApplication::restoreOverrideCursor ();
    Q_EMIT (sceneChanged ());
    statusBar()->showMessage
            (QString ("Coplanar face(s) merged"), DELAY_STATUSMSG);
}

bool is_external(CDT::Face_handle fh)
{
    return fh->info().is_external;
}

int number_of_existing_edge(CDT::Face_handle fh)
{
    unsigned res=0;
    for (int i=0; i<3; ++i)
        if (fh->info().exist_edge[i]) ++res;
    return res;
}

int get_free_edge(CDT::Face_handle fh)
{
    CGAL_assertion( number_of_existing_edge(fh)==2 );
    for (int i=0; i<3; ++i)
        if (!fh->info().exist_edge[i]) return i;

    CGAL_assertion(false);
    return -1;
}

void constrained_delaunay_triangulation(LCC &lcc, Dart_handle d1)
{
    Vector_3 normal = CGAL::compute_normal_of_cell_2(lcc,d1);
    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    //inserting the constraints edge by edge
    LCC::Dart_of_orbit_range<1>::iterator
            it(lcc.darts_of_orbit<1>(d1).begin());

    CDT::Vertex_handle previous=LCC::null_handle, first=LCC::null_handle,
            vh=LCC::null_handle;

    for (LCC::Dart_of_orbit_range<1>::iterator
         itend(lcc.darts_of_orbit<1>(d1).end()); it!=itend; ++it)
    {
        vh = cdt.insert(lcc.point(it));
        vh->info().dh=it;
        if( first==nullptr )
        {
            first=vh;
        }
        if( previous!=nullptr)
        {
            CGAL_assertion( previous !=vh );
            cdt.insert_constraint(previous,vh);
        }

        previous=vh;
    }
    cdt.insert_constraint(previous,first);
    CGAL_assertion(cdt.is_valid());

    // sets mark is_external
    for( CDT::All_faces_iterator fit = cdt.all_faces_begin(),
         fitend = cdt.all_faces_end(); fit != fitend; ++fit)
    {
        fit->info().is_external = true;
        fit->info().is_process = false;
        fit->info().exist_edge[0]=false;
        fit->info().exist_edge[1]=false;
        fit->info().exist_edge[2]=false;
    }

    std::queue<CDT::Face_handle> face_queue;
    CDT::Face_handle face_internal = nullptr;

    face_queue.push(cdt.infinite_vertex()->face());
    while(! face_queue.empty() )
    {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(!fh->info().is_process)
        {
            fh->info().is_process = true;
            for(int i = 0; i <3; ++i)
            {
                if(!cdt.is_constrained(std::make_pair(fh, i)))
                {
                    face_queue.push(fh->neighbor(i));
                }
                else if (face_internal==nullptr)
                {
                    face_internal = fh->neighbor(i);
                }
            }
        }
    }
    if ( face_internal!=nullptr )
        face_queue.push(face_internal);

    while(! face_queue.empty() )
    {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(!fh->info().is_process)
        {
            fh->info().is_process = true;
            fh->info().is_external = false;
            for(int i = 0; i <3; ++i)
            {
                if(!cdt.is_constrained(std::make_pair(fh, i)))
                {
                    face_queue.push(fh->neighbor(i));
                }
            }
        }
    }

    for( CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(),
         eitend = cdt.finite_edges_end(); eit != eitend; ++eit)
    {
        CDT::Face_handle fh = eit->first;
        int index = eit->second;
        CDT::Face_handle opposite_fh = fh->neighbor(index);
        if(cdt.is_constrained(std::make_pair(fh, index)))
        {
            fh->info().exist_edge[index]=true;
            opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

            if ( !fh->info().is_external && number_of_existing_edge(fh)==2 )
                face_queue.push(fh);
            if ( !opposite_fh->info().is_external &&
                 number_of_existing_edge(opposite_fh)==2 )
                face_queue.push(opposite_fh);
        }
    }

    while( !face_queue.empty() )
    {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        CGAL_assertion( number_of_existing_edge(fh)>=2 ); // i.e. ==2 or ==3
        CGAL_assertion( !fh->info().is_external );

        if (number_of_existing_edge(fh)==2)
        {
            int index = get_free_edge(fh);
            CDT::Face_handle opposite_fh = fh->neighbor(index);

            CGAL_assertion( !fh->info().exist_edge[index] );
            CGAL_assertion( !opposite_fh->info().
                            exist_edge[cdt.mirror_index(fh,index)] );
            // triangle is (vc, vb, va)
            const CDT::Vertex_handle va = fh->vertex(cdt. cw(index));
            const CDT::Vertex_handle vb = fh->vertex(cdt.ccw(index));
            const CDT::Vertex_handle vc = fh->vertex(index);

            Dart_handle dd1 = nullptr;
            for (LCC::Dart_of_cell_range<0, 2>::iterator it(lcc.darts_of_cell<0, 2>(va->info().dh).begin());
                 dd1==nullptr && it.cont(); ++it)
            {
                if (lcc.point(lcc.beta<1>(it))==vc->point())
                    dd1=it;
            }

            Dart_handle dd2 = nullptr;
            for (LCC::Dart_of_cell_range<0, 2>::iterator it(lcc.darts_of_cell<0, 2>(vb->info().dh).begin());
                 dd2==nullptr && it.cont(); ++it)
            {
                if (lcc.point(lcc.beta<0>(it))==vc->point())
                    dd2=it;
            }

            //       assert(((lcc.beta<0,0>(dd1)==dd2) || lcc.beta<1,1>(dd1)==dd2));

            Dart_handle ndart=lcc.insert_cell_1_in_cell_2(dd1, dd2);
            va->info().dh=lcc.beta<2>(ndart);

            fh->info().exist_edge[index]=true;
            opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

            if ( !opposite_fh->info().is_external &&
                 number_of_existing_edge(opposite_fh)==2 )
                face_queue.push(opposite_fh);
        }
    }
}

void MainWindow::on_actionTriangulate_all_facets_triggered()
{
    QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    std::vector<LCC::Dart_handle> v;
    for (LCC::One_dart_per_cell_range<2>::iterator
         it(scene.lcc->one_dart_per_cell<2>().begin()); it.cont(); ++it)
    {
        if ( scene.lcc->info<3>(it).is_filled_and_visible() ||
             (!scene.lcc->is_free<3>(it) &&
              scene.lcc->info<3>(scene.lcc->beta<3>(it)).is_filled_and_visible()) )
            v.push_back(it);
    }

    for (std::vector<LCC::Dart_handle>::iterator itv(v.begin());
         itv!=v.end(); ++itv)
        constrained_delaunay_triangulation(*scene.lcc, *itv);

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to triangulate all filled faces: "
            <<timer.time()<<" seconds."<<std::endl;
#endif

    recreate_whole_volume_list();

    QApplication::restoreOverrideCursor ();
    Q_EMIT (sceneChanged ());
    statusBar()->showMessage
            (QString ("All visible and filled faces were triangulated"), DELAY_STATUSMSG);
}

bool MainWindow::is_volume_in_list(LCC::Attribute_handle<3>::type ah)
{
    for(int row=0; row < volumeList->rowCount(); ++row)
    {
        LCC::Attribute_type<3>::type* ptr=
                reinterpret_cast<LCC::Attribute_type<3>::type*>
                ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

        if(ptr==&(scene.lcc->get_attribute<3>(ah))) return true;
    }

    return false;
}

void MainWindow::update_volume_list_add(LCC::Attribute_handle<3>::type ah)
{
    // CGAL_assertion( !is_volume_in_list(ah) );

    volumeList->disconnect(this);

    int newRow = volumeList->rowCount();
    volumeList->setRowCount(newRow+1);

    QTableWidgetItem* volumeLabel = new QTableWidgetItem
            (QString((scene.lcc->get_attribute<3>(ah).info().label().c_str())));
    volumeLabel->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    volumeLabel->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    volumeList->setItem(newRow,0,volumeLabel);

    QTableWidgetItem* fillCB = new QTableWidgetItem;
    fillCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
    fillCB->setTextAlignment(Qt::AlignCenter);
    if ( scene.lcc->get_attribute<3>(ah).info().is_filled() )
        fillCB->setCheckState(Qt::Checked);
    else
        fillCB->setCheckState(Qt::Unchecked);
    volumeList->setItem(newRow,1, fillCB);

    QTableWidgetItem* hiddenCB = new QTableWidgetItem();
    hiddenCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
    hiddenCB->setTextAlignment(Qt::AlignCenter);
    if ( scene.lcc->get_attribute<3>(ah).info().is_visible() )
        hiddenCB->setCheckState(Qt::Unchecked);
    else
        hiddenCB->setCheckState(Qt::Checked);
    volumeList->setItem(newRow,2,hiddenCB);

    QTableWidgetItem* attribHandle = new QTableWidgetItem;
    attribHandle->setData
            (Qt::UserRole,
             reinterpret_cast<quintptr>(&scene.lcc->get_attribute<3>(ah)));

    volumeList->setItem(newRow,3,attribHandle);

    connectVolumeListHandlers();
}

void MainWindow::update_volume_list_remove(int i)
{
    CGAL_assertion(i<volumeList->rowCount());
    volumeList->removeRow(i);
}

void MainWindow::update_volume_list_remove(LCC::Attribute_handle<3>::type ah)
{
    for(int row=0; row < volumeList->rowCount(); ++row)
    {
        LCC::Attribute_type<3>::type* ptr=
                reinterpret_cast<LCC::Attribute_type<3>::type*>
                ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

        if(ptr==&scene.lcc->get_attribute<3>(ah))
        {
            update_volume_list_remove(row);
            return;
        }
    }
}

void MainWindow::update_volume_list_all_ckeckstates()
{
    volumeList->disconnect(this);

    for(int row=0; row < volumeList->rowCount(); ++row)
    {
        LCC::Attribute_type<3>::type* ptr=
                reinterpret_cast<LCC::Attribute_type<3>::type*>
                ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

        if ( ptr->info().is_filled() )
            volumeList->item(row,1)->setCheckState(Qt::Checked);
        else
            volumeList->item(row,1)->setCheckState(Qt::Unchecked);

        if ( !ptr->info().is_visible() )
            volumeList->item(row,2)->setCheckState(Qt::Checked);
        else
            volumeList->item(row,2)->setCheckState(Qt::Unchecked);
    }

    connectVolumeListHandlers();
}

void MainWindow::recreate_whole_volume_list()
{
    volumeList->clearContents();
    volumeList->setRowCount(0);

    for (LCC::Attribute_range<3>::type::iterator
         it=scene.lcc->attributes<3>().begin(),
         itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
        update_volume_list_add(it);
}

void MainWindow::onCellChanged(int row, int col)
{
    volumeList->disconnect(this);

    LCC::Attribute_type<3>::type* ptr=
            reinterpret_cast<LCC::Attribute_type<3>::type*>
            ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

    if ( col==1 )
    {
        ptr->info().negate_filled();
    }
    else if ( col==2 )
    {
        ptr->info().negate_visible();
        if ( !ptr->info().is_visible() )
            volumeList->item(row,1)->setFlags
                    (volumeList->item(row,1)->flags()^Qt::ItemIsEnabled);
        else
            volumeList->item(row,1)->setFlags
                    (volumeList->item(row,1)->flags()|Qt::ItemIsEnabled);
    }

    connectVolumeListHandlers();
    Q_EMIT( sceneChanged());
}

void MainWindow::onHeaderClicked(int col)
{
    if(col != 0)
    {
        volumeList->disconnect(this);

        for(int i = 0; i < volumeList->rowCount(); ++i)
        {
            LCC::Attribute_type<3>::type* ptr=
                    reinterpret_cast<LCC::Attribute_type<3>::type*>
                    ( volumeList->item(i,3)->data(Qt::UserRole).value<quintptr>() );

            switch(qApp->keyboardModifiers())
            {
            case(Qt::ShiftModifier):
                if (col==1)
                    ptr->info().set_filled(false);
                else if (col==2)
                {
                    ptr->info().set_visible(true);
                    volumeList->item(i,1)->setFlags
                            (volumeList->item(i,1)->flags()|Qt::ItemIsEnabled);
                }
                volumeList->item(i,col)->setCheckState(Qt::Unchecked);
                break;
            case(Qt::ControlModifier):
                if (col==1)
                    ptr->info().negate_filled();
                else if (col==2)
                {
                    ptr->info().negate_visible();
                    if ( !ptr->info().is_visible() )
                        volumeList->item(i,1)->setFlags
                                (volumeList->item(i,1)->flags()^Qt::ItemIsEnabled);
                    else
                        volumeList->item(i,1)->setFlags
                                (volumeList->item(i,1)->flags()|Qt::ItemIsEnabled);
                }
                volumeList->item(i,col)->
                        setCheckState(volumeList->item(i,col)->checkState() ?
                                          Qt::Unchecked: Qt::Checked);
                break;
            default:
                if (col==1)
                    ptr->info().set_filled(true);
                else if (col==2)
                {
                    if ( ptr->info().is_visible() )
                    {
                        ptr->info().set_visible(false);
                        volumeList->item(i,1)->setFlags
                                (volumeList->item(i,1)->flags()^Qt::ItemIsEnabled);
                    }
                }
                volumeList->item(i,col)->setCheckState(Qt::Checked);
                break;
            }
        }

        connectVolumeListHandlers();
        Q_EMIT( sceneChanged());
    }
}

void MainWindow::mark_all_filled_and_visible_volumes(LCC::size_type amark)
{
    for (LCC::Attribute_range<3>::type::iterator
         it=scene.lcc->attributes<3>().begin(),
         itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
    {
        if ( scene.lcc->get_attribute<3>(it).info().is_filled_and_visible() &&
             !scene.lcc->is_marked(it->dart(), amark) )
            CGAL::mark_cell<LCC,3>(*scene.lcc,
                                   scene.lcc->get_attribute<3>(it).dart(), amark);
    }
}

void MainWindow::on_actionGenerate_IndoorGML_triggered()
{

#ifdef CGAL_PROFILE_LCC_DEMO
    CGAL::Timer timer;
    timer.start();
#endif

    std::cout << "Generating IndoorGML data..." << std::endl;
    std::size_t e = generateNetworkFromLCC(*scene.lcc, this);

    if(e==0){
        QMessageBox msgBox;
        msgBox.setText("No Edge could be reconstructed. There may be an issue with the IfcRelSpaceBoundary relationships in your IFC file.");
        msgBox.exec();
    }

    statusBar ()->showMessage (QString ("Generated IndoorGML Data."),
                                   DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());

#ifdef CGAL_PROFILE_LCC_DEMO
    timer.stop();
    std::cout<<"Time to generate IndoorGML data: "
            <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionReset_scene_view_triggered()
{
    viewer->showEntireScene();
    Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionExport_IndoorGML_v1_triggered()
{
    lcc2igmlHandler.write_IndoorGML1(*scene.lcc, this);
}

void MainWindow::on_actionExport_IndoorGML_v2_triggered()
{
    lcc2igmlHandler.write_IndoorGML2(*scene.lcc, this);
}
#undef DELAY_STATUSMSG


