#ifndef IFC_LOADER_UI_H
#define IFC_LOADER_UI_H

#include <QDialog>

#include "ui_IFC_loader.h"
#include "ui_load_spaces.h"
#include "IFC_Importer.h"

class IFC_Loader : public QDialog
{
    Q_OBJECT

public:
    explicit IFC_Loader(QWidget *parent=0,
                        std::vector<std::wstring> *storey_names = new std::vector<std::wstring>());
    ~IFC_Loader();
    bool done_set,
         load_full,
         load_all_storeys,
         do_simplif,
         load_spaces, correct_spaces,
         on_all_objects,
         no_sew2,
         only_closed_meshes,
         openings, structurals, columns, roof_slabs, furniture;

    std::vector<bool> checked_storeys;
    std::vector<const char*> selected_elem;


private Q_SLOTS:

    void on_LoadFullModel_toggled(bool checked);

    void on_LoadAllStoreys_toggled(bool checked);

    void on_SelectStrorey_toggled(bool checked);

    void on_LoadAllObjects_toggled(bool checked);

    void on_LoadSelectedObjects_toggled(bool checked);

    void on_StoreyList_itemDoubleClicked(QListWidgetItem *item);

    void on_buttonBox_accepted();

    void on_SimplifyObjects_toggled(bool checked);

    void on_LoadSpaces_toggled(bool checked);

    void on_Nosew2_toggled(bool checked);

    void on_LoadOpenings_toggled(bool checked);

    void on_LoadStructurals_toggled(bool checked);

    void on_LoadColumn_toggled(bool checked);

    void on_LoadRoofs_toggled(bool checked);

    void on_LoadFurniture_toggled(bool checked);

    void on_OnlyClosedMeshes_toggled(bool checked);

    void on_SpaceCorrection_toggled(bool checked);

private:
    Ui::IFC_Loader *IFC_ui;

};


class DialogLoadSpaces
#ifdef IFCPP_ON
        : public QDialog
#endif
{
  Q_OBJECT

private Q_SLOTS:
#ifdef IFCPP_ON
    // For the spaces
    void on_check_uncheck_spaces_toggled(bool);
    void on_inner_space_elem_toggled(bool);
    void on_simplify_elem_toggled(bool);
    void on_aabb_simplif_toggled(bool);
#endif

public:

  std::map< std::string, CGAL::Space, LCCtools::cmp_string> space_set_uid;
  bool initialized, spaces_from_ifc;

  ~DialogLoadSpaces();
  Ui::LoadSpaces *IFCSpaceLoader_ui;


  // Methods specific to IFC++
#ifdef IFCPP_ON
  DialogLoadSpaces(QWidget *parent=0 );

  std::vector< shared_ptr<IfcSpace> > space_set;
  std::vector<bool> checked_spaces;
  bool simplif_elem;
  int simplif_method;

  void Set_space_uid_map(std::string uid, Dart_handle& d, shared_ptr<IfcSpace>& space)
  {
      CGAL::Space new_space;
      new_space.uid = uid;
      new_space.dart = d;
      new_space.to_IfcSpace = space;
      space_set_uid[uid] = new_space;
  }

  void Set_space_list(std::vector< std::vector< shared_ptr<IfcSpace> > >& spaces_of_storeys, std::vector<bool>& checked_storeys)
  {
      assert( checked_storeys.size() == spaces_of_storeys.size() );
      for(uint i=0; i<checked_storeys.size(); i++)
          if (checked_storeys[i])
              space_set.insert( space_set.end(), spaces_of_storeys[i].begin(), spaces_of_storeys[i].end() );

      for(int j=0; j<space_set.size(); j++)
      {
          QString long_name, short_name;
          if ( space_set[j]->m_LongName )
          {
              long_name = QString::fromWCharArray((space_set[j]->m_LongName->m_value).c_str());
              long_name = long_name + " ";
          }
          if ( space_set[j]->m_Name )
              short_name = QString::fromWCharArray((space_set[j]->m_Name->m_value).c_str());

          IFCSpaceLoader_ui->list_of_spaces->addItem( long_name + short_name );
      }

      // Make the list checkable
      QListWidgetItem* item = 0;
      for(uint i = 0; i < IFCSpaceLoader_ui->list_of_spaces->count(); i++)
      {
          item = IFCSpaceLoader_ui->list_of_spaces->item(i);
          item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
          item->setCheckState(Qt::Checked);
          checked_spaces.push_back( false );
      }
      if (IFCSpaceLoader_ui->list_of_spaces->count() > 0)
          IFCSpaceLoader_ui->check_uncheck_spaces->setCheckState(Qt::Checked);
  }

  // Initialize the UI sheet
  void init()
  {
      space_set.clear();
      simplif_method = 2;
      simplif_elem = false;
      initialized = true;
      spaces_from_ifc = true;
      checked_spaces.clear();

      IFCSpaceLoader_ui->list_of_spaces->clear();
      IFCSpaceLoader_ui->obb_simplif->setChecked(true);
      IFCSpaceLoader_ui->aabb_simplif->setChecked(false);
      IFCSpaceLoader_ui->simplify_elem->setChecked(false);
      IFCSpaceLoader_ui->simplif_method_frame->setEnabled(false);
      IFCSpaceLoader_ui->inner_space_elem->setCheckState( Qt::Unchecked );
      IFCSpaceLoader_ui->surrounding_space_elem->setCheckState( Qt::Unchecked );
      IFCSpaceLoader_ui->void_of_openings->setCheckState( Qt::Checked );
  }

#else
  DialogLoadSpaces()
  {
      initialized = false;
      spaces_from_ifc = false;
  }
#endif

  // Initialize the DialogLoadSpace able to handle other than IFC files
  void init(LCC& alcc)
  {
      std::string label;
      space_set_uid.clear();
      typename LCC::Base::One_dart_per_cell_range<3>::iterator
              it = alcc.one_dart_per_cell<3>().begin(),
              itend = alcc.one_dart_per_cell<3>().end();
      for(; it!=itend; it++)
      {
          assert( alcc.attribute<3>(it) != LCC::null_handle );
          label = alcc.info<3>(it).label();
          if(label.length() >= 5 && label.substr( label.length() - 5 ) == "Space")
          {
              CGAL::Space new_space;
              new_space.uid = alcc.info<3>(it).id();
              new_space.dart = it;
              space_set_uid[ new_space.uid ] = new_space;
          }
      }

      initialized = true;
      spaces_from_ifc = false;

      std::cout << "Initialized " << space_set_uid.size() << " spaces." << std::endl;
  }

};

#endif // IFC_LOADER_H
