#include "IFC_loader_ui.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                 UI for loading IFC models                                     ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
IFC_Loader::IFC_Loader(QWidget *parent, std::vector<std::wstring> *storey_names) :
    QDialog(parent),
    IFC_ui(new Ui::IFC_Loader)
{
    IFC_ui->setupUi(this);

//    IFC_ui->LoadFullModel->setEnabled();

    IFC_ui->StoreyList->setEnabled(true);
    IFC_ui->SpecificObjectsFrame->setEnabled(true);


    for(uint i=0; i<storey_names->size(); i++)
    {
        const QString one_name = QString::fromWCharArray( (*storey_names)[i].c_str() );
        IFC_ui->StoreyList->addItem(one_name);
    }

    checked_storeys.clear();

    // Make the list checkable
    QListWidgetItem* item = 0;
    for(uint i = 0; i < IFC_ui->StoreyList->count(); i++)
    {
        item = IFC_ui->StoreyList->item(i);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Checked);
        checked_storeys.push_back( true );
    }

    done_set = false;
    load_full = false;
    load_all_storeys = true;
    do_simplif = false;
    load_spaces = true;
    correct_spaces = true;
    on_all_objects = false;
    no_sew2 = false;
    only_closed_meshes = false;
    openings = false;
    structurals = false;
    columns = false;
    roof_slabs = false;
    furniture = false;
}

IFC_Loader::~IFC_Loader()
{
    delete IFC_ui;
}

void IFC_Loader::on_LoadFullModel_toggled(bool checked)
{
    if (checked)
    {
        load_full = true;
        IFC_ui->SelectStoreysFrame->setDisabled(true);
        IFC_ui->SelectObjectsFrame->setDisabled(true);
    }
    else
    {
        load_full = false;
        IFC_ui->SelectStoreysFrame->setEnabled(true);
        IFC_ui->SelectObjectsFrame->setEnabled(true);
    }
}

void IFC_Loader::on_LoadAllStoreys_toggled(bool checked)
{
    if (checked)
    {
        load_all_storeys = true;
        IFC_ui->StoreyList->setDisabled(true);
        IFC_ui->SelectStorey->setChecked(false);
    }
}

void IFC_Loader::on_SelectStrorey_toggled(bool checked)
{
    if (checked)
    {
        load_all_storeys = false;
        IFC_ui->StoreyList->setEnabled(true);
        IFC_ui->LoadAllStoreys->setChecked(false);
    }
}

void IFC_Loader::on_LoadAllObjects_toggled(bool checked)
{
    if (checked)
    {
        on_all_objects = true;
        IFC_ui->SpecificObjectsFrame->setDisabled(true);
        IFC_ui->LoadSelectedObjects->setChecked(false);

        QListWidgetItem* item = 0;
        for(uint i = 0; i < IFC_ui->StoreyList->count(); i++)
        {
            item = IFC_ui->StoreyList->item(i);
            item->setCheckState(Qt::Unchecked);
            checked_storeys[i] = false;
        }
    }
}

void IFC_Loader::on_LoadSelectedObjects_toggled(bool checked)
{
    if (checked)
    {
        on_all_objects = false;
        IFC_ui->SpecificObjectsFrame->setEnabled(true);
        IFC_ui->LoadAllObjects->setChecked(false);
    }
}

void IFC_Loader::on_StoreyList_itemDoubleClicked(QListWidgetItem *item)
{
    if (item->checkState() == Qt::Unchecked)
        item->setCheckState(Qt::Checked);
    else
        item->setCheckState(Qt::Unchecked);
}

void IFC_Loader::on_buttonBox_accepted()
{
    // Detect selected storeys to load
    if (load_all_storeys)
    {
        // Set all the story as checked
        for(uint i = 0; i < IFC_ui->StoreyList->count(); i++)
            checked_storeys[i] = true;
    }
    else
    {
        QListWidgetItem* item = 0;
        for(uint i = 0; i < IFC_ui->StoreyList->count(); i++)
        {
            item = IFC_ui->StoreyList->item(i);
            if (item->checkState() == Qt::Checked)
                checked_storeys[i] = true;
        }
    }

    done_set = true;

    if (openings)
    {
        selected_elem.push_back( "IfcDoor" );
        selected_elem.push_back( "IfcWindow" );
    }

    if (structurals)
    {
        selected_elem.push_back( "IfcWall" );
        selected_elem.push_back( "IfcSlab" );
        selected_elem.push_back( "IfcCovering" );
        selected_elem.push_back( "IfcWallStandardCase" );
    }

    if (columns)
    {
        selected_elem.push_back( "IfcColumn" );
        selected_elem.push_back( "IfcBeam" );
    }

    if (roof_slabs)
        selected_elem.push_back( "IfcRoof" );

    if (furniture)
        selected_elem.push_back( "IfcFurnishingElement" );
}

void IFC_Loader::on_SimplifyObjects_toggled(bool checked)
{
    do_simplif = checked;
}

void IFC_Loader::on_LoadSpaces_toggled(bool checked)
{
    load_spaces = checked;
    IFC_ui->SpaceCorrection->setEnabled(checked);
}

void IFC_Loader::on_Nosew2_toggled(bool checked)
{
    no_sew2 = checked;
}

void IFC_Loader::on_LoadOpenings_toggled(bool checked)
{
    openings = checked;
}

void IFC_Loader::on_LoadStructurals_toggled(bool checked)
{
    structurals = checked;
}

void IFC_Loader::on_LoadColumn_toggled(bool checked)
{
    columns = checked;
}

void IFC_Loader::on_LoadRoofs_toggled(bool checked)
{
    roof_slabs = checked;
}

void IFC_Loader::on_LoadFurniture_toggled(bool checked)
{
    furniture = checked;
}

void IFC_Loader::on_OnlyClosedMeshes_toggled(bool checked)
{
    only_closed_meshes = checked;
}

void IFC_Loader::on_SpaceCorrection_toggled(bool checked)
{
    correct_spaces = checked;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                 UI for loading IFC Spaces                                     ///
/////////////////////////////////////////////////////////////////////////////////////////////////////



#ifdef IFCPP_ON

DialogLoadSpaces::DialogLoadSpaces(QWidget *parent):
    QDialog(parent),
    IFCSpaceLoader_ui(new Ui::LoadSpaces)
{
    IFCSpaceLoader_ui->setupUi (this);
    initialized = false;
    spaces_from_ifc = false;
}

DialogLoadSpaces::~DialogLoadSpaces()
{
    delete IFCSpaceLoader_ui;
}

void DialogLoadSpaces::on_check_uncheck_spaces_toggled(bool checked)
{
    if (checked)
        for(uint i=0; i<IFCSpaceLoader_ui->list_of_spaces->count(); i++)
            IFCSpaceLoader_ui->list_of_spaces->item(i)->setCheckState( Qt::Checked );
    else
        for(uint i=0; i<IFCSpaceLoader_ui->list_of_spaces->count(); i++)
            IFCSpaceLoader_ui->list_of_spaces->item(i)->setCheckState( Qt::Unchecked );
}

void DialogLoadSpaces::on_inner_space_elem_toggled(bool checked)
{
    IFCSpaceLoader_ui->simplify_elem->setEnabled(checked);

    if (checked && IFCSpaceLoader_ui->simplify_elem->isChecked())
        IFCSpaceLoader_ui->simplif_method_frame->setEnabled(checked);
    else if (!checked)
    {
        simplif_elem = false;
        IFCSpaceLoader_ui->simplif_method_frame->setEnabled(checked);
    }
}

void DialogLoadSpaces::on_simplify_elem_toggled(bool checked)
{
    IFCSpaceLoader_ui->simplif_method_frame->setEnabled(checked);
    simplif_elem = checked;
}

void DialogLoadSpaces::on_aabb_simplif_toggled(bool checked)
{
    IFCSpaceLoader_ui->obb_simplif->setChecked( !checked );
    if (checked)
        simplif_method = 1;
    else
        simplif_method = 2;
}
#endif
