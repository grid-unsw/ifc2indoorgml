#include "addthematiclayer.h"
#include "ui_addthematiclayer.h"

AddThematicLayer::AddThematicLayer(QWidget *parent, std::string* id,
                                   std::vector<std::string>* existingLayers) :
    QDialog(parent),
    ui(new Ui::AddThematicLayer)
{
    ui->setupUi(this);

    done = false;
    ui->IL_id->setText( QString( (*id).c_str() ) );
//    ui->ILCtableWidget = new QTableWidget(0,2,this);
    ui->ILCtableWidget->setColumnCount(2);

    ui->ILCtableWidget->verticalHeader()->hide();
    ui->ILCtableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    QStringList labels(QString(tr("Current Thematic Layer(s)")));
    labels.append(QString(tr("typeOfTopoExpression")));
    ui->ILCtableWidget->setHorizontalHeaderLabels(labels);

    if ( existingLayers->size() > 0 )
    {
        if ( ui->ILCtableWidget->rowCount() == 1
             && strcmp(ui->ILCtableWidget->item(0,0)->text().toStdString().c_str(), "No exisiting layer") != 0 )
        {
            ui->ILCtableWidget->removeRow(0);
        }

        // Adding the already existing ThematicLayers in the list + topological relationship options
//        int rowCount = ui->ILCtableWidget->rowCount();

//        std::cout << "\nNum of layers: " << existingLayers->size() << std::endl;

        for( int i=0; i<existingLayers->size(); i++ )
        {
//            std::cout << "Layer " << i << ": " << existingLayers->operator[](i) << std::endl;
            QTableWidgetItem* existingLayer = new QTableWidgetItem ( ((*existingLayers)[i]).c_str() );
            // Set the existing layer selectable
            existingLayer->setFlags(existingLayer->flags() | Qt::ItemIsUserCheckable);
            existingLayer->setCheckState(Qt::Unchecked);
            ui->ILCtableWidget->setRowCount(i+1);
            ui->ILCtableWidget->setItem(i, 0, existingLayer);

            // All possible InterLayer topological relationships in a drop-down menu
            QComboBox *topoExpression = new QComboBox( parent );
            QStringList tExp = { "Contains", "Crosses", "Equals", "Overlaps", "Within", "Other"};
            topoExpression->addItems( tExp );
            ui->ILCtableWidget->setCellWidget(i, 1, topoExpression);
//            rowCount++;
        }
    }
    else
    {
        ui->ILCtableWidget->setRowCount( ui->ILCtableWidget->rowCount()+1 );
        ui->ILCtableWidget->setItem(0, 0, new QTableWidgetItem ("No exisiting layer"));
    }

//    ui->ILCtableWidget->resizeColumnsToContents();
}

AddThematicLayer::~AddThematicLayer()
{
    delete ui;
}

void AddThematicLayer::on_buttonBox_accepted()
{
    layerID = ui->IL_id->text().toStdString();
    if ( strcmp(ui->ILCtableWidget->item(0,0)->text().toStdString().c_str(), "No exisiting layer") != 0 )
    {
        int i;
        for( i=0; i<ui->ILCtableWidget->rowCount(); i++)
        {
            QComboBox *topoExpression = new QComboBox (ui->ILCtableWidget->cellWidget(i,1));
            selectedLayers[ ui->ILCtableWidget->item(i,0)->text().toStdString() ]
                    = topoExpression->currentText().toStdString();

//            std::cout << "\tAdded layer -- " << ui->ILCtableWidget->item(i,0)->text().toStdString() << std::endl;
        }

        done = true;
    }
    else
    {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "Create independent layer", "Do you want to create an independent ThemanticLayer? (with currently visible & filled entities in the scene)",
                                    QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            done = true;
        }
    }
}

