#ifndef ADDTHEMATICLAYER_H
#define ADDTHEMATICLAYER_H

#include <QDialog>
#include <QComboBox>
#include <QTableWidget>
#include <QMessageBox>
#include "typedefs.h"
//#include <IndoorGML2LCC.h>

namespace Ui {
class AddThematicLayer;
}

class AddThematicLayer : public QDialog
{
    Q_OBJECT

public:
    explicit AddThematicLayer(QWidget *parent = nullptr, std::string* id = new std::string(),
                              std::vector<std::string> *existingLayers = new std::vector<std::string>());
//                              IndoorGML2::IndoorFeatures* InFt = new IndoorGML2::IndoorFeatures());
    ~AddThematicLayer();

    bool done;
    std::string layerID;
    std::map<std::string, std::string> selectedLayers;


private Q_SLOTS:
    void on_buttonBox_accepted();

private:
    Ui::AddThematicLayer *ui;
};

#endif // ADDTHEMATICLAYER_H
