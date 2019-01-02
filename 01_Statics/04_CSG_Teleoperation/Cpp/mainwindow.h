#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QKeyEvent>
#include "InverseKinematicsCSG.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow{
    Q_OBJECT

private:
    Ui::MainWindow* ui;
    QtCharts::QChart charts[4];
    QtCharts::QLineSeries legs[4][6];
    QtCharts::QLineSeries ee[4];

    Vector3d pE = 0.4*Vector3d::UnitZ();
    Matrix3d RE = Matrix3d::Identity();
    Vector6d L = inverseKinematicsCSG(pE,RE);

    const RowVector3d ortho_x = (Ry(pi/6)*Rz(pi*5/12)).row(1);
    const RowVector3d ortho_y = (Ry(pi/6)*Rz(pi*5/12)).row(2);

public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private:
    void initializeCharts();
    void refreshCharts();
    void keyPressEvent(QKeyEvent* event);
    void resetKinematicsSolver();
};

#endif // MAINWINDOW_H
