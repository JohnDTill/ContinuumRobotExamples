#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow){
    ui->setupUi(this);
    initializeCharts();
}

MainWindow::~MainWindow(){
    delete ui;
}

void MainWindow::initializeCharts(){
    QtCharts::QChartView* view_front = new QtCharts::QChartView(&charts[0]);
    QtCharts::QChartView* view_top = new QtCharts::QChartView(&charts[1]);
    QtCharts::QChartView* view_side = new QtCharts::QChartView(&charts[2]);
    QtCharts::QChartView* view_ortho = new QtCharts::QChartView(&charts[3]);
    charts[0].setTitle("<b>Front</b>");
    charts[1].setTitle("<b>Top</b>");
    charts[2].setTitle("<b>Side</b>");
    charts[3].setTitle("<b>Orthographic</b>");
    ui->leftPane->addWidget(view_top);
    ui->leftPane->addWidget(view_front);
    ui->rightPane->addWidget(view_ortho);
    ui->rightPane->addWidget(view_side);

    QPen leg_pen(QColor::fromRgb(88,89,91)), ee_pen(QColor::fromRgb(255,130,0));
    leg_pen.setWidthF(1.5);
    ee_pen.setWidthF(2);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 6; j++){
            legs[i][j].setPen(leg_pen);
            charts[i].addSeries(&legs[i][j]);
        }
        ee[i].setPen(ee_pen);
        charts[i].addSeries(&ee[i]);

        charts[i].createDefaultAxes();
        charts[i].legend()->hide();
    }

    charts[0].axisX()->setRange(-0.25, 0.25);
    charts[0].axisY()->setRange(0, 0.6);
    charts[0].axisX()->setTitleText("x (m)");
    charts[0].axisY()->setTitleText("z (m)");

    charts[1].axisX()->setRange(-0.25, 0.25);
    charts[1].axisY()->setRange(-0.25, 0.25);
    charts[1].axisX()->setTitleText("x (m)");
    charts[1].axisY()->setTitleText("y (m)");

    charts[2].axisX()->setRange(-0.25, 0.25);
    charts[2].axisY()->setRange(0, 0.6);
    charts[2].axisX()->setTitleText("y (m)");
    charts[2].axisY()->setTitleText("z (m)");

    charts[3].axisX()->setRange(-0.2, 0.2);
    charts[3].axisY()->setRange(0, 0.6);
    charts[3].axisX()->hide();
    charts[3].axisY()->hide();

    refreshCharts();
}

void MainWindow::refreshCharts(){
    QList<QPointF> ee_points[4];

    for(int i = 0; i < 6; i++){
        QList<QPointF> leg_points[4];
        Matrix3Xd leg_data = getLegCenterline(i);
        int N = static_cast<int>(leg_data.cols());

        for(int j = 0; j < N; j++){
            Vector3d p = leg_data.col(j);

            leg_points[0] << QPointF(p.x(), p.z());
            leg_points[1] << QPointF(p.x(), p.y());
            leg_points[2] << QPointF(p.y(), p.z());
            leg_points[3] << QPointF(ortho_x*p, ortho_y*p+0.1);
        }
        for(int j = 0; j < 4; j++){
            legs[j][i].replace(leg_points[j]);
            ee_points[j] << leg_points[j].last();
        }
    }

    for(int i = 0; i < 4; i++){
        ee_points[i] << ee_points[i].first();
        ee[i].replace(ee_points[i]);
    }
}

void MainWindow::keyPressEvent(QKeyEvent* event){
    const double v_incr = 1e-4;
    const double w_incr = 1e-3;
    const int steps = 10;

    Vector3d vE = Vector3d::Zero();
    Vector3d wE = Vector3d::Zero();

    switch(event->key()){
        case Qt::Key_A: vE.x() -= v_incr; break;
        case Qt::Key_D: vE.x() += v_incr; break;
        case Qt::Key_S: vE.y() -= v_incr; break;
        case Qt::Key_W: vE.y() += v_incr; break;
        case Qt::Key_Q: vE.z() -= v_incr; break;
        case Qt::Key_E: vE.z() += v_incr; break;
        case Qt::Key_I: wE.x() -= w_incr; break;
        case Qt::Key_K: wE.x() += w_incr; break;
        case Qt::Key_J: wE.y() -= w_incr; break;
        case Qt::Key_L: wE.y() += w_incr; break;
        case Qt::Key_U: wE.z() -= w_incr; break;
        case Qt::Key_O: wE.z() += w_incr; break;
        case Qt::Key_Space: resetKinematicsSolver();
    }

    try{
        for(int i = 0; i < steps; i++){
            pE += vE;
            RE += hat(wE)*RE;
            L = inverseKinematicsCSG(pE,RE);
        }
    }catch(...){
        resetKinematicsSolver();
    }
    refreshCharts();
}

void MainWindow::resetKinematicsSolver(){
    resetCSG();
    pE = 0.4*Vector3d::UnitZ();
    RE = Matrix3d::Identity();
    L = inverseKinematicsCSG(pE,RE);
}
