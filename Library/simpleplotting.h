//Must include "QT += core charts widgets" in .pro file

#ifdef QT_CORE_LIB
#ifndef SIMPLEPLOTTING_H
#define SIMPLEPLOTTING_H

#include "eigenmatrixtypes.h"

#include <QtWidgets/QApplication>
#include <QtCharts/QChartView> //Qt's plotting library to visualize results
#include <QtCharts/QLineSeries>
#include <QMessageBox>

namespace ContinuumRobotLibrary{

/*! Visualizes an x and y series. The thread of execution is paused while the window is open. */
static void plot(VectorXd x, VectorXd y, QString title = "", QString x_title = "", QString y_title = ""){
    if(x.size()!=y.size()){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Plotting Error");
        msgBox.setText("x and y vectors must be the same size.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    if(x.size()==0){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Plotting Error");
        msgBox.setText("The plotted series are empty.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();
        return;
    }

    int argc = 0;
    char* argv[1];
    char null = '\0';
    argv[0] = &null;
    QApplication a(argc, argv);

    QtCharts::QChart* chart = new QtCharts::QChart();
    chart->setTitle(title);

    QtCharts::QLineSeries* p = new QtCharts::QLineSeries(chart);
    for (int i = 0; i < x.size(); i++)
        p->append(QPointF(x(i), y(i)));
    chart->addSeries(p);

    chart->createDefaultAxes();
    chart->axisX()->setTitleText(x_title);
    chart->axisY()->setTitleText(y_title);
    chart->legend()->hide();

    QtCharts::QChartView* chartView = new QtCharts::QChartView(chart);
    chartView->resize(400,600);
    chartView->show();

    a.exec();
}

}

#endif // SIMPLEPLOTTING_H
#endif // QT_CORE_LIB
