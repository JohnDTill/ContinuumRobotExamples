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
    chart->axes(Qt::Horizontal).back()->setTitleText(x_title);
    chart->axes(Qt::Vertical).back()->setTitleText(y_title);
    chart->legend()->hide();

    QtCharts::QChartView* chartView = new QtCharts::QChartView(chart);
    chartView->resize(400,600);
    chartView->show();

    a.exec();
}

/*! Visualizes multiple x and y series. The thread of execution is paused while the window is open. */
static void plot(std::vector<VectorXd> x, std::vector<VectorXd> y, std::vector<Vector3d> rgb, QString title = "", QString x_title = "", QString y_title = ""){
    if(x.size()!=y.size()){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Plotting Error");
        msgBox.setText("There must be an equal number of x and y series.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    if(rgb.size()!=x.size()){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Plotting Error");
        msgBox.setText("Each series must have an rgb color specification.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    if(x.size()==0){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Plotting Error");
        msgBox.setText("No data series provided.");
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

    for(unsigned long long i = 0; i < x.size(); i++){
        VectorXd& x_i = x[i];
        VectorXd& y_i = y[i];

        if(x_i.size()!=y_i.size()){
            QMessageBox msgBox;
            msgBox.setWindowTitle("Plotting Error");
            msgBox.setText("x and y series at index " + QString::number(i) + " must be the same size.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.setDefaultButton(QMessageBox::Ok);
            msgBox.exec();
            return;
        }
        if(x_i.size()==0){
            QMessageBox msgBox;
            msgBox.setWindowTitle("Plotting Error");
            msgBox.setText("The plotted series at index " + QString::number(i) + " is empty.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.setDefaultButton(QMessageBox::Ok);
            msgBox.exec();
            return;
        }

        QtCharts::QLineSeries* p = new QtCharts::QLineSeries(chart);
        for (int j = 0; j < x_i.size(); j++)
            p->append(QPointF(x_i(j), y_i(j)));
        p->setColor(QColor::fromRgb(rgb[i](0), rgb[i](1), rgb[i](2)));
        chart->addSeries(p);
    }

    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).back()->setTitleText(x_title);
    chart->axes(Qt::Vertical).back()->setTitleText(y_title);
    chart->legend()->hide();

    QtCharts::QChartView* chartView = new QtCharts::QChartView(chart);
    chartView->resize(400,600);
    chartView->show();

    a.exec();
}

}

#endif // SIMPLEPLOTTING_H
#endif // QT_CORE_LIB
