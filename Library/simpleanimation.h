//Must include "QT += core charts widgets" in .pro file

#ifdef QT_CORE_LIB
#ifndef SIMPLEANIMATION_H
#define SIMPLEANIMATION_H

#include "eigenmatrixtypes.h"

#include <QtWidgets/QApplication>
#include <QtCharts/QChartView> //Qt's plotting library to visualize results
#include <QtCharts/QLineSeries>
#include <QMessageBox>
#include <thread>
#include <QTime>

namespace ContinuumRobotLibrary{

//There must be separate threads to update the animation
//and to assume responsibility for the plotting window.
//A thread is created to run this function.
static void workerAnimationFunc(MatrixXd x,
                              MatrixXd y,
                              double dt,
                              double t0,
                              QString title,
                              QtCharts::QLineSeries* p,
                              QtCharts::QChart* c){
    QTime plot_time = QTime::currentTime();
    double t = t0;
    int cols = static_cast<int>(x.cols());
    for(int i = 0; i < x.rows(); i++){
        //Approximately do not exceed the time step
        while(plot_time.msecsTo(QTime::currentTime()) < 1000*dt);
        plot_time = QTime::currentTime();
        t += dt;

        //Plot the series at this time
        QList<QPointF> points;
        points.reserve(cols);
        for(int j = 0; j < cols; j++)
            points << QPointF(x(i,j),y(i,j));
        p->replace(points);
        c->setTitle("<b>" + title + "</b><br>t = " + QString::number(t) + 's');
    }
}

/*! Shows an animated plot of the evolution of x and y series over time.
    Each row of x and y represents the data series at an instant in time.
    The thread of execution is paused while the window is open. */
inline void playAnimation(MatrixXd x,
                          MatrixXd y,
                          double dt = 1,
                          QString title = "",
                          QString x_title = "",
                          QString y_title = "",
                          double t0 = 0){

    if(x.cols()!=y.cols() || x.rows()!=y.rows()){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Animation Error");
        msgBox.setText("x and y series must be the same size.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    if(x.cols()==0 || x.rows()==0){
        QMessageBox msgBox;
        msgBox.setWindowTitle("Animation Error");
        msgBox.setText("The animated series are empty.");
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
    chart->setTitle("<b>" + title + "</b><br>t = " + QString::number(t0) + 's');

    QtCharts::QLineSeries* p = new QtCharts::QLineSeries(chart);
    for (int i = 0; i < x.cols(); i++)
        p->append(QPointF(x(0,i), y(0,i)));
    chart->addSeries(p);

    chart->createDefaultAxes();
    chart->axisX()->setTitleText(x_title);
    chart->axisY()->setTitleText(y_title);
    chart->legend()->hide();

    double x_min = x.minCoeff();
    double x_max = x.maxCoeff();
    double y_min = y.minCoeff();
    double y_max = y.maxCoeff();
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    chart->axisX()->setRange(x_min - 0.05*x_range, x_max + 0.05*x_range);
    chart->axisY()->setRange(y_min - 0.05*y_range, y_max + 0.05*y_range);

    QtCharts::QChartView* chartView = new QtCharts::QChartView(chart);
    chartView->resize(600,600);
    chartView->show();

    std::thread t(workerAnimationFunc, x, y, dt, t0, title, p, chart);
    a.exec();
}

}

#endif // SIMPLEANIMATION_H
#endif // QT_CORE_LIB
