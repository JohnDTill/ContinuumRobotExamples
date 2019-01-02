QT       += core gui charts widgets

INCLUDEPATH += "../../../Library"

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = CSG_Teleoperation
TEMPLATE = app

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp

HEADERS += \
        mainwindow.h \
        InverseKinematicsCSG.h

FORMS += \
        mainwindow.ui
		
RC_ICONS += CSG.ico

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
