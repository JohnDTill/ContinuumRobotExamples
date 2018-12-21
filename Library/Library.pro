#TEMPLATE = aux
TEMPLATE = app
CONFIG += console c++11
QT += core charts widgets

HEADERS += \
    commonmath.h \
    convexoptimization.h \
    eigenmatrixtypes.h \
    numericaldifferentiation.h \
    numericalintegration.h \
    numericalpde.h \
    simpleplotting.h \
    timemanagement.h \
    simpleanimation.h
