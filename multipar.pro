TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$system(root-config --incdir)
LIBS += $$system(root-config --libs) -lMinuit -lSpectrum -lMathCore

SOURCES += \
        main.cpp
