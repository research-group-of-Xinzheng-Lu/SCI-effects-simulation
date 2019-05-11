TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES +=  createSAM.cpp \
	Building.cpp \
	FloorParam.cpp \
	HazusSAM_Generator.cpp \
	InterstoryParam.cpp
    csvparser.c

HEADERS += \
    Building.h \
	FloorParam.h \
	HazusSAM_Generator.h \
	InterstoryParam.h

win32 {
INCLUDEPATH += $$PWD\..\..\include
LIBS += -L$$PWD\..\..\lib -llibjansson-4
}
