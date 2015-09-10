#-------------------------------------------------
#
# Project created by QtCreator 2015-09-09T20:48:54
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = blocknufft
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    blocknufft1d.cpp blocknufft2d.cpp blocknufft3d.cpp

QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

HEADERS += \
    blocknufft1d.h blocknufft2d.h blocknufft3d.h
