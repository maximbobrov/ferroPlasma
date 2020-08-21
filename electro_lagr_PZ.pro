TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    electronlagrangian.cpp \
    globals.cpp \
    Main.cpp \
    efieldlagrangian.cpp \
    lagrangiansolver.cpp \
    pzsolver.cpp \
    multisolver.cpp

HEADERS += \
    electronlagrangian.h \
    globals.h \
    efieldlagrangian.h \
    lagrangiansolver.h \
    pzsolver.h \
    multisolver.h
#QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11

#LIBS+=  -lGL -lGLU -lglut -lm
LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32

