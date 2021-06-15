TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    globals.cpp \
    Main.cpp \
    efieldlagrangian.cpp \
    pzsolver.cpp \
    multisolver.cpp

HEADERS += \
    globals.h \
    efieldlagrangian.h \
    pzsolver.h \
    multisolver.h
#QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11

#QMAKE_CXXFLAGS_RELEASE += -O3   -msse -std=c++11

#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11
#QMAKE_LFLAGS += -O3  -msse -std=c++11
unix{
LIBS+=  -lGL -lGLU -lglut -lm
}
win32{
LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32
}
