TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    globals.cpp \
    Main.cpp \
    phi_mult.cpp \
    sse_sum.cpp

HEADERS += \
    phi_mult.h \
    globals.h \
    sse_sum.h
#QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11

#LIBS += -lGL -lGLU -lglut -lpthread
LIBS += -L$$PWD/my_lib -lglut32

