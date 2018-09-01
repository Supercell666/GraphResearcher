TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    bomce.h \
    clusterisator.h \
    dna_alg.h \
    dsu.h \
    equal_cluster_metric.h \
    graph.h \
    graphutility.h \
    scd.h \
    strchange.h \
    wg_test.h \
    list.h

QMAKE_LFLAGS += -pthread
