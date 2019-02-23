TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_LFLAGS += -pthread

SOURCES += \
        main.cpp

HEADERS += \
    bomce.h \
    bomce_utils.h \
    dna_alg.h \
    equal_cluster_metric.h \
    forward_list.h \
    graph.h \
    graphutility.h \
    list.h \
    scd.h \
    strchange.h \
    wg_test.h
