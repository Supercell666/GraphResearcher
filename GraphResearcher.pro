TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

<<<<<<< HEAD
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
=======
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
>>>>>>> a2c79c1ed8cb811c9c12aec25dbae1cb0a5d6b14
