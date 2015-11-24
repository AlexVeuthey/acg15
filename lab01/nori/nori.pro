TARGET = nori

SOURCES += src/common.cpp \
           src/ao.cpp \
           src/bitmap.cpp \
           src/block.cpp \
           src/gui.cpp \
           src/independent.cpp \
           src/kdtree.cpp \
           src/main.cpp \
           src/medium.cpp \
           src/mesh.cpp \
           src/obj.cpp \
           src/object.cpp \
           src/parser.cpp \
           src/perspective.cpp \
           src/proplist.cpp \
           src/random.cpp \
           src/rfilter.cpp \
           src/scene.cpp

HEADERS += $$PWD/include/nori/*.h \
           include/nori/gui.h

# Exercise 1
SOURCES += $$PWD/src/ex1/depth*.cpp
SOURCES += src/diffuse.cpp  # Only in Ex1.


RESOURCES += data/resources.qrc

INCLUDEPATH += include

DEPENDPATH += include data
OBJECTS_DIR = build
RCC_DIR = build
MOC_DIR = build
UI_DIR = build
DESTDIR = bin

QT += xml xmlpatterns opengl

macx {
        OBJECTIVE_SOURCES += src/support_osx.m
        QMAKE_LFLAGS += -framework Cocoa -lobjc -mmacosx-version-min=10.7 -stdlib=libc++
        QMAKE_CXXFLAGS_X86_64 += -mmacosx-version-min=10.7
        QMAKE_CXXFLAGS += -stdlib=libc++
        # Install OpenEXR with homebrew
        QMAKE_LIBDIR += /usr/local/lib
        INCLUDEPATH += /usr/local/include/OpenEXR
        # Homebrew version of boost somehow doesn't work.
        INCLUDEPATH += $$PWD/include/boost1.49_min
}

unix {
        QMAKE_CXXFLAGS_RELEASE -= -O3
        QMAKE_CXXFLAGS += -O2 -march=nocona -msse2 -mfpmath=sse -fstrict-aliasing
        LIBS += -lIlmImf -lIex
        INCLUDEPATH += $$PWD/include/OpenEXR
        QMAKE_RPATHDIR = += $$PWD/lib
        LIBS += -L$$PWD/lib
}

win32 {
        # You will have to update the following three lines based on where you have installed OpenEXR
        QMAKE_LIBDIR += ./openexr/lib/x64
        INCLUDEPATH += ./openexr/include
        LIBPATH += C:\your\path\to\nori\openexr\lib\i386

        INCLUDEPATH += ./include/boost1.49_min

        QMAKE_CXXFLAGS += /O2 /fp:fast /GS- /GL /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS
        QMAKE_LDFLAGS += /LTCG
        SOURCES += src/support_win32.cpp

        LIBS += IlmImf.lib Iex.lib IlmThread.lib Imath.lib Half.lib
}
