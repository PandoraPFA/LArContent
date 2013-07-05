ifndef PROJECT_DIR
    PROJECT_DIR = $(PANDORA_DIR)/LArContent
    PROJECT_LIBRARY_DIR = $(PANDORA_DIR)/lib
else
    PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib
endif

ifdef MONITORING
    DEFINES = -DMONITORING=1
endif

INCLUDES  = -I$(PROJECT_DIR)/include
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
ifdef MONITORING
    INCLUDES += -I$(PANDORA_DIR)/PandoraMonitoring/include
endif

CC = g++
CFLAGS = -c -Wall -g -w -fPIC -O2
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES  = $(wildcard $(PROJECT_DIR)/src/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/ClusterAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Clustering/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/ClusterSeedAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/ClusterSplitting/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Helpers/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Monitoring/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Objects/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Reclustering/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/ThreeDSeed/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/TwoDSeed/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Utility/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/Vertex/*.cc)

OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

LIBS = -L$(PANDORA_DIR)/lib -lPandoraSDK

ifdef MONITORING
    LIBS += -lPandoraMonitoring
endif

ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS = $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libLArContent.so

all: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

-include $(DEPENDS)

%.o:%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

install:
ifdef INCLUDE_TARGET
	rsync -r --exclude=.svn $(PROJECT_DIR)/include/ ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PROJECT_LIBRARY_DIR)/libLArContent.so ${LIB_TARGET}
endif

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(LIBRARY)
