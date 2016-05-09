ifndef PROJECT_DIR
    PROJECT_DIR = $(PANDORA_DIR)/LArContent
    PROJECT_LIBRARY_DIR = $(PANDORA_DIR)/lib
else
    PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib
endif

CC = g++
CFLAGS = -c -g -fPIC -O2 -Wall -Wextra -pedantic -Wshadow -Werror -std=c++11
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

LIBS = -L$(PANDORA_DIR)/lib -lPandoraSDK
ifdef MONITORING
    LIBS += -lPandoraMonitoring
endif
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)/include/
PROJECT_LIBRARY = $(PROJECT_LIBRARY_DIR)/libLArContent.so

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
ifdef MONITORING
    INCLUDES += -I$(PANDORA_DIR)/PandoraMonitoring/include
endif

ifdef MONITORING
    DEFINES = -DMONITORING=1
endif

SOURCES  = $(wildcard $(PROJECT_DIR)/src/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArCheating/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArCustomParticles/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArHelpers/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArMonitoring/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArObjects/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArPersistency/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArPlugins/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArStitching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArEventBuilding/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArHitCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArLongitudinalTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArPfoMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArShowerFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArShowerMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArTrackFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArTransverseTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDReco/LArThreeDBase/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArClusterAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArClusterCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArClusterMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArClusterSplitting/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDReco/LArSeedFinding/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArUtility/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArVertex/*.cc)
OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

all: library

library: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(PROJECT_LIBRARY)

-include $(DEPENDS)

%.o:%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(PROJECT_LIBRARY)

install:
ifdef INCLUDE_TARGET
	rsync -r --exclude=.svn $(PROJECT_INCLUDE_DIR) ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PROJECT_LIBRARY) ${LIB_TARGET}
endif
