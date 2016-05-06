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

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)
PROJECT_LIBRARY = $(PROJECT_LIBRARY_DIR)/libLArContent.so

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
ifdef MONITORING
    INCLUDES += -I$(PANDORA_DIR)/PandoraMonitoring/include
endif

ifdef MONITORING
    DEFINES = -DMONITORING=1
endif

SOURCES  = $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArCheating/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArCustomParticles/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArHelpers/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArMonitoring/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArObjects/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArPersistency/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArPlugins/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArStitching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArLongitudinalTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArPfoMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArShowerFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArTrackFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArThreeDReco/LArThreeDBase/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArClusterCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArClusterMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArTwoDReco/LArSeedFinding/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArUtility/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArContent/LArVertex/*.cc)
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
	rsync -r --include '*/' --include '*.h' --exclude '*' --prune-empty-dirs $(PROJECT_INCLUDE_DIR)/ ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PROJECT_LIBRARY) ${LIB_TARGET}
endif
