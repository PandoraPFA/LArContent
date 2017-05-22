ifndef PROJECT_DIR
    PROJECT_DIR = $(PANDORA_DIR)/LArContent
    PROJECT_LIBRARY_DIR = $(PANDORA_DIR)/lib
else
    PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib
endif

CC = g++
CFLAGS = -c -g -fPIC -O2 -Wall -Wextra -pedantic -Werror -std=c++11
#-Wshadow removed for Eigen
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
INCLUDES += -I$(EIGEN_INC)/

ifdef MONITORING
    DEFINES = -DMONITORING=1
endif

SOURCES  = $(wildcard $(PROJECT_DIR)/larpandoracontent/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArCheating/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArCustomParticles/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArHelpers/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArMonitoring/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArObjects/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArPersistency/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArPlugins/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArStitching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArEventBuilding/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArHitCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArPfoMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArPfoRecovery/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArShowerFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArShowerMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArTrackFragments/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArThreeDReco/LArThreeDBase/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTrackShowerId/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/LArClusterAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/LArClusterCreation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/LArClusterMopUp/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/LArClusterSplitting/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArTwoDReco/LArCosmicRay/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArUtility/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/larpandoracontent/LArVertex/*.cc)
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
