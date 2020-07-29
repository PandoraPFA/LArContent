/**
 *  @file   larpandoracontent/LArDLContentRegistrations.h
 *
 *  @brief  Macro definitions for deep learning algorithm and tool registration
 *
 *  $Log: $
 */

#ifndef LAR_DL_CONTENT_REGISTRATIONS_H
#define LAR_DL_CONTENT_REGISTRATIONS_H 1

#include "larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.h"

#define LAR_DL_ALGORITHM_LIST(d)                                                                                                   \
    d("LArDeepLearningTrackShowerId",           DeepLearningTrackShowerIdAlgorithm)

#endif // #ifndef LAR_DL_CONTENT_REGISTRATIONS_H

