/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleIdTool.cc
 *
 *  @brief  Implementation of the beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

BeamParticleIdTool::BeamParticleIdTool() :
    m_selectAllBeamParticles(true),
    m_selectOnlyFirstSliceBeamParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::SelectOutputPfos(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (beamSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // First, simple approach
    if (m_selectAllBeamParticles || m_selectOnlyFirstSliceBeamParticles)
    {
        for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        {
            const PfoList &sliceOutput((m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex))) ?
                beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));

            selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
        }

        return;
    }

    // Now start to examine topology of beam slice hypotheses
    for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        bool usebeamHypothesis(false);
        const PfoList &pfoListBeam(beamSliceHypotheses.at(sliceIndex));

        try
        {
            PfoList allConnectedPfoList;
            LArPfoHelper::GetAllConnectedPfos(pfoListBeam, allConnectedPfoList);

            CaloHitList caloHitList3D;
            LArPfoHelper::GetCaloHits(allConnectedPfoList, TPC_3D, caloHitList3D);

            if (!caloHitList3D.empty())
            {
                CartesianVector centroid(0.f, 0.f, 0.f);
                LArPcaHelper::EigenVectors eigenVecs;
                LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
                LArPcaHelper::RunPca(caloHitList3D, centroid, eigenValues, eigenVecs);

                // Major axis of slice
                //const CartesianVector &fitDirection(eigenVecs.at(0));
            }            
        }
        catch (const StatusCodeException &) {}

        const PfoList &sliceOutput(usebeamHypothesis ? beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectAllBeamParticles", m_selectAllBeamParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectOnlyFirstSliceBeamParticles", m_selectOnlyFirstSliceBeamParticles));

    if (m_selectAllBeamParticles == m_selectOnlyFirstSliceBeamParticles)
    {
        std::cout << "BeamParticleIdTool::ReadSettings - exactly one of SelectAllBeamParticles and SelectOnlyFirstSliceBeamParticles must be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
