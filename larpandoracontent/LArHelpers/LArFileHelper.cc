/**
 *  @file   larpandoracontent/LArHelpers/LArFileHelper.cc
 *
 *  @brief  Implementation of the file helper class.
 *
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"

#include <cstdlib>
#include <sys/stat.h>

using namespace pandora;

namespace lar_content
{

std::string LArFileHelper::FindFileInPath(const std::string &unqualifiedFileName, const std::string &environmentVariable, const std::string &delimiter)
{
    StringVector filePaths;
    const char *const pFilePathList(std::getenv(environmentVariable.c_str()));

    if (pFilePathList)
        XmlHelper::TokenizeString(pFilePathList, filePaths, delimiter);

    // Always test unqualified file name too
    filePaths.push_back("");

    for (const std::string &filePath : filePaths)
    {
        const std::string qualifiedFileNameAttempt(filePath + "/" + unqualifiedFileName);
        struct stat fileInfo;

        if (0 == stat(qualifiedFileNameAttempt.c_str(), &fileInfo))
            return qualifiedFileNameAttempt;
    }

    std::cout << "Unable to find file  " << unqualifiedFileName << " in any path specified by environment variable " << environmentVariable
              << ", delimiter " << delimiter << std::endl;
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

} // namespace lar_content
