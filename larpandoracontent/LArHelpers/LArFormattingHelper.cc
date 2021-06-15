/**
 *  @file   larpandoracontent/LArHelpers/LArFormattingHelper.cc
 *
 *  @brief  Implementation of the lar formatting helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

#include <iomanip>

using namespace pandora;

namespace lar_content
{

void LArFormattingHelper::SetStyle(const LArFormattingHelper::Style style, std::ostream &stream)
{
    LArFormattingHelper::PrintFormatCharacter(static_cast<unsigned int>(style), stream);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::SetColor(const LArFormattingHelper::Color color, std::ostream &stream)
{
    LArFormattingHelper::PrintFormatCharacter(static_cast<unsigned int>(color), stream);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::ResetStyle(std::ostream &stream)
{
    LArFormattingHelper::SetStyle(REGULAR, stream);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::ResetColor(std::ostream &stream)
{
    LArFormattingHelper::SetColor(DEFAULT, stream);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Reset(std::ostream &stream)
{
    LArFormattingHelper::ResetColor(stream);
    LArFormattingHelper::ResetStyle(stream);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::PrintFormatCharacter(const unsigned int code, std::ostream &stream)
{
    if (!(stream << LArFormattingHelper::GetFormatCharacter(code)))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//--------------------------------------------------------------------------------------------------------------------------------------

std::string LArFormattingHelper::GetFormatCharacter(const unsigned int code)
{
    const std::string startFormattingCode("\x1B[");
    const std::string endFormattingCode("m");
    return (startFormattingCode + std::to_string(code) + endFormattingCode);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::PrintHeader(const std::string &title, const unsigned int width)
{
    LArFormattingHelper::SetStyle(BOLD);
    std::cout << std::endl << std::string(width, '-') << "\r-" << title << std::endl;
    LArFormattingHelper::Reset();
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::PrintRule(const unsigned int width)
{
    LArFormattingHelper::PrintHeader("", width);
}

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

LArFormattingHelper::Table::Table(const StringVector &columnTitles, const unsigned int precision) :
    m_columnTitles(columnTitles),
    m_precision(precision)
{
    m_stringstream.precision(m_precision);

    for (const std::string &title : m_columnTitles)
        m_widths.push_back(title.length());
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::CheckAndSetSeparatorColumn()
{
    if (m_columnTitles[m_elements.size() % m_columnTitles.size()] == "")
    {
        m_format.push_back("");
        m_elements.push_back("");
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

bool LArFormattingHelper::Table::IsSeparatorColumn(const unsigned int column) const
{
    if (column >= m_columnTitles.size())
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return (m_columnTitles[column] == "");
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::UpdateColumnWidth()
{
    if (m_widths.empty() || m_elements.empty())
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    const unsigned int currentElementIndex(m_elements.size() - 1);
    const unsigned int column(currentElementIndex % m_widths.size());
    m_widths[column] = std::max(static_cast<unsigned int>(m_elements[currentElementIndex].length()), m_widths[column]);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::Print() const
{
    if (m_columnTitles.empty())
        return;

    this->PrintColumnTitles();
    this->PrintHorizontalLine();
    this->PrintTableElements();
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::PrintColumnTitles() const
{
    for (unsigned int i = 0, nColumns = m_columnTitles.size(); i < nColumns; ++i)
        this->PrintTableCell(m_columnTitles[i], LArFormattingHelper::GetFormatCharacter(LArFormattingHelper::BOLD), i);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::PrintHorizontalLine() const
{
    for (unsigned int i = 0, nColumns = m_columnTitles.size(); i < nColumns; ++i)
        this->PrintTableCell(std::string(m_widths[i], '-'), LArFormattingHelper::GetFormatCharacter(LArFormattingHelper::REGULAR), i);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::PrintTableElements() const
{
    if (m_columnTitles.empty())
        return;

    const unsigned int nRows((m_elements.size() + m_columnTitles.size() - 1) / m_columnTitles.size());
    const unsigned int nColumns(m_columnTitles.size());
    if (nRows * nColumns != m_elements.size())
    {
        std::cout << "LArFormattingHelper::Table::PrintTableElements - Error: Number of table elements added doesn't fill a whole row" << std::endl;
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    for (unsigned int i = 0; i < nRows * nColumns; i++)
        this->PrintTableCell(m_elements[i], m_format[i], i);
}

//--------------------------------------------------------------------------------------------------------------------------------------

void LArFormattingHelper::Table::PrintTableCell(const std::string &value, const std::string &format, const unsigned int index) const
{
    if (m_columnTitles.empty())
        return;

    const unsigned int column(index % m_widths.size());

    const std::string separator(this->IsSeparatorColumn(column) ? "" : " ");
    std::cout << "|" << format << separator << std::setw(m_widths[column]) << std::internal << value << separator
              << std::resetiosflags(std::ios::adjustfield);
    LArFormattingHelper::ResetStyle();

    if (column == m_columnTitles.size() - 1)
        std::cout << "|" << std::endl;
}

} // namespace lar_content
