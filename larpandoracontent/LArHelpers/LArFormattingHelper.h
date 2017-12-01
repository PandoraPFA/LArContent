/**
 *  @file   larpandoracontent/LArHelpers/LArFormatting.h
 *
 *  @brief  Header file for the formatting helper class
 *
 *  $Log: $
 */
#ifndef LAR_FORMATTING_HELPER_H
#define LAR_FORMATTING_HELPER_H 1

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace lar_content
{

/**
 *  @brief  LArFormattingHelper class
 */
class LArFormattingHelper
{
public:

    enum Colour : int;
    enum Style : int;

    /**
     *  @brief  Set the format style (to standard output stream)
     *
     *  @param  style the style of choice
     */
    static void SetStyle(const Style &style, std::ostream &stream=std::cout);
    
    /**
     *  @brief  Set the text colour (of standard output stream)
     *
     *  @param  colour the colour of choice
     */
    static void SetColour(const Colour &colour, std::ostream &stream=std::cout);
    
    /**
     *  @brief  Reset the style of the standard output stream
     */
    static void ResetStyle(std::ostream &stream=std::cout);

    /**
     *  @brief  Reset the text colour of the standard output stream
     */
    static void ResetColour(std::ostream &stream=std::cout);

    /**
     *  @brief  Reset the formatting and text colour of the standard output stream
     */
    static void Reset(std::ostream &stream=std::cout);

    /**
     *  @brief  Print a formatting character to the standard output stream
     *
     *  @param  code the formatting code to output
     */
    static void PrintFormatCharacter(const unsigned int &code, std::ostream &stream=std::cout);

    /**
     *  @brief  Get a formatting character
     *
     *  @param  code the formatting code to output
     */
    static std::string GetFormatCharacter(const unsigned int &code);

    /**
     *  @brief  Print a header line of a given width
     *
     *  @param  title the title of the header
     *  @param  width the width of the header line
     */
    static void PrintHeader(const std::string &title = "", const unsigned int &width = 140);

    /**
     *  @brief  Print a horizontal rule
     *
     *  @param  width the width of the rule line
     */
    static void PrintRule(const unsigned int &width = 140);

    enum Style : int
    {
        REGULAR = 0,
        BOLD = 1,
        DIM = 2,
        UNDERLINED = 4,
        INVERTED = 7
    };
    
    enum Colour : int
    {
        DEFAULT = 39,
        BLACK = 30,
        RED = 31,
        GREEN = 32,
        YELLOW = 33,
        BLUE = 34,
        MAGENTA = 35,
        CYAN = 36,
        LIGHT_GRAY = 37,
        DARK_GRAY = 90,
        LIGHT_RED = 91,
        LIGHT_GREEN = 92,
        LIGHT_YELLOW = 93,
        LIGHT_BLUE = 94,
        LIGHT_MAGENTA = 95,
        LIGHT_CYAN = 96,
        WHITE = 97
    };

    //--------------------------------------------------------------------------------------------------------------------------------------

    class Table 
    {
    public:

        /**
         * @brief  Table constructor
         *
         * @param  columnTitles the string columns titles use empty string for a separator column
         * @param  precision the number of significant figures to display for number type table elements
         */
        Table(const std::vector<std::string> &columnTitles, const unsigned int &precision=3);

        /**
         * @brief  Add an element to the table into the next (non-separator) column
         *
         * @param  value the element to add to the table
         * @param  style the style of the element
         * @param  colour the colour of the element
         */
        template<class T>
        void AddElement(const T &value, const Style &style=REGULAR, const Colour &colour=DEFAULT);

        /**
         * @brief  Print the table
         */
        void Print();
    private:
        const std::vector<std::string>  m_columnTitles;     ///< The vector of columns titles in the table
        const unsigned int              m_precision;        ///< The number of significant figures to use when displaying number types
        std::vector<std::string>        m_elements;         ///< The vector of flattened table elements
        std::vector<std::string>        m_format;           ///< The formatting of each table element
        std::vector<unsigned int>       m_widths;           ///< The widths of each column (in units of number of characters)
        std::stringstream               m_stringstream;     ///< The stringstream to print objects to

        /**
         *  @brief  Check if the next table cell is in a separator column, if so add a blank element
         */
        void CheckForSeparatorColumn();

        /**
         *  @brief  If the supplied column is a separator (vertical rule)
         *
         *  @param  column the column index to check
         *
         *  @return true/false
         */
        bool IsSeparatorColumn(const unsigned int &column);
    
        /**
         *  @brief  Update the width of the last column in which an element was added
         */
        void UpdateColumnWidth();

        /**
         *  @brief  Print the column titles
         *
         *  @param  widths a vector of the number of characters to include in each column 
         */
        void PrintColumnTitles();

        /**
         *  @brief  Print a horizontal line
         *
         *  @param  widths a vector of the number of characters to include in each column 
         */
        void PrintHorizontalLine();

        /**
         *  @brief  Print the table elements
         *
         *  @param  widths a vector of the number of characters to include in each column 
         */
        void PrintTableElements();

        /**
         *  @brief  Print a table cell
         *
         *  @param  value the string to print in the cell
         *  @param  format a formatting string (see GetFormatCharacter(...))
         *  @param  index the index of the table cell (0 --> number of elements) used to find the column
         *  @param  widths a vector of the number of characters to include in each column 
         */
        void PrintTableCell(const std::string &value, const std::string &format, const unsigned int &index);
    };
    
};

template<class T>
inline void LArFormattingHelper::Table::AddElement(const T &value, const Style &style, const Colour &colour)
{
    this->CheckForSeparatorColumn();

    m_stringstream << value;
    m_elements.push_back(static_cast<std::string>(m_stringstream.str()));
    m_format.push_back(LArFormattingHelper::GetFormatCharacter(style) + LArFormattingHelper::GetFormatCharacter(colour));
    m_stringstream.str("");

    this->UpdateColumnWidth();
    this->CheckForSeparatorColumn();
}

} // namespace lar_content

#endif // #ifndef LAR_FORMATTING_HELPER_H
