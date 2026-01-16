/**
 *  @file   larpandoradlcontent/LArHelpers/LArCanvasHelper.cc
 *
 *  @brief  Implementation of the lar deep learning helper helper class.
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArHelpers/LArCanvasHelper.h"

namespace lar_dl_content
{

void LArCanvasHelper::DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight)
{
    // Set the starting position for each circle bounding the ring
    int c1{inner}, r1{0}, c2{outer}, r2{0};
    int inner2{inner * inner}, outer2{outer * outer};
    while (c2 >= r2)
    {
        // Set the output pixel location
        int rp2{r2}, cp2{c2};
        // We're still within the octant for the inner ring, so use the inner pixel location (see Update comment below)
        // Note also that the inner row is always the same as the outer row, so no need to define rp1
        int cp1{c1};
        if (c1 <= r1)
        { // We've completed the arc of the inner ring already, so just move radially out from here (see Update comment below)
            cp1 = r2;
        }
        // Fill the pixels from inner to outer in the current row and their mirror pixels in the other octants
        for (int c = cp1; c <= cp2; ++c)
        {
            canvas[row + rp2][col + c] += weight;
            if (rp2 != c)
                canvas[row + c][col + rp2] += weight;
            if (rp2 != 0 && cp2 != 0)
            {
                canvas[row - rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col - rp2] += weight;
            }
            if (rp2 != 0)
            {
                canvas[row - rp2][col + c] += weight;
                if (rp2 != c)
                    canvas[row + c][col - rp2] += weight;
            }
            if (cp2 != 0)
            {
                canvas[row + rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col + rp2] += weight;
            }
        }
        // Only update the inner location while it remains in the octant (outer ring also remains in the octant of course, but the logic of
        // the update means that the inner ring can leave its octant before the outer ring is complete, so we need to stop that)
        if (c1 > r1)
            LArCanvasHelper::Update(inner2, c1, r1);
        // Update the outer location - increase the row position with every step, decrease the column position if conditions are met
        LArCanvasHelper::Update(outer2, c2, r2);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArCanvasHelper::Update(const int radius2, int &col, int &row)
{
    // Bresenham midpoint circle algorithm to determine if we should update the column position
    // This obscure looking block of code uses bit shifts and integer arithmetic to perform this check as efficiently as possible
    const int a{1 - (col << 2)};
    const int b{col * col + row * row - radius2 + (row << 2) + 1};
    const int c{(a << 2) * b + a * a};
    if (c < 0)
    {
        --col;
        ++row;
    }
    else
        ++row;
}

} // namespace lar_dl_content
