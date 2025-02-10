/**
 *  @file   larpandoradlcontent/LArHelpers/LArCanvasHelper.h
 *
 *  @brief  Header file for the lar deep learning helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CANVAS_HELPER_H
#define LAR_CANVAS_HELPER_H 1

namespace lar_dl_content
{

/**
 *  @brief  LArCanvasHelper class
 */
class LArCanvasHelper
{
public:
    /**
     *  @brief  Add a filled ring to the specified canvas.
     *          The ring has an inner radius based on the minimum predicted distance to the vertex and an outer radius based on the maximum
     *          predicted distance to the vertex. The centre of the ring is the location of the hit used to predict the distance to the
     *          vertex. Each pixel to be filled is augmented by the specified weight. In this way, once all hits have been considered, a
     *          consensus view emerges of the likely vertex location based on the overlap of various rings centred at different locations.
     *
     *          The underlying implementation is a variant of the Bresenham midpoint circle algorithm and therefore only computes pixel
     *          coordinates for one octant of each circle (one of radius 'inner', one of radius 'outer') and interpolates the fill between
     *          points using integer arithmetic, guaranteeing each pixel of the ring is filled once and only once, and then mirrored to the
     *          remaining seven octants.
     *
     *  @param  networkOutput The TorchOutput object populated by the network inference step
     *  @param  pixelVector The vector of populated pixels
     *  @param  thresholds The fractional distance thresholds representing the classes predicted by the network
     *  @param  columnOffset The output column offset for the canvas
     *  @param  rowOffset The output row offset for the canvas
     *  @param  width The output width for the canvas
     *  @param  height The output height for the canvas
     */
    static void DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight);

    /**
     *  @brief  Update the coordinates along the loci of a circle.
     *          When drawing the ring we need an efficient means to determine the next pixel defining the inner and outer loci of the ring.
     *          This update function uses the Bresenham midpoint circle update function to determine this location. The row position is
     *          always incremented by 1 pixel, the column position is left unchanged, or decremented by 1 pixel to best follow the arc of
     *          the true underlying circle.
     *
     *  @param  radius2 The squared radius of the circle under consideration
     *  @param  col The input/output column position to (potentially) update
     *  @param  row The input/output row position to update
     */
    static void Update(const int radius, int &col, int &row);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_CANVAS_HELPER_H

