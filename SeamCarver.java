/* 
 * The MIT License
 *
 * Copyright 2018 Manish Joshi.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

import edu.princeton.cs.algs4.Picture;
import java.util.Stack;

public class SeamCarver {

    private Picture picture;    // stores defensive copy of argument to the constructor
    private double[][] energy;  // stores energy of pixels
    private double[] distTo;    // vertex indexed array of shortest distances
    private int[] edgeTo;       // vertex indexed array of optimal incident edges
    private boolean transposed; // is Transposed ?
    private int height;
    private int width;

    /**
     * Creates a SeamCraver object based on the given Picture.
     *
     * @param picture
     */
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException("Null passed to constructor");
        }

        // make defensive copy of provided picture object
        this.picture = new Picture(picture);

        this.height = picture.height();
        this.width = picture.width();

        // setup energy array for first time, the function does its job
        energy = null;
        setEnergyMatrix(null);          // don't bother with null argument yet

        // not yet transposed dude, hang on
        this.transposed = false;
    }

    /**
     * Returns the current Picture.
     *
     * @return the current Picture
     */
    public Picture picture() {
        if (isTransposed()) {
            transpose();
        }
        return new Picture(picture);
    }

    /**
     * Returns width of the current Picture.
     *
     * @return width of the current Picture
     */
    public int width() {
        return width;
    }

    /**
     * Returns height of the current Picture.
     *
     * @return height of the current Picture
     */
    public int height() {
        return height;
    }

    /**
     * Returns the energy of the pixel at column x and row y.
     *
     * @param x the column index
     * @param y the row index
     * @return the energy of the pixel at column x and row y
     */
    public double energy(int x, int y) {
        if (isTransposed()) {
            transpose();
        }
        return energyAt(x, y);
    }

    // same as energy() but does not check for transposition
    private double energyAt(int x, int y) {
        validateCoordinates(x, y);
        if (isEdgePixel(x, y)) {
            return 1000;
        }
        // Calulate the dual gradient energy function
        return Math.sqrt(xGradientSqr(x, y) + yGradientSqr(x, y));
    }

    /**
     * Returns the sequence of indices for vertical seam.
     *
     * @return the sequence of indices for vertical seam
     */
    public int[] findVerticalSeam() {
        // we don't want transpose for finding vertical seam
        if (isTransposed()) {
            transpose();    // straighten if we have a transpose
        }
        return findSeam();      // find the seam and return
    }

    /**
     * Returns the sequence of indices for horizontal seam.
     *
     * @return the sequence of indices for horizontal seam
     */
    public int[] findHorizontalSeam() {
        // we run findSeam() for the transpose, which actually gives verticalSeam
        // but in this case, for the transpose. 
        // In Short, we want a transposed picture
        if (!isTransposed()) {
            transpose();        // so we get a transpose if needed
        }
        return findSeam();      // done! Get the fuck out of here.
    }

    /**
     * Removes the provided horizontal seam from the Picture.
     *
     * @param seam the seam to be removed
     */
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null || seam.length != width()) {
            throw new java.lang.IllegalArgumentException();
        }
        if (!isTransposed()) {
            transpose();
        }
        removeSeam(seam);
    }

    /**
     * Removes the provided vertical seam from the Picture.
     *
     * @param seam the seam to be removed
     */
    public void removeVerticalSeam(int[] seam) {
        if (seam == null || seam.length != height()) {
            throw new java.lang.IllegalArgumentException();
        }
        if (isTransposed()) {
            transpose();
        }
        removeSeam(seam);       // should be done !
    }

    /**
     * *************************************************************************
     * General Helper functions.
     * *************************************************************************
     */
    // is Picture transposed currently
    private boolean isTransposed() {
        return transposed;
    }

    // transpose the picture
    private void transpose() {
        width = picture.width();            // just make sure width and height
        height = picture.height();          // are proper becuase we messed with them right ?

        Picture tPicture = new Picture(height(), width());

        // remember Energy matrix is oriented differently
        double[][] tEnergy = new double[width()][height()];

        // iterate through all pixels in picture to make a transpose
        for (int col = 0; col < width(); ++col) {
            for (int row = 0; row < height(); ++row) {
                // copy the color values keeping transposition in mind
                tPicture.setRGB(row, col, picture.getRGB(col, row));

                // copy the energy matrix as transpose
                tEnergy[col][row] = energy[row][col];
            }
        }
        picture = tPicture;         // set Picture to Transpose of itself
        energy = tEnergy;           // set energy matrix to transpose
        width = picture.width();    // set new width and height
        height = picture.height();
        transposed = !transposed;
    }

    /* Fixes or Creates energy matrix, Don't call this method if your Energy matrix
    doesn't need to be fixed
     */
    private void setEnergyMatrix(int[] seam) {

        // set up energy array if it's the first time
        if (energy == null) {
            energy = new double[height()][width()];
            for (int col = 0; col < width(); ++col) {
                for (int row = 0; row < height(); ++row) {
                    energy[row][col] = energyAt(col, row);
                }
            }
        } // fix the removed pixel entries if it's not the first time
        else {

            double[][] newEnergy = new double[height()][width()];
            // Now here's the deal, we reuse the energy matrix 
            // using System.arrayCopy()
            for (int row = 0; row < height(); ++row) {
                double[] originalRow = energy[row];
                double[] newRow = new double[width()];

                System.arraycopy(originalRow, 0, newRow, 0, seam[row]);
                System.arraycopy(originalRow, seam[row] + 1, newRow, seam[row], width() - seam[row]);

                // recalculate energies of pixels along the seam both left and right
                if (seam[row] > 0) {
                    newRow[seam[row] - 1] = energyAt(seam[row] - 1, row);
                }
                if (seam[row] < width()) {
                    newRow[seam[row]] = energyAt(seam[row], row);
                }

                newEnergy[row] = newRow;
            }
            energy = newEnergy;         // Done!
        }
    }

    // initialise and setup the data structures for search
    private void setupSeamSearch() {

        // init distTo
        distTo = new double[height() * width()];

        // initialize edgeTo
        edgeTo = new int[width() * height()];

        // iterate the pixels, do some stuff
        for (int row = 0; row < height(); ++row) {
            for (int col = 0; col < width(); ++col) {

                int v = vertexAt(col, row);

                edgeTo[v] = -1;         // set dummy value

                // set distances to top row of pixels
                if (row == 0) {
                    distTo[v] = energy[row][col];
                } // set distances INFINITY for all other pixels
                else {
                    distTo[v] = Double.POSITIVE_INFINITY;
                }
                // works like charm ! Done dude!
            }
        }
    }

    // actually returns a vertical seam crafted using Graph Processing functions
    private int[] findSeam() {

        width = picture.width();            // just make sure width and height
        height = picture.height();          // are proper becuase we messed with them right ?

        // Setup the data structures required for the search
        setupSeamSearch();

        // Kickstart relaxation
        initRelaxation();

        int[] vertexShortestPath = getShortestPath();
        int[] seam = new int[vertexShortestPath.length];

        // pull off corresponding column indices from vertexShortestPath
        for (int i = 0; i < vertexShortestPath.length; ++i) {
            int col = vertexShortestPath[i] % width();
            seam[i] = col;
        }

        // check if findSeam() was called for a transposed picture. If So ? do some stuff
        if (isTransposed()) {
            // we are not transposing back yet,
            // Just taking a lazy approach and fix the width and height only
            // we will transpose thoroughly when required
            width = picture.height();
            height = picture.width();
        }
        return seam;
    }

    // actually deletes a vertical seam 
    private void removeSeam(int[] seam) {
        width = picture.width();            // just make sure width and height
        height = picture.height();          // are proper becuase we messed with them right ?

        // deal with the new picture, set it up discarding the seam pixels
        Picture newPicture = new Picture(width() - 1, height());
        for (int row = 0; row < height(); ++row) {
            for (int col = 0; col < seam[row]; ++col) {
                newPicture.setRGB(col, row, picture.getRGB(col, row));
            }
            for (int col = seam[row] + 1; col < width(); ++col) {
                newPicture.setRGB(col - 1, row, picture.getRGB(col, row));
            }
        }
        width = width() - 1;    // decrease width informally
        picture = newPicture;       // done, eh ?
        setEnergyMatrix(seam);  // remove holes from energy matrix

        // check if findSeam() was called for a transposed picture. If So ? do some stuff
        if (isTransposed()) {
            // we are not transposing back yet,
            // Just taking a lazy approach and fix the width and height only
            // we will transpose thoroughly when required
            width = picture.height();
            height = picture.width();
        }
    }

    // if the pixel is located on edge
    private boolean isEdgePixel(int col, int row) {
        return row == 0 || row == height() - 1 || col == width() - 1 || col == 0;
    }

    // returns the X-Gradient of a point at column x and row y
    private double xGradientSqr(int x, int y) {
        int rgb1 = picture.getRGB(x - 1, y);
        int rgb2 = picture.getRGB(x + 1, y);
        return rgbDiffSquare(rgb1, rgb2);
    }

    // returns the X-Gradient of a point at column x and row y
    private double yGradientSqr(int x, int y) {
        int rgb1 = picture.getRGB(x, y - 1);
        int rgb2 = picture.getRGB(x, y + 1);
        return rgbDiffSquare(rgb1, rgb2);
    }

    private double rgbDiffSquare(int rgb1, int rgb2) {
        int red1 = (rgb1 >> 16) & 0xFF;
        int green1 = (rgb1 >> 8) & 0xFF;
        int blue1 = rgb1 & 0xFF;

        int red2 = (rgb2 >> 16) & 0xFF;
        int green2 = (rgb2 >> 8) & 0xFF;
        int blue2 = rgb2 & 0xFF;

        int dR = red2 - red1;
        int dG = green2 - green1;
        int dB = blue2 - blue1;

        return (dR * dR) + (dG * dG) + (dB * dB);
    }

    /**
     * *************************************************************************
     * Graph helper functions.
     * *************************************************************************
     */
    // convert coordinates to Graph vertex
    private int vertexAt(int col, int row) {
        return row * width() + col;
    }

    // I know the name suffers from poor grammer
    private boolean isValidCoordinates(int x, int y) {
        return x >= 0 && x < width() && y >= 0 && y < height();
    }

    // validates coordinates
    private void validateCoordinates(int col, int row) {
        if (!isValidCoordinates(col, row)) {
            throw new IllegalArgumentException("col and row exceed the boundaries, col: " + col + ", row : " + row);
        }
    }

    // adjacency list of vertex at col, row
    private Iterable<Integer> adj(int col, int row) {
        Stack<Integer> adj = new Stack<>();

        // bottom
        int x = col;
        int y = row + 1;
        if (isValidCoordinates(x, y)) {
            adj.push(vertexAt(x, y));
        }

        // bottom left
        x = col - 1;
        y = row + 1;
        if (isValidCoordinates(x, y)) {
            adj.push(vertexAt(x, y));
        }

        // bottim right
        x = col + 1;
        y = row + 1;
        if (isValidCoordinates(x, y)) {
            adj.add(vertexAt(x, y));
        }

        return adj;
    }

    private void relax(int col, int row) {
        int v = vertexAt(col, row);
        for (int w : adj(col, row)) {
            int colW = w % width();
            int rowW = w / width();
            if (distTo[v] + energy[rowW][colW] < distTo[w]) {
                distTo[w] = distTo[v] + energy[rowW][colW];
                edgeTo[w] = v;
            }
        }
    }

    private void initRelaxation() {
        for (int row = 0; row < height(); ++row) {
            for (int col = 0; col < width(); ++col) {
                relax(col, row);
            }
        }
    }

    private int[] getShortestPath() {

        // calculate the end vertex of shortest path
        int nearestBottomVertex = vertexAt(width() - 1, height() - 1);
        for (int col = 0, row = height() - 1; col < width(); ++col) {
            if (distTo[nearestBottomVertex] > distTo[vertexAt(col, row)]) {
                nearestBottomVertex = vertexAt(col, row);
            }
        }

        int[] shortestPath = new int[height()];

        for (int x = nearestBottomVertex, row = height() - 1; row >= 0 && x != -1; x = edgeTo[x]) {
            shortestPath[row] = x;
            row--;
        }
        return shortestPath;
    }
}
