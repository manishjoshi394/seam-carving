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

/**
 * The class provides an API for Seam Carving.
 * <p>
 * It has methods to find optimal seams and remove horizontal and vertical seams
 * from the image.
 *
 * @author Manish Joshi
 */
public class SeamCarver {

    private Picture picture;    // stores defensive copy of argument to the constructor
    private boolean transposed; // is Transposed ?
    private int height;         // stores height 
    private int width;          // stores width of picture

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
        validateSeam(seam);
        if (seam.length != width() || height() <= 1) {
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
        validateSeam(seam);
        if (seam.length != height() || width() <= 1) {
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

        // iterate through all pixels in picture to make a transpose
        for (int col = 0; col < width(); ++col) {
            for (int row = 0; row < height(); ++row) {
                // copy the color values keeping transposition in mind
                tPicture.setRGB(row, col, picture.getRGB(col, row));
            }
        }
        picture = tPicture;         // set Picture to Transpose of itself
        width = picture.width();    // set new width and height
        height = picture.height();
        transposed = !transposed;
    }

    // initialise and setup the data structures for search
    private void setupSeamSearch(double[] distTo, int[] edgeTo, double[][] energy) {
        // iterate the pixels, do some stuff
        for (int row = 0; row < height(); ++row) {
            for (int col = 0; col < width(); ++col) {

                int v = vertexAt(col, row);

                edgeTo[v] = -1;         // set dummy value

                // set distances to top row of pixels
                if (row == 0) {
                    distTo[v] = energy[col][row];
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

        double[] distTo;    // vertex indexed array of shortest distances
        int[] edgeTo;       // vertex indexed array of optimal incident edges
        double[][] energy;

        // init distTo
        distTo = new double[height() * width()];

        // initialize edgeTo
        edgeTo = new int[width() * height()];

        // setup local energy matrix
        energy = new double[width()][height()];
        for (int col = 0; col < width(); ++col) {
            for (int row = 0; row < height(); ++row) {
                energy[col][row] = energyAt(col, row);
            }
        }

        // Setup the data structures required for the search
        setupSeamSearch(distTo, edgeTo, energy);

        // Kickstart relaxation
        initRelaxation(distTo, edgeTo, energy);

        int[] vertexShortestPath = getShortestPath(distTo, edgeTo);
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

        // check if findSeam() was called for a transposed picture. If So ? do some stuff
        if (isTransposed()) {
            // we are not transposing back yet,
            // Just taking a lazy approach and fix the width and height only
            // we will transpose thoroughly when required
            width = picture.height();
            height = picture.width();
        }
    }

    // validate seam
    private void validateSeam(int[] seam) {
        if (seam == null) {
            throw new java.lang.IllegalArgumentException();
        }
        for (int i = 0; i < seam.length - 1; ++i) {
            if (Math.abs(seam[i] - seam[i + 1]) > 1) {
                throw new java.lang.IllegalArgumentException("Adjacent pixels in seam have distances more than 1");
            }
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

    // relax all edges connected to the vertex at col, row
    private void relax(int col, int row, double[] distTo, int[] edgeTo, double[][] energy) {
        int v = vertexAt(col, row);
        for (int w : adj(col, row)) {
            int colW = w % width();
            int rowW = w / width();
            if (distTo[v] + energy[colW][rowW] < distTo[w]) {
                distTo[w] = distTo[v] + energy[colW][rowW];
                edgeTo[w] = v;
            }
        }
    }

    // kickstart relaxation
    private void initRelaxation(double[] distTo, int[] edgeTo, double[][] energy) {
        for (int row = 0; row < height(); ++row) {
            for (int col = 0; col < width(); ++col) {
                relax(col, row, distTo, edgeTo, energy);
            }
        }
    }

    // returns an array of vertices in the shortest path
    private int[] getShortestPath(double[] distTo, int[] edgeTo) {

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
