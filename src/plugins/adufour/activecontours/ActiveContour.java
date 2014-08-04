package plugins.adufour.activecontours;

import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.system.SystemUtil;
import icy.system.thread.Processor;

import java.awt.Color;

import javax.media.j3d.BoundingBox;
import javax.media.j3d.BoundingSphere;
import javax.vecmath.Point3d;

import plugins.adufour.activecontours.ActiveContours.ROIType;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.nchenouard.spot.Detection;

/**
 * A generic active contour
 * 
 * @author Alexandre Dufour
 */
public abstract class ActiveContour extends Detection implements Iterable<Point3d>
{
    protected final Processor      processor      = new Processor(SystemUtil.getAvailableProcessors() * 2);
    
    protected final SlidingWindow  convergence;
    
    protected final EzVarDouble    resolution;
    
    protected final BoundingSphere boundingSphere = new BoundingSphere();
    
    protected final BoundingBox    boundingBox    = new BoundingBox();
    
    protected ActiveContour(EzVarDouble contour_resolution, SlidingWindow convergenceWindow)
    {
        super(0, 0, 0, 0);
        
        this.resolution = contour_resolution;
        this.convergence = convergenceWindow;
        
        setColor(Color.getHSBColor((float) Math.random(), 0.8f, 0.9f));
        
        processor.setDefaultThreadName(getClass().getSimpleName());
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    protected ActiveContour(ActiveContour contour)
    {
        this(contour.resolution, new SlidingWindow(contour.convergence.getSize()));
        
        setX(contour.x);
        setY(contour.y);
        setZ(contour.z);
        
        setColor(contour.getColor());
        
        initFrom(contour);
    }
    
    @Override
    public abstract ActiveContour clone();
    
    /**
     * Adds the specified point to the contour. It is up to the implementing classes to determine
     * whether the points should be added in a specific order or based on gemoetrical rules
     * 
     * @param p
     *            the point to add
     * @return the index indicating where the point has been inserted in the internal list
     */
    protected abstract int addPoint(Point3d p);
    
    /**
     * Checks whether the contour is self-intersecting. Depending on the given parameters, a
     * self-intersection can be considered as a loop or as a contour division.
     * 
     * @param minDistance
     *            the distance threshold between non-neighboring points to detect self-intersection
     * @return null if either no self-intersection is detected or if one of the new contours is too
     *         small, or an array of contours with 0 elements if both contours are too small, and 2
     *         elements if both contours are viable
     */
    protected abstract ActiveContour[] checkSelfIntersection(double minDistance);
    
    protected abstract void initFrom(ActiveContour contour);
    
    /**
     * Update the axis constraint force, which adjusts the takes the final forces and normalize them
     * to keep the contour shape along its principal axis <br>
     * WARNING: this method directly update the array of final forces used to displace the contour
     * points. It should be used among the last to keep it most effective
     * 
     * @param weight
     */
    abstract void computeAxisForces(double weight);
    
    abstract void computeBalloonForces(double weight);
    
    /**
     * Update edge term of the contour evolution according to the image gradient
     * 
     * @param weight
     * @param edgeData
     *            a sequence containing the edge information (one channel per edge direction)
     */
    abstract void computeEdgeForces(Sequence edgeData, int channel, double weight);
    
    /**
     * Update region term of the contour evolution according to the Chan-Vese-Mumford-Shah
     * functional
     * 
     * @param imageData
     *            the image data (must a double-type image of range [0-1])
     * @param weight
     *            the weight of the data attachment term
     * @param cin
     *            the intensity mean inside the contour
     * @param cout
     *            the intensity mean outside the contour
     * @param sensitivity
     *            set 1 for default, greater than 1 for high SNRs and vice-versa
     */
    abstract void computeRegionForces(Sequence imageData, int channel, double weight, double sensitivity, double cin, double cout);
    
    abstract void computeInternalForces(double weight);
    
    abstract void computeVolumeConstraint(double targetVolume);
    
    /**
     * Computes the feedback forces yielded by the penetration of the current contour into the
     * target contour
     * 
     * @param target
     *            the contour that is being penetrated
     * @return the number of actual point-mesh intersection tests
     */
    abstract int computeFeedbackForces(ActiveContour target);
    
    /**
     * Compute the average image intensity inside the contour on the specified image data, and fill
     * out the mask buffer to allow the global exterior mean to be computed
     * 
     * @param imageData
     *            the data on which the average intensity should be computed
     * @param channel
     *            the channel on which the average intensity should be computed
     * @param buffer
     *            the binary mask buffer where this contour should be rasterised
     * @return the average intensity inside the contour
     * @throws TopologyException
     *             if the contour becomes extremely thin to the point where it contains no pixel to
     *             measure intensity
     */
    public abstract double computeAverageIntensity(Sequence imageData, int channel, Sequence buffer) throws TopologyException;
    
    /**
     * Tests whether the given point is inside the contour, and if so returns the penetration depth
     * of this point. <br>
     * 
     * @param p
     *            a point to test
     * @return <ul>
     *         <li/>if <code>p</code> is outside: <code>0</code>
     *         <li/>if <code>p</code> is inside: the distance from <code>p</code> to the contour
     *         edge
     *         </ul>
     */
    public abstract double contains(Point3d p);
    
    public BoundingBox getBoundingBox()
    {
        BoundingBox bbox = new BoundingBox();
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double minZ = Double.MAX_VALUE;
        double maxX = 0.0;
        double maxY = 0.0;
        double maxZ = 0.0;
        
        for (Point3d p : this)
        {
            if (p.x < minX) minX = p.x;
            if (p.y < minY) minY = p.y;
            if (p.z < minZ) minZ = p.z;
            if (p.x > maxX) maxX = p.x;
            if (p.y > maxY) maxY = p.y;
            if (p.z > maxZ) maxZ = p.z;
        }
        
        bbox.setLower(minX, minY, minZ);
        bbox.setUpper(maxX, maxY, maxZ);
        
        return bbox;
    }
    
    /**
     * @param order
     *            the dimension (a.k.a. norm) to compute:<br/>
     *            <ul>
     *            <li>0: number of points,</li>
     *            <li>1: perimeter (2D) or surface area (3D),</li>
     *            <li>2: surface (2D) or volume (3D)</li>
     *            </ul>
     * @return the dimension for the specified order
     */
    public abstract double getDimension(int order);
    
    abstract void move(ROI field, double timeStep);
    
    /**
     * Re-samples the Contour according to an 'average distance between points' criterion. This
     * method ensures that the distance between two consecutive points is strictly comprised between
     * a minimum value and a maximum value. In order to avoid oscillatory behavior, 'max' and 'min'
     * should verify the following relations: min < 1, max > 1, 2*min <= max.
     * 
     * @param minFactor
     *            the minimum distance between two points.
     * @param maxFactor
     *            the maximum distance between two points.
     */
    public abstract void reSample(double minFactor, double maxFactor) throws TopologyException;
    
    /**
     * @return a ROI representing the contour
     * @deprecated use {@link #toROI(ROIType)} instead
     */
    public abstract ROI toROI();
    
    /**
     * @param type
     *            the type of ROI to export
     * @param sequence
     *            a sequence which can be used to retrieve the pixel size (may be needed if the
     *            contour is defined in real units rather than pixel units)
     * @return a ROI representing the contour
     * @throws UnsupportedOperationException
     *             if the contour cannot be exported in the requested type
     */
    public abstract ROI toROI(ROIType type, Sequence sequence) throws UnsupportedOperationException;
    
    /**
     * Paints the contour onto the specified sequence with the specified value
     * 
     * @param output
     * @param value
     */
    public abstract void toSequence(Sequence output, double value);
    
    /**
     * Updates the contour's meta-data (i.e. data that can be computed directly from the actual
     * contour data, but stored locally for fast repetitive access). This includes:<br/>
     * <ul>
     * <li>bounding box</li>
     * <li>bounding sphere</li>
     * <li>contour normals</li>
     * <li>center of mass</li>
     * </ul>
     */
    protected void updateMetaData()
    {
        // center of mass
        Point3d center = new Point3d();
        
        // bounding sphere
        double radius = 0;
        
        // bounding box
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double minZ = Double.MAX_VALUE;
        double maxX = 0.0;
        double maxY = 0.0;
        double maxZ = 0.0;
        
        // center calculation
        double nbPts = 0;
        for (Point3d p : this)
        {
            nbPts++;
            
            center.add(p);
            if (p.x < minX) minX = p.x;
            if (p.y < minY) minY = p.y;
            if (p.z < minZ) minZ = p.z;
            if (p.x > maxX) maxX = p.x;
            if (p.y > maxY) maxY = p.y;
            if (p.z > maxZ) maxZ = p.z;
        }
        
        center.scale(1.0 / nbPts);
        setX(center.x);
        setY(center.y);
        setZ(center.z);
        
        boundingSphere.setCenter(center);
        
        // radius calculation
        for (Point3d p : this)
        {
            double d = p.distance(center);
            
            if (d > radius) radius = d;
        }
        
        boundingSphere.setRadius(radius);
        
        boundingBox.setLower(minX, minY, minZ);
        boundingBox.setUpper(maxX, maxY, maxZ);
        
        updateNormals();
    }
    
    protected abstract void updateNormals();
    
    /**
     * Optional method to clean intermediate resources before destroying the contours (e.g. VTK
     * overlays)
     */
    protected abstract void clean();
}
