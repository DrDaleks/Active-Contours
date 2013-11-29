package plugins.adufour.activecontours;

import icy.roi.ROI;
import icy.sequence.Sequence;

import java.awt.Color;

import javax.media.j3d.BoundingSphere;
import javax.vecmath.Point3d;

import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.nchenouard.spot.Detection;

/**
 * A generic active contour
 * 
 * @author Alexandre Dufour
 * 
 * @param <C>
 *            the type of contour (see {@link Polygon2D})
 * @param <R>
 *            the type of ROI defining the field bounds for this contour
 */
public abstract class ActiveContour extends Detection implements Iterable<Point3d>
{
    final ActiveContours           owner;
    
    SlidingWindow                  convergence;
    
    protected boolean              normalsNeedUpdate;
    
    protected final BoundingSphere boundingSphere = new BoundingSphere();
    
    protected boolean              boundingSphereNeedsUpdate;
    
    protected final EzVarDouble    contour_resolution;
    
    protected final EzVarInteger   contour_minArea;
    
    protected ActiveContour(ActiveContours owner, EzVarDouble contour_resolution, EzVarInteger contour_minArea, SlidingWindow convergenceWindow)
    {
        super(0, 0, 0, 0);
        
        this.owner = owner;
        this.contour_resolution = contour_resolution;
        this.contour_minArea = contour_minArea;
        this.convergence = convergenceWindow;
        
        setColor(Color.getHSBColor((float) Math.random(), 0.8f, 0.9f));
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    protected ActiveContour(ActiveContour contour)
    {
        this(contour.owner, contour.contour_resolution, contour.contour_minArea, new SlidingWindow(contour.convergence.size));
        // TODO later on, clone the resolution variables as well if needed
        
        setX(contour.x);
        setY(contour.y);
        setZ(contour.z);
        
        setColor(contour.getColor());
        
        initFrom(contour);
    }
    
    @Override
    public abstract ActiveContour clone();
    
    protected abstract void initFrom(ActiveContour contour);
    
    public void setConvergenceWindow(SlidingWindow window)
    {
        convergence = window;
    }
    
    protected abstract void addPoint(Point3d p);
    
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
    
    abstract void move(ROI field, double timeStep);
    
    protected abstract void updateNormalsIfNeeded();
    
    public BoundingSphere getBoundingSphere()
    {
        if (!boundingSphereNeedsUpdate) return boundingSphere;
        
        Point3d center = new Point3d();
        double radius = 0;
        
        // center calculation
        {
            double nbPts = 0;
            for (Point3d p : this)
            {
                center.add(p);
                nbPts++;
            }
            center.scale(1.0 / nbPts);
        }
        boundingSphere.setCenter(center);
        
        // radius calculation
        {
            for (Point3d p : this)
            {
                double d = p.distance(center);
                
                if (d > radius) radius = d;
            }
        }
        boundingSphere.setRadius(radius);
        boundingSphereNeedsUpdate = false;
        return boundingSphere;
    }
    
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
    
    /**
     * Computes the feedback forces yielded by the penetration of the current contour into the
     * target contour
     * 
     * @param target
     *            the contour that is being penetrated
     * @return the number of actual point-mesh intersection tests
     */
    abstract int computeFeedbackForces(ActiveContour target);
    
    public abstract double computeAverageIntensity(Sequence imageData, int channel, Sequence buffer);
    
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
    
    /**
     * 
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
    
    /**
     * 
     * @return a ROI representing the contour
     */
    public abstract ROI toROI();
}
