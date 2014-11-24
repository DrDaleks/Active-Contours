package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.canvas.IcyCanvas2D;
import icy.gui.frame.progress.AnnounceFrame;
import icy.main.Icy;
import icy.roi.BooleanMask2D;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.system.IcyHandledException;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import plugins.adufour.activecontours.ActiveContours.ROIType;
import plugins.adufour.activecontours.SlidingWindow.Operation;
import plugins.adufour.morphology.FillHolesInROI;
import plugins.adufour.vars.lang.Var;
import plugins.adufour.vars.lang.VarDouble;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DEllipse;
import plugins.kernel.roi.roi2d.ROI2DPolygon;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.kernel.roi.roi2d.ROI2DShape;

public class Polygon2D extends ActiveContour
{
    private static final class Segment implements Iterable<Point3d>
    {
        final ArrayList<Point3d> points;
        
        Segment(Point3d head, Point3d tail)
        {
            points = new ArrayList<Point3d>(2);
            points.add(head);
            points.add(tail);
        }
        
        final Point3d getHead()
        {
            return points.get(0);
        }
        
        final Point3d getTail()
        {
            return points.get(points.size() - 1);
        }
        
        final void addHead(Point3d p)
        {
            points.add(0, p);
        }
        
        final void addHead(Segment s)
        {
            for (int i = 0; i < s.points.size(); i++)
                points.add(i, s.points.get(i));
        }
        
        final void addTail(Point3d p)
        {
            points.add(p);
        }
        
        public Iterator<Point3d> iterator()
        {
            return points.iterator();
        }
    }
    
    private final FillHolesInROI holeFiller = new FillHolesInROI();
    
    final ArrayList<Point3d>     points     = new ArrayList<Point3d>();
    
    Path2D.Double                path       = new Path2D.Double();
    
    double                       cout       = 0.0;
    
    private Vector3d[]           modelForces;
    
    private Vector3d[]           contourNormals;
    
    private Vector3d[]           feedbackForces;
    
    private boolean              counterClockWise;
    
    protected Polygon2D(Var<Double> sampling, SlidingWindow convergenceWindow)
    {
        super(sampling, convergenceWindow);
        
        setColor(Color.getHSBColor((float) Math.random(), 0.8f, 0.9f));
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    public Polygon2D(Polygon2D contour)
    {
        this(contour.sampling, new SlidingWindow(contour.convergence.getSize()));
        
        setColor(contour.getColor());
        int n = contour.points.size();
        
        points.ensureCapacity(n);
        contourNormals = new Vector3d[n];
        
        for (int i = 0; i < n; i++)
        {
            contourNormals[i] = new Vector3d();
            addPoint(new Point3d(contour.points.get(i)));
        }
        
        updateNormals();
        
        // // artifically grow the contour by a tiny bit
        // for (int i = 0; i < n; i++)
        // {
        // Point3d p = points.get(i);
        // p.scaleAdd(2.0, contourNormals[i], p);
        // }
        
        updatePath();
        counterClockWise = contour.counterClockWise;
    }
    
    public Polygon2D(Var<Double> sampling, SlidingWindow convergenceWindow, ROI2D roi)
    {
        this(sampling, convergenceWindow);
        
        if (!(roi instanceof ROI2DEllipse) && !(roi instanceof ROI2DRectangle) && !(roi instanceof ROI2DPolygon) && !(roi instanceof ROI2DArea))
        {
            throw new IcyHandledException("Active contours: invalid ROI. Only Rectangle, Ellipse, Polygon and Area are supported");
        }
        
        if (roi instanceof ROI2DArea)
        {
            // dilate ROI
            roi = dilateROI(roi, 1, 1);
            
            // fill holes first
            BooleanMask2D mask = roi.getBooleanMask(true);
            if (holeFiller.fillHoles(mask)) ((ROI2DArea) roi).setAsBooleanMask(mask);
            
            try
            {
                triangulate((ROI2DArea) roi, sampling.getValue());
            }
            catch (TopologyException e)
            {
                int z = roi.getZ();
                int t = roi.getT();
                roi = new ROI2DEllipse(roi.getBounds2D());
                roi.setZ(z);
                roi.setT(t);
            }
            
        }
        
        if (points.size() <= 4)
        {
            // replace by ellipse
            int z = roi.getZ();
            int t = roi.getT();
            roi = new ROI2DEllipse(roi.getBounds2D());
            roi.setZ(z);
            roi.setT(t);
            
            points.clear();
            
            // convert the ROI into a linked list of points
            
            double[] segment = new double[6];
            
            PathIterator pathIterator = ((ROI2DShape) roi).getPathIterator(null, 0.1);
            
            // first segment is necessarily a "move to" operation
            
            pathIterator.currentSegment(segment);
            addPoint(new Point3d(segment[0], segment[1], z));
            
            while (!pathIterator.isDone())
            {
                if (pathIterator.currentSegment(segment) == PathIterator.SEG_LINETO)
                {
                    addPoint(new Point3d(segment[0], segment[1], z));
                }
                pathIterator.next();
                
                // the final one should be a "close" operation, do nothing
            }
            
            final int nbPoints = points.size();
            if (modelForces == null || modelForces.length != nbPoints)
            {
                modelForces = new Vector3d[nbPoints];
                contourNormals = new Vector3d[nbPoints];
                feedbackForces = new Vector3d[nbPoints];
                
                for (int i = 0; i < nbPoints; i++)
                {
                    modelForces[i] = new Vector3d();
                    contourNormals[i] = new Vector3d();
                    feedbackForces[i] = new Vector3d();
                }
            }
            
            counterClockWise = (getAlgebraicInterior() > 0);
        }
        
        updateMetaData();
    }
    
    protected void addPoint(Point3d p)
    {
        points.add(p);
    }
    
    /**
     * Checks whether the contour is self-intersecting. Depending on the given parameters, a
     * self-intersection can be considered as a loop or as a contour division.
     * 
     * @param minDistance
     *            the distance threshold between non-neighboring points to detect self-intersection
     * @return <code>null</code> if either no self-intersection is detected or if one of the new
     *         contours is too small, or an array of contours with 0 elements if both contours are
     *         too small, and 2 elements if both contours are viable
     */
    protected Polygon2D[] checkSelfIntersection(double minDistance)
    {
        int i = 0, j = 0, n = points.size();
        Point3d p_i = null, p_j = null;
        
        double divisionDistQ = boundingSphere.getRadius() * 0.1;
        divisionDistQ *= divisionDistQ;
        
        double minDistanceQ = minDistance;
        minDistanceQ *= minDistanceQ;
        
        boolean selfIntersection = false;
        
        loop:
        for (i = 0; i < n; i++)
        {
            p_i = points.get(i);
            Vector3d n_i = contourNormals[i];
            
            for (j = i + 2; j < n - 1; j++)
            {
                p_j = points.get(j);
                Vector3d n_j = contourNormals[j];
                
                double distQ = p_i.distanceSquared(p_j);
                
                if (distQ < minDistanceQ)
                {
                    // local self-intersection
                    // deal with the special case that i and j are 2 points away
                    
                    if (i == 0 && j == n - 2)
                    {
                        n--;
                        points.remove(n);
                        continue;
                    }
                    else if (i == 1 && j == n - 1)
                    {
                        points.remove(0);
                        n--;
                        continue;
                    }
                    else if (j == i + 2)
                    {
                        points.remove(i + 1);
                        n--;
                        continue;
                    }
                    
                    // a real self-intersection is happening
                    selfIntersection = true;
                    break loop;
                }
                else
                {
                    // self-intersection always involved opposite normals
                    if (n_i.dot(n_j) > -0.8) continue;
                    
                    // look for division
                    
                    // are points sufficiently close?
                    // => use the bounding radius / 2
                    if (distQ < divisionDistQ)
                    {
                        // are points located "in front" of each other?
                        if ((j - i) > (2 * n / 5) && (j - i) < (3 * n / 5))
                        {
                            Vector3d vi1 = new Vector3d(points.get((i + n - 4) % n));
                            Vector3d vi2 = new Vector3d(points.get((i + 4) % n));
                            vi1.sub(p_i);
                            vi2.sub(p_i);
                            vi1.normalize();
                            vi2.normalize();
                            vi1.cross(vi1, vi2);
                            double vi_lQ = vi1.lengthSquared();
                            // System.out.println(vi1.lengthSquared());
                            Vector3d vj1 = new Vector3d(points.get((j + n - 4) % n));
                            Vector3d vj2 = new Vector3d(points.get((j + 4) % n));
                            vj1.sub(p_j);
                            vj2.sub(p_j);
                            vj1.normalize();
                            vj2.normalize();
                            vj1.cross(vj1, vj2);
                            double vj_lQ = vj1.lengthSquared();
                            // System.out.println(vj1.lengthSquared());
                            
                            // curvature has to be at a relative minimum
                            if (vi_lQ < 0.5 && vj_lQ < 0.5) continue;
                            
                            // System.out.println(vi_lQ + "  |  " + vj_lQ);
                            
                            // a real self-intersection is happening
                            selfIntersection = true;
                            break loop;
                        }
                    }
                }
                
            }
        }
        
        if (!selfIntersection) return null;
        
        Point3d center = new Point3d();
        
        int nPoints = j - i;
        Polygon2D child1 = new Polygon2D(sampling, new SlidingWindow(this.convergence.getSize()));
        for (int p = 0; p < nPoints; p++)
        {
            Point3d pp = points.get(p + i);
            center.add(pp);
            child1.addPoint(pp);
        }
        center.scale(1.0 / nPoints);
        child1.setX(center.x);
        child1.setY(center.y);
        child1.setZ(center.z);
        child1.setT(getT());
        
        center.set(0, 0, 0);
        
        nPoints = i + n - j;
        Polygon2D child2 = new Polygon2D(sampling, new SlidingWindow(this.convergence.getSize()));
        for (int p = 0, pj = p + j; p < nPoints; p++, pj++)
        {
            Point3d pp = points.get(pj % n);
            center.add(pp);
            child2.addPoint(pp);
        }
        center.scale(1.0 / nPoints);
        child2.setX(center.x);
        child2.setY(center.y);
        child2.setZ(center.z);
        child2.setT(getT());
        
        // determine whether the intersection is a loop or a division
        // rationale: check the normal of the two colliding points (i & j)
        // if they point away from the junction => division
        // if they point towards the junction => loop
        
        Vector3d n_i = contourNormals[i];
        Vector3d i_j = new Vector3d(p_j.x - p_i.x, p_j.y - p_i.y, 0);
        
        if (n_i.dot(i_j) < 0)
        {
            // division => keep c1 and c2 if their size is ok
            
            double c1area = child1.getDimension(2), c2area = child2.getDimension(2);
            
            // if only one of the two children has a size lower than minArea, then the division
            // should be considered as an artifact loop, the other child thus is the new contour
            
            if (child1.points.size() < 10 || c1area < c2area / 5)
            {
                // remove c1 (too small)
                points.clear();
                points.addAll(child2.points);
                return null;
            }
            
            if (child2.points.size() < 10 || c2area < c1area / 5)
            {
                // remove c2 (too small)
                points.clear();
                points.addAll(child1.points);
                return null;
            }
            
            // keep both then...
            return new Polygon2D[] { child1, child2 };
        }
        else
        {
            // loop => keep only the contour with correct orientation
            // => the contour with a positive algebraic area
            
            if (child1.getAlgebraicInterior() < 0)
            {
                // c1 is the outer loop => keep it
                points.clear();
                points.addAll(child1.points);
                return null;
            }
            else
            {
                // c1 is the inner loop => keep c2
                points.clear();
                points.addAll(child2.points);
                return null;
            }
        }
    }
    
    @Override
    protected void clean()
    {
        // nothing to clean (everything should be garbage-collected)
    }
    
    @Override
    public Polygon2D clone()
    {
        return new Polygon2D(this);
    }
    
    /**
     * Update the axis constraint force, which adjusts the takes the final forces and normalize them
     * to keep the contour shape along its principal axis <br>
     * WARNING: this method directly update the array of final forces used to displace the contour
     * points. It should be used among the last to keep it most effective
     * 
     * @param weight
     */
    @Override
    void computeAxisForces(double weight)
    {
        Vector3d axis = new Vector3d();
        int s = (int) getDimension(0);
        
        // Compute the object axis as the vector between the two most distant
        // contour points
        // TODO this is not optimal, geometric moments should be used
        {
            double maxDistSq = 0;
            Vector3d vec = new Vector3d();
            
            for (int i = 0; i < s; i++)
            {
                Point3d vi = points.get(i);
                
                for (int j = i + 1; j < s; j++)
                {
                    Point3d vj = points.get(j);
                    
                    vec.sub(vi, vj);
                    double dSq = vec.lengthSquared();
                    
                    if (dSq > maxDistSq)
                    {
                        maxDistSq = dSq;
                        axis.set(vec);
                    }
                }
            }
            axis.normalize();
        }
        
        // To drive the contour along the main object axis, each displacement
        // vector is scaled by the scalar product between its normal and the main axis.
        {
            for (int i = 0; i < s; i++)
            {
                Vector3d normal = contourNormals[i];
                
                // dot product between normalized vectors ranges from -1 to 1
                double colinearity = Math.abs(normal.dot(axis)); // now from 0 to 1
                
                // goal: adjust the minimum using the weight, but keep max to 1
                double threshold = Math.max(colinearity, 1 - weight);
                
                if (normal != null) modelForces[i].scale(threshold);
            }
        }
    }
    
    @Override
    void computeBalloonForces(double weight)
    {
        int n = points.size();
        
        for (int i = 0; i < n; i++)
        {
            Vector3d f = modelForces[i];
            
            f.x += weight * contourNormals[i].x;
            f.y += weight * contourNormals[i].y;
        }
    }
    
    /**
     * Update edge term of the contour evolution according to the image gradient
     * 
     * @param weight
     * @param edge_data
     */
    @Override
    void computeEdgeForces(Sequence edgeData, int channel, double weight)
    {
        int n = points.size();
        
        Vector3d grad = new Vector3d();
        
        int width = edgeData.getWidth();
        int height = edgeData.getHeight();
        float[] data = edgeData.getDataXYAsFloat(0, (int) Math.round(getZ()), channel);
        
        for (int i = 0; i < n; i++)
        {
            Point3d p = points.get(i);
            Vector3d f = modelForces[i];
            
            // compute the gradient (2nd order)
            double nextX = getPixelValue(data, width, height, p.x + 0.5, p.y);
            if (nextX == 0) continue;
            double prevX = getPixelValue(data, width, height, p.x - 0.5, p.y);
            if (prevX == 0) continue;
            double nextY = getPixelValue(data, width, height, p.x, p.y + 0.5);
            if (nextY == 0) continue;
            double prevY = getPixelValue(data, width, height, p.x, p.y - 0.5);
            if (prevY == 0) continue;
            grad.set(nextX - prevX, nextY - prevY, 0.0);
            
            grad.scale(weight);
            
            f.add(grad);
        }
    }
    
    @Override
    void computeRegionForces(Sequence imageData, int channel, double weight, double sensitivity, double cin, double cout)
    {
        // sensitivity should be high for dim objects, low for bright objects...
        // ... but none of the following options work properly
        // sensitivity *= 1/(1+cin);
        // sensitivity = sensitivity * cin / (0.01 + cout);
        // sensitivity = sensitivity / (2 * Math.max(cout, cin));
        // sensitivity = sensitivity / (Math.log10(cin / cout));
        
        Point3d p;
        Vector3d f, norm, cvms = new Vector3d();
        double val, inDiff, outDiff, sum;
        int n = points.size();
        
        int width = imageData.getWidth();
        int height = imageData.getHeight();
        
        int myZ = (int) Math.round(getZ());
        float[] _data = imageData.getDataXYAsFloat(0, myZ, channel);
        if (_data == null) throw new IllegalArgumentException("Contour.getZ() = " + getZ() + "; Stack size = " + imageData.getSizeZ());
        
        for (int i = 0; i < n; i++)
        {
            p = points.get(i);
            f = modelForces[i];
            norm = contourNormals[i];
            
            // bounds check
            // if (p.x <= 1 || p.y <= 1 || p.x >= width - 2 || p.y >= height - 2) continue;
            
            val = getPixelValue(_data, width, height, p.x, p.y);
            
            inDiff = val - cin;
            inDiff *= inDiff;
            
            outDiff = val - cout;
            outDiff *= outDiff;
            
            sum = weight * sampling.getValue() * (sensitivity * outDiff) - (inDiff / sensitivity);
            
            cvms.scale(counterClockWise ? -sum : sum, norm);
            
            // if (cvms.dot(norm) < 0) cvms.scale(1.25);
            
            f.add(cvms);
        }
        
    }
    
    @Override
    void computeInternalForces(double weight)
    {
        if (feedbackForces == null) return;
        
        int n = points.size();
        
        if (n < 3) return;
        
        Vector3d f;
        Point3d prev, curr, next;
        
        weight /= sampling.getValue();
        
        // first point
        prev = points.get(n - 1);
        curr = points.get(0);
        next = points.get(1);
        
        f = feedbackForces[0];
        
        f.x += weight * (prev.x - 2 * curr.x + next.x);
        f.y += weight * (prev.y - 2 * curr.y + next.y);
        
        // middle points
        for (int i = 1; i < n - 1; i++)
        {
            f = feedbackForces[i];
            prev = points.get(i - 1);
            curr = points.get(i);
            next = points.get(i + 1);
            
            f.x += weight * (prev.x - 2 * curr.x + next.x);
            f.y += weight * (prev.y - 2 * curr.y + next.y);
        }
        
        // last point
        f = feedbackForces[n - 1];
        prev = points.get(n - 2);
        curr = points.get(n - 1);
        next = points.get(0);
        
        f.x += weight * (prev.x - 2 * curr.x + next.x);
        f.y += weight * (prev.y - 2 * curr.y + next.y);
    }
    
    void computeVolumeConstraint(double targetVolume)
    {
        // 1) compute the difference between target and current volume
        double volumeDiff = targetVolume - getDimension(2);
        // if (volumeDiff > 0): contour too small, should no longer shrink
        // if (volumeDiff < 0): contour too big, should no longer grow
        
        int n = points.size();
        
        for (int i = 0; i < n; i++)
        {
            Vector3d mf = modelForces[i];
            
            // 2) check whether the final force has same direction as the outer normal
            double forceNorm = mf.dot(contourNormals[i]);
            // if forces have same direction (forceNorm > 0): contour is growing
            // if forces have opposite direction (forceNorm < 0): contour is shrinking
            
            if (forceNorm * volumeDiff < 0)
            {
                // forceNorm and volumeDiff have opposite signs because:
                // - contour too small (volumeDiff > 0) and shrinking (forceNorm < 0)
                // or
                // - contour too large (volumeDiff < 0) and growing (forceNorm > 0)
                // => in both cases, constrain the final force accordingly
                mf.scale(1.0 / (1.0 + Math.abs(volumeDiff) / 10));
            }
        }
    }
    
    /**
     * Computes the feedback forces yielded by the penetration of the current contour into the
     * target contour
     * 
     * @param target
     *            the contour that is being penetrated
     * @return the number of actual point-mesh intersection tests
     */
    @Override
    int computeFeedbackForces(ActiveContour target)
    {
        Point3d targetCenter = new Point3d();
        target.boundingSphere.getCenter(targetCenter);
        
        double targetRadiusSq = target.boundingSphere.getRadius();
        targetRadiusSq *= targetRadiusSq;
        
        double penetration = 0;
        
        int tests = 0;
        int index = 0;
        
        for (Point3d p : points)
        {
            double distanceSq = p.distanceSquared(targetCenter);
            
            if (distanceSq < targetRadiusSq)
            {
                tests++;
                
                if ((penetration = target.getDistanceToEdge(p)) > 0)
                {
                    Vector3d feedbackForce = feedbackForces[index];
                    
                    // feedbackForce.scale(-penetration, contourNormals[index]);
                    feedbackForce.scaleAdd(-penetration * 0.5, contourNormals[index], feedbackForce);
                    
                    modelForces[index].scale(0.05);
                }
            }
            index++;
        }
        
        return tests;
    }
    
    private static void createEdge(ArrayList<Segment> segments, double xStart, double yStart, double xEnd, double yEnd)
    {
        double EPSILON = 0.00001;
        
        Point3d head = new Point3d(xStart, yStart, 0);
        Point3d tail = new Point3d(xEnd, yEnd, 0);
        
        if (segments.size() == 0)
        {
            segments.add(new Segment(head, tail));
            return;
        }
        
        int insertAtTailOf = -1, insertAtHeadOf = -1;
        
        for (int i = 0; i < segments.size(); i++)
        {
            if (tail.distance(segments.get(i).getHead()) <= EPSILON) insertAtHeadOf = i;
            else if (head.distance(segments.get(i).getTail()) <= EPSILON) insertAtTailOf = i;
        }
        
        if (insertAtTailOf >= 0)
        {
            if (insertAtHeadOf >= 0)
            {
                segments.get(insertAtHeadOf).addHead(segments.get(insertAtTailOf));
                segments.remove(insertAtTailOf);
            }
            else
            {
                segments.get(insertAtTailOf).addTail(tail);
            }
        }
        else if (insertAtHeadOf >= 0)
        {
            segments.get(insertAtHeadOf).addHead(head);
        }
        else
        {
            segments.add(new Segment(head, tail));
        }
    }
    
    /**
     * Perform a morphological dilation on the specified ROI by the given radius in each dimension
     * 
     * @param roi
     *            the ROI to dilate
     * @param xRadius
     *            the radius in pixels along X
     * @param yRadius
     *            the radius in pixels along X
     * @param zRadius
     *            the radius in pixels along Z (not used if <code>roi</code> is 2D)
     * @return a new, dilated ROI of type "area"
     */
    private static ROI2D dilateROI(ROI2D roi, int xRadius, int yRadius)
    {
        int rx = xRadius, rrx = rx * rx;
        int ry = yRadius, rry = ry * ry;
        
        BooleanMask2D m2 = roi.getBooleanMask(true);
        ROI2DArea r2 = new ROI2DArea(m2);
        
        r2.setT(roi.getT());
        r2.setZ(roi.getZ());
        r2.setC(roi.getC());
        
        r2.beginUpdate();
        
        for (Point p : m2.getContourPoints())
        {
            // Brute force
            for (int y = -ry; y <= ry; y++)
                for (int x = -rx; x <= rx; x++)
                    if (x * x / rrx + y * y / rry <= 1.0)
                    {
                        if (!m2.contains(p.x + x, p.y + y)) r2.addPoint(p.x + x, p.y + y);
                    }
        }
        r2.endUpdate();
        
        return r2;
    }
    
    /**
     * Calculates the 2D image value at the given real coordinates by bilinear interpolation
     * 
     * @param imageFloat
     *            the image to sample (must be of type {@link DataType#DOUBLE})
     * @param x
     *            the X-coordinate of the point
     * @param y
     *            the Y-coordinate of the point
     * @return the interpolated image value at the given coordinates
     */
    private float getPixelValue(float[] data, int width, int height, double x, double y)
    {
        int i = (int) Math.round(x);
        int j = (int) Math.round(y);
        
        // if (i <= 0 || i >= width - 1 || j <= 0 || j >= height - 1) return 0f;
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i > width - 2) i = width - 2;
        if (j > height - 2) j = height - 2;
        
        float value = 0;
        
        final int offset = i + j * width;
        final int offset_plus_1 = offset + 1; // saves 1 addition
        
        x -= i;
        y -= j;
        
        final double mx = 1 - x;
        final double my = 1 - y;
        
        value += mx * my * data[offset];
        value += x * my * data[offset_plus_1];
        value += mx * y * data[offset + width];
        value += x * y * data[offset_plus_1 + width];
        
        return value;
    }
    
    /**
     * Computes the algebraic area of the current contour. The returned value is negative if the
     * contour points are order clockwise and positive if ordered counter-clockwise. The contour's
     * surface is just the absolute value of this algebraic surface
     * 
     * @return
     */
    protected double getAlgebraicInterior()
    {
        int nm1 = points.size() - 1;
        double area = 0;
        
        // all points but the last
        for (int i = 0; i < nm1; i++)
        {
            Point3d p1 = points.get(i);
            Point3d p2 = points.get(i + 1);
            area += (p2.x * p1.y - p1.x * p2.y) * 0.5;
        }
        
        // last point
        Point3d p1 = points.get(nm1);
        Point3d p2 = points.get(0);
        area += (p2.x * p1.y - p1.x * p2.y) * 0.5;
        
        return area;
    }
    
    public double getDimension(int order)
    {
        if (points.size() <= 1) return 0;
        
        switch (order)
        {
        
        case 0: // number of points
        {
            return points.size();
        }
        
        case 1: // perimeter
        {
            int size = points.size();
            
            Point3d p1 = points.get(size - 1);
            Point3d p2 = points.get(0);
            
            double perimeter = p1.distance(p2);
            
            for (int i = 0; i < size - 1; i++)
            {
                // shift pair of points by one index
                p1 = p2;
                p2 = points.get(i + 1);
                perimeter += p1.distance(p2);
            }
            
            return perimeter;
        }
        case 2: // area
        {
            return Math.abs(getAlgebraicInterior());
        }
        default:
            throw new UnsupportedOperationException("Dimension " + order + " not implemented");
        }
    }
    
    /**
     * Tests whether the given point is inside the contour, and if so returns the penetration depth
     * of this point. <br>
     * This methods computes the number of intersections between the contour and a semi-infinite
     * line starting from the contour center and passing through the given point. The point is thus
     * considered inside if the number of intersections is odd (Jordan curve theorem).<br>
     * Implementation note: the AWT Line2D class only provides a "segment to segment" intersection
     * test instead of a "semi-infinite line to segment" test, meaning that one must "fake" a
     * semi-infinite line using a big segment. This is done by building a segment originating from
     * the given point and leaving in the opposite direction of the contour center. The full segment
     * can be written in algebraic coordinates as
     * 
     * <pre>
     * [PQ] where Q = P + n * CP
     * </pre>
     * 
     * , where n is chosen arbitrarily large.
     * 
     * @param c
     *            a contour
     * @param p
     *            a point to test
     * @return true if the point is inside the contour
     */
    public double getDistanceToEdge(Point3d p)
    {
        Point3d q = new Point3d(p.x + 10000 * (p.x - x), p.y + 10000 * (p.y - y), 0);
        
        int nb = 0;
        int nbPtsM1 = points.size() - 1;
        double dist = 0, minDist = Double.MAX_VALUE;
        
        // all points but the last
        for (int i = 0; i < nbPtsM1; i++)
        {
            Point3d p1 = points.get(i);
            Point3d p2 = points.get(i + 1);
            
            if (Line2D.linesIntersect(p1.x, p1.y, p2.x, p2.y, p.x, p.y, q.x, q.y))
            {
                nb++;
                dist = Line2D.ptLineDist(p1.x, p1.y, p2.x, p2.y, p.x, p.y);
                if (dist < minDist) minDist = dist;
            }
        }
        
        // last point
        Point3d p1 = points.get(nbPtsM1);
        Point3d p2 = points.get(0);
        if (Line2D.linesIntersect(p1.x, p1.y, p2.x, p2.y, p.x, p.y, q.x, q.y))
        {
            nb++;
            dist = Line2D.ptLineDist(p1.x, p1.y, p2.x, p2.y, p.x, p.y);
            if (dist < minDist) minDist = dist;
        }
        
        // return (nb % 2) == 0;
        return (nb % 2 == 1) ? minDist : 0.0;
        
    }
    
    @Override
    public Iterator<Point3d> iterator()
    {
        return points.iterator();
    }
    
    void move(ROI field, double timeStep)
    {
        Vector3d force = new Vector3d();
        double maxDisp = sampling.getValue() * timeStep;
        
        int n = points.size();
        
        if (modelForces == null || modelForces.length != n) return;
        
        for (int index = 0; index < n; index++)
        {
            Point3d p = points.get(index);
            
            // apply model forces if p lies within the area of interest
            if (field != null && field.contains(p.x, p.y, 0, 0, 0) && modelForces[index] != null)
            {
                if (p.x < 1 || p.x > field.getBounds5D().getSizeX() - 2) modelForces[index].scale(0.1);
                if (p.y < 1 || p.y > field.getBounds5D().getSizeY() - 2) modelForces[index].scale(0.1);
                force.set(modelForces[index]);
            }
            else
            {
                feedbackForces[index].scale(0.25);
            }
            
            // apply feedback forces all the time
            force.add(feedbackForces[index]);
            
            force.scale(timeStep);
            
            double disp = force.length();
            
            if (disp > maxDisp) force.scale(maxDisp / disp);
            
            p.add(force);
            
            force.set(0, 0, 0);
            
            modelForces[index].set(0, 0, 0);
            feedbackForces[index].set(0, 0, 0);
        }
        
        updateMetaData();
        
        // compute some convergence criterion
        
        if (convergence == null) return;
        
        convergence.push(getDimension(2));
        
    }
    
    @Override
    public boolean hasConverged(Operation operation, double epsilon)
    {
        Double value = convergence.computeCriterion(operation);
        return value != null && value <= epsilon / 100;
    }
    
    AnnounceFrame f = null; // new AnnounceFrame("ready");
                            
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        // only paint detections on the current frame
        if (getT() != canvas.getPositionT()) return;
        
        if (g != null)
        {
            // 2D viewer
            Rectangle2D viewBounds = ((IcyCanvas2D) canvas).canvasToImage(canvas.getBounds());
            g.clip(viewBounds);
            
            float fontSize = (float) canvas.canvasToImageLogDeltaX(30);
            g.setFont(new Font("Trebuchet MS", Font.BOLD, 10).deriveFont(fontSize));
            
            double stroke = Math.max(canvas.canvasToImageLogDeltaX(3), canvas.canvasToImageLogDeltaY(3));
            
            g.setColor(getColor());
            
            g.setStroke(new BasicStroke((float) stroke, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
            
            g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
            
            synchronized (path)
            {
                g.draw(path);
            }
            // this line displays the average intensity inside the object
            // g.drawString(StringUtil.toString(cin, 2), (float) getX() + 3, (float) getY());
        }
        else
        {
            // TODO other viewers
        }
    }
    
    int cpt = 0;
    
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
    @Override
    public void reSample(double minFactor, double maxFactor) throws TopologyException
    {
        if (getDimension(0) < 15) throw new TopologyException(this, new Polygon2D[] {});
        
        double minLength = sampling.getValue() * minFactor;
        double maxLength = sampling.getValue() * maxFactor;
        
        // optimization to avoid multiple points.size() calls (WARNING: n must
        // be updated manually whenever points is changed)
        int n = points.size();
        
        if (contourNormals == null)
        {
            // first pass: update normals once
            contourNormals = new Vector3d[n];
            for (int i = 0; i < n; i++)
                contourNormals[i] = new Vector3d();
            
            updateNormals();
        }
        
        Polygon2D[] children = cpt % 2 == 0 ? null : checkSelfIntersection(sampling.getValue());
        cpt++;
        
        if (children != null) throw new TopologyException(this, children);
        
        // update the number of total points
        n = points.size();
        boolean noChange = false;
        
        while (noChange == false)
        {
            noChange = true;
            
            // all points but the last
            for (int i = 0; i < n - 1; i++)
            {
                if (n < 4) throw new TopologyException(this, new Polygon2D[] {});
                
                Point3d pt1 = points.get(i);
                Point3d pt2 = points.get(i + 1);
                
                double distance = pt1.distance(pt2);
                
                if (distance < minLength)
                {
                    noChange = false;
                    pt2.set((pt1.x + pt2.x) * 0.5, (pt1.y + pt2.y) * 0.5, (pt1.z + pt2.z) * 0.5);
                    points.remove(i);
                    i--; // comes down to i-1+1 when looping
                    n--;
                }
                else if (distance > maxLength)
                {
                    noChange = false;
                    
                    points.add(i + 1, new Point3d((pt1.x + pt2.x) * 0.5, (pt1.y + pt2.y) * 0.5, (pt1.z + pt2.z) * 0.5));
                    i++; // comes down to i+=2 when looping
                    n++;
                }
            }
            
            // last point
            Point3d pt1 = points.get(n - 1);
            Point3d pt2 = points.get(0);
            
            if (pt1.distance(pt2) < minLength)
            {
                noChange = false;
                pt2.set((pt1.x + pt2.x) * 0.5, (pt1.y + pt2.y) * 0.5, (pt1.z + pt2.z) * 0.5);
                points.remove(n - 1);
                n--;
            }
            else if (pt1.distance(pt2) > maxLength)
            {
                noChange = false;
                points.add(new Point3d((pt1.x + pt2.x) * 0.5, (pt1.y + pt2.y) * 0.5, (pt1.z + pt2.z) * 0.5));
                n++;
            }
        }
        
        // re-sampling is done => update internal structures
        
        final int nbPoints = n;
        if (modelForces == null || modelForces.length != nbPoints)
        {
            modelForces = new Vector3d[nbPoints];
            contourNormals = new Vector3d[nbPoints];
            feedbackForces = new Vector3d[nbPoints];
            
            for (int i = 0; i < nbPoints; i++)
            {
                modelForces[i] = new Vector3d();
                contourNormals[i] = new Vector3d();
                feedbackForces[i] = new Vector3d();
            }
        }
        
        updateMetaData();
    }
    
    private void triangulate(ROI2DArea roi, double resolution) throws TopologyException
    {
        ArrayList<Segment> segments = new ArrayList<Segment>();
        
        Rectangle bounds = roi.getBounds();
        
        int grid = 1;// Math.max(1, (int) Math.round(resolution));
        double halfgrid = 0.5 * grid;
        
        int cubeWidth = grid;
        int cubeHeight = grid * bounds.width;
        int cubeDiag = cubeWidth + cubeHeight;
        
        boolean[] mask = roi.getBooleanMask(roi.getBounds(), true);
        // erase first line and first row to ensure closed contours
        java.util.Arrays.fill(mask, 0, bounds.width - 1, false);
        for (int o = 0; o < mask.length; o += bounds.width)
            mask[o] = false;
        
        for (int j = 0; j < bounds.height; j += grid)
            for (int i = 0, index = j * bounds.width; i < bounds.width; i += grid, index += grid)
            {
                // The image is divided into square cells containing two
                // triangles each:
                //
                // a---b---
                // |../|../
                // |./.|./.
                // |/..|/..
                // c---d---
                //
                // By convention I choose to turn around the object in a
                // clockwise fashion
                // Warning: to ensure connectivity, the objects must NOT touch
                // the image border, strange behavior may occur otherwise
                
                boolean a = mask[index];
                boolean b = (i + grid < bounds.width) && mask[index + cubeWidth];
                boolean c = (j + grid < bounds.height) && mask[index + cubeHeight];
                boolean d = (i + grid < bounds.width) && (j + grid < bounds.height) && mask[index + cubeDiag];
                
                // For each triangle, check for difference between image values
                // to determine the contour location
                // => there are 6 possible combinations in each triangle, that
                // is 12 per cube
                
                if (a != b)
                {
                    if (b == c) // diagonal edge
                    {
                        if (a == false) // b,c are inside
                        {
                            createEdge(segments, i, j + 0.5, i + halfgrid, j);
                            
                        }
                        else
                        // b,c are outside
                        {
                            createEdge(segments, i + halfgrid, j, i, j + halfgrid);
                            
                        }
                    }
                    else
                    // a = c -> vertical edge
                    {
                        if (a == false) // a,c are outside
                        {
                            createEdge(segments, i + halfgrid, j + halfgrid, i + halfgrid, j);
                            
                        }
                        else
                        // a,c are inside
                        {
                            createEdge(segments, i + halfgrid, j, i + halfgrid, j + halfgrid);
                            
                        }
                    }
                }
                else // a = b -> horizontal edge only if c is different
                if (a != c)
                {
                    if (a == false) // a,b are outside
                    {
                        createEdge(segments, i, j + halfgrid, i + halfgrid, j + halfgrid);
                        
                    }
                    else
                    // a,b are inside
                    {
                        createEdge(segments, i + halfgrid, j + halfgrid, i, j + halfgrid);
                        
                    }
                }
                
                if (c != d)
                {
                    if (b == c) // diagonal edge
                    {
                        if (c == false) // b,c are outside
                        {
                            createEdge(segments, i + halfgrid, j + grid, i + grid, j + halfgrid);
                            
                        }
                        else
                        // b,c are inside
                        {
                            createEdge(segments, i + grid, j + halfgrid, i + halfgrid, j + grid);
                            
                        }
                    }
                    else
                    // b = d -> vertical edge
                    {
                        if (c == false) // b,d are inside
                        {
                            createEdge(segments, i + halfgrid, j + grid, i + halfgrid, j + halfgrid);
                            
                        }
                        else
                        // b,d are outside
                        {
                            createEdge(segments, i + halfgrid, j + halfgrid, i + halfgrid, j + grid);
                            
                        }
                    }
                }
                else // c = d -> horizontal edge only if b is different
                if (b != c)
                {
                    if (b == false) // c,d are inside
                    {
                        createEdge(segments, i + halfgrid, j + halfgrid, i + grid, j + halfgrid);
                        
                    }
                    else
                    // c,d are outside
                    {
                        createEdge(segments, i + grid, j + halfgrid, i + halfgrid, j + halfgrid);
                        
                    }
                }
            }
        
        if (segments.size() == 0) return;
        
        for (Point3d p : segments.get(0))
        {
            p.x += bounds.x;
            p.y += bounds.y;
            p.z = roi.getZ();
            addPoint(p);
        }
        
        // at this point the triangulated contour has an actual resolution of halfgrid
        // if 2*resolution < desired_resolution, resample() will loop and destroy the contour
        // decimate the contour by a factor 2 recursively until 2*resolution >= desired_resolution
        
        double current_resolution_doubled = halfgrid * 2;
        while (current_resolution_doubled < resolution * 0.7)
        {
            for (int i = 0; i < points.size(); i++)
                points.remove(i);
            current_resolution_doubled *= 2;
        }
        
        reSample(0.8, 1.4);
    }
    
    protected void updateNormals()
    {
        int n = points.size();
        
        // first point
        {
            Point3d p1 = points.get(n - 1);
            Point3d p2 = points.get(1);
            contourNormals[0].normalize(new Vector3d(p2.y - p1.y, p1.x - p2.x, 0));
        }
        
        // middle points
        for (int i = 1; i < n - 1; i++)
        {
            Point3d p1 = points.get(i - 1);
            Point3d p2 = points.get(i + 1);
            contourNormals[i].normalize(new Vector3d(p2.y - p1.y, p1.x - p2.x, 0));
        }
        
        // last point
        {
            Point3d p1 = points.get(n - 2);
            Point3d p2 = points.get(0);
            contourNormals[n - 1].normalize(new Vector3d(p2.y - p1.y, p1.x - p2.x, 0));
        }
    }
    
    @Override
    protected void updateMetaData()
    {
        super.updateMetaData();
        updatePath();
    }
    
    private void updatePath()
    {
        if (Icy.getMainInterface().isHeadLess()) return;
        
        synchronized (path)
        {
            path.reset();
            
            int nbPoints = points.size();
            
            Point3d p = points.get(0);
            path.moveTo(p.x, p.y);
            
            for (int i = 1; i < nbPoints; i++)
            {
                p = points.get(i);
                path.lineTo(p.x, p.y);
            }
            
            path.closePath();
        }
    }
    
    @Override
    public ROI2D toROI()
    {
        return toROI(ROIType.AREA, null);
    }
    
    @Override
    public ROI2D toROI(ROIType type, Sequence sequence)
    {
        ROI2D roi;
        
        switch (type)
        {
        case AREA:
            roi = new ROI2DArea();
            ((ROI2DArea) roi).addShape(path);
            break;
        
        case POLYGON:
            List<Point2D> p2d = new ArrayList<Point2D>(points.size());
            for (Point3d p : this)
                p2d.add(new Point2D.Double(p.x, p.y));
            roi = new ROI2DPolygon(p2d);
            break;
        
        default:
            throw new IllegalArgumentException("ROI of type " + type + " cannot be exported yet");
        }
        
        roi.setT(t);
        return roi;
    }
    
    @Override
    public double computeAverageIntensity(Sequence imageData_float, int channel, BooleanMask3D mask)
    {
        int myZ = (int) Math.round(getZ());
        
        if (myZ == -1 && imageData_float.getSizeZ() == 1) myZ = 0;
        
        float[] _data = imageData_float.getDataXYAsFloat(0, myZ, channel);
        if (_data == null) throw new IllegalArgumentException("Contour.getZ() = " + getZ() + "; Stack size = " + imageData_float.getSizeZ());
        
        boolean[] _mask = (mask == null ? null : mask.mask.get(myZ).mask);
        
        int sizeX = imageData_float.getWidth();
        int sizeY = imageData_float.getHeight();
        
        // compute the interior mean intensity
        double inSum = 0, inCpt = 0;
        // double outSum = 0, outCpt = 0;
        Point3d minBounds = new Point3d();
        Point3d maxBounds = new Point3d();
        boundingBox.getLower(minBounds);
        boundingBox.getUpper(maxBounds);
        
        int minX = Math.max((int) minBounds.x - 10, -1);
        int maxX = Math.min((int) maxBounds.x + 10, sizeX);
        int minY = Math.max((int) minBounds.y - 10, 0);
        int maxY = Math.min((int) maxBounds.y + 10, sizeY);
        
        int n = points.size();
        Point3d p1 = null, p2 = null;
        TreeSet<Integer> crosses = new TreeSet<Integer>();
        
        for (int j = minY; j < maxY; j++)
        {
            int offset = j * sizeX + minX;
            
            crosses.clear();
            crosses.add(minX);
            
            for (int p = 0; p < n - 1; p++)
            {
                p1 = points.get(p);
                p2 = points.get(p + 1);
                
                if (j > Math.min(p1.y, p2.y) && j < Math.max(p1.y, p2.y)) crosses.add((int) Math.round((p1.x + p2.x) * 0.5));
            }
            p1 = points.get(0);
            
            if (j > Math.min(p1.y, p2.y) && j < Math.max(p1.y, p2.y)) crosses.add((int) Math.round((p1.x + p2.x) * 0.5));
            
            crosses.add(maxX);
            
            int nC = crosses.size();
            
            if (nC == 2)
            {
                // for (int i = minX + 1; i < maxX - 1; i++, offset++)
                // {
                // outSum += _data[offset];
                // outCpt++;
                // }
            }
            else
            {
                boolean in = false;
                
                Iterator<Integer> it = crosses.iterator();
                int start = it.next();
                if (start < 0) start = 0;
                
                while (it.hasNext())
                {
                    int end = it.next();
                    if (end > sizeX) end = sizeX;
                    
                    if (in)
                    {
                        for (int i = start; i < end; i++, offset++)
                        {
                            if (offset < 0 || offset >= _data.length) continue;
                            _mask[offset] = true;
                            inSum += _data[offset];
                            inCpt++;
                        }
                    }
                    else
                    {
                        // for (int i = start; i < end; i++, offset++)
                        // {
                        // outSum += _data[offset];
                        // outCpt++;
                        // }
                        offset += end - start + 1;
                    }
                    
                    start = end;
                    in = !in;
                }
            }
        }
        // cout = outSum / outCpt;
        return inSum / inCpt;
    }
    
    @Override
    public void toSequence(Sequence output, double value)
    {
        int myZ = (int) Math.round(getZ());
        int myT = (int) Math.round(getT());
        
        Object _mask = output.getDataXY(myT, myZ, 0);
        
        int sizeX = output.getWidth();
        int sizeY = output.getHeight();
        
        // compute the interior mean intensity
        Point3d minBounds = new Point3d();
        Point3d maxBounds = new Point3d();
        boundingBox.getLower(minBounds);
        boundingBox.getUpper(maxBounds);
        
        int minX = Math.max((int) minBounds.x - 2, -10);
        int maxX = Math.min((int) maxBounds.x + 2, sizeX + 10);
        int minY = Math.max((int) minBounds.y - 1, 0);
        int maxY = Math.min((int) maxBounds.y + 1, sizeY);
        
        int n = points.size();
        Point3d p1 = null, p2 = null;
        TreeSet<Integer> crosses = new TreeSet<Integer>();
        
        for (int j = minY; j < maxY; j++)
        {
            int lineOffset = j * sizeX;
            int offset = lineOffset + (minX < 0 ? 0 : minX);
            
            crosses.clear();
            crosses.add(minX);
            
            for (int p = 0; p < n - 1; p++)
            {
                p1 = points.get(p);
                p2 = points.get(p + 1);
                
                if (j >= Math.min(p1.y, p2.y) && j <= Math.max(p1.y, p2.y))
                {
                    // crosses.add((int) Math.round((p1.x + p2.x) * 0.5));
                    int cross = (int) Math.round(p1.x + ((j - p1.y) * (p2.x - p1.x) / (p2.y - p1.y)));
                    if (crosses.contains(cross))
                    {
                        crosses.remove(cross);
                        Array1DUtil.setValue(_mask, lineOffset + cross, value);
                    }
                    else
                    {
                        crosses.add(cross);
                    }
                }
            }
            p1 = points.get(0);
            
            if (j >= Math.min(p1.y, p2.y) && j <= Math.max(p1.y, p2.y))
            {
                // crosses.add((int) Math.round((p1.x + p2.x) * 0.5));
                int cross = (int) Math.round(p1.x + ((j - p1.y) * (p2.x - p1.x) / (p2.y - p1.y)));
                if (crosses.contains(cross))
                {
                    crosses.remove(cross);
                    Array1DUtil.setValue(_mask, lineOffset + cross, value);
                }
                else
                {
                    crosses.add(cross);
                }
            }
            
            crosses.add(maxX);
            
            int nC = crosses.size();
            
            if (nC > 2)
            {
                boolean in = false;
                
                Iterator<Integer> it = crosses.iterator();
                int start = it.next();
                if (start < 0) start = 0;
                
                while (it.hasNext())
                {
                    int end = it.next();
                    if (end > sizeX) end = sizeX;
                    
                    if (in)
                    {
                        for (int i = start; i < end; i++, offset++)
                        {
                            // if (start < 0 || end >= sizeY) continue;
                            Array1DUtil.setValue(_mask, offset, value);
                        }
                    }
                    else offset += end - start;
                    
                    start = end;
                    in = !in;
                }
            }
        }
    }
    
    /**
     * Perform a raster scan on the contour, and calculate various information based on the
     * specified parameters. The raster scan algorithm is a 2D version of the line-scan algorithm
     * developed in: <i>Dufour et al., 3D active meshes: fast discrete deformable models for cell
     * tracking in 3D time-lapse microscopy. IEEE Transactions on Image Processing 20, 2011</i>
     * 
     * @param updateLocalMask
     *            indicates whether the local boolean mask should be updated (this is used to get a
     *            {@link BooleanMask3D} version of this mesh)
     * @param imageData
     *            (set to <code>null</code> if not needed) a sequence that will be used to compute
     *            the average intensity inside the mesh (note that the T and C have to be different
     *            than <code>-1</code> if the sequence has more than 1 time point and 1 channel)
     * @param averageIntensity
     *            (only used if <code>imageData</code> is provided) a variable that will be hold the
     *            average image intensity inside the contour after the scan is complete
     * @param bufferData
     *            (set to <code>null</code> if not needed) a byte buffer data in the form [z][xy] of
     *            same dimensions as the image data that will be filled with value <code>1</code>
     *            inside the contour
     */
    public void rasterScan(final boolean updateLocalMask, final Sequence imageData, VarDouble averageIntensity, final BooleanMask3D imageMask)
    {
        
    }
}
