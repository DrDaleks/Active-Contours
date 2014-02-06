package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROI3D;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.rectangle.Rectangle3D;
import icy.util.StringUtil;
import icy.vtk.VtkUtil;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Stack;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import javax.media.j3d.BoundingSphere;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import plugins.adufour.activecontours.Mesh3D.Face;
import plugins.adufour.activecontours.Mesh3D.Vertex;
import plugins.adufour.activemeshes.mesh.ContourSplittingException;
import plugins.adufour.activemeshes.mesh.Mesh;
import plugins.adufour.activemeshes.mesh.MeshException;
import plugins.adufour.activemeshes.mesh.MeshSplittingException;
import plugins.adufour.ezplug.EzException;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DEllipse;
import plugins.kernel.roi.roi2d.ROI2DPolygon;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.kernel.roi.roi2d.ROI2DShape;
import plugins.kernel.roi.roi3d.ROI3DArea;
import vtk.vtkActor;
import vtk.vtkDoubleArray;
import vtk.vtkPoints;
import vtk.vtkPolyData;
import vtk.vtkPolyDataMapper;
import vtk.vtkRenderer;

public class Mesh3D extends ActiveContour
{
    final ArrayList<Face>   faces    = new ArrayList<Face>();
    
    final ArrayList<Vertex> vertices = new ArrayList<Vertex>();
    
    private VTKMesh         vtkMesh;
    
    protected Mesh3D(ActiveContours owner, EzVarDouble contour_resolution, EzVarInteger contour_minArea, SlidingWindow convergenceWindow)
    {
        super(owner, contour_resolution, contour_minArea, convergenceWindow);
        
        setColor(Color.getHSBColor((float) Math.random(), 0.8f, 0.9f));
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    public Mesh3D(Mesh3D contour)
    {
        this(contour.owner, contour.contour_resolution, contour.contour_minArea, new SlidingWindow(contour.convergence.size));
        // TODO later on, clone the resolution variables as well if needed
        
        setX(contour.x);
        setY(contour.y);
        setZ(contour.z);
        
        setColor(contour.getColor());
        points.ensureCapacity(contour.points.size());
        
        for (Point3d p : contour.points)
            addPoint(new Point3d(p));
        
        // in any case, don't forget to close the path
        updatePath();
        
        counterClockWise = (getAlgebraicInterior() > 0);
    }
    
    public Mesh3D(ActiveContours owner, EzVarDouble contour_resolution, EzVarInteger contour_minArea, SlidingWindow convergenceWindow, ROI2D roi)
    {
        this(owner, contour_resolution, contour_minArea, convergenceWindow);
        
        if (!(roi instanceof ROI2DEllipse) && !(roi instanceof ROI2DRectangle) && !(roi instanceof ROI2DPolygon) && !(roi instanceof ROI2DArea))
        {
            throw new EzException("Wrong ROI type. Only Rectangle, Ellipse, Polygon and Area are supported", true);
        }
        
        setZ(roi.getZ());
        
        boolean contourOK = false;
        double area = 0;
        
        if (roi instanceof ROI2DArea)
        {
            try
            {
                triangulate((ROI2DArea) roi, contour_resolution.getValue());
            }
            catch (TopologyException e)
            {
                roi = new ROI2DRectangle(roi.getBounds2D());
            }
            
            if (points.size() == 0)
            {
                // replace by ellipse
                roi = new ROI2DEllipse(roi.getBounds2D());
            }
            else
            {
                area = getAlgebraicInterior();
                if (Math.abs(area) < contour_minArea.getValue())
                {
                    roi = new ROI2DEllipse(roi.getBounds2D());
                }
                else
                {
                    contourOK = true;
                }
            }
        }
        
        if (!contourOK)
        {
            points.clear();
            
            // convert the ROI into a linked list of points
            
            double[] segment = new double[6];
            
            PathIterator pathIterator = ((ROI2DShape) roi).getPathIterator(null, 0.1);
            
            // first segment is necessarily a "move to" operation
            
            pathIterator.currentSegment(segment);
            addPoint(new Point3d(segment[0], segment[1], 0));
            
            while (!pathIterator.isDone())
            {
                if (pathIterator.currentSegment(segment) == PathIterator.SEG_LINETO)
                {
                    addPoint(new Point3d(segment[0], segment[1], 0));
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
            
            updatePath();
            
            area = getAlgebraicInterior();
            
        }
        
        counterClockWise = (area > 0);
    }
    
    protected void addPoint(Point3d p)
    {
        points.add(p);
    }
    
    @Override
    public Mesh3D clone()
    {
        return new Mesh3D(this);
    }
    
    /**
     * Checks whether the contour is self-intersecting. Depending on the given parameters, a
     * self-intersection can be considered as a loop or as a contour division.
     * 
     * @param minSpacing
     *            the distance threshold between non-neighboring points to detect self-intersection
     * @param minArea
     *            if a self-intersection is detected, this value specifies if the new contours are
     *            kept or eliminated
     * @return null if either no self-intersection is detected or if one of the new contours is too
     *         small, or an array of Contour2Ds with 0 elements if both contours are too small, and
     *         2 elements if both contours are viable
     */
    private Polygon2D[] checkSelfIntersection(double minSpacing, double minArea)
    {
        int i = 0, j = 0, n = points.size();
        Point3d p_i = null, p_j = null;
        
        boolean intersection = false;
        
        for (i = 0; i < n; i++)
        {
            p_i = points.get(i);
            
            for (j = i + 2; j < n - 1; j++)
            {
                p_j = points.get(j);
                
                intersection = (p_i.distance(p_j) < minSpacing);
                
                if (intersection)
                {
                    // deal with the special case that i and j are 2 points away
                    
                    if (i == 0 && j == n - 2)
                    {
                        n--;
                        points.remove(n);
                        intersection = false;
                    }
                    else if (i == 1 && j == n - 1)
                    {
                        points.remove(0);
                        n--;
                        intersection = false;
                    }
                    else if (j == i + 2)
                    {
                        points.remove(i + 1);
                        n--;
                        intersection = false;
                    }
                }
                
                if (intersection) break;
            }
            if (intersection) break;
        }
        
        if (!intersection) return null;
        
        Point3d center = new Point3d();
        
        int nPoints = j - i;
        Polygon2D c1 = new Polygon2D(this.owner, contour_resolution, contour_minArea, new SlidingWindow(this.convergence.size));
        for (int p = 0; p < nPoints; p++)
        {
            Point3d pp = points.get(p + i);
            center.add(pp);
            c1.addPoint(pp);
        }
        center.scale(1.0 / nPoints);
        c1.setX(center.x);
        c1.setY(center.y);
        c1.setZ(center.z);
        c1.setT(getT());
        
        center.set(0, 0, 0);
        
        nPoints = i + n - j;
        Polygon2D c2 = new Polygon2D(this.owner, contour_resolution, contour_minArea, new SlidingWindow(this.convergence.size));
        for (int p = 0, pj = p + j; p < nPoints; p++, pj++)
        {
            Point3d pp = points.get(pj < n ? pj : pj - n);
            center.add(pp);
            c2.addPoint(pp);
        }
        center.scale(1.0 / nPoints);
        c2.setX(center.x);
        c2.setY(center.y);
        c2.setZ(center.z);
        c2.setT(getT());
        
        // determine whether the intersection is a loop or a division
        // rationale: check the normal of the two colliding points (i & j)
        // if they point away from the junction => division
        // if they point towards the junction => loop
        
        Vector3d n_i = contourNormals[i];
        Vector3d i_j = new Vector3d(p_j.x - p_i.x, p_j.y - p_i.y, 0);
        
        if (n_i.dot(i_j) < 0)
        {
            // division => keep c1 and c2 if their size is ok
            
            double c1area = c1.getDimension(2), c2area = c2.getDimension(2);
            
            // if only one of the two children has a size lower than minArea, then the division
            // should be considered as an artifact loop, the other child thus is the new contour
            
            if (c1area > minArea)
            {
                if (c2area > minArea) return new Polygon2D[] { c1, c2 };
                
                points.clear();
                points.addAll(c1.points);
                
                return null;
            }
            else
            {
                if (c2area < minArea) return new Polygon2D[] {};
                
                points.clear();
                points.addAll(c2.points);
                
                return null;
            }
        }
        else
        {
            // loop => keep only the contour with correct orientation
            // => the contour with a positive algebraic area
            
            if (c1.getAlgebraicInterior() < 0)
            {
                // c1 is the outer loop => keep it
                points.clear();
                points.addAll(c1.points);
                return null;
            }
            else
            {
                // c1 is the inner loop => keep c2
                points.clear();
                points.addAll(c2.points);
                return null;
            }
        }
    }
    
    private double cin = 0, cout = 0;
    
    @Override
    public double computeAverageIntensity(Sequence imageData_float, int channel, Sequence buffer_data)
    {
        IcyBufferedImage data = imageData_float.getImage(0, 0);
        float[] _data = data.getDataXYAsFloat(channel);
        
        int sizeX = data.getWidth();
        int sizeY = data.getHeight();
        
        // if (bufferGraphics == null) bufferGraphics = buffer.createGraphics();
        
        // fill the contour contents in the buffer
        // bufferGraphics.fill(path);
        
        // compute the interior mean intensity
        double inSum = 0, inCpt = 0, outSum = 0, outCpt = 0;
        Rectangle bounds = path.getBounds();
        
        int minX = Math.max(bounds.x - bounds.width, 0);
        int maxX = Math.min(bounds.x + bounds.width * 2, sizeX);
        int minY = Math.max(bounds.y - bounds.height, 0);
        int maxY = Math.min(bounds.y + bounds.height * 2, sizeY);
        
        for (int j = minY; j < maxY; j++)
        {
            int offset = j * sizeX + minX;
            
            for (int i = minX; i < maxX; i++, offset++)
            {
                if (bounds.contains(i, j) && path.contains(i, j))
                {
                    double value = _data[offset];
                    inSum += value;
                    inCpt++;
                }
                else
                {
                    double value = _data[offset];
                    outSum += value;
                    outCpt++;
                }
            }
        }
        
        cin = inSum / inCpt;
        cout = outSum / outCpt;
        return cin;
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
            updateNormalsIfNeeded();
            
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
        updateNormalsIfNeeded();
        
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
        updateNormalsIfNeeded();
        
        Vector3d grad = new Vector3d();
        
        int width = edgeData.getWidth();
        int height = edgeData.getHeight();
        float[] data = edgeData.getDataXYAsFloat(0, (int) Math.round(getZ()), channel);
        
        for (int i = 0; i < n; i++)
        {
            Point3d p = points.get(i);
            Vector3d f = modelForces[i];
            // compute the gradient (2nd order)
            double nextX = getPixelValue(data, width, height, p.x + 1, p.y);
            double prevX = getPixelValue(data, width, height, p.x - 1, p.y);
            double nextY = getPixelValue(data, width, height, p.x, p.y + 1);
            double prevY = getPixelValue(data, width, height, p.x, p.y - 1);
            grad.set(nextX - prevX, nextY - prevY, 0.0);
            grad.scale(weight);
            f.add(grad);
        }
    }
    
    @Override
    void computeRegionForces(Sequence imageData, int channel, double weight, double sensitivity, double cin, double cout)
    {
        // uncomment these 2 lines for mean-based information
        sensitivity = sensitivity / Math.max(cout, cin);
        
        updateNormalsIfNeeded();
        
        Point3d p;
        Vector3d f, norm;
        double val, inDiff, outDiff, sum;
        int n = points.size();
        
        int width = imageData.getWidth();
        int height = imageData.getHeight();
        float[] data = imageData.getDataXYAsFloat(0, (int) Math.round(getZ()), channel);
        
        for (int i = 0; i < n; i++)
        {
            p = points.get(i);
            f = modelForces[i];
            
            norm = contourNormals[i];
            
            // bounds check
            if (p.x < 1 || p.y < 1 || p.x >= width - 1 || p.y >= height - 1) continue;
            
            // invert the following lines for mean-based information
            val = getPixelValue(data, width, height, p.x, p.y);
            inDiff = val - cin;
            outDiff = val - cout;
            sum = weight * contour_resolution.getValue() * (sensitivity * (outDiff * outDiff) - (inDiff * inDiff));
            
            if (counterClockWise)
            {
                f.x -= sum * norm.x;
                f.y -= sum * norm.y;
            }
            else
            {
                f.x += sum * norm.x;
                f.y += sum * norm.y;
            }
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
        
        weight /= contour_resolution.getValue();
        
        // first point
        prev = points.get(n - 1);
        curr = points.get(0);
        next = points.get(1);
        
        f = feedbackForces[0];
        f.x += weight * (prev.x + next.x - 2 * curr.x);
        f.y += weight * (prev.y + next.y - 2 * curr.y);
        
        // middle points
        for (int i = 1; i < n - 1; i++)
        {
            f = feedbackForces[i];
            prev = points.get(i - 1);
            curr = points.get(i);
            next = points.get(i + 1);
            
            f.x += weight * (prev.x + next.x - 2 * curr.x);
            f.y += weight * (prev.y + next.y - 2 * curr.y);
        }
        
        // last point
        f = feedbackForces[n - 1];
        prev = points.get(n - 2);
        curr = points.get(n - 1);
        next = points.get(0);
        
        f.x += weight * (prev.x - 2 * curr.x + next.x);
        f.y += weight * (prev.y - 2 * curr.y + next.y);
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
        updateNormalsIfNeeded();
        
        BoundingSphere targetSphere = target.getBoundingSphere();
        Point3d targetCenter = new Point3d();
        targetSphere.getCenter(targetCenter);
        double targetRadius = targetSphere.getRadius();
        
        double penetration = 0;
        
        int tests = 0;
        
        int index = 0;
        for (Point3d p : points)
        {
            Vector3d feedbackForce = feedbackForces[index];
            
            double distance = p.distance(targetCenter);
            
            if (distance < targetRadius)
            {
                tests++;
                
                if ((penetration = target.contains(p)) > 0)
                {
                    feedbackForce.scaleAdd(-penetration, contourNormals[index], feedbackForce);
                }
            }
            index++;
        }
        
        return tests;
    }
    
    public double contains(Point3d p)
    {
        // FIXME this is *THE* hot-spot
        
        double epsilon = 1.0e-12;
        
        Vector3d edge1 = new Vector3d(), edge2 = new Vector3d();
        Vector3d vp = new Vector3d(), vt = new Vector3d(), vq = new Vector3d();
        Vector3d ray = new Vector3d();
        
        // the input point belongs to another mesh
        // 1) trace a ray from that point outwards (away from my center)
        // 2) count the intersections with my boundary
        // 3) if the number is odd, the point is inside
        // 4) if the point is inside, measure the distance to the edge
        
        // ray.negate(vTest.normal); // was supper buggy upon close contacts !
        
        // create a ray from the mass center to the given point
        ray.sub(p, getMassCenter());
        
        double penetration = Double.MAX_VALUE;
        int crossCount = 0;
        
        double det, u, v, distance;
        
        for (Face f : faces)
        {
            Point3d v1 = vertices.get(f.v1).position;
            Point3d v2 = vertices.get(f.v2).position;
            Point3d v3 = vertices.get(f.v3).position;
            
            edge1.sub(v3, v1);
            edge2.sub(v2, v1);
            
            vp.cross(ray, edge2);
            
            det = edge1.dot(vp);
            
            if (det < epsilon) continue;
            
            vt.sub(p, v1);
            
            u = vt.dot(vp);
            
            if (u < 0 || u > det) continue;
            
            vq.cross(vt, edge1);
            
            v = ray.dot(vq);
            
            if (v < 0 || u + v > det) continue;
            
            distance = edge2.dot(vq) / det;
            
            if (distance < 0) continue;
            
            if (penetration > distance) penetration = distance;
            
            crossCount++;
        }
        
        return crossCount % 2 == 1 ? penetration : 0;
    }
    
    public double getDimension(int order)
    {
        switch (order)
        {
            case 0: // number of points
            {
                int nbPts = 0;
                
                for (Vertex v : vertices)
                    if (v != null) nbPts++;
                
                return nbPts;
            }
            case 1:
            case 2: // surface and volume are both computed as algebraic sums
                // over each face
            {
                double surface = 0, volume = 0;
                
                Vector3d v12 = new Vector3d();
                Vector3d v13 = new Vector3d();
                Vector3d cross = new Vector3d();
                
                Vector3d v1 = new Vector3d();
                Vector3d v2 = new Vector3d();
                Vector3d v3 = new Vector3d();
                
                for (Face f : faces)
                {
                    v1.set(vertices.get(f.v1).position);
                    v2.set(vertices.get(f.v2).position);
                    v3.set(vertices.get(f.v3).position);
                    
                    v12.sub(v2, v1);
                    v13.sub(v3, v1);
                    
                    cross.cross(v12, v13);
                    
                    double surf = cross.length() * 0.5f;
                    
                    if (order == 1)
                    {
                        surface += surf;
                    }
                    else
                    {
                        cross.normalize();
                        
                        // from the divergence theorem:
                        // V = sum_i( surface_of_i * normal_of_i 'dot' any_point_in_i ) / 3
                        
                        volume += surf * cross.dot(v1) / 3.0;
                    }
                    
                }
                
                double value = (order == 1) ? surface : volume;
                return value;
            }
            default:
                throw new UnsupportedOperationException("Dimension " + order + " not implemented");
        }
    }
    
    public Point3d getMassCenter()
    {
        return new Point3d(getX(), getY(), getZ());
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
        final int i = (int) Math.round(x);
        final int j = (int) Math.round(y);
        
        if (i < 0 || i >= width - 1 || j < 0 || j >= height - 1) return 0;
        
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
    
    @Override
    protected void initFrom(ActiveContour contour)
    {
        int n = 0;
        for (Point3d p : contour)
        {
            addPoint(p);
            x += p.x;
            y += p.y;
            n++;
        }
        x /= n;
        y /= n;
    }
    
    @Override
    public Iterator<Point3d> iterator()
    {
        // return a "tweaked" iterator that will skip null entries automatically
        
        return new Iterator<Point3d>()
        {
            Iterator<Vertex> vertexIterator       = vertices.iterator();
            
            Vertex           next;
            
            boolean          hasNext;
            
            boolean          hasNextWasCalledOnce = false;
            
            @Override
            public void remove()
            {
                vertexIterator.remove();
            }
            
            @Override
            public Point3d next()
            {
                hasNextWasCalledOnce = false;
                return next.position;
            }
            
            @Override
            public boolean hasNext()
            {
                if (hasNextWasCalledOnce) return hasNext;
                
                hasNextWasCalledOnce = true;
                
                if (!vertexIterator.hasNext()) return false;
                
                do
                {
                    next = vertexIterator.next();
                } while (next == null && vertexIterator.hasNext());
                
                return next != null;
            }
        };
    }
    
    void move(ROI field, double timeStep)
    {
        double maxDisp = contour_resolution.getValue() * timeStep;
        
        // prepare to compute the new mass center
        double xSum = 0, ySum = 0, zSum = 0;
        int nbPoints = 0;
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            // apply model forces if p lies within the area of interest
            if (field != null && !field.contains(v.position.x, v.position.y, v.position.z, 0, 0))
            {
                // vertex is out of bounds => don't drive it anywhere
                v.drivingForces.set(0, 0, 0);
            }
            
            // add feedback forces
            v.drivingForces.add(v.feedbackForces);
            
            // compute the final force
            v.drivingForces.scale(timeStep);
            
            // compute the length of the final force
            double disp = v.drivingForces.length();
            
            // adjust the force if too strong (~ CFL condition)
            if (disp > maxDisp) v.drivingForces.scale(maxDisp / disp);
            
            // move the vertex
            v.position.add(v.drivingForces);
            
            // accumulate coordinates
            xSum += v.position.x;
            ySum += v.position.y;
            zSum += v.position.z;
            nbPoints++;
            
            // reset forces
            v.drivingForces.set(0, 0, 0);
            v.feedbackForces.set(0, 0, 0);
        }
        
        setX(xSum / nbPoints);
        setY(ySum / nbPoints);
        setZ(zSum / nbPoints);
        
        normalsNeedUpdate = true;
        boundingSphereNeedsUpdate = true;
        
        // compute some convergence criterion
        
        if (convergence == null) return;
        
        if (cout != 0.0)
        {
            convergence.add(cin);
        }
        else
        {
            convergence.add(getDimension(2));
        }
    }
    
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        if (g != null)
        {
            // 2D viewer
            
            if (getT() != canvas.getPositionT() || !enabled || g == null) return;
            
            float fontSize = (float) canvas.canvasToImageLogDeltaX(30);
            g.setFont(new Font("Trebuchet MS", Font.BOLD, 10).deriveFont(fontSize));
            
            double stroke = Math.max(canvas.canvasToImageLogDeltaX(3), canvas.canvasToImageLogDeltaY(3));
            
            g.setStroke(new BasicStroke((float) stroke, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
            
            g.setColor(getColor());
            
            // FIXME draw something in 2D (points? raster mesh?)
        }
        else
        {
            // 3D viewer
            
            // nothing to do here (the VTK viewer should update itself automatically)
        }
    }
    
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
        if (getDimension(2) < contour_minArea.getValue()) throw new TopologyException(this, new Polygon2D[] {});
        
        double minLength = contour_resolution.getValue() * minFactor;
        double maxLength = contour_resolution.getValue() * maxFactor;
        
        // if there are 2 faces only in the mesh, it should be destroyed
        
        if (faces.size() == 2)
        {
            throw new TopologyException(this, null);
        }
        
        // FIXME proper self-intersection (not just division)
        
        boolean noChange = false;
        
        int cpt = -1;
        
        while (noChange == false)
        {
            noChange = true;
            
            cpt++;
            
            // we are looking for 2 faces f1 = a-b-c1 and f2 = b-a-c2
            // such that they share an edge a-b that is either
            // - shorter than the low-threshold (resolution * min)
            // or
            // - longer than the high-threshold (resolution * max)
            
            for (int i = 0; i < faces.size(); i++)
            {
                boolean split = false, merge = false;
                
                Face f1 = faces.get(i);
                Integer v1 = 0, v2 = 0, f1v3 = 0;
                
                // Check all edges of this face
                
                Integer[] f1v123 = { f1.v1, f1.v2, f1.v3 };
                Integer[] f1v231 = { f1.v2, f1.v3, f1.v1 };
                Integer[] f1v312 = { f1.v3, f1.v1, f1.v2 };
                
                for (int v = 0; v < 3; v++)
                {
                    v1 = f1v123[v];
                    v2 = f1v231[v];
                    f1v3 = f1v312[v];
                    
                    double edgeLength = vertices.get(v1).position.distance(vertices.get(v2).position);
                    
                    if (edgeLength < minLength)
                    {
                        merge = true;
                        break;
                    }
                    
                    if (edgeLength > maxLength)
                    {
                        split = true;
                        break;
                    }
                }
                
                if (split == merge) // both are false
                {
                    continue; // to the next face
                }
                
                // If the code runs here, a change will occur
                noChange = false;
                
                // => we need the second associated face for that edge
                
                Face f2 = null;
                Integer f2v3 = -1;
                
                // start from the current face => O(N)
                for (int j = i + 1; j < faces.size(); j++)
                {
                    f2 = faces.get(j);
                    
                    // check if f2 contains [v1-v2]
                    if (v1.compareTo(f2.v1) == 0 && v2.compareTo(f2.v3) == 0)
                    {
                        f2v3 = f2.v2;
                        break;
                    }
                    else if (v1.compareTo(f2.v2) == 0 && v2.compareTo(f2.v1) == 0)
                    {
                        f2v3 = f2.v3;
                        break;
                    }
                    else if (v1.compareTo(f2.v3) == 0 && v2.compareTo(f2.v2) == 0)
                    {
                        f2v3 = f2.v1;
                        break;
                    }
                }
                
                // CASE 0: THE MESH IS INCONSISTENT //
                
                if (f2v3.compareTo(0) < 0)
                {
                    // should never happen:
                    // if f2 does not exist, then [v1-v2] only belongs to a single face (f1)
                    // => this means the mesh is inconsistent (not closed)
                    
                    System.err.println("Problem in face " + i + ":");
                    System.err.print("  " + f1.v1.intValue() + " : ");
                    for (Integer nn : vertices.get(v1).neighbors)
                        System.err.print(nn.intValue() + "  ");
                    System.err.println();
                    System.err.print("  " + f1.v2.intValue() + " : ");
                    for (Integer nn : vertices.get(v2).neighbors)
                        System.err.print(nn.intValue() + "  ");
                    System.err.println();
                    System.err.print("  " + f1.v3.intValue() + " : ");
                    for (Integer nn : vertices.get(f1.v3).neighbors)
                        System.err.print(nn.intValue() + "  ");
                    System.err.println();
                    
                    System.err.println("The mesh is removed from further computations");
                    
                    throw new TopologyException(this, null);
                }
                
                // CASE 1: MERGE //
                
                else if (merge)
                {
                    // Check first if the edge to merge is at the base of a tetrahedron
                    // if so, delete the whole tetrahedron
                    
                    if (vertices.get(f1v3).neighbors.size() == 3)
                    {
                        deleteTetrahedron(f1v3, v1, v2);
                        return;
                    }
                    if (vertices.get(f2v3).neighbors.size() == 3)
                    {
                        deleteTetrahedron(f2v3, v2, v1);
                        return;
                    }
                    
                    Vertex vx1 = vertices.get(v1);
                    Vertex vx2 = vertices.get(v2);
                    
                    for (Integer n : vx1.neighbors)
                        if (vx2.neighbors.contains(n) && n.compareTo(f1v3) != 0 && n.compareTo(f2v3) != 0)
                        {
                            splitContourAtVertices(v1, v2, n);
                            return;
                        }
                    
                    // Here, the normal merge operation can be implemented
                    
                    // remove the 2 faces
                    faces.remove(f1);
                    faces.remove(f2);
                    
                    // move v1 to the middle of v1-v2
                    vx1.position.interpolate(vx2.position, 0.5);
                    
                    // remove v2 from its neighborhood...
                    for (Integer n : vx2.neighbors)
                    {
                        Vertex vxn = vertices.get(n);
                        if (vxn == null) continue;
                        vxn.neighbors.remove(v2);
                        
                        // ...and add v2's neighbors to v1
                        // except for f1v3 and f2v3 and ... v1 !
                        if (n.compareTo(f1v3) != 0 && n.compareTo(f2v3) != 0 && n.compareTo(v1) != 0)
                        {
                            vx1.neighbors.add(n);
                            vxn.neighbors.add(v1);
                        }
                    }
                    
                    // all the faces pointing to v2 must now point to v1
                    for (Face f : faces)
                    {
                        if (f.v1.compareTo(v2) == 0)
                        {
                            f.v1 = v1;
                        }
                        else if (f.v2.compareTo(v2) == 0)
                        {
                            f.v2 = v1;
                        }
                        else if (f.v3.compareTo(v2) == 0)
                        {
                            f.v3 = v1;
                        }
                    }
                    
                    // delete everything
                    vertices.set(v2, null);
                    
                }
                
                // CASE 2: SPLIT //
                
                else
                {
                    // 1) remove the old faces
                    faces.remove(f1);
                    faces.remove(f2);
                    
                    // 2) remove the vertices from their mutual neighbor list
                    vertices.get(v1).neighbors.remove(v2);
                    vertices.get(v2).neighbors.remove(v1);
                    
                    // check if the edge should just be inverted
                    
                    Vertex vx1 = vertices.get(f1v3);
                    Vertex vx2 = vertices.get(f2v3);
                    
                    if (!vx1.neighbors.contains(f2v3) && vx1.distanceTo(vx2) < maxLength)
                    {
                        // just invert the edge
                        
                        // 3) create the two new faces
                        addFace(f1v3, v1, f2v3);
                        addFace(f2v3, v2, f1v3);
                    }
                    else
                    {
                        // split the edge and its faces
                        
                        // 3) create a vertex in the middle of the edge
                        Integer c = addVertexBetween(v1, v2);
                        
                        // 4) create 4 new faces around the new vertex
                        addFace(v1, c, f1v3);
                        addFace(f1v3, c, v2);
                        addFace(v1, f2v3, c);
                        addFace(c, f2v3, v2);
                    }
                }
                
                break;
            }
            
            // prevent infinite loop
            if (cpt > vertices.size()) noChange = true;
        }
        
        updatePath();
        
        normalsNeedUpdate = true;
        boundingSphereNeedsUpdate = true;
    }
    
    public void addFace(Integer v1, Integer v2, Integer v3) throws TopologyException
    {
        if (v1.compareTo(v2) == 0 || v1.compareTo(v3) == 0 || v2.compareTo(v3) == 0) throw new TopologyException(this, null);
        
        if (v1 < 0) v1 = -v1 - 1;
        if (v2 < 0) v2 = -v2 - 1;
        if (v3 < 0) v3 = -v3 - 1;
        
        {
            Vertex v = vertices.get(v1);
            if (!v.neighbors.contains(v2)) v.neighbors.add(v2);
            if (!v.neighbors.contains(v3)) v.neighbors.add(v3);
        }
        {
            Vertex v = vertices.get(v2);
            if (!v.neighbors.contains(v1)) v.neighbors.add(v1);
            if (!v.neighbors.contains(v3)) v.neighbors.add(v3);
        }
        {
            Vertex v = vertices.get(v3);
            if (!v.neighbors.contains(v1)) v.neighbors.add(v1);
            if (!v.neighbors.contains(v2)) v.neighbors.add(v2);
        }
        
        faces.add(new Face(v1, v2, v3));
    }
    
    /**
     * If it doens't exist, adds a new vertex in the center of the specified vertices to the vertex
     * list
     * 
     * @param point
     *            the point to add
     * @return A signed integer i defined as follows: <br>
     *         i>=0 : vertex is new and stored at position i <br>
     *         i<0 : vertex already existed at position -(i+1)
     */
    public Integer addVertex(Point3d point)
    {
        Integer index, nullIndex = -1;
        
        for (index = 0; index < vertices.size(); index++)
        {
            Vertex v = vertices.get(index);
            
            if (v == null)
            {
                // The current position in the vertex list is null.
                // To avoid growing the data array, this position
                // can be reused to store a new vertex if needed
                nullIndex = index;
                continue;
            }
            
            if (v.position.distanceSquared(point) < contour_resolution.getValue() * 0.00001) return -index - 1;
        }
        
        // if code runs until here, the vertex must be created
        Vertex v = new Vertex(new Point3d(point));
        
        // if there is a free spot in the ArrayList, use it
        if (nullIndex >= 0)
        {
            index = nullIndex;
            vertices.set(index, v);
        }
        else
        {
            vertices.add(v);
        }
        
        return index;
    }
    
    /**
     * If it doens't exist, adds a new vertex in the center of the specified vertices to the vertex
     * list
     * 
     * @param v1
     *            the first vertex index
     * @param v2
     *            the second vertex index
     * @return A signed integer i defined as follows: <br>
     *         i>=0 : vertex is new and stored at position i <br>
     *         i<0 : vertex already existed at position -(i+1)
     */
    private Integer addVertexBetween(Integer v1, Integer v2)
    {
        Point3d newPosition = new Point3d(vertices.get(v1).position);
        newPosition.interpolate(vertices.get(v2).position, 0.5);
        
        return addVertex(newPosition);
    }
    
    /**
     * Deletes a tetrahedron from the mesh, and fill the hole with a new face
     * 
     * @param topv
     *            the vertex at the top of the tetrahedron
     * @param v1
     *            one of the three vertices at the base of the tetrahedron
     * @param v2
     *            another of the vertices at the base of the tetrahedron
     */
    private void deleteTetrahedron(Integer topv, Integer v1, Integer v2)
    {
        // find the third bottom vertex
        Integer v3 = -1;
        for (Integer n : vertices.get(topv).neighbors)
        {
            if (n.compareTo(v1) == 0 || n.compareTo(v2) == 0) continue;
            v3 = n;
            break;
        }
        
        // remove the top vertex from the neighborhood
        vertices.get(v1).neighbors.remove(topv);
        vertices.get(v2).neighbors.remove(topv);
        vertices.get(v3).neighbors.remove(topv);
        
        // delete the top vertex
        vertices.set(topv, null);
        
        // find the three faces and delete them
        for (int i = 0; i < faces.size(); i++)
        {
            Face f = faces.get(i);
            if (f.v1.compareTo(topv) == 0 || f.v2.compareTo(topv) == 0 || f.v3.compareTo(topv) == 0) faces.remove(i--);
        }
        
        // create the new face to close the hole left by the tetrahedron
        faces.add(new Face(v1, v2, v3));
    }
    
    private void extractVertices(Integer seedIndex, ArrayList<Integer> visitedIndices, ArrayList<Vertex> oldPoints, ArrayList<Vertex> newPoints)
    {
        Stack<Integer> seeds = new Stack<Integer>();
        seeds.add(seedIndex);
        
        while (!seeds.isEmpty())
        {
            Integer seed = seeds.pop();
            
            if (visitedIndices.contains(seed)) continue;
            
            visitedIndices.add(seed);
            Vertex v = oldPoints.get(seed);
            newPoints.set(seed, v);
            
            for (int i = 0; i < v.neighbors.size(); i++)
            {
                Integer n = v.neighbors.get(i);
                if (oldPoints.get(n) == null)
                {
                    v.neighbors.remove(n);
                    i--;
                    continue;
                }
                seeds.push(n);
            }
        }
    }
    
    private void extractFaces(ArrayList<Vertex> pointsList, ArrayList<Face> sourceFacesList, ArrayList<Face> targetFacesList)
    {
        for (int i = 0; i < sourceFacesList.size(); i++)
        {
            Face f = sourceFacesList.get(i);
            
            if (pointsList.get(f.v1) != null || pointsList.get(f.v2) != null || pointsList.get(f.v3) != null)
            {
                targetFacesList.add(f);
                sourceFacesList.remove(i--);
            }
        }
    }
    
    /**
     * Splits the current contour using the 'cutting' face defined by the given vertices. <br>
     * 
     * <pre>
     * How this works:
     *  - separate all vertices on each side of the cutting face (without considering the vertices of the cutting face), 
     *  - separate all faces touching at least one vertex of each group (will include faces touching the cutting face),
     *  - create a contour with each group of vertices and faces,
     *  - add the cutting face and its vertices to each created contour
     * </pre>
     * 
     * @param v1
     * @param v2
     * @param v3
     * @throws TopologyException
     * 
     */
    private void splitContourAtVertices(Integer v1, Integer v2, Integer v3) throws TopologyException
    {
        Mesh3D[] children = new Mesh3D[2];
        
        ArrayList<Integer> visitedIndexes = new ArrayList<Integer>(vertices.size());
        visitedIndexes.add(v1);
        visitedIndexes.add(v2);
        visitedIndexes.add(v3);
        
        int seed;
        
        for (int child = 0; child < 2; child++)
        {
            // pick any non-null and non-visited vertex as seed
            for (seed = 0; seed < vertices.size(); seed++)
                if (vertices.get(seed) != null && !visitedIndexes.contains(seed)) break;
            
            if (seed == vertices.size())
            {
                System.err.println("Mesh splitting error (pass " + (child + 1) + "): no valid seed found");
                throw new TopologyException(this, null);
            }
            
            ArrayList<Face> newFaces = new ArrayList<Face>();
            ArrayList<Vertex> newPoints = new ArrayList<Vertex>(vertices.size());
            for (int i = 0; i < vertices.size(); i++)
                newPoints.add(null);
            
            extractVertices(seed, visitedIndexes, vertices, newPoints);
            extractFaces(newPoints, faces, newFaces);
            
            // Add the vertices of the cutting face
            for (Integer v : new Integer[] { v1, v2, v3 })
            {
                Vertex vx = vertices.get(v);
                
                if (vx == null) System.err.println(v.intValue() + " is null");
                
                // create a clone for each vertex (position and neighbors)
                Vertex newV = new Vertex(vertices.get(v));
                
                // check the neighborhood to remove the neighbors that belong to
                // the other mesh
                // (these neighbors will point to null in the current point
                // list)
                for (int i = 0; i < newV.neighbors.size(); i++)
                {
                    Integer n = newV.neighbors.get(i);
                    if (n.compareTo(v1) != 0 && n.compareTo(v2) != 0 && n.compareTo(v3) != 0 && newPoints.get(n) == null)
                    {
                        newV.neighbors.remove(n);
                        i--;
                    }
                }
                newPoints.set(v, newV);
            }
            
            for (Face f : newFaces)
                if (f.contains(v1) && f.contains(v2))
                {
                    // if the edge v1-v2 appears counter-clockwisely in f,
                    // the new face must be clockwise and vice-versa
                    newFaces.add(f.isCounterClockwise(v1, v2) ? new Face(v1, v3, v2) : new Face(v1, v2, v3));
                    break;
                }
            
            Mesh3D newContour = new Mesh(newPoints, newFaces, mesh.getResolution(), mesh.getVTKMesh() != null);
            
            newContour.setT(getT());
            newContour.updateMassCenter();
            
            if (newContour.getDimension(2) >= contour_minArea.getValue()) children[child] = newContour;
        }
        
        if (children[0] == null)
        {
            if (children[1] == null) throw new TopologyException(this, null);
        }
        else
        {
            if (children[1] != null) throw new TopologyException(this, children);
        }
    }
    
    protected void updateNormalsIfNeeded()
    {
        if (!normalsNeedUpdate) return;
        
        Vector3d v31 = new Vector3d();
        Vector3d v12 = new Vector3d();
        Vector3d v23 = new Vector3d();
        
        for (Face f : faces)
        {
            // Accumulate face normals in each vertex
            
            Vertex v1 = vertices.get(f.v1);
            Vertex v2 = vertices.get(f.v2);
            Vertex v3 = vertices.get(f.v3);
            
            v31.sub(v1.position, v3.position);
            v12.sub(v2.position, v1.position);
            v23.sub(v3.position, v2.position);
            
            // normal at v1 = [v1 v2] ^ [v1 v3] = [v3 v1] ^ [v1 v2]
            v1.normal.x += v31.y * v12.z - v31.z * v12.y;
            v1.normal.y += v31.z * v12.x - v31.x * v12.z;
            v1.normal.z += v31.x * v12.y - v31.y * v12.x;
            
            // normal at v2 = [v2 v3] ^ [v2 v1] = [v1 v2] ^ [v2 v3]
            v2.normal.x += v12.y * v23.z - v12.z * v23.y;
            v2.normal.y += v12.z * v23.x - v12.x * v23.z;
            v2.normal.z += v12.x * v23.y - v12.y * v23.x;
            
            // normal at v3 = [v3 v1] ^ [v3 v2] = [v2 v3] ^ [v3 v1]
            v3.normal.x += v23.y * v31.z - v23.z * v31.y;
            v3.normal.y += v23.z * v31.x - v23.x * v31.z;
            v3.normal.z += v23.x * v31.y - v23.y * v31.x;
        }
        
        // Normalize the accumulated normals
        for (Vertex v : vertices)
            if (v != null) v.normal.normalize();
        
        normalsNeedUpdate = false;
    }
    
    private void updatePath()
    {
        double x, y;
        synchronized (path)
        {
            path.reset();
            
            int nbPoints = points.size();
            
            Point3d p = points.get(0);
            path.moveTo(p.x, p.y);
            x = p.x;
            y = p.y;
            
            for (int i = 1; i < nbPoints; i++)
            {
                p = points.get(i);
                path.lineTo(p.x, p.y);
                x += p.x;
                y += p.y;
                // TODO
                // Point3d pp = points.get(i-1);
                // path.quadTo((pp.x + p.x)/2, (pp.y + p.y)/2, p.x, p.y);
            }
            path.closePath();
            
            setX(x / nbPoints);
            setY(y / nbPoints);
        }
    }
    
    /**
     * Returns the coordinates of the Contour's "boundary box". The boundary box is defined as the
     * smallest rectangle (2D or 3D) area which entirely contains the Contour.
     * 
     * @param min
     *            a point that will represent the top left-hand corner of the box
     * @param max
     *            a point that will represent the bottom right-hand corner of the box
     */
    public Rectangle3D computeBounds3D()
    {
        double x = Double.MAX_VALUE, y = Double.MAX_VALUE, z = Double.MAX_VALUE;
        double sizeX = 0, sizeY = 0, sizeZ = 0;
        
        for (Point3d p : this)
        {
            if (p == null) continue;
            if (p.x < x) x = p.x;
            if (p.x > sizeX) sizeX = p.x;
            if (p.y < y) y = p.y;
            if (p.y > sizeY) sizeY = p.y;
            if (p.z < z) z = p.z;
            if (p.z > sizeZ) sizeZ = p.z;
        }
        
        return new Rectangle3D.Double(x, y, z, sizeX, sizeY, sizeZ);
    }
    
    @Override
    public ROI3D toROI()
    {
        ROI3DArea area = new ROI3DArea();
        
        Rectangle3D bounds = computeBounds3D();
        
        final int scanLine = (int) bounds.getSizeZ() + 1;
        
        final int minX = Math.max(0, (int) Math.floor(boxMin.x / resolution.x) - 1);
        final int minY = Math.max(0, (int) Math.floor(boxMin.y / resolution.y) - 1);
        final int minZ = Math.max(0, (int) Math.floor(boxMin.z / resolution.z) - 1);
        
        final int maxY = Math.min(dimensions.y - 1, (int) Math.ceil(boxMax.y / resolution.y));
        final int maxZ = Math.min(dimensions.z - 1, (int) Math.ceil(boxMax.z / resolution.z));
        
        if (maxZ < minZ)
        {
            String message = "Cannot rasterize sequence\n";
            message += "Make sure the following values seem correct (and report this error if not):\n";
            message += "Mesh bounding box: " + boxMin.x + " x " + boxMin.y + " x " + boxMin.z + "\n";
            message += "Image resolution: " + StringUtil.toString(resolution.x, 3) + " x " + StringUtil.toString(resolution.y, 3) + " x " + StringUtil.toString(resolution.z, 3) + "\n";
            throw new IllegalArgumentException(message);
        }
        
        final double epsilon = 1.0e-12;
        
        final Vector3d direction = new Vector3d(1, 0, 0);
        
        short[][] mask_Z_XY = unsignedShortSequence.getDataXYZAsShort(t, 0);
        
        ArrayList<Future<?>> tasks = new ArrayList<Future<?>>(maxZ - minZ + 1);
        
        for (int k = minZ; k <= maxZ; k++) // '<' changed to '<=' from single-thread version
        {
            final short[] maskSlice = (mask_Z_XY != null) ? mask_Z_XY[k] : null;
            
            final int slice = k;
            
            tasks.add(multiThreadService.submit(new Runnable()
            {
                @Override
                public void run()
                {
                    Vector3d edge1 = new Vector3d(), edge2 = new Vector3d();
                    Vector3d vp = new Vector3d(), vt = new Vector3d(), vq = new Vector3d();
                    
                    ArrayList<Integer> crossDistancesList = new ArrayList<Integer>(4);
                    
                    Point3d origin = new Point3d(minX * resolution.x, minY * resolution.y, slice * resolution.z);
                    
                    for (int j = minY; j < maxY; j++, origin.y += resolution.y)
                    {
                        int lineOffset = j * scanLine;
                        
                        crossDistancesList.clear();
                        int crosses = 0;
                        
                        for (Face f : faces)
                        {
                            Point3d v1 = vertices.get(f.v1).position;
                            Point3d v2 = vertices.get(f.v2).position;
                            Point3d v3 = vertices.get(f.v3).position;
                            
                            if (origin.y < v1.y && origin.y < v2.y && origin.y < v3.y) continue;
                            if (origin.z < v1.z && origin.z < v2.z && origin.z < v3.z) continue;
                            if (origin.y > v1.y && origin.y > v2.y && origin.y > v3.y) continue;
                            if (origin.z > v1.z && origin.z > v2.z && origin.z > v3.z) continue;
                            
                            edge1.sub(v2, v1);
                            edge2.sub(v3, v1);
                            
                            vp.cross(direction, edge2);
                            
                            double det = edge1.dot(vp);
                            
                            if (Math.abs(det) < epsilon) continue;
                            
                            double inv_det = 1.0 / det;
                            
                            vt.sub(origin, v1);
                            double u = vt.dot(vp) * inv_det;
                            if (u < 0 || u > 1.0) continue;
                            
                            vq.cross(vt, edge1);
                            double v = direction.dot(vq) * inv_det;
                            if (v < 0.0 || u + v > 1.0) continue;
                            
                            double distance = edge2.dot(vq) * inv_det;
                            
                            Integer distPx = minX + (int) Math.round(distance / resolution.x);
                            
                            if (distPx < 0) continue;
                            
                            if (!crossDistancesList.contains(distPx))
                            {
                                crossDistancesList.add(distPx);
                                crosses++;
                            }
                            else
                            {
                                // if distPx already exists, then the mesh "slices" the same voxel
                                // twice
                                // (round-off error thus gives the same distance for the 2 crosses)
                                // => don't add the (same) cross distance
                                // instead, remove the existing one to discard all crosses for that
                                // voxel
                                crossDistancesList.remove(distPx);
                                crosses--;
                            }
                        }
                        
                        // ignore the following cases:
                        // - crosses = 0 --> the ray does not cross the mesh
                        // - crosses = 1 --> the ray touches the mesh on an edge (only the first
                        // distance
                        // was recorded)
                        if (crosses < 2) continue;
                        
                        if (crosses == 2)
                        {
                            // optimization for the most frequent case: a ray crosses
                            // twice (in & out)
                            
                            int cross1 = crossDistancesList.get(0);
                            int cross2 = crossDistancesList.get(1);
                            
                            // the ray touches the mesh on an edge
                            
                            if (cross1 < cross2) // crosses are ordered
                            {
                                for (int i = cross1; i < cross2; i++)
                                    if (maskSlice != null) maskSlice[i + lineOffset] = id;
                            }
                            else
                            // invert the order
                            {
                                for (int i = cross2; i < cross1; i++)
                                    if (maskSlice != null) maskSlice[i + lineOffset] = id;
                            }
                        }
                        else
                        {
                            int[] crossDistances = new int[crosses];
                            
                            for (int item = 0; item < crosses; item++)
                                crossDistances[item] = crossDistancesList.get(item).intValue();
                            
                            java.util.Arrays.sort(crossDistances);
                            
                            int nbSegments = crossDistances.length / 2;
                            int crossOffset, start, stop;
                            
                            for (int segment = 0; segment < nbSegments; segment++)
                            {
                                crossOffset = segment << 1;
                                start = crossDistances[crossOffset];
                                stop = crossDistances[crossOffset + 1];
                                
                                for (int i = start; i < stop; i++)
                                    if (maskSlice != null) maskSlice[i + lineOffset] = id;
                            }
                        }
                    }
                }
            }));
        }
        
        try
        {
            for (Future<?> task : tasks)
                task.get();
        }
        catch (InterruptedException e)
        {
            e.printStackTrace();
        }
        catch (ExecutionException e)
        {
            e.printStackTrace();
        }
        
        area.setT(t);
        return area;
    }
    
    /**
     * Structural element of triangular mesh
     * 
     * @author Alexandre Dufour
     * 
     */
    public static class Vertex
    {
        public final Point3d            position;
        
        public final ArrayList<Integer> neighbors;
        
        public final Vector3d           drivingForces;
        
        public final Vector3d           feedbackForces;
        
        public final Vector3d           normal;
        
        public double                   distanceToCenter = 0;
        
        public Vertex(Vertex v)
        {
            this(v.position, v.neighbors);
        }
        
        public Vertex(Point3d position)
        {
            this(position, 0);
        }
        
        public Vertex(Point3d position, int nbNeighbors)
        {
            this.position = new Point3d(position);
            if (Double.isNaN(position.x))
            {
                System.out.println("new vertex is NaN");
            }
            this.neighbors = new ArrayList<Integer>(nbNeighbors);
            this.drivingForces = new Vector3d();
            this.feedbackForces = new Vector3d();
            this.normal = new Vector3d();
        }
        
        public Vertex(Point3d position, ArrayList<Integer> neighbors)
        {
            this(position, neighbors.size());
            for (Integer i : neighbors)
                this.neighbors.add(new Integer(i));
        }
        
        public double distanceTo(Vertex v)
        {
            return distanceTo(v.position);
        }
        
        public double distanceTo(Point3d p)
        {
            return position.distance(p);
        }
        
        public String toString()
        {
            return "Vertex {" + position.x + "," + position.y + "," + position.z + "}, " + neighbors.size() + " neighbors";
        }
    }
    
    /**
     * Structural element of triangular mesh
     * 
     * @author Alexandre Dufour
     * 
     */
    public static class Face
    {
        Integer v1, v2, v3;
        
        /**
         * Constructs a new mesh face with given vertices indices. Note that vertices must be given
         * in counter-clockwise order.
         * 
         * @param v1
         *            the first vertex index
         * @param v2
         *            the second vertex index
         * @param v3
         *            the third vertex index
         */
        Face(Integer v1, Integer v2, Integer v3)
        {
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }
        
        /**
         * Returns true if the specified vertex index is referred to by this face
         * 
         * @param v
         *            the vertex index to look for
         * @return true if the index is contained in this face, false otherwise
         */
        public boolean contains(Integer v)
        {
            return (v.compareTo(v1) == 0 || v.compareTo(v2) == 0 || v.compareTo(v3) == 0);
        }
        
        /**
         * Calculates the area of this face.
         * 
         * @param points
         *            - vertex list
         * @return the area of the face
         */
        public double getArea(ArrayList<Vertex> points)
        {
            Vector3d a = new Vector3d(points.get(v1).position);
            Vector3d b = new Vector3d(points.get(v2).position);
            
            a.sub(points.get(v3).position);
            b.sub(points.get(v3).position);
            
            a.cross(a, b);
            return 0.5 * a.length();
        }
        
        /**
         * Calculates the area of this face with unit sphere pre-scaling.
         * 
         * @author Michael Reiner
         * @param points
         *            - vertex list
         * @return the area of the face
         */
        public double getArea(ArrayList<Vertex> points, Point3d center)
        {
            Point3d p1 = new Point3d(points.get(v1).position);
            Point3d p2 = new Point3d(points.get(v2).position);
            Point3d p3 = new Point3d(points.get(v3).position);
            
            // if the vertices are not on the unit-sphere, project the points onto it (for correct
            // integral)
            p1.scale(1. / p1.distance(center));
            p2.scale(1. / p2.distance(center));
            p3.scale(1. / p3.distance(center));
            
            Vector3d a = new Vector3d(p1);
            Vector3d b = new Vector3d(p2);
            
            a.sub(p3);
            b.sub(p3);
            
            a.cross(a, b);
            return 0.5 * a.length();
        }
        
        /**
         * Returns true if the specified vertex indices are ordered counter-clockwisely in the
         * current face. An exception is raised if the indices do not belong to this face
         * 
         * @param v1
         *            the first vertex
         * @param v2
         *            the second vertex
         * @return true if the edge v1-v2 is counter-clockwise for the current face, false otherwise
         * @throws IllegalArgumentException
         *             thrown if the specified vertex indices do not belong to this face
         */
        public boolean isCounterClockwise(Integer v1, Integer v2) throws IllegalArgumentException
        {
            if (v1.compareTo(this.v1) == 0)
            {
                if (v2.compareTo(this.v2) == 0) return true;
                
                if (v2.compareTo(this.v3) == 0) return false;
                
                throw new IllegalArgumentException("Vertex index " + v2 + " does not belong to this face");
            }
            
            if (v1.compareTo(this.v2) == 0)
            {
                if (v2.compareTo(this.v3) == 0) return true;
                
                if (v2.compareTo(this.v1) == 0) return false;
                
                throw new IllegalArgumentException("Vertex index " + v2 + " does not belong to this face");
            }
            
            if (v1.compareTo(this.v3) == 0)
            {
                if (v2.compareTo(this.v1) == 0) return true;
                
                if (v2.compareTo(this.v2) == 0) return false;
                
                throw new IllegalArgumentException("Vertex index " + v2 + " does not belong to this face");
            }
            
            throw new IllegalArgumentException("Vertex index " + v1 + " does not belong to this face");
        }
        
        public Point3d[] getCoords(ArrayList<Vertex> points, Point3d[] out)
        {
            out[0] = points.get(v1).position;
            out[1] = points.get(v2).position;
            out[2] = points.get(v3).position;
            
            return out;
        }
        
        public String toString()
        {
            return "Face [" + v1 + "," + v2 + "," + v3 + "]";
        }
    }
    
    private class VTKMesh
    {
        private final vtkDoubleArray    vtkVerticesCoords = new vtkDoubleArray();
        private final vtkPoints         vtkVerticesInfo   = new vtkPoints();
        
        public final vtkPolyData        polyData          = new vtkPolyData();
        private final vtkPolyDataMapper polyDataMapper    = new vtkPolyDataMapper();
        public final vtkActor           actor             = new vtkActor();
        private vtkRenderer             renderer;
        
        VTKMesh(Mesh3D mesh)
        {
            this.vtkVerticesCoords.SetNumberOfComponents(3);
            this.vtkVerticesInfo.SetData(this.vtkVerticesCoords);
            this.polyData.SetPoints(this.vtkVerticesInfo);
            
            this.polyDataMapper.SetInput(this.polyData);
            
            this.actor.SetMapper(this.polyDataMapper);
        }
        
        public void update()
        {
            
            if (topology.updating) return;
            
            synchronized (topology)
            {
                int nFaces = faces.size();
                
                if (nFaces == 0) return;
                
                topology.lock = true;
                
                double[] vertices = new double[this.mesh.vertices.size() * 3];
                
                int cIndex = 0;
                for (Vertex vertex : Mesh3D.this.vertices)
                {
                    if (vertex == null)
                    {
                        cIndex += 3;
                    }
                    else
                    {
                        vertices[(cIndex++)] = vertex.position.x;
                        vertices[(cIndex++)] = vertex.position.y;
                        vertices[(cIndex++)] = vertex.position.z;
                    }
                }
                this.vtkVerticesCoords.SetJavaArray(vertices);
                
                int[] faces = new int[nFaces * 4];
                
                int vIndex = 0;
                for (Face face : Mesh3D.this.faces)
                {
                    faces[(vIndex++)] = 3;
                    faces[(vIndex++)] = face.v1.intValue();
                    faces[(vIndex++)] = face.v2.intValue();
                    faces[(vIndex++)] = face.v3.intValue();
                }
                
                this.polyData.SetPolys(VtkUtil.getCells(nFaces, faces));
                
                topology.lock = false;
            }
        }
        
        public void clean()
        {
            if ((this.renderer != null) && (this.actor != null)) this.renderer.RemoveActor(this.actor);
        }
    }
}
