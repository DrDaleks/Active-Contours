package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.Overlay;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI3D;
import icy.sequence.Sequence;
import icy.type.DataType;

import java.awt.Graphics2D;
import java.util.HashMap;
import java.util.Iterator;

import javax.vecmath.Point3d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

import org.w3c.dom.Node;

import plugins.adufour.activecontours.ActiveContours.ROIType;
import plugins.adufour.activecontours.SlidingWindow.Operation;
import plugins.adufour.roi.mesh.MeshTopologyException;
import plugins.adufour.roi.mesh.Vertex3D;
import plugins.adufour.roi.mesh.polygon.ROI3DTriangularMesh;
import plugins.adufour.vars.lang.Var;
import plugins.adufour.vars.lang.VarDouble;
import plugins.kernel.roi.roi3d.ROI3DArea;

public class Mesh3D extends ActiveContour
{
    /**
     * an active vertex is a vertex that carries motion information
     * 
     * @author Alexandre Dufour
     */
    private static class ActiveVertex extends Vertex3D
    {
        public final Vector3d imageForces      = new Vector3d();
        
        public final Vector3d internalForces   = new Vector3d();
        
        public final Vector3d feedbackForces   = new Vector3d();
        
        public final Vector3d volumeConstraint = new Vector3d();
        
        public ActiveVertex(ActiveVertex v)
        {
            super(v.position, v.neighbors);
        }
        
        public ActiveVertex(Point3d position)
        {
            this(position, 0);
        }
        
        public ActiveVertex(Point3d position, int nbNeighbors)
        {
            super(position, nbNeighbors);
        }
        
        @Override
        public Vertex3D clone()
        {
            return new ActiveVertex(this);
        }
    }
    
    /**
     * An ActiveMesh is defined here as a 3D surface mesh with {@link Triangle} as elementary face
     * and {@link ActiveVertex} as elementary vertex
     * 
     * @author Alexandre Dufour
     */
    public static class ActiveMesh extends ROI3DTriangularMesh
    {
        public ActiveMesh()
        {
            
        }
        
        public ActiveMesh(double sampling, ROI3D roi, Tuple3d pixelSize)
        {
            super(sampling, roi, pixelSize);
        }
        
        @Override
        public ActiveVertex createVertex(Point3d position)
        {
            return new ActiveVertex(position);
        }
    }
    
    final ActiveMesh mesh;
    
    /**
     * DO NOT USE! This constructor is for XML loading purposes only
     */
    public Mesh3D()
    {
        mesh = new ActiveMesh();
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    public Mesh3D(Mesh3D contour)
    {
        super(contour.sampling, new SlidingWindow(contour.convergence.getSize()));
        
        mesh = (ActiveMesh) contour.mesh.clone();
        setColor(contour.getColor());
        mesh.setColor(getColor());
    }
    
    public Mesh3D(Var<Double> sampling, Tuple3d pixelSize, ROI3D roi, SlidingWindow convergenceWindow)
    {
        super(sampling, convergenceWindow);
        
        mesh = new ActiveMesh(sampling.getValue(), roi, pixelSize);
        mesh.setColor(getColor());
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
        Vector3d axis = getMajorAxis();
        
        // To drive the contour along the main object axis, each displacement
        // vector is scaled by the scalar product between its normal and the main axis.
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            // dot product between normalized vectors ranges from -1 to 1
            double colinearity = Math.abs(v.normal.dot(axis)); // now from 0 to 1
            
            // goal: adjust the minimum using the weight, but keep max to 1
            double threshold = Math.max(colinearity, 1 - weight);
            
            ((ActiveVertex) v).imageForces.scale(threshold);
        }
    }
    
    @Override
    void computeBalloonForces(double weight)
    {
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            ((ActiveVertex) v).imageForces.scaleAdd(weight, v.normal, ((ActiveVertex) v).imageForces);
        }
    }
    
    /**
     * Update edge term of the contour evolution according to the image gradient
     * 
     * @param weight
     * @param edgeData
     */
    @Override
    void computeEdgeForces(Sequence edgeData, int channel, double weight)
    {
        Vector3d grad = new Vector3d();
        Point3d prev = new Point3d();
        
        Point3d p = new Point3d();
        double pixelSizeX = edgeData.getPixelSizeX();
        double pixelSizeY = edgeData.getPixelSizeY();
        double pixelSizeZ = edgeData.getPixelSizeZ();
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            // convert from metric to image space
            p.set(v.position.x / pixelSizeX, v.position.y / pixelSizeY, v.position.z / pixelSizeZ);
            
            // compute the gradient (2nd order)
            
            grad.x = getPixelValue(edgeData, p.x + 0.5, p.y, p.z);
            grad.y = getPixelValue(edgeData, p.x, p.y + 0.5, p.z);
            grad.z = getPixelValue(edgeData, p.x, p.y, p.z + 0.5);
            
            prev.x = getPixelValue(edgeData, p.x - 0.5, p.y, p.z);
            prev.y = getPixelValue(edgeData, p.x, p.y - 0.5, p.z);
            prev.z = getPixelValue(edgeData, p.x, p.y, p.z - 0.5);
            
            grad.sub(prev);
            grad.scale(weight);
            ((ActiveVertex) v).imageForces.add(grad);
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
        
        Vector3d regionForce = new Vector3d();
        
        double pixelSizeX = imageData.getPixelSizeX();
        double pixelSizeY = imageData.getPixelSizeY();
        double pixelSizeZ = imageData.getPixelSizeZ();
        
        double val, inDiff, outDiff;
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            Point3d p = v.position;
            
            regionForce.set(v.normal);
            
            val = getPixelValue(imageData, p.x / pixelSizeX, p.y / pixelSizeY, p.z / pixelSizeZ);
            
            inDiff = val - cin;
            inDiff *= inDiff;
            
            outDiff = val - cout;
            outDiff *= outDiff;
            
            regionForce.scale(weight * sampling.getValue() * (sensitivity * outDiff) - (inDiff / sensitivity));
            
            ((ActiveVertex) v).imageForces.add(regionForce);
        }
        
    }
    
    @Override
    void computeInternalForces(double weight)
    {
        Vector3d internalForce = new Vector3d();
        
        weight /= sampling.getValue();
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            internalForce.scale(-v.neighbors.size(), v.position);
            
            for (Integer nn : v.neighbors)
                internalForce.add(mesh.getVertex(nn).position);
            
            internalForce.scale(weight);
            
            ((ActiveVertex) v).internalForces.add(internalForce);
        }
    }
    
    void computeVolumeConstraint(double targetVolume)
    {
        // 1) compute the difference between target and current volume
        double volumeDiff = targetVolume - getDimension(2);
        // if (volumeDiff > 0): contour too small, should no longer shrink
        // if (volumeDiff < 0): contour too big, should no longer grow
        
        Vector3d avgFeedback = new Vector3d();
        int cpt = 0;
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            Vector3d feedbackForce = ((ActiveVertex) v).feedbackForces;
            
            // 2) check whether the final force has same direction as the outer normal
            double forceNorm = ((ActiveVertex) v).imageForces.dot(v.normal);
            
            // if forces have same direction (forceNorm > 0): contour is growing
            // if forces have opposite direction (forceNorm < 0): contour is shrinking
            
            // if (forceNorm * volumeDiff < 0)
            // {
            // // forceNorm and volumeDiff have opposite signs because:
            // // - contour too small (volumeDiff > 0) and shrinking (forceNorm < 0)
            // // or
            // // - contour too large (volumeDiff < 0) and growing (forceNorm > 0)
            // // => in both cases, constrain the final force accordingly
            // ((ActiveVertex) v).drivingForces.scale(1.0 / (1.0 + Math.abs(volumeDiff)));
            // }
            
            // estimate an average feedback
            if (forceNorm > 0 && volumeDiff < 0)
            {
                avgFeedback.add(feedbackForce);
                cpt++;
            }
        }
        
        if (avgFeedback.length() > 0)
        {
            avgFeedback.scale(1.0 / cpt);
            avgFeedback.scale(Math.abs(volumeDiff / targetVolume) / 1.5);
            
            // move the entire mesh (ugly, but amazingly efficient!!)
            for (Vertex3D v : mesh.getVertices())
            {
                if (v != null) ((ActiveVertex) v).volumeConstraint.add(avgFeedback);
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
        double targetRadius = target.boundingSphere.getRadius();
        
        Point3d lowerBounds = new Point3d();
        Point3d upperBounds = new Point3d();
        target.boundingBox.getLower(lowerBounds);
        target.boundingBox.getUpper(upperBounds);
        
        double feedback = 0;
        Vector3d feedbackForce = new Vector3d();
        
        int tests = 0;
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            if (v.position.x < lowerBounds.x || v.position.y < lowerBounds.y || v.position.z < lowerBounds.z) continue;
            
            if (v.position.x > upperBounds.x || v.position.y > upperBounds.y || v.position.z > upperBounds.z) continue;
            
            double distance = v.position.distance(targetCenter);
            
            if (distance < targetRadius)
            {
                tests++;
                
                if ((feedback = target.getDistanceToEdge(v.position)) > 0)
                {
                    feedbackForce.set(v.normal);
                    feedbackForce.scale(-feedback * 10);
                    ((ActiveVertex) v).feedbackForces.add(feedbackForce);
                }
            }
        }
        
        return tests;
    }
    
    public double getCurvature(Point3d pt)
    {
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            if (!v.position.equals(pt)) continue;
            
            Vector3d sum = new Vector3d();
            Vector3d diff = new Vector3d();
            for (Integer n : v.neighbors)
            {
                Point3d neighbor = mesh.getVertex(n).position;
                diff.sub(neighbor, pt);
                sum.add(diff);
            }
            sum.scale(1.0 / v.neighbors.size());
            
            return sum.length() * Math.signum(sum.dot(v.normal));
        }
        
        return 0;
    }
    
    public double getDimension(int order)
    {
        switch (order)
        {
        case 0:
            return mesh.getNumberOfVertices(true);
        case 1:
            return mesh.getNumberOfContourPoints();
        case 2:
            return mesh.getNumberOfPoints();
        }
        return Double.NaN;
    }
    
    Point3d getMassCenter(boolean convertToImageSpace)
    {
        return mesh.getMassCenter(convertToImageSpace);
    }
    
    /**
     * @return The major axis of this contour, i.e. an unnormalized vector formed of the two most
     *         distant contour points
     */
    public Vector3d getMajorAxis()
    {
        return mesh.getMajorAxis();
    }
    
    /**
     * Calculates the 3D image value at the given real coordinates by tri-linear interpolation
     * 
     * @param imageFloat
     *            the image to sample (must be of type {@link DataType#FLOAT})
     * @param x
     *            the X-coordinate of the point
     * @param y
     *            the Y-coordinate of the point
     * @param z
     *            the Z-coordinate of the point
     * @return the interpolated image value at the given coordinates
     */
    private float getPixelValue(Sequence data, double x, double y, double z)
    {
        // "center" the coordinates to the center of the pixel
        x -= 0.5;
        y -= 0.5;
        z -= 0.5;
        
        int width = data.getSizeX();
        int height = data.getSizeY();
        int depth = data.getSizeZ();
        
        final int i = (int) Math.floor(x);
        final int j = (int) Math.floor(y);
        final int k = (int) Math.floor(z);
        
        if (i < 0 || i >= width - 1) return 0;
        if (j < 0 || j >= height - 1) return 0;
        if (k < 0 || k >= depth - 1) return 0;
        
        final int pixel = i + j * width;
        final int east = pixel + 1; // saves 3 additions
        final int south = pixel + width; // saves 1 addition
        final int southeast = south + 1; // saves 1 addition
        
        float[] currSlice = data.getDataXYAsFloat(0, k, 0);
        float[] nextSlice = data.getDataXYAsFloat(0, k + 1, 0);
        
        float value = 0;
        
        x -= i;
        y -= j;
        z -= k;
        
        final double mx = 1 - x;
        final double my = 1 - y;
        final double mz = 1 - z;
        
        value += mx * my * mz * currSlice[pixel];
        value += x * my * mz * currSlice[east];
        value += mx * y * mz * currSlice[south];
        value += x * y * mz * currSlice[southeast];
        value += mx * my * z * nextSlice[pixel];
        value += x * my * z * nextSlice[east];
        value += mx * y * z * nextSlice[south];
        value += x * y * z * nextSlice[southeast];
        
        return value;
    }
    
    @Override
    public double getX()
    {
        // get this in pixel units
        return mesh.getMassCenter(true).x;
    }
    
    @Override
    public double getY()
    {
        // get this in pixel units
        return mesh.getMassCenter(true).y;
    }
    
    @Override
    public double getZ()
    {
        // get this in pixel units
        return mesh.getMassCenter(true).z;
    }
    
    @Override
    public void setT(int t)
    {
        super.setT(t);
        mesh.setT(t);
    }
    
    @Override
    public Iterator<Point3d> iterator()
    {
        // return a "tweaked" iterator that will skip null entries automatically
        
        return new Iterator<Point3d>()
        {
            Iterator<Vertex3D> vertexIterator       = mesh.getVertices().iterator();
            
            Vertex3D           next;
            
            boolean            hasNext;
            
            boolean            hasNextWasCalledOnce = false;
            
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
                }
                while (next == null && vertexIterator.hasNext());
                
                return next != null;
            }
        };
    }
    
    void move(ROI field, double timeStep)
    {
        Vector3d force = new Vector3d();
        double maxDisp = sampling.getValue() * timeStep;
        
        Point3d p = new Point3d();
        Tuple3d pixelSize = mesh.getPixelSize();
        
        for (Vertex3D v : mesh.getVertices())
        {
            if (v == null) continue;
            
            ActiveVertex av = (ActiveVertex) v;
            
            // get the vertex location in image space
            p.set(v.position.x / pixelSize.x, v.position.y / pixelSize.y, v.position.z / pixelSize.z);
            
            // apply model forces if p lies within the area of interest
            if (field != null && p.z >= 0 && field.contains(p.x, p.y, p.z, 0, 0))
            {
                if (av.volumeConstraint.length() > 0) av.position.add(av.volumeConstraint);
                
                force.set(av.imageForces);
                
                force.add(av.internalForces);
                
                force.add(av.feedbackForces);
            }
            else
            {
                // force.set(av.imageForces);
                force.set(av.internalForces);
                force.scale(0.1);
            }
            
            force.scale(timeStep);
            
            double disp = force.length();
            // forces cannot be larger than the max authorized displacement (stability condition)
            if (disp > maxDisp) force.scale(maxDisp / disp);
            
            // move the vertex
            v.position.add(force);
            
            // reset forces
            av.imageForces.set(0, 0, 0);
            av.internalForces.set(0, 0, 0);
            av.feedbackForces.set(0, 0, 0);
            av.volumeConstraint.set(0, 0, 0);
        }
        
        updateMetaData();
        
        // compute some convergence criterion
        
        if (convergence == null) return;
        
        convergence.push(mesh.getNumberOfPoints());
    }
    
    @Override
    public boolean hasConverged(Operation operation, double epsilon)
    {
        Double value = convergence.computeCriterion(operation);
        return value != null && value <= epsilon / 10;
    }
    
    @Override
    protected void updateMetaData()
    {
        super.updateMetaData();
        
        mesh.roiChanged();
    }
    
    @Override
    public boolean saveToXML(Node node)
    {
        if (!super.saveToXML(node)) return false;
        
        return mesh.saveToXML(node);
    }
    
    @Override
    public boolean loadFromXML(Node node)
    {
        if (!super.loadFromXML(node)) return false;
        
        boolean success = mesh.loadFromXML(node);
        
        mesh.setT(getT());
        
        return success;
    }
    
    @Override
    public void reSample(double minFactor, double maxFactor) throws TopologyException
    {
        try
        {
            mesh.reSampleToAverageDistance(sampling.getValue(), 0.4);
        }
        catch (MeshTopologyException e)
        {
            if (e.children == null) throw new TopologyException(this, null);
            
            Mesh3D[] children = new Mesh3D[e.children.length];
            
            for (int i = 0; i < children.length; i++)
            {
                children[i] = new Mesh3D(sampling, mesh.getPixelSize(), e.children[i], convergence);
            }
            
            throw new TopologyException(this, children);
        }
    }
    
    @Override
    public ActiveContour clone()
    {
        return new Mesh3D(this);
    }
    
    @Override
    protected void addPoint(Point3d p)
    {
        mesh.addVertex(mesh.createVertex(p));
    }
    
    @Override
    protected ActiveContour[] checkSelfIntersection(double minDistance)
    {
        // TODO
        return null;
    }
    
    @Override
    public double computeAverageIntensity(Sequence imageData, int channel, BooleanMask3D mask) throws TopologyException
    {
        VarDouble avg = new VarDouble("avg", 0.0);
        mesh.rasterScan(imageData, avg, mask);
        return avg.getValue();
    }
    
    @Override
    public double getDistanceToEdge(Point3d p)
    {
        return mesh.getDistanceToMesh(p);
    }
    
    @Override
    public ROI toROI()
    {
        return toROI(ROIType.POLYGON, null);
    }
    
    @Override
    public ROI toROI(ROIType type, Sequence sequence) throws UnsupportedOperationException
    {
        ROI3D roi;
        
        switch (type)
        {
        case POLYGON: {
            roi = mesh.clone();
            break;
        }
        case AREA: {
            roi = new ROI3DArea(mesh.getBooleanMask(true));
            break;
        }
        default:
            throw new UnsupportedOperationException("Cannot export a ROI of type: " + type);
        }
        
        roi.setT(mesh.getT());
        roi.setColor(getColor());
        return roi;
    }
    
    @Override
    public void toSequence(Sequence output, double value)
    {
        // TODO
    }
    
    @Override
    protected void updateNormals()
    {
        mesh.updateNormals();
    }
    
    @Override
    protected void clean()
    {
        for (Overlay overlay : overlays.values())
            overlay.remove();
        overlays.clear();
    }
    
    HashMap<IcyCanvas, Overlay> overlays = new HashMap<IcyCanvas, Overlay>();
    
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        if (overlays.containsKey(canvas)) return;
        
        Overlay overlay = mesh.getOverlay();
        overlay.setCanBeRemoved(true);
        overlay.setName("Active Mesh overlay");
        overlays.put(canvas, overlay);
        sequence.addOverlay(overlay);
    }
}
