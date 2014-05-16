package plugins.adufour.activecontours;

import icy.canvas.Canvas3D;
import icy.canvas.IcyCanvas;
import icy.roi.ROI;
import icy.roi.ROI3D;
import icy.sequence.Sequence;
import icy.system.SystemUtil;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import icy.type.point.Point3D;
import icy.type.rectangle.Rectangle3D;
import icy.vtk.VtkUtil;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.media.j3d.BoundingBox;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import plugins.adufour.activecontours.ActiveContours.ROIType;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.quickhull.QuickHull3D;
import vtk.vtkActor;
import vtk.vtkDoubleArray;
import vtk.vtkPoints;
import vtk.vtkPolyData;
import vtk.vtkPolyDataMapper;
import vtk.vtkRenderer;

public class Mesh3D extends ActiveContour
{
    static final ExecutorService threadPool = Executors.newFixedThreadPool(SystemUtil.getAvailableProcessors());
    
    /**
     * Structural element of triangular mesh
     * 
     * @author Alexandre Dufour
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
    
    final ArrayList<Face> faces = new ArrayList<Face>();
    
    /**
     * Structural element of triangular mesh
     * 
     * @author Alexandre Dufour
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
    
    final ArrayList<Vertex> vertices = new ArrayList<Vertex>();
    
    protected Mesh3D(ActiveContours owner, EzVarDouble contour_resolution, SlidingWindow convergenceWindow)
    {
        super(owner, contour_resolution, convergenceWindow);
        
        setColor(Color.getHSBColor((float) Math.random(), 0.8f, 0.9f));
    }
    
    /**
     * Creates a clone of the specified contour
     * 
     * @param contour
     */
    public Mesh3D(Mesh3D contour)
    {
        this(contour.owner, contour.resolution, new SlidingWindow(contour.convergence.getSize()));
        
        setX(contour.x);
        setY(contour.y);
        setZ(contour.z);
        
        setColor(contour.getColor());
        
        vertices.ensureCapacity(contour.vertices.size());
        for (Vertex v : contour.vertices)
            vertices.add(v == null ? null : new Vertex(v));
        
        faces.ensureCapacity(contour.faces.size());
        for (Face f : contour.faces)
            faces.add(new Face(f.v1, f.v2, f.v3));
    }
    
    public Mesh3D(ActiveContours owner, double xyRes, double zRes, EzVarDouble contour_resolution, SlidingWindow convergenceWindow, ROI3D roi)
    {
        this(owner, contour_resolution, convergenceWindow);
        
        // just get the convex envelope, it's good enough
        
        // Point3D.Integer[] points = roi.getBooleanMask(true).getContourPoints();
        // Point3d[] doublePts = new Point3d[points.length];
        // for (int i = 0; i < points.length; i++)
        // doublePts[i] = new Point3d(points[i].x * xyRes, points[i].y * xyRes, points[i].z * zRes);
        //
        // QuickHull3D q3 = new QuickHull3D(doublePts);
        //
        // q3.triangulate();
        
        try
        {
            double PHI = (1.0 + Math.sqrt(5.0)) / 2.0;
            
            // Declaration of the vertices
            int zA = addPoint(new Point3d(PHI, 1, 0));
            int zB = addPoint(new Point3d(-PHI, 1, 0));
            int zC = addPoint(new Point3d(-PHI, -1, 0));
            int zD = addPoint(new Point3d(PHI, -1, 0));
            int yA = addPoint(new Point3d(1, 0, PHI));
            int yB = addPoint(new Point3d(1, 0, -PHI));
            int yC = addPoint(new Point3d(-1, 0, -PHI));
            int yD = addPoint(new Point3d(-1, 0, PHI));
            int xA = addPoint(new Point3d(0, PHI, 1));
            int xB = addPoint(new Point3d(0, -PHI, 1));
            int xC = addPoint(new Point3d(0, -PHI, -1));
            int xD = addPoint(new Point3d(0, PHI, -1));
            
            // Declaration of the faces
            addFace(yA, xA, yD);
            addFace(yA, yD, xB);
            addFace(yB, yC, xD);
            addFace(yB, xC, yC);
            
            addFace(zA, yA, zD);
            addFace(zA, zD, yB);
            addFace(zC, yD, zB);
            addFace(zC, zB, yC);
            
            addFace(xA, zA, xD);
            addFace(xA, xD, zB);
            addFace(xB, xC, zD);
            addFace(xB, zC, xC);
            
            addFace(xA, yA, zA);
            addFace(xD, zA, yB);
            addFace(yA, xB, zD);
            addFace(yB, zD, xC);
            addFace(yD, xA, zB);
            addFace(yC, zB, xD);
            addFace(yD, zC, xB);
            addFace(yC, xC, zC);
            
            Rectangle3D r3 = roi.computeBounds3D();
            for(Point3d pt : this)
            {
                pt.x += r3.getCenterX() * xyRes;
                pt.y += r3.getCenterY() * xyRes;
                pt.z += r3.getCenterZ() * zRes;
            }
            
//            reSample(0.6, 1.4);
            
            updateMetaData();
        }
        catch (TopologyException e)
        {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
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
    protected int addPoint(Point3d point)
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
            
            if (v.position.distanceSquared(point) < resolution.getValue() * 0.00001) return -index - 1;
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
        
        return addPoint(newPosition);
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
    protected Mesh3D[] checkSelfIntersection(double minSpacing)
    {
        // TODO
        return null;
    }
    
    @Override
    public double computeAverageIntensity(Sequence imageData_float, int channel, Sequence buffer_data)
    {
        double inSum = 0, inCpt = 0;
        
        // direction of scan: along X
        final Vector3d direction = new Vector3d(1, 0, 0);
        
        final int scanLine = imageData_float.getSizeX();
        final int height = imageData_float.getSizeY();
        final int depth = imageData_float.getSizeZ();
        
        final double resX = imageData_float.getPixelSizeX();
        final double resY = imageData_float.getPixelSizeY();
        final double resZ = imageData_float.getPixelSizeZ();
        
        Point3d boxMin = new Point3d(), boxMax = new Point3d();
        BoundingBox bbox = getBoundingBox();
        bbox.getLower(boxMin);
        bbox.getUpper(boxMax);
        
        final int minX = Math.max(0, (int) Math.floor(boxMin.x / resX) - 1);
        final int minY = Math.max(0, (int) Math.floor(boxMin.y / resY) - 1);
        final int minZ = Math.max(0, (int) Math.floor(boxMin.z / resZ) - 1);
        
        final int maxY = Math.min(height - 1, (int) Math.ceil(boxMax.y / resY) + 1);
        final int maxZ = Math.min(depth - 1, (int) Math.ceil(boxMax.z / resZ) + 1);
        
        // epsilon used for the line-triangle intersection test
        final double epsilon = 1.0e-13;
        
        ArrayList<Future<double[]>> sliceTasks = new ArrayList<Future<double[]>>(maxZ - minZ + 1);
        
        for (int k = minZ; k <= maxZ; k++)
        {
            final short[] maskSlice = buffer_data.getDataXYAsShort(0, k, 0);
            final float[] dataSlice = imageData_float.getDataXYAsFloat(0, k, 0);
            
            final int slice = k;
            
            Callable<double[]> sliceTask = new Callable<double[]>()
            {
                @Override
                public double[] call() throws Exception
                {
                    Point3d origin = new Point3d(minX * resX, minY * resY, slice * resZ);
                    
                    Vector3d edge1 = new Vector3d();
                    Vector3d edge2 = new Vector3d();
                    Vector3d vp = new Vector3d();
                    Vector3d vt = new Vector3d();
                    Vector3d vq = new Vector3d();
                    
                    ArrayList<Integer> intersections = new ArrayList<Integer>(4);
                    
                    double localInSum = 0;
                    int localInCpt = 0;
                    
                    for (int j = minY; j < maxY; j++, origin.y += resY)
                    {
                        int lineOffset = j * scanLine;
                        
                        intersections.clear();
                        int nbCrosses = 0;
                        
                        for (Face f : faces)
                        {
                            Point3d v1 = vertices.get(f.v1).position;
                            Point3d v2 = vertices.get(f.v2).position;
                            Point3d v3 = vertices.get(f.v3).position;
                            
                            // raw intersection test with the face bounds
                            
                            if (origin.y < v1.y && origin.y < v2.y && origin.y < v3.y) continue;
                            if (origin.z < v1.z && origin.z < v2.z && origin.z < v3.z) continue;
                            if (origin.y > v1.y && origin.y > v2.y && origin.y > v3.y) continue;
                            if (origin.z > v1.z && origin.z > v2.z && origin.z > v3.z) continue;
                            
                            // more precise ray-triangle intersection test
                            
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
                            
                            // intersection found, compute the distance
                            
                            double distance = edge2.dot(vq) * inv_det;
                            Integer distanceInVoxels = minX + (int) Math.round(distance / resX);
                            
                            if (distanceInVoxels < 0) continue;
                            
                            if (!intersections.contains(distanceInVoxels))
                            {
                                intersections.add(distanceInVoxels);
                                nbCrosses++;
                            }
                            else
                            {
                                // if a distance appears twice, then the mesh "slices" the same
                                // voxel twice (rounding off to the same integer distance)
                                // => it's a double crossing
                                // ==> remove the first occurrence (forget about it)
                                intersections.remove(distanceInVoxels);
                                nbCrosses--;
                            }
                        }
                        
                        // if nbCrosses = 0, the ray never crosses the mesh
                        // if nbCrosses = 1, the ray is tangent to the mesh
                        if (nbCrosses < 2) continue;
                        
                        Collections.sort(intersections);
                        
                        if (nbCrosses == 3)
                        {
                            intersections.remove(1);
                            nbCrosses--;
                        }
                        
                        for (int cross = 0; cross < nbCrosses; cross += 2)
                        {
                            int start = intersections.get(cross);
                            int stop = intersections.get(cross + 1);
                            
                            for (int i = start; i < stop; i++)
                            {
                                localInCpt++;
                                int offset = i + lineOffset;
                                localInSum += dataSlice[offset];
                                if (maskSlice != null) maskSlice[offset] = 1;
                            }
                        }
                    }
                    
                    return (localInCpt == 0) ? new double[2] : new double[] { localInSum, localInCpt };
                }
            };
            
            sliceTasks.add(threadPool.submit(sliceTask));
        }
        
        try
        {
            for (Future<double[]> sliceTask : sliceTasks)
            {
                double[] sliceStats = sliceTask.get();
                inSum += sliceStats[0];
                inCpt += sliceStats[1];
            }
        }
        catch (InterruptedException e)
        {
            // preserve the interrupted state
            Thread.currentThread().interrupt();
        }
        catch (ExecutionException e)
        {
            e.printStackTrace();
        }
        
        if (inCpt == 0)
        {
            // no voxel inside the mesh => mesh is becoming extremely flat
            System.err.println("Warning: flat mesh detected (probably on the volume edge)");
        }
        
        return inSum / inCpt;
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
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            // dot product between normalized vectors ranges from -1 to 1
            double colinearity = Math.abs(v.normal.dot(axis)); // now from 0 to 1
            
            // goal: adjust the minimum using the weight, but keep max to 1
            double threshold = Math.max(colinearity, 1 - weight);
            
            v.drivingForces.scale(threshold);
        }
    }
    
    @Override
    void computeBalloonForces(double weight)
    {
        for (Vertex v : vertices)
        {
            v.drivingForces.scaleAdd(weight, v.normal, v.drivingForces);
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
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            // compute the gradient (2nd order)
            
            grad.x = getPixelValue(edgeData, v.position.x + 0.5, v.position.y, v.position.z);
            grad.y = getPixelValue(edgeData, v.position.x, v.position.y + 0.5, v.position.z);
            grad.z = getPixelValue(edgeData, v.position.x, v.position.y, v.position.z + 0.5);
            
            prev.x = getPixelValue(edgeData, v.position.x - 0.5, v.position.y, v.position.z);
            prev.y = getPixelValue(edgeData, v.position.x, v.position.y - 0.5, v.position.z);
            prev.z = getPixelValue(edgeData, v.position.x, v.position.y, v.position.z - 0.5);
            
            grad.sub(prev);
            grad.scale(weight);
            v.drivingForces.add(grad);
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
        
        double val, inDiff, outDiff, sum;
        
        int width = imageData.getSizeX();
        int height = imageData.getSizeY();
        int depth = imageData.getSizeZ();
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            Point3d p = v.position;
            
            // bounds check
            if (p.x < 1 || p.x >= width - 1) continue;
            if (p.y < 1 || p.y >= height - 1) continue;
            if (p.z < 1 || p.z >= depth - 1) continue;
            
            val = getPixelValue(imageData, p.x, p.y, p.z);
            
            inDiff = val - cin;
            inDiff *= inDiff;
            
            outDiff = val - cout;
            outDiff *= outDiff;
            
            sum = weight * resolution.getValue() * (sensitivity * outDiff) - (inDiff / sensitivity);
            
            v.drivingForces.scaleAdd(sum, v.normal, v.drivingForces);
        }
        
    }
    
    @Override
    void computeInternalForces(double weight)
    {
        Vector3d internalForce = new Vector3d();
        
        weight /= resolution.getValue();
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            internalForce.scale(-v.neighbors.size(), v.position);
            
            for (Integer nn : v.neighbors)
                internalForce.add(vertices.get(nn).position);
            
            internalForce.scale(weight);
            
            v.drivingForces.add(internalForce);
        }
    }
    
    void computeVolumeConstraint(double targetVolume)
    {
        // 1) compute the difference between target and current volume
        double volumeDiff = targetVolume - getDimension(2);
        // if (volumeDiff > 0): contour too small, should no longer shrink
        // if (volumeDiff < 0): contour too big, should no longer grow
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            // 2) check whether the final force has same direction as the outer normal
            double forceNorm = v.drivingForces.dot(v.normal);
            
            // if forces have same direction (forceNorm > 0): contour is growing
            // if forces have opposite direction (forceNorm < 0): contour is shrinking
            
            if (forceNorm * volumeDiff < 0)
            {
                // forceNorm and volumeDiff have opposite signs because:
                // - contour too small (volumeDiff > 0) and shrinking (forceNorm < 0)
                // or
                // - contour too large (volumeDiff < 0) and growing (forceNorm > 0)
                // => in both cases, constrain the final force accordingly
                v.drivingForces.scale(1.0 / (1.0 + Math.abs(volumeDiff)));
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
        
        double penetration = 0;
        
        int tests = 0;
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            double distance = v.position.distance(targetCenter);
            
            if (distance < targetRadius)
            {
                tests++;
                
                if ((penetration = target.contains(v.position)) > 0)
                {
                    v.feedbackForces.scaleAdd(-penetration, v.normal, v.feedbackForces);
                }
            }
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
     * @return The major axis of this contour, i.e. an unnormalized vector formed of the two most
     *         distant contour points
     */
    public Vector3d getMajorAxis()
    {
        Vector3d axis = new Vector3d();
        int nbPoints = (int) getDimension(0);
        
        // TODO this is not optimal, geometric moments should be used
        
        double maxDistSq = 0;
        Vector3d vec = new Vector3d();
        
        for (int i = 0; i < nbPoints; i++)
        {
            Vertex v1 = vertices.get(i);
            if (v1 == null) continue;
            
            for (int j = i + 1; j < nbPoints; j++)
            {
                Vertex v2 = vertices.get(j);
                if (v2 == null) continue;
                
                vec.sub(v1.position, v2.position);
                
                double dSq = vec.lengthSquared();
                
                if (dSq > maxDistSq)
                {
                    maxDistSq = dSq;
                    axis.set(vec);
                }
            }
        }
        
        return axis;
    }
    
    /**
     * Calculates the 3D image value at the given real coordinates by tri-linear interpolation
     * 
     * @param imageFloat
     *            the image to sample (must be of type {@link DataType#DOUBLE})
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
        int width = data.getSizeX();
        int height = data.getSizeY();
        int depth = data.getSizeZ();
        
        final int i = (int) Math.round(x);
        final int j = (int) Math.round(y);
        final int k = (int) Math.round(z);
        
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
                }
                while (next == null && vertexIterator.hasNext());
                
                return next != null;
            }
        };
    }
    
    void move(ROI field, double timeStep)
    {
        Vector3d force = new Vector3d();
        double maxDisp = resolution.getValue() * timeStep;
        
        Point3d center = new Point3d();
        
        for (Vertex v : vertices)
        {
            if (v == null) continue;
            
            // apply model forces if p lies within the area of interest
            if (field != null && field.contains(v.position.x, v.position.y, 0, 0, 0)) force.set(v.drivingForces);
            
            // apply feedback forces all the time
            force.add(v.feedbackForces);
            
            force.scale(timeStep);
            
            double disp = force.length();
            
            if (disp > maxDisp) force.scale(maxDisp / disp);
            
            v.position.add(force);
            
            center.add(v.position);
            
            // reset forces
            v.drivingForces.set(0, 0, 0);
            v.feedbackForces.set(0, 0, 0);
        }
        
        updateMetaData();
        
        // compute some convergence criterion
        
        if (convergence == null) return;
        
        convergence.push(getDimension(2));
        
    }
    
    VTKMesh mesh = null;
    
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        // only paint detections on the current frame
        if (getT() != canvas.getPositionT()) return;
        
        if (g != null)
        {
            // 2D viewer
            
            float fontSize = (float) canvas.canvasToImageLogDeltaX(30);
            g.setFont(new Font("Trebuchet MS", Font.BOLD, 10).deriveFont(fontSize));
            
            double stroke = Math.max(canvas.canvasToImageLogDeltaX(3), canvas.canvasToImageLogDeltaY(3));
            
            g.setStroke(new BasicStroke((float) stroke, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
            
            g.setColor(getColor());
            
            g.drawString("[3D mesh]", (float) x, (float) y);
            
            // TODO draw something more in 2D (points? raster mesh?)
        }
        else if (canvas instanceof Canvas3D)
        {
            // 3D viewer
            
            if (mesh == null)
            {
                mesh = new VTKMesh();
                mesh.update();
                ((Canvas3D) canvas).getRenderer().AddActor(mesh.actor);
            }
        }
        else
        {
            // other viewers
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
        double minLength = resolution.getValue() * minFactor;
        double maxLength = resolution.getValue() * maxFactor;
        
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
            // if (cpt > vertices.size()) noChange = true;
        }
        
        updateMetaData();
    }
    
    @Override
    protected void updateMetaData()
    {
        super.updateMetaData();
        
        if (mesh != null) mesh.update();
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
            ArrayList<Vertex> newVertices = new ArrayList<Vertex>(vertices.size());
            for (int i = 0; i < vertices.size(); i++)
                newVertices.add(null);
            
            extractVertices(seed, visitedIndexes, vertices, newVertices);
            extractFaces(newVertices, faces, newFaces);
            
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
                    if (n.compareTo(v1) != 0 && n.compareTo(v2) != 0 && n.compareTo(v3) != 0 && newVertices.get(n) == null)
                    {
                        newV.neighbors.remove(n);
                        i--;
                    }
                }
                newVertices.set(v, newV);
            }
            
            for (Face f : newFaces)
                if (f.contains(v1) && f.contains(v2))
                {
                    // if the edge v1-v2 appears counter-clockwisely in f,
                    // the new face must be clockwise and vice-versa
                    newFaces.add(f.isCounterClockwise(v1, v2) ? new Face(v1, v3, v2) : new Face(v1, v2, v3));
                    break;
                }
            
            Mesh3D newMesh = new Mesh3D(owner, resolution, new SlidingWindow(convergence.getSize()));
            newMesh.setMeshData(newVertices, newFaces);
            newMesh.setT(getT());
            
            if (newMesh.getDimension(0) > 4) children[child] = newMesh;
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
    
    private void setMeshData(ArrayList<Vertex> newVertices, ArrayList<Face> newFaces)
    {
        synchronized (vertices)
        {
            this.vertices.clear();
            this.vertices.addAll(newVertices);
            this.faces.clear();
            this.faces.addAll(newFaces);
        }
        updateMetaData();
    }
    
    protected void updateNormals()
    {
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
        return toROI(ROIType.AREA);
    }
    
    @Override
    public ROI3D toROI(ROIType type)
    {
        // TODO Auto-generated method stub
        return null;
    }
    
    private class VTKMesh
    {
        private final vtkDoubleArray    vtkVerticesCoords = new vtkDoubleArray();
        private final vtkPoints         vtkVerticesInfo   = new vtkPoints();
        
        public final vtkPolyData        polyData          = new vtkPolyData();
        private final vtkPolyDataMapper polyDataMapper    = new vtkPolyDataMapper();
        public final vtkActor           actor             = new vtkActor();
        private vtkRenderer             renderer;
        
        VTKMesh()
        {
            this.vtkVerticesCoords.SetNumberOfComponents(3);
            this.vtkVerticesInfo.SetData(this.vtkVerticesCoords);
            this.polyData.SetPoints(this.vtkVerticesInfo);
            
            this.polyDataMapper.SetInput(this.polyData);
            
            this.actor.SetMapper(this.polyDataMapper);
        }
        
        public void update()
        {
            synchronized (vertices)
            {
                int nFaces = faces.size();
                
                if (nFaces == 0) return;
                
                double[] vertices = new double[Mesh3D.this.vertices.size() * 3];
                
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
            }
        }
        
        public void clean()
        {
            if ((this.renderer != null) && (this.actor != null)) this.renderer.RemoveActor(this.actor);
        }
    }
    
    @Override
    public void toSequence(Sequence output, double value)
    {
        // direction of scan: along X
        final Vector3d direction = new Vector3d(1, 0, 0);
        
        final int scanLine = output.getSizeX();
        final int height = output.getSizeY();
        final int depth = output.getSizeZ();
        
        final double resX = 1;// output.getPixelSizeX();
        final double resY = 1;// output.getPixelSizeY();
        final double resZ = 1;// output.getPixelSizeZ();
        
        Point3d boxMin = new Point3d(), boxMax = new Point3d();
        BoundingBox bbox = getBoundingBox();
        bbox.getLower(boxMin);
        bbox.getUpper(boxMax);
        
        final int minX = Math.max(0, (int) Math.floor(boxMin.x / resX) - 1);
        final int minY = Math.max(0, (int) Math.floor(boxMin.y / resY) - 1);
        final int minZ = Math.max(0, (int) Math.floor(boxMin.z / resZ) - 1);
        
        final int maxX = Math.min(scanLine - 1, (int) Math.ceil(boxMax.x / resX) + 1);
        final int maxY = Math.min(height - 1, (int) Math.ceil(boxMax.y / resY) + 1);
        final int maxZ = Math.min(depth - 1, (int) Math.ceil(boxMax.z / resZ) + 1);
        
        // epsilon used for the line-triangle intersection test
        final double epsilon = 1.0e-12;
        
        ArrayList<Future<?>> sliceTasks = new ArrayList<Future<?>>(maxZ - minZ + 1);
        
        for (int k = minZ; k <= maxZ; k++)
        {
            final Object maskSlice = output.getDataXY(0, k, 0);
            
            final int slice = k;
            
            Runnable sliceTask = new Runnable()
            {
                @Override
                public void run()
                {
                    for (int j = minY; j < maxY; j++)
                        for (int i = minX; i < maxX; i++)
                        {
                            if (contains(new Point3d(i, j, slice)) > 0)
                            {
                                Array1DUtil.setValue(maskSlice, i + scanLine * j, 3000);
                            }
                        }
                }
            };
            
            sliceTasks.add(threadPool.submit(sliceTask));
        }
        
        try
        {
            for (Future<?> sliceTask : sliceTasks)
            {
                sliceTask.get();
            }
        }
        catch (InterruptedException e)
        {
            // preserve the interrupted state
            Thread.currentThread().interrupt();
        }
        catch (ExecutionException e)
        {
            e.getCause().printStackTrace();
            e.printStackTrace();
        }
        
    }
}
