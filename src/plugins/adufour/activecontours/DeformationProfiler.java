package plugins.adufour.activecontours;

import icy.file.FileUtil;
import icy.gui.dialog.SaveDialog;
import icy.gui.frame.progress.AnnounceFrame;
import icy.gui.util.GuiUtil;
import icy.math.ArrayMath;
import icy.plugin.interface_.PluginBundled;
import icy.preferences.XMLPreferences;
import icy.util.XLSUtil;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.vecmath.Point3d;
import javax.vecmath.SingularMatrixException;
import javax.vecmath.Vector3d;

import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import plugins.adufour.quickhull.QuickHull3D;
import plugins.adufour.roi3d.mesh.Vertex3D;
import plugins.adufour.roi3d.mesh.polygon.ROI3DTriangularMesh;
import plugins.fab.trackmanager.PluginTrackManagerProcessor;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * WARNING: this is an experimental class and may change in the future, use at your own risks!
 * 
 * @author Alexandre Dufour
 */
public class DeformationProfiler extends PluginTrackManagerProcessor implements PluginBundled
{
    private static XMLPreferences preferences = null;
    
    private enum Descriptors
    {
        None("None"),
        Norm1("Surface area"),
        Norm2("Volume"),
        Roundness("Roundness (min radius / max radius)"),
        Convexity("Solidity (area / convex area)"),
        Smoothness("Contour smoothness"),
        Sphericity("Sphericity"),
        Elongation("Elongation factor"),
        Flatness("Flatness");
        
        private Descriptors(String name)
        {
            this.name = name;
        }
        
        final String name;
        
        @Override
        public String toString()
        {
            return name;
        }
    }
    
    private JComboBox jComboDescriptors     = new JComboBox(Descriptors.values());
    
    private JButton   jButtonSaveToVTK      = new JButton("Export meshes to VTK files");
    
    private JButton   jButtonSaveShapeToXLS = new JButton("Export shape measures to XLS");
    
    private JPanel    chartPanel            = new JPanel();
    
    private double    tScale;
    
    public DeformationProfiler()
    {
        if (preferences == null) preferences = getPreferencesRoot();
        
        super.getDescriptor().setDescription("Monitor the 3D deformation over time");
        super.setName("3D Deformation Profiler");
        
        panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
        
        panel.add(Box.createVerticalStrut(5));
        
        panel.add(new JLabel("Warning: this plug-in is in development and may change at any time", JLabel.CENTER));
        
        panel.add(Box.createVerticalStrut(5));
        
        JPanel panelExport = new JPanel();
        panelExport.setLayout(new BoxLayout(panelExport, BoxLayout.X_AXIS));
        panelExport.add(jButtonSaveToVTK);
        panelExport.add(Box.createHorizontalStrut(5));
        panelExport.add(jButtonSaveShapeToXLS);
        panel.add(panelExport);
        
        panel.add(Box.createVerticalStrut(5));
        
        panel.add(GuiUtil.createLineBoxPanel(Box.createHorizontalStrut(10), new JLabel("Plot descriptor:"), Box.createHorizontalStrut(10), jComboDescriptors,
                Box.createHorizontalStrut(10)));
        
        jComboDescriptors.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                Compute();
            }
        });
        
        jButtonSaveToVTK.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                new Thread()
                {
                    public void run()
                    {
                        // Restore last used folder
                        String vtkFolder = preferences.get("vtkFolder", null);
                        
                        String path = SaveDialog.chooseFile("Export shape information", vtkFolder, "Mesh");
                        
                        if (path == null) return;
                        
                        // store the folder in the preferences
                        preferences.put("vtkFolder", FileUtil.getDirectory(path));
                        
                        AnnounceFrame message = new AnnounceFrame("Saving VTK files...", 0);
                        try
                        {
                            exportMeshToVTK(path);
                        }
                        finally
                        {
                            message.close();
                        }
                    }
                }.start();
            }
        });
        
        jButtonSaveShapeToXLS.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                new Thread()
                {
                    public void run()
                    {
                        // Restore last used folder
                        String xlsFolder = preferences.get("xlsFolder", null);
                        
                        String path = SaveDialog.chooseFile("Export shape information", xlsFolder, "Shape", ".xls");
                        
                        if (path == null) return;
                        
                        // store the folder in the preferences
                        preferences.put("xlsFolder", FileUtil.getDirectory(path));
                        
                        AnnounceFrame message = new AnnounceFrame("Saving shape information...", 0);
                        try
                        {
                            exportShapeToXLS(path);
                        }
                        finally
                        {
                            message.close();
                        }
                    }
                }.start();
            }
        });
        
        panel.add(chartPanel);
    }
    
    @Override
    public void Close()
    {
    }
    
    @Override
    public synchronized void Compute()
    {
        chartPanel.removeAll();
        
        if (!super.isEnabled()) return;
        
        if (trackPool.getDisplaySequence() == null) return;
        
        tScale = trackPool.getDisplaySequence().getTimeInterval();
        
        ChartPanel chart = null;
        
        Descriptors descriptor = (Descriptors) jComboDescriptors.getSelectedItem();
        
        switch (descriptor)
        {
        case None:
            return;
            
        case Norm1:
            double[][] dim1 = computeDimension(1);
            chart = createChartPanel(dim1, descriptor.toString(), "Time (sec.)", "\u03BC" + "m\u00B2");
            break;
        
        case Norm2:
            double[][] dim2 = computeDimension(2);
            chart = createChartPanel(dim2, descriptor.toString(), "Time (sec.)", "\u03BC" + "m\u00B3");
            break;
        
        case Roundness:
            double[][] roundness = computeRoundness();
            chart = createChartPanel(roundness, descriptor.toString(), "Time (sec.)", "%");
            break;
        
        case Sphericity:
            double[][] radiivar = computeSphericity();
            chart = createChartPanel(radiivar, descriptor.toString(), "Time (sec.)", "%");
            break;
        
        case Convexity:
            double[][] convexHullDiff = computeConvexity();
            chart = createChartPanel(convexHullDiff, descriptor.toString(), "Time (sec.)", "%");
            break;
        
        case Smoothness:
            double[][] roughness = computeSmoothness();
            chart = createChartPanel(roughness, descriptor.toString(), "Time (sec.)", "%");
            break;
        
        case Elongation:
            double[][] elongation = computeElongation();
            chart = createChartPanel(elongation, descriptor.toString(), "Time (sec.)", "A.U.");
            break;
        
        case Flatness:
            double[][] flatness = computeFlatness();
            chart = createChartPanel(flatness, descriptor.toString(), "Time (sec.)", "%");
            break;
        
        default:
            throw new UnsupportedOperationException("\"" + descriptor.toString() + "\" measure is not yet implemented");
        }
        jComboDescriptors.setSelectedItem(Descriptors.None);
        
        if (chart != null)
        {
            // replace default chart colors by detection colors (taken from t=0)
            XYItemRenderer renderer = ((XYPlot) chart.getChart().getPlot()).getRenderer();
            for (TrackSegment ts : trackPool.getTrackSegmentList())
                renderer.setSeriesPaint(trackPool.getTrackIndex(ts), ts.getFirstDetection().getColor());
            
            chartPanel.add(chart);
        }
        
        super.panel.updateUI();
    }
    
    private void exportShapeToXLS(String xlsPath)
    {
        try
        {
            WritableWorkbook wb = XLSUtil.createWorkbook(xlsPath);
            
            saveToXls(wb, "Surface area", computeDimension(1));
            saveToXls(wb, "Volume", computeDimension(2));
            saveToXls(wb, "Roundness", computeRoundness());
            saveToXls(wb, "Convexity", computeConvexity());
            saveToXls(wb, "Sphericity", computeSphericity());
            saveToXls(wb, "Smoothness", computeSmoothness());
            
            wb.write();
            wb.close();
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
        catch (WriteException e)
        {
            e.printStackTrace();
        }
    }
    
    private void exportMeshToVTK(String vtkPath)
    {
        int cpt = 0;
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            String filePrefix = "_#" + (cpt < 1000 ? "0" : "");
            filePrefix += (cpt < 100 ? "0" : "");
            filePrefix += (cpt < 10 ? "0" : "") + cpt;
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                int t = det.getT();
                String fileName = filePrefix + "_T" + (t < 10 ? "000" : t < 100 ? "00" : t < 1000 ? "0" : "") + t + ".vtk";
                try
                {
                    if (det.getDetectionType() != Detection.DETECTIONTYPE_VIRTUAL_DETECTION)
                    {
                        BufferedWriter writer = new BufferedWriter(new FileWriter(vtkPath + fileName));
                        ((Mesh3D) det).mesh.saveToVTK(writer);
                    }
                }
                catch (FileNotFoundException e)
                {
                    e.printStackTrace();
                }
                catch (IOException e)
                {
                    e.printStackTrace();
                }
            }
            
            cpt++;
        }
    }
    
    private ChartPanel createChartPanel(double[][] array, String title, String xLabel, String yLabel)
    {
        XYSeriesCollection plot = new XYSeriesCollection();
        
        for (int track = 0; track < array.length; track++)
        {
            XYSeries series = new XYSeries(track);
            
            double[] row = array[track];
            
            for (int frame = 0; frame < row.length; frame++)
            {
                double value = row[frame];
                
                if (value != 0) series.add(frame * tScale, value);
            }
            
            plot.addSeries(series);
        }
        
        JFreeChart chart = ChartFactory.createXYLineChart(title, xLabel, yLabel, plot, PlotOrientation.VERTICAL, // orientation
                true, // include legend
                true, // tooltips
                false // urls
                );
        
        return new ChartPanel(chart, 500, 300, 500, 300, 500, 300, false, false, true, true, true, true);
    }
    
    private void saveToXls(WritableWorkbook book, String pageTitle, double[][] array)
    {
        WritableSheet sheet = XLSUtil.createNewPage(book, pageTitle);
        
        XLSUtil.setCellString(sheet, 0, 0, "Track \\ Time");
        
        // time labels
        for (int t = 0; t < trackPool.getDisplaySequence().getSizeT(); t++)
            XLSUtil.setCellNumber(sheet, t + 1, 0, t * tScale);
        
        // indexes + data
        for (int i = 0; i < trackPool.getTrackSegmentList().size(); i++)
        {
            XLSUtil.setCellNumber(sheet, 0, i + 1, i);
            double[] row = array[i];
            for (int j = 0; j < row.length; j++)
                if (row[j] != 0) XLSUtil.setCellNumber(sheet, j + 1, i + 1, row[j]);
        }
    }
    
    /**
     * Measures the convexity of all objects over time. The convexity is computed as the difference
     * between the object and its convex hull, expressed as a percentage of the convex hull. Hence,
     * the value is 1 if the object is convex, and lower than 1 otherwise
     * 
     * @return
     */
    private double[][] computeConvexity()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Mesh3D contour = (Mesh3D) det;
                
                Point3d[] points = new Point3d[(int) contour.getDimension(0)];
                int i = 0;
                for (Point3d p : contour)
                    if (p != null) points[i++] = p;
                
                QuickHull3D q3d = new QuickHull3D(points);
                int[][] hullFaces = q3d.getFaces();
                Point3d[] hullPoints = q3d.getVertices();
                double hullVolume = 0;
                
                Vector3d v12 = new Vector3d();
                Vector3d v13 = new Vector3d();
                Vector3d cross = new Vector3d();
                
                for (int[] face : hullFaces)
                {
                    Point3d v1 = hullPoints[face[0]];
                    Point3d v2 = hullPoints[face[1]];
                    Point3d v3 = hullPoints[face[2]];
                    
                    v12.sub(v2, v1);
                    v13.sub(v3, v1);
                    cross.cross(v12, v13);
                    
                    double surf = cross.length() * 0.5f;
                    
                    cross.normalize();
                    hullVolume += surf * cross.dot(new Vector3d(v1)) / 3.0;
                }
                
                result[trackIndex][det.getT()] = 100 * contour.getDimension(2) / hullVolume;
            }
        }
        
        return result;
    }
    
    private double[][] computeDimension(int order)
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                double value = ((Mesh3D) det).getDimension(order);
                result[trackIndex][det.getT()] = value;
            }
        }
        
        return result;
    }
    
    private double[][] computeElongation()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Point3d radii = new Point3d();
                
                computeEllipse((Mesh3D) det, null, radii, null, null);
                
                result[trackIndex][det.getT()] = radii.x / radii.y;
            }
        }
        
        return result;
    }
    
    private double[][] computeFlatness()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Point3d radii = new Point3d();
                
                computeEllipse((Mesh3D) det, null, radii, null, null);
                
                result[trackIndex][det.getT()] = 100 * (1 - radii.y / radii.x);
            }
        }
        
        return result;
    }
    
    /**
     * <b>Description: </b>The smoothness provides a global dimension-less measure of the local
     * contour smoothness, i.e., regardless of its global shape. It is essentially given as the
     * inverse of the so-called <i>average roughness</i> measure (following the ISO 4287 standard),
     * which measures how the contour deviates from a hypothetical "smoothed" version. This
     * descriptor is local, such that the roughness is first computed locally for each surface point
     * according to a given sampling distance, while the global roughness measure is given as the
     * average roughness over the entire surface.<br/>
     * <b>Calculation: </b>For each surface point, the local roughness is given by the standard
     * deviation of the distances from the considered point and its neighborhood to the object's
     * mass center. For a regularly sampled mesh, the neighborhood is defined here by the first two
     * rings of neighbors around the considered point. Finally, the roughness measures for all
     * points of the surface are averaged, and the smoothness measure is calculated as
     * <code>smoothness = 1 / (1 + roughness)</code><br/>
     * <b>Interpretation: </b>The smoothness measure calculated here is simpler to interpret than
     * the raw roughness: its value ranges from 0 to 1, where values close to 1 indicate a smooth
     * object, while values close to 0 indicate objects with a rough surface (e.g. spikes)
     * 
     * @return
     */
    private double[][] computeSmoothness()
    {
        ExecutorService service = Executors.newCachedThreadPool();
        
        final double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        ArrayList<Future<?>> results = new ArrayList<Future<?>>();
        
        long i = System.nanoTime();
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            final int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                final ROI3DTriangularMesh mesh = ((Mesh3D) det).mesh;
                // final double samplingDistance = contour.sampling.getValue() * 3;
                
                results.add(service.submit(new Runnable()
                {
                    @Override
                    public void run()
                    {
                        // Point3d center = contour.getMassCenter();
                        
                        double[] localRoughness = new double[(int) mesh.getNorm(0)];
                        
                        int n = 0;
                        
                        // compute local roughness for each mesh point
                        for (Vertex3D v1 : mesh.vertices)
                        {
                            if (v1 == null) continue;
                            
                            // 1) Build a neighborhood for each vertex w.r.t. the sampling distance
                            ArrayList<Vertex3D> neighborhood = new ArrayList<Vertex3D>();
                            // neighborhood.add(v1);
                            
                            // add the first ring of neighbors
                            for (Integer n1 : v1.neighbors)
                            {
                                Vertex3D vn1 = mesh.vertices.get(n1);
                                if (vn1 == null || vn1 == v1 || neighborhood.contains(vn1)) continue;
                                
                                neighborhood.add(vn1);
                                
                                // add the second ring of neighbors
                                for (Integer n2 : vn1.neighbors)
                                {
                                    Vertex3D vn2 = mesh.vertices.get(n2);
                                    if (vn2 == null || vn2 == v1 || vn2 == vn1 || neighborhood.contains(vn2)) continue;
                                    
                                    neighborhood.add(vn2);
                                }
                            }
                            
                            // 2) Compute local roughness for the current neighborhood
                            
                            double[] distances = new double[neighborhood.size()];
                            for (int i = 0; i < distances.length; i++)
                            {
                                distances[i] = neighborhood.get(i).normal.angle(v1.normal);
                            }
                            
                            localRoughness[n++] = ArrayMath.std(distances, false) / ArrayMath.mean(distances);
                        }
                        
                        result[trackIndex][mesh.getT()] = 100 / (1 + ArrayMath.mean(localRoughness));
                    }
                }));
            }
        }
        for (Future<?> futureResult : results)
        {
            try
            {
                futureResult.get();
            }
            catch (Exception e)
            {
                e.printStackTrace();
            }
        }
        long j = System.nanoTime();
        System.out.println("Execution time: " + (j - i) / 1000000 + "ms");
        service.shutdown();
        return result;
    }
    
    /**
     * <b>Description: </b>The roundness (following the ISO 1101 standard) is a measure based on the
     * smallest inscribed and largest circumscribed spheres needed to best fit the interior and
     * exterior of the object.<br/>
     * <b>Calculation: </b>The roundness is approximated here by the ratio between the smallest and
     * largest distance from the mass center to any point on the surface.<br/>
     * <b>Interpretation: </b>A value of 1 corresponds to a perfect sphere, while a value close to 0
     * describes flat objects (e.g. star-shaped objects with flat branches).<br/>
     * The figure below summarises the differences one may expect:<br/>
     * <br/>
     * <img alt="" src=
     * "data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/4gfMSUNDX1BST0ZJTEUAAQEAAAe8YXBwbAIgAABtbnRyR1JBWVhZWiAH0AACAA4ADAAAAABhY3NwQVBQTAAAAABub25lAAAAAAAAAAAAAAAAAAAAAAAA9tYAAQAAAADTLWFwcGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAVkZXNjAAAAwAAAAG9kc2NtAAABMAAABi5jcHJ0AAAHYAAAADh3dHB0AAAHmAAAABRrVFJDAAAHrAAAAA5kZXNjAAAAAAAAABVHZW5lcmljIEdyYXkgUHJvZmlsZQAAAAAAAAAAAAAAFUdlbmVyaWMgR3JheSBQcm9maWxlAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAbWx1YwAAAAAAAAAeAAAADHNrU0sAAAAqAAABeGRhREsAAAA0AAABomNhRVMAAAAsAAAB1nB0QlIAAAAqAAACAnVrVUEAAAAsAAACLGZyRlUAAAAqAAACWGh1SFUAAAAuAAACgnpoVFcAAAAQAAACsG5iTk8AAAAsAAACwGNzQ1oAAAAkAAAC7GhlSUwAAAAgAAADEGl0SVQAAAAuAAADMHJvUk8AAAAkAAADXmRlREUAAAA6AAADgmtvS1IAAAAYAAADvHN2U0UAAAAuAAAD1HpoQ04AAAAQAAAEAmphSlAAAAAWAAAEEmVsR1IAAAAkAAAEKHB0UE8AAAA4AAAETG5sTkwAAAAqAAAEhGVzRVMAAAAoAAAErnRoVEgAAAAkAAAE1nRyVFIAAAAiAAAE+mZpRkkAAAAsAAAFHGhySFIAAAA6AAAFSHBsUEwAAAA2AAAFgnJ1UlUAAAAmAAAFuGFyRUcAAAAoAAAF3mVuVVMAAAAoAAAGBgBWAWEAZQBvAGIAZQBjAG4A/QAgAHMAaQB2AP0AIABwAHIAbwBmAGkAbABHAGUAbgBlAHIAZQBsACAAZwByAOUAdABvAG4AZQBiAGUAcwBrAHIAaQB2AGUAbABzAGUAUABlAHIAZgBpAGwAIABkAGUAIABnAHIAaQBzACAAZwBlAG4A6AByAGkAYwBQAGUAcgBmAGkAbAAgAEMAaQBuAHoAYQAgAEcAZQBuAOkAcgBpAGMAbwQXBDAEMwQwBDsETAQ9BDgEOQAgBD8EQAQ+BEQEMAQ5BDsAIABHAHIAYQB5AFAAcgBvAGYAaQBsACAAZwDpAG4A6QByAGkAcQB1AGUAIABnAHIAaQBzAMEAbAB0AGEAbADhAG4AbwBzACAAcwB6APwAcgBrAGUAIABwAHIAbwBmAGkAbJAadShwcJaOgnJfaWPPj/AARwBlAG4AZQByAGkAcwBrACAAZwByAOUAdABvAG4AZQBwAHIAbwBmAGkAbABPAGIAZQBjAG4A/QAgAWEAZQBkAP0AIABwAHIAbwBmAGkAbAXkBegF1QXkBdkF3AAgAEcAcgBhAHkAIAXbBdwF3AXZAFAAcgBvAGYAaQBsAG8AIABnAHIAaQBnAGkAbwAgAGcAZQBuAGUAcgBpAGMAbwBQAHIAbwBmAGkAbAAgAGcAcgBpACAAZwBlAG4AZQByAGkAYwBBAGwAbABnAGUAbQBlAGkAbgBlAHMAIABHAHIAYQB1AHMAdAB1AGYAZQBuAC0AUAByAG8AZgBpAGzHfLwYACAARwByAGEAeQAg1QS4XNMMx3wARwBlAG4AZQByAGkAcwBrACAAZwByAOUAcwBrAGEAbABlAHAAcgBvAGYAaQBsZm6QGnBwXqZjz4/wZYdO9k4AgiwwsDDsMKQw1zDtMNUwoTCkMOsDkwO1A70DuQO6A8wAIAPAA8EDvwPGA68DuwAgA7MDugPBA7kAUABlAHIAZgBpAGwAIABnAGUAbgDpAHIAaQBjAG8AIABkAGUAIABjAGkAbgB6AGUAbgB0AG8AcwBBAGwAZwBlAG0AZQBlAG4AIABnAHIAaQBqAHMAcAByAG8AZgBpAGUAbABQAGUAcgBmAGkAbAAgAGcAcgBpAHMAIABnAGUAbgDpAHIAaQBjAG8OQg4bDiMORA4fDiUOTA4qDjUOQA4XDjIOFw4xDkgOJw5EDhsARwBlAG4AZQBsACAARwByAGkAIABQAHIAbwBmAGkAbABpAFkAbABlAGkAbgBlAG4AIABoAGEAcgBtAGEAYQBwAHIAbwBmAGkAaQBsAGkARwBlAG4AZQByAGkBDQBrAGkAIABwAHIAbwBmAGkAbAAgAHMAaQB2AGkAaAAgAHQAbwBuAG8AdgBhAFUAbgBpAHcAZQByAHMAYQBsAG4AeQAgAHAAcgBvAGYAaQBsACAAcwB6AGEAcgBvAVsAYwBpBB4EMQRJBDgEOQAgBEEENQRABEsEOQAgBD8EQAQ+BEQEOAQ7BEwGRQZEBkEAIAYqBjkGMQZKBkEAIABHAHIAYQB5ACAGJwZEBjkGJwZFAEcAZQBuAGUAcgBpAGMAIABHAHIAYQB5ACAAUAByAG8AZgBpAGwAZQAAdGV4dAAAAABDb3B5cmlnaHQgMjAwNyBBcHBsZSBJbmMuLCBhbGwgcmlnaHRzIHJlc2VydmVkLgBYWVogAAAAAAAA81EAAQAAAAEWzGN1cnYAAAAAAAAAAQHNAAD/4QD2RXhpZgAATU0AKgAAAAgABwESAAMAAAABAAEAAAEaAAUAAAABAAAAYgEbAAUAAAABAAAAagEoAAMAAAABAAIAAAExAAIAAAAeAAAAcgEyAAIAAAAUAAAAkIdpAAQAAAABAAAApAAAAAAAAABgAAAAAQAAAGAAAAABQWRvYmUgUGhvdG9zaG9wIENTMiBNYWNpbnRvc2gAMjAxMzowMjoyOCAxNzo1MDowOQAABJAEAAIAAAAUAAAA2qABAAMAAAAB//8AAKACAAQAAAABAAABlKADAAQAAAABAAABSAAAAAAyMDEzOjAyOjI4IDE3OjMzOjIwAP/hAfZodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIj4KICAgICAgICAgPHhtcDpDcmVhdGVEYXRlPjIwMTMtMDItMjhUMTc6MzM6MjA8L3htcDpDcmVhdGVEYXRlPgogICAgICAgICA8eG1wOkNyZWF0b3JUb29sPkFkb2JlIFBob3Rvc2hvcCBDUzIgTWFjaW50b3NoPC94bXA6Q3JlYXRvclRvb2w+CiAgICAgICAgIDx4bXA6TW9kaWZ5RGF0ZT4yMDEzLTAyLTI4VDE3OjUwOjA5PC94bXA6TW9kaWZ5RGF0ZT4KICAgICAgPC9yZGY6RGVzY3JpcHRpb24+CiAgIDwvcmRmOlJERj4KPC94OnhtcG1ldGE+Cv/bAEMAAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQICAgICAgICAgICAwMDAwMDAwMDA//AAAsIAUgBlAEBEQD/xAAfAAABBQEBAQEBAQAAAAAAAAAAAQIDBAUGBwgJCgv/xAC1EAACAQMDAgQDBQUEBAAAAX0BAgMABBEFEiExQQYTUWEHInEUMoGRoQgjQrHBFVLR8CQzYnKCCQoWFxgZGiUmJygpKjQ1Njc4OTpDREVGR0hJSlNUVVZXWFlaY2RlZmdoaWpzdHV2d3h5eoOEhYaHiImKkpOUlZaXmJmaoqOkpaanqKmqsrO0tba3uLm6wsPExcbHyMnK0tPU1dbX2Nna4eLj5OXm5+jp6vHy8/T19vf4+fr/2gAIAQEAAD8A/vwoooooooooooooooooooooooooooooooooooooooooooooooooooooooooor5K/bz+KfjX4IfsWftTfF/4caxFoHxA+G/wL+JHjDwXrk+nadq8Wk+KND8M6he6HqEmlavY6lpWoLbajFG5huYJbeTGJF2k18tfsufCr9rz4yfs4fs+/GDxd/wUg+PSa78V/gv8NPibrtt4f+CP7FVhp1vqXj/wRofimex0v7Z+zbq7R6fp0+qskRcyvIqgs3OK+j0/Zs/aCXdv/wCChf7TEmen/FsP2JEx+X7KJzUy/s4/H4fe/wCCgf7Sb/Lj/kmX7Fg+f+/x+yx+lUv+Gbf2j/s0kf8Aw8M/aH+0H/V3I+Ef7GWE/wB6E/s0lJPzWib9m79o9otsH/BQz9oWKX/nrL8Iv2NJu/P7tf2a4ByPeli/Zv8A2kFi2y/8FCv2gpZf+eo+EP7G8ef+AD9m8iqcv7NP7ThGLf8A4KL/ALQKepuPgx+xpOT/AN+v2crUV8m/tt6B+2b+zF+yp8dvj74X/wCChXxO1fW/hf4Dn8S6RpXiP4A/spTaTd31rqFlCY9QOnfCHS7lorpLkodsiCMndzjB+nIP2Zf2uUaMz/8ABSb4zzhUjEip+z7+yJCJHSLZK6/8WVkaMSy/PjLY6Diph+zL+1jjDf8ABSP42k785X4Cfsfqdv8Ac+b4GOPxq1J+zX+1MxTy/wDgov8AHCIK2ZM/A39kFzIv90k/AgbT7gVF/wAM0ftWDaR/wUd+Nxw+59/wI/ZAO9Of3fy/AqPaBnr1p9x+zV+1U8OLb/goz8bILjH+tk+BX7IlxAffyB8DYZP/ACLWdD+zT+2KhzN/wUi+Jsw5+Ufs4fstR/qPhsx61D/wzV+2iFXH/BSPx4XCbSX/AGY/2ZzGz/39ieDY2B9g+K+YP2v9F/bx/Zo+Ams/FnQ/+CguteJtW0/4gfA7wdHY63+y18A/sIt/i58dfhz8I7y+kXStLsJ1XRLLxy94Pmbcbcg8N8v0Bd/sxft9yXSzWf8AwU7122gEcitZy/sf/s8XMbSNny5TL5EFwBH6BsNUH/DMf/BQcpID/wAFPr0O/wDq2T9jL4CgRd+Ve/ff+JFQXP7MX/BRJ4nFr/wVDMEzH5JJ/wBij4GXMaZ/6ZJr1sz/APfdWP8Ahmr/AIKI793/AA820/Z/c/4Ym+D+P++v+Ez3frTpP2af+ChjRYj/AOCm0Ec3/PQ/sV/Bh4/+/R8UBv8Ax+mJ+zR/wUQCusn/AAU4tXLMhSRf2KPg0kiKM7058VvG2/I5K8Ypzfs4f8FGVJ+z/wDBS3ww4x01D9hz4a3P/pF8UdLNfN3xL1X/AIKR/C39qP8AY/8AgEv7bnwb8TWP7RsX7QE+ra7qn7F0Fi2if8KZ+H+ieLbFE0vSfjmBqn9uXWqmOX/iY6WYUXdH5p+QfZ0Xwm/4KBIcv+2j8Aph6P8AsRa6v6xftcR0zTPhN/wUFsZd95+2n8ANaj/54337D+vW4P8AwPS/2vNPfJ+tbkHgH9vu0s1jP7UX7K+qXYPz3OofsV/E63Dj1EenftwwIuPTaazv+EA/4KIR3iFP2pv2SZ7D/TvMFz+xP8VVulLj/QAnkftxxRyLEfvksv8AwKs+98A/8FJfMWOw/ag/Y4e3/iuLr9i/4vQ3fviKP9tq6tj+lVB4K/4KbxR71/aN/YkurjzM+TP+yL8bYLby+6+ZB+2Q8+f85rWl8Gf8FIJJrcp+0b+xXBCP+PgJ+xv8bnZ/9wS/txv/AOhLXzR+1/8AFT/gpH+yh+zz8Qvj3H8X/wBin4hjwNP4TRfC9x+yn8b/AAsL9fF/xB8KeCIidYg/bN137ONKj8QG4P8Ao7eYVKHaMOPdZ/h9/wAFPdTlDSftSfsY+GYu8Wh/scfFzWSfffrv7YiHNZA+FP8AwVLSfeP2zP2RZoP+eM37EfxATP4w/teB/wDx6tAfDz/gqLaXEjwftT/sX6tA5+WHVf2OPi3ZeTx/A+l/tgh35/vE1U/4V7/wVT2f8nSfsS+Z5uf+TPvjBt8r0z/w1znd7Y/GqV58PP8AgrBI872f7Uf7EMCv/qIn/ZI+LzLD0/jb9qaV2/HdVK48Bf8ABXR7u1Nv+0p+wVFYx/8AHyG/ZV+NzXU/03ftOvGn4Faq3Pgb/gsC7iO3/aI/YEig73Dfs1fHE3X/AH6b9oKS3P5ivnv9oL4m/wDBVH9miD4K/EP4jfG/9jHxJ4H8SftJfs6fCHxj4O8A/AX4qeH9d1vRfiv8VPDPgrWv7D1rxH8UvG32G/TTdTnkUtEAgBcMGVVP7eUUUV8Z/tk/t3/An9hzQvh1e/FlPiH4r8Z/GXxk3gD4M/Br4L/D7X/iv8Z/iz4rtrFtW1bTvAvgDw3FJf6jB4e0ZDd6hdzPb2dpEUV5RLNBHL8o/tZ/8FmvgD+yh8d/HvwAb4DfthftFeIvgh8PfDvxW/ad8V/swfBGD4o+BP2XfAHiyxv9Z0DXPjVrdx4u8OXejS3/AIZ06bWxa2NrqFx/YkM10FPkvHX6m+AfHfhD4peBfBfxN+H2v2PirwF8RfCfh3x14J8UaW0j6b4j8I+LdIs9f8Oa9p7zRxStZavo+oQ3ERdFbZIMgHiutoor4C/4Kqxxy/8ABNn9uSOV/LR/2Y/i8Gfj5f8AikNSIPPH3q9N/YL5/YZ/YwICqD+yf+zqQqJ5SKP+FQeDsKkfWNR2XsOK+sKKKKKK/OD/AIK8Osf/AATX/a8d2kRV+FzEtF/rB/xUeggFffJr9H6KKKKKK/Nj/grfcLafsK+PbpovtHkfF/8AZClS3882/nun7YvwEZY/NBGM4zjviv0noooooor8zv2n45Jf+Cjv/BLcquY4LX9uK5kf+4w+CnhK1Qf8CN2fyr9MaKKKKKK/NL/gr+7J/wAE8vjuU+82rfAyMf8AbX9oj4TRH9Hr9LaKKKKKK/NH/gqg2fgX8Erfr9t/bo/YbtcKcTtv/aW+Hz4te5uPk49s1+l1FFFfjR/wUc8K/Hf4WftZ/sM/8FAvhH+zl4+/a48J/s0eFf2oPhL8Wfgv8I7rwzN8Z9F0D9ovS/hZ/ZHxX+E/hnxbq/h7SvF+p+H9S+GTaZq2nRX0N7PYaqjxDyobh0+I/gXq/wC1R+yx+y5+2x+178VP2DP2gvid+2d/wVR/au+IuseAP2YPhp4Z0j4h+K/hZ8Prj4YN8Ov2XPBX7TXjK31G08J/Df4deDPDnhOR9e128mNvpa64kUlslyz2yftX/wAE8f2fPFX7KH7Cn7I37NnjrUrbVvG/wT/Z8+Fvw88ZXthdNe6YfFfh3wnptn4itdIvGWNrrRdP1hZreykKqXtYoztXOB9kUUV8Sf8ABSvTF1b/AIJ4/txWjnCxfsofH3Uv+BaN8MfEusIPxexArsP2Erea0/Yg/Y2tblt1xbfsqfs8wTt13TQ/CPwhHKfxdTX1XXk/x0+M/gb9nb4P/Eb44fEm/k0/wT8MvCup+K9ektxC99dxWEJNpo2kQXM9rDea/wCINRkhsdPtzJGbm9uIogQXFfy2w/8AB2J4CtbxrDxN+w78QtDuopBFeRW/xp8Mas2nseD9sx4EsvKx+Nem6V/wdV/swv5U/iL9mH4+6fYGfyZ73w/q3w/8Spb/AO0yXGs+H9//AH0K/pa+E3xV8B/HH4Z+Bvi/8MPENl4q+H/xG8NaV4s8Ka9YPmK+0nVrdLiETRk+bZahaMzQXdrKFuLS6ikgmVZY3Ueh1+aX/BYeKOf/AIJqftYQyxGaOXwPoMbRj+Lf8QPCAH5Mc1+ltFFFFFFfnF/wVdsotS/Yw1/TpJYYZtR+PX7GVhZyT25uUF3d/tk/AW3iHldyRI3PpxX6O0UUUUUV+bH7R1g99/wUl/4JoSrEWTSPA/7d+sSSBsCJR4B+DOhrlf49769gelfpPRRRUccsco3RSJIoJUlGDAMOq5UkZBqSivzN/wCCwYB/4J6/GsPKsIPin9nf960P2hUJ/ab+DeCYTw/PHt1r9MqKKZJIsSNJI6xxorPI7sFRUUFnd3YgIqAZJPAFJFLHPGssMsc0TjcksTrJG45GVkQlGGR2qSivzG/4Kv3bWP7Pfwhu0TdLB+25+xFJE/GYZF/aT+H+2UZyMjp+NfpzRRRX40/8FGfE/wAevir+15+wZ+wJ8Iv2g/H/AOyr4L/aI0L9pj40fG34wfCQaDp/xi1rwf8As36b8LE0f4T/AAs8U+JdK13T/C2oeKvEHxRjvdXv4bSW9t9O05VQGKadH+DfFPgL4+f8FNf2mf8AgqN8YdE/bz/ay/ZS+Cv7A3jG5/ZZ/ZQ0n9nH4wXHww+HX/C/fhB8KrTx38fviT8efDdrpdzp/wAXNP0P4jeJrLS/sN7IkP8AYkd5ZuwZkeL9uv8AgnT8fvGH7VH7Bv7H/wC0b8Q7OKx8ffGf9nb4UfEDxtHbWken2Vz4s8QeD9Lu/EWqaZYRJHFZaTrOrNNeWkKjbFbTxqCQAT9nV8y+E/2x/wBmTx18ada/Z58KfF7w5q/xf0LUPFWj3vhOK31q2jufEHgP7OfHHhjQ/El7pVt4T8T+LfBaXKyavpOm393qWmxJI9zBGkUrJ6T8YvjV8K/2fvAOrfFD4zeN9F+H/gTRZtPtb7xBrkk4hN9q17Dp2laXYWVnDd6lq+r6pfXCRW1naQz3U8jBURjXyh+1X8V/hz8e/wDgmr+2J8Qfg74z0H4geDvEX7JH7TtppOv6BcfbLRr+0+EfjjTr/TL22dYrzTtX0vUY3t7uxuo4bq2mRo5Y0cFR7V+xdC9v+x1+ydbyDEkH7NPwJhcejxfC7wsjD/vpa+lq/j9/4OeP25LWxX4UfsEeENdhiTVGt/jX+0RJp+pi11Gz0LT4dSX4T+C3eC+ZRHqesWl3rmpW9xbrJEtno80blZyD/EFd+Jzq2r6lP9ruzZ5vJ/3+Lkjnk5570XOuX09zbXHnwhJjacxQfZzzn2Ga/pd/4Im/8FqIf2E5bT9nT4/rrfiH9lXxfr2p67Z+K9Osra/1b4AeItXurb+0Nc+wWWL3Vvh1r15MJ9Ss4w8trIGv9OSXzLi3uf8AQM0PXNG8T6Jo/iTw7qlhrnh7xFpWn65oOtaVdw32maxo2rWcOoaZqmm31u8lveWGoWVxHLDKjMkkbhgSDX5v/wDBY9tn/BNf9p5jbSXijRfh75ltHjzJYj8X/h8JgueMiMk/hX6bUVxfxF+Ifgv4S+BPF/xO+I/iLT/CXgPwH4e1XxV4s8SarI0djo+h6NayXuoXkwjSSeZlhiIjiiSSaeQrHGjyMqn/ADa/2nv+C2X7cHxr/ac1b9oP4dfH34rfA7wjoHiWVfg98LPCniO5tvBfhbwfbarHPpWleMvC9nEfDnxC8Q+ILCJZNVOq6ZqMd8QqRokMccSfpR8Ff+DqT9pbQNd8HWnx++B3wY+IXga2tbHTvGOt/DaHxj4G8f6lLFCBe+JNPXVfEXirwe995AFzNYRWdvbPISkc0CMuz+3TwJ438NfEvwP4N+I/gzUk1nwd8QPCfh3xv4T1iJHji1Xw14r0iz13QtSjSVUkjS+0u/ilCsAwDciusr8/v+CmvkN+yvFBcRvJFdftLfsRW7BOqh/2z/gI24+w2Y+pr9AaK5bxr418I/Djwj4k8e+P/Eui+DvBfhDR77xB4o8U+I9QttK0PQdF02B7m+1PU9Qu5Ire0tLaGMszMw9BzxX8FH/BTP8A4OK/j78bPiNqPgb9ifx94o+AXwI8PX0mn6b4w0fTtP0z4n/FO+sp1lXxNf61qCTat4J8NX4ULpWl2KpeXClm1IhXW3h5j9h//g5V/ax/Z7v7nw7+01cTftcfDPZ5cU2vXej+Gvix4duka6d7rRfH1lo1tB4i0e5kud8tvrtjeXuI40juLGMMjf02fsNf8F7/ANiv9uD4neGPgdpVj8S/gz8YvGNpO3hLw18UtJ8ProPjPV7S1kvr3QPB/i3wv4i1+3u9VtrGJ5Fi1O20iScoUiV5SiP+3dfnX8eUlf8A4KS/8E8Smdkfwb/b2lmx02/Yv2YYlz6/vJhX6KUV/MZ/wV9/4L++D/2Wl8S/s+fsdaz4f8eftD6fLc6R44+JUlraa78PvgrqEb3NtPotuLwNoni/4l2s0BL2cnm6bp5Gyf7Tch7NP4Zvjl+1X8e/jz4lPi34zfGn4jfFnXpNtrZ6l468beINdeztYyxXR7SHUtWWHRLZS5ISMBVJ4qx8Df2l/jJ+zH8RPD/xd+BvxL8S/Dn4geG9Ts7u31Xw7rV0trqcGn6gp1XR9etLRZLbxB4f8Qbcapp2oq2lyLkMCOK/1Q/2Cv2uPDH7cn7Jvwd/aW8M2q6TJ4+8NhfFvhvcm/wp4/0K5m0Pxt4eCC7vZ10618Q2M0mmyTus91pU1rcsiecFH2DX5mf8FhGjX/gnx8ZjKJCg8Zfs27hDD9pk5/aj+CwG2H+Pk8+g5r9M6K+Kv+Cgf7bPw+/4J/fsveO/2jPH9u2sS6O1n4b8A+DYLhLW98ffEnxCLiLwr4RtLmT5baO4e2mu72bDNbabZ3MypI0axv8A5nn7aH/BS/8Aav8A23/HWq+JPj58Uta1LQ4tRj1Hw58MNAvNR0T4ReCnigS1jt9A8EPrMmmQ3pgjVZLq4aXWbgjdJK7ZY/OvwD/a0+Of7MvjWz+IfwO+K/j34X+L7BzcHUvCfiu9sLa9t2cO2ka7olmH0HX9AdlG/T9RVkYDBWv9IX/gin/wVKT/AIKYfs86xfeOLHSdD/aH+DNzoeg/GLTPD9tcWnh3XYPEMWqP4R+IHh+1uWmOnW3ieLQryO908SyNYahaTDCQy261+ztfmX/wVdfH7OfwwXy5JPM/bS/Yii2x9Tn9pz4bHJ9VBH51+mlFFFfkv/wUa+HX7H/7Sfxx/Yl/ZA/aC0j46aH8cfi9qH7QnxI/ZX+Pf7P3jO8+Ffjv4Ga38CfA/hXVPihf2nxX8PeK9F8XeF38Y+FfGVtZR2dvp2s2GpPFm7ige2tZ1+Vf+GNP+Cefib9hH9oj9jz4O+Pv2tvht+zl+wh8aviwn7X3hj4O/E3WfBHxF/aU8ZaV8I7X4tfGL4d/GT4j+IbWbxF8V/B/xS8KfE6zfVRaajo6XkyQ2cd1a2tqIl/Yn9kHxv8AC74m/smfsv8AxJ+B/g+7+HnwW+IX7O/wU8c/CHwBf2GnaVfeB/hh4t+G3hrX/APhC90zR9S1jSNOvPDXhXULSylgtby6t4nhKRzSoFdvomvwy+Df7Pn7Q+iaJ+xj+x3qnwX8XeH9F/Y6/aYu/jN4w/aiu9d8FRfDn4ieA9APxtn8Kah4Gt9E8bzePdQ+IfxeuPHEFv4j02/0mGLS4bu/a6lnWa3ab7F/aS0DUf2nPgjrt4/wQ/aM0Hx18Af2nNK8SfCe08JXnwZ8MfFS+8T/AAX+IUWmaL8dvhQvxd1jUvhH4m8Kah4Zvr7U9MsPFCRRazZ77Z7eKWSGUfIfhv4NfF74Hf8ABNH/AIKm+JfjPN4quPFXxtsP22/2g9L0/wCI1z8NoPiDp3hXxJ+z9a+H9AtPiGnwjt4/hdofi3UY/BLX9/baIZbCze82B9yui/pv+x5u/wCGR/2Wd/3/APhnL4Ib/wDe/wCFZ+GN3619G1/kmf8ABSHUvjXrP7dH7WEvx/uILr4rS/Hb4kaf4sCXT3toq6Hri6LoGl6A0jNJp3hqDw9pdmmjgnjSQg6CvheC0sEmkjPneZdHjPHbHU4Nb+kx27m+sU8kST/6me+gwFBzj68ivqr9kD9nPxx+1p8evg/+zb4Aljj8XfGDxaPD8ep3dnqVxaeHNEa1mvfGXirU4oUeeTR/DPhy2vNRYICdtmcV/rKfCv4ceG/g58Mfhz8I/BqXcXhD4XeBfCfw88LR6hOl3qI8PeC9B0/w5o39oXUcNul1ff2fp0fnSiNBJLltozivg3/gsPDDcf8ABOD9pS3uJfIhn0/4awyTf88lk+Mvw6Tf2+6TX6YUV/mtf8Fov+Cm37W37T/7Qvxm/Zz+Imq6f8O/gn8CfjJ8SPh1ovwf8D3WsW+ha5rHw+8caz4XtPGXj3WbiK2vvH/iaCHS1uYRPFFpOnhPMsrHzpJZJPxKvNK1Oy03SrybUNImi123vb+Cy0vUP7WvNN/s2/bTGOp6cOdMJI79RXPI032u5g8uaSO1m/10PSf6ZxX+ih/wbw/8FEvBv7Sn7MXhf9kTXLTWNF+Of7KHw+0vw9Ob0G90X4gfCfRNQh8PeFvGHh/U4YhHaXfhy1vNP0rVNNnPm28zQTRPKk7x239FdfCP/BR+N5f2YrdYxM7/APDTn7CJ8qA4lmB/bm/Z0Xyl5B+bd69RX3dXnfxc+Jfh74LfCn4nfGLxdFqc/hT4T/D3xp8S/E8Oi2aX+sy+HvAfhzUvFOtRaTYvNape6nJpulSC3iaSMSS7VLLnNf5n/wDwUQ/4LAftU/8ABQjV9R8PePdfPgD4ILrC6j4Y+BPgq6Nl4Us5NLWIWM/jXUXRdd+IervJG0omuxFZwyu/2eCBG21+Q/2vVvD19bPc3dpqWPtf2fyZzd+RjrgnHXNV9SurWC3sp/ImzeT3v2j9x9kFwB1PTIPNfo5/wSK0+/1r/gpx+wvp2ipcecn7Q/gbWnWKYssGnaDqFxrmrrIOcx/2NpN1uH92v9Wivz7+Ndm13/wUe/YJfJ26Z8Bv289UOHxk/bP2SdIAZc/OP+JyfpX6CV4/+0HYfEnVPgJ8b9M+Dd4+n/F7UfhB8S7D4VahHPFayWPxIvPBmtW/ga8S5uHjt7eS28USWrh3YIhXJIANf4/PjrT/ABBaaxrOneIVvdP1uy1G/wBP1jSr+AjWbXW9M1ArqenanvG4asNX+05B5Brzy4SN4pZDJ5eMZzwbjPPvVfT5k3yf62WPJP74kA+2evev9Ab/AINPPi0viL9kT9on4LXMwn1D4V/HyDxrblXBW28N/F3wfptvpVls/hKa58NtVnz3+0e1f1V1+an/AAV5iSf9gH4uwyfcm8dfszRN9JP2qPgov8zX6V0V/nv/APB0J+2Nqfxi/a70T9mXw1qzJ8P/ANlfQ0tNatYLjzYNU+MHxI0bRtZ17WGFmgkuIfDnhHUbDSljlZjb3S34XaXcH+WKQPdXBkmkm9uc59D65zVaBL+Z5ERPM8q37f8ALD8vav6zf+DSifxAn7aH7QFoJpE8OP8Asp6pd31qn2jyZtTHxi+GyeHbqdZcp9pisZdSVcc7HNf381+ZH/BVdv8AixvwIt1lkjlvf28v2FrOBY/vTSyftKeAXEX0KRsfwr9N6KKK/C3/AIKiH9oD4Z/t2f8ABLn9rP4P/sl/HP8Aa08Ifs+aD+3bofxN8L/AW08J3finQpvjV8Pvgv4T8EXM6+MfEvhfSkt76+0m9k5uNxjs5cAttDfmr8L/ANor9t3wL4N/4Kp+G9Q/4I7f8FCLy6/b1+Ovxn+Knw6ns9G+C/2fwbo/xK/Zq+FvwT0rTvGvn/FaKQalYa94EuLu4+wrdxfY5o9jtJuRf6Lv+CcPgDxl8KP+CeX7Bvwt+I3h7UPCPxC+Gv7GP7LvgDx34U1ZI49V8MeMvB3wP8DeHfE/h7U44pJYk1DRdb06e2mCuyiSIgEjmvs6iivjT/goxeW9h/wT5/bpvLpxHBB+x5+0u8jH3+DHjRVX6u7AD3Nd7+x0Qf2Rv2WCHWQH9nH4HkOn3XB+GXhjDr/st1HtX0dX+dN/wc8fs9an8Lf+CjJ+MGlE/wBk/tI/C/wV44hNrALf7N4o8AabH8K/EdiSPknljsPDekag7gZL35z8xJP831u9w8MtpvhiB3Gc9uvzZ5znFb+m6TO93HHFmRBBeefPz67gT97jJr+5z/g2E/YJHg/wJ4z/AG9viBpFt/bPxMtLv4afAKK7tIWuNK8AaNq1wnxA8a2bywytAfFvim0Gk2EsbxTw6bpU8bAxXeK/rfr8y/8AgsZKkP8AwTe/aSkli8+L7L8LUmh/56Qy/G34bRSp+Mbmv00or/Mc/wCC9nwF8S/BP/gqJ+0UdR0q/tfCvxd1bSPjh4J1a6tlEPiHRfiBomkSeKdR07a7AroHxFtNa07JILfYsnGa/ILTZMeYEEVtLLOTcXucXX2P/nxGAe/51r6Tq1vPdyjUIbS1udJ+a3/0f/j45/4l1hqIIBA5r+mf/g1gg8St+3X8X7yy0+4uPDsf7KnjCy8SasrldJgv7n4tfCW60CGNRgSavK1jdAr/AM8d57V/fFXwZ/wUflij/Zx8PtNCZ0P7Vn7CmYhkb8ftq/ANsZ/h+7X3nXN+MfC+l+OPCPinwXrttb3uieMPDet+FtYs7qAXNrd6V4g0y60rUba6tmdFuLee0u3V0LAMpIzzX+Rf+0X+z58R/wBlr42+NPgj8cvCmo+FPiH8OdXu9JvNL1IIttqETbTpPifRbmEvBrngvxHoxF/pckZZJEYEHBr5hkNpbrG/l3fmRcf67PA9s443Veu47XUrn+0pILs+VA2f35tQbs8jsa/qa/4Ng/2GvF3xD/amvf21vE2jarZ/Cj9nzw34l8NeANbu7dY9L8VfGLxrpGp+ErvTdKeVGXUY/BPw71y9k1GRDuhu77Tj3BH99VfA/wAY7eGf/gov+wzI5xLZ/s/ft3XUI55Laz+x/Zye3CXdffFFfwj/APB09+yT8NPht8a/gF+0R8PvB6eGtZ/aF0n4o6Z8W59DtY7TRNc8a+AZvAtzo3i28tIlWC38Y+JNL8YTw384AOqRaahkDTiWST+Re+sJLfMMckTRmH/XYx+uT0qC1tvssNy/7oSH/lvz1J7fXNf2of8ABoz/AGn/AG7+3J9kYHQP7A/Z4/tYQ4+zDXftvxi/sbYQOWGli7x7V/avX5uf8FbkeX9gn4sRxp5ryePP2Z4wg5J3/tT/AAVU/kpzX6R0V/k//wDBXLTPEGm/8FKf267bXVuFuD+038Tr+JZ/9bPomtaxHrOgCAnrAdA1az2f7OK/MryJ55fLHnfZhPd/Z/O7fzUHFUrSxukaV0TgT+55JIHJ75Nf6FX/AAa4fsReIfgZ+zP8Qv2rPH2lT6T4h/aguPDlh8PtN1C0iiv7X4Q+A7jX5tP8RxyxXTyJZfEHxB4hnnhglhjYWem2s6s6XC4/qYr8u/8Agq/Yx6l8HP2bLKXlbj/goL+w5GcjP3vj34YHcHrmv1Eooor4W/bf/ZYf9ozQPCmv3X7cP7Wn7E3h74SWXjbXPE3iP9mP4s/D74TaV4k0jVbXQrm61H4qav8AED4dePbCXSvBNr4blmsplewjs0vrxpmkV1MX8v7/AAl/aw/bR8fSfDX/AIJF/wDBTr/grV8W/AHhvxN/Y/xR/wCCg37Q/wC0L4E0v9izwwdKvhD4i8NfBvSvDP7PfhP4h/tU+OIUgngz4evtL8O2dwbeSXVZLecSL/Y98F/BHiP4ZfB34T/Dbxh8Q9f+Lvi74ffDTwJ4I8U/FjxVCtv4o+J/iPwn4X0rQdb+IfiSBLq+SDX/ABrqdhLqV4gmmC3Fy4DvjcfS6KK/PD/grTrCaF/wTO/bjvpDGqSfs2/E3SiZc+WP7d0G40IZx3zqXy+9e9/sZZ/4Y+/ZRy24/wDDNfwLy394/wDCr/C2W/E19KV+DX/Bcz/glr8Y/wDgpR4P/Z9b4E6/8MdC8b/BjxD4/a+j+J+t+JvDml6p4a8e2HhY3MNlqvhbwr4rvRe2mreDLRlSSBE2SuQ+eD/Lz4i/4NjP+Coel+dPpeg/s8eJcNxbaJ8aruzuLkdMGXXvBXhuHj/acZqb4Uf8G4//AAU/vPib4E8M+O/hl4N+Hvw21fX9H0z4h/Eez+K/w18Wt4Y8K6hqCDxDrNl4fsfF8GueINR0vSpZhDAtkDLI4CnBLD/Qy+HXgDwp8KPh/wCBvhf4D0qLQ/BHw48IeG/Ang/RIXeSLSfDHhLRrPQdB02OSVnllWy0uwij3OS7bckkmuzr8wv+Cyxtx/wTZ/aR+1P5du0fwkSSTrsEnx2+GCK//AXYGv09or8mP+CxP/BPnTv29v2UPE+ieDvAvh7xF+0t8PorbxD8Ate1LUrDw5qFrqqa5ot14m8JP4lvkNpDonjHw9YTwPbXrLp73y200rR+UJF/kc/Zo/4Nvf8AgoN8bNfSL4ueHPCf7LngJbof2j4p+IOtaL4s8WXWnrcpFND4U+HvgbxD4gkvb2MZkWLWbzRLB0XKs2QD/Rj+yh/wbVf8E/fgF9m134wWPib9q7x5ES7aj8Sbqbwz8PrWcSh1l074X+ENQtdMuwYxsdNbvdcRskgKcY/e7wV4F8EfDbw3pvgz4deDfCvgHwfo0bQ6R4U8FeHtJ8LeG9JikdpXi03QtCtLHS7GN5HLMsUSgk5NdXXwP/wUn1OPTf2atFEkRm+2/tTfsM2yr5Rmw0H7aPwF1Qt5f8WF04j8a++KK+ZP2lP2NP2W/wBsHw7D4Z/aU+B3gL4s2Npa3VnpWoeI9IEfinw/b3oP2uPwz410uTT/ABf4ZFy2Gk+wX1v5jAFgSBj+TX9sn/g1V8Y3HxOm8QfsKfErwHp/wt1qGGe68CfHfxb4qtPFHgvVvuT2/hnxTovgDxlb+IPD0rL5omv44tVgDeV5lxtMz/Wv7Fv/AAa1/s6/Dmy0fxd+2j4+1v48/EAfZrjVPhx8PtW1XwV8D7NojOsuk3V+lnpfxG8eWzAxOl1NPojKQyG2YYev6dPhp8L/AIc/BnwR4f8Ahp8JvA/hX4cfD7wrZDT/AA54O8GaJYeHvD2j2gYu0VjpmmwW9rE00zmSR9u+WRmdyzEk95Xw18UrRLj/AIKH/sbTldz6f+zJ+3ZdKd+Cnm+OP2IbIttz8wYXWPavuWivj39uX9ib4Pft+/s/6/8As/8Axlh1O20q8v7XxJ4S8WaBOLfxJ4C8caVa6hZaN4t0N5d9tcT21nqt1a3FtOrQ3ljdzwNt8wOn8jHij/g0r/aEj8Z3Vl4T/as+C+sfDg3cbWmv+KPCnjjQPGotzGvnXV54M0iDxF4Zm1SORmCsNXVXAySudo67wl/waOfERtavLPxv+2b4DsvCTwsbfUPDfwd1/X/EUs6j92k+lav448P6bbxMerreSkf3Wr+nX/gnh/wTq+B3/BN/4W6p8K/g5ZanrN/4nbQNf+I/xX8Tz6e/i34keLLG0vdPZbuw02xs7PQPDXh2DJ0jToTJDZxX0kYaSQSzzfoNX5uf8Fbknk/YN+KKWqeZcN8Rf2XfLQTfZyzD9q/4IE/vv4PlH49O9fpHRX8kX/Bwz/wR8+KH7Q/iu0/bU/ZT8B3/AI/+IH/CLW/hT48fDXwrDLqPjfxVaeGrFbfwb8QvCfh+4mc+LNU0zSLKHRr7StP8rUZIbfT57WOV4rhk/h91X4WePrPxCnhbVvB/jbT/ABHBdG1PhO/8H+ILLxKL/tpx8NNpS60GA7YzX9D/APwSf/4N9/j3+01428O/FL9rbwN4j+C37LujX+m6tJ4d8Y2ep+E/ij8XYbKRLpPDOg+Fr1ofEPhfwhrFrIItS1vUPLmngYDS93zyRf6F+kaRpPh7SdL0HQdM0/RNC0TTrLSNG0bSbO307StI0nTbaOz07TNM0+zjhtLDT7CzgSKGGJUjijQKoCgCtGvy9/4KoyTr8P8A9kGC3Uk3v/BSL9h60kxbfaiI5PjHpzsfJPX5oxz+Hev1Cooor8J/+C8X7F/wx/at+AXw98XftBf8FBNT/YU/Zz/Z/wBf1vxn8TE1TRPDPi34T/F3WdUvPBz+ANO+Jfgfxlrdj4Z8fz+FtX8PTpo2h3VhrA1O71eSFLSSQojfkz+zT+1d4c+KOufDz4BfAL/g6E17whPqsmleBvhB4G13/glT+zh+z14Q1pYkt9N0Dwv8Mpfid+z58OfBd/A7tFZ6fY6bMfMlMcMEZZ41b+wj4aeH/FnhL4cfD/wr498e3nxV8deGfBPhTw/40+KGo+H9D8J6h8SPFmjaFYad4j8e33hbwxBbeGvDV54w1i2m1GXT9OjjsbN7gwwKsSIB21FFflr/AMFsJRD/AMErP215dnmbfhIfk8kz7s+KfDq48r+LrX1v+xnF5H7H/wCylB/zx/Zs+BkX/fv4YeF0/pX0nRRRRRX5i/8ABZOze+/4JtftJWkcC3LSwfCgiBujiP45fDOVs/7qoT+Ffp1RRRRRRX55f8FOFhk/Zw8IW88csiXf7XH7DlvmHAliZ/2vvgu3nR8H51VCPoa/Q2iiiiiivh/4jzXK/wDBQ/8AZMhIT7K/7Kf7b0kbD/WGWP4lfsPLdBhnhAJINvqSa+4KKKKKKK/O3/gqtb/av2JPHdv9tbTvN+LP7Ja/bV+9B/xl58CjvHXk4x+NfolRRRRRRX5g/wDBUMxjw7+w2kqI6Sf8FO/2HY8ONyhj8T2MZIPBPmKMe+K/T6iiiv57/wDgrv8AEv8AZu8Ef8FA/wDgjjZ/tl+Ovhn4O/Zh0vxD+278YvEMfxp1jR9K+Fd58Xfhb8J/hjoXwU1fxHH4jmTw/e6r4Y1f4n6jcaQLhXeLUHR48MCG0P8AgoL/AMFGP+CLP7VH7F/7SXwU8eftp/sU/FG08U/Bv4iJ4Y8MXHxh+HurapH47s/COr3PgfVvCKpqr32n+MtJ8TR202l3NoUuobtUMbA8H9L/APgm7448S/E3/gnd+wT8SfGepTax4w+IX7F37LfjjxXq9wzvcar4l8WfA7wLr2u6lO8jySPNfapfyysWZiWY5JNfaNFFfmz/AMFhbD+0P+CX37ccPlRzfZ/2ffG2qeXLko39i20WsHOOcj7Dke9cp+zr+yn8TvFf7PX7O3iO3/b8/bW8K2F98EPhDqdn4S8Px/si2miaRFceAfDdxFpsL6n+yPqPiS4tLRf3YF7fXUpXPmO7ZNezT/sa/FC4TY//AAUV/brT/aguP2P7d/8AvuH9kFDWRF+w38RIv+cjn/BQGX/rr4m/Zab8OP2UVGKcn7DvxETr/wAFGv8AgoA/z7/n8TfsuH/gPH7Ka/L7VJZfsQ/EbTyTa/8ABRn9vzJ/5/fEP7K2rAfRdY/ZP1Bf0q1/wxf8UP3p/wCHjn7eWZW35+1/se/uz6RD/hj7aq+2CKbY/sW/FDT23wf8FHf29JT6X95+x7qa/wDfOo/sfXQr4g/4Kc/s3/E7wZ+wp8dtX1D9uX9rnx3ptu3wwjl8N+MNN/ZFl0aWO8+NPw6tXuL278NfsteDPExg0wT/AGg7NWh4i+Yuu5G+75f2P/ixND5Lf8FGP25F+bf5sdv+xZDN/u74/wBjRPl9sVUh/Yz+KsMTxJ/wUe/bvYP/ABzN+xlPKv8AuST/ALG0jrRJ+xx8XDLaSQ/8FH/254UtyfNiNv8AsXzLcr2Enm/sckhh68/SgfsZ/FQXv20/8FHf27yf3n+imT9jT7F+86/6OP2OAuU/h54rR/4ZH+Lv2j7V/wAPGP23vMznb/Zn7D32f/wE/wCGK/sv/jlEn7JfxkktLq1H/BRr9tZfPTbFcDQf2HBc2xII3JMP2LlZmyfaubt/2Nfj5blnH/BTL9s+aQjj7R4X/YxmiB9RC/7KBH618af8FA/2aPj5o/wS+H0uq/8ABRD9pDXbC4/aw/ZC0maPxP8ADf8AZNgsrLUPE/7TXws0PQPFKS/D39n34eawmpeBvEmoW2r2KNe/YPPs4lkgJVZE+5L39lz9pqa58+z/AOCk/wC0tZxj/l3f4P8A7Ed1H+Of2XYcmsyb9lP9q2bBH/BTn9pGIj/nl8EP2IlX6lf+GZiT+dXLf9lr9qWITCX/AIKX/tIz+YSYi3wV/YiUwn2x+zEd6+3FWD+y5+03j5f+Clf7TmfVvg7+ws36L+yglWB+zD+0p5Sof+CkP7TRkH3pf+FRfsN/MMc5j/4ZS2jJ9KQfsw/tLcbv+CkP7SzYPP8AxZ/9h8Z/81XPNSyfsz/tLGHy4f8Agoz+0Ykn/PeX4PfsVSv75RP2Y4Yj+VfGnjr9m39qhf24P2c9PP8AwUU+Ljahd/s3fteX+k67c/AL9l2bWNL0rTviT+ybb61oVvHY/C/T/Dcg1xvEekTtdXelXJhbQAqENeEL9RyfsvftmnGz/gpt8XE/eAnP7Nn7I75Tun/JJVwfeq7/ALLX7a20CL/gp98WA2zBaf8AZk/ZHmy/98CL4VW2PpzUH/DLf7cmbbH/AAVA+Iu2I/6Tn9lb9lYvdD03DwIFhyf7qmp5P2XP22mOU/4KefFCP5cY/wCGYP2T3G7+/wDN8Nc59ulL/wAMxftv+ZBj/gpp48EKf65P+GWf2YTNP9Jv+EOEcX4RmqF3+zD+3xLMj2n/AAU98RWsS/fhk/ZC/Z0uvM69ZPsMLrn2pD+zP/wUD3sR/wAFNpwh+6h/Y1+BhK/8D/tUbvyr4s/4KF/Af9tjw7+yx4m1Txj/AMFBf+Ez8Pj4q/sz2Uvh5/2T/hVoPnX2r/tM/B7RfDt2uqaDrY1VW0PxNqFrqLRRKWvRa/ZgFWUmvtC6/Zs/4KH3L3bR/wDBTHTbBJmzapZfsT/CMi0X+5u1DxtqDXH1bFU4P2ZP+CjEabZv+CoNpcv/AH2/Yi+DEX/jsXi4CpW/Zn/4KKFNq/8ABTqxV/N3bz+xJ8HG/df88tv/AAmIG7/a/Ss1f2YP+ClAZyf+CpekOp+4G/YT+DoKfVl8fLv/ACFK37M//BTEuGX/AIKheD1TvGf2Cfhq+f8Agf8Awt9DVj/hnD/gpeY41/4eZ/D4OszySS/8MD+CcywtnZCV/wCF67UMf94cnuKbN+zl/wAFN3RVh/4KY/C6BxLuaQf8E/vC0gaL/nkUf9oggH/aBH0r4i/bH+D/AO2v4K1n9hPUPj1+2h4G+PXgm5/4KSfsXWcng7Q/2UtA+DF6uq2/j+TU7XWo/FVh8V/G1wfLm0t1NqLZFcXBO/5Ajf0G0UVIEBAPPIz1/wDrV+cX/BQf4TftZfE6D4Tr+y38AP8Agnb8dZNFm8cN44T9v2L4gyweGE1BPCQ8ON8Kh4D+GvxDYS6y1lfDW/tX2MYtbHy/N/eeX/OX+1N+0F+23+y38SfB/wAAn/4J6/8ABBT9oT9qj4gTWbeEf2Tv2YvBnx1+KXx8u9FupLfz/GPiPQLr4FeGfC3w08CabaXH2qfWvFes6FpxtY5JIpZRFIF/rw/Z+i8cW/wF+CMHxO8BeDvhV8SofhF8NoviF8L/AIdyWs3w/wDhv44j8G6KnizwF4FmsZrixl8HeD9fFxp2mNDI8TWVvEUZlwT67RRXwd/wVHsTqP8AwTb/AG8bZY/NKfskfH69EecZGmfDLxJqRPUfdFrn8K9w/ZMkM37K/wCzPMRgy/s/fBqQj0Mnw58NsR+Zr6Booooor8xf+CyDyx/8E4P2iXhuJLWVZfg2Unih8+WM/wDC/fhaMrFn5yQce2c1+nVFFFFFFfnj/wAFL7uCx+CHwjurlITBD+2z+wtJNNOcR2ca/tW/CpmvD7w4/Wv0Ooooooor4s8YXkj/APBRH9nrTvPuEjg/Yx/bE1AwLj7NK03xt/YdtdzjBJkjEWR9a+06KKKKKK/PX/gqRKI/2OPEIIkYy/HP9jWFREMuTJ+2X8Albb6HZmv0Koooooor8yP+Cmqs1n+wSF2f8pN/2QGbf3VNc8TudvP38LxX6b0UVOOg+gr8cP8AgsVb/wDBW3xF8OPhb8Pv+CV+geEPM8a6r4stv2i/iBc/EDwT8OPiz4N8F2Q8KL4f0j4MeKPiHb694W8NeJPGltqGtwza2+hazdaM1pbzW0cUzK9fA/7Gvgv/AIKNfsJeEtZ8P/AX/gg/8J4/FPjS9/tr4sfGnx5/wVj+Hvj748/GzxVK73F74s+LnxY8Qfs7T+KPGGr3t/PNciFpYtNtJriX7Ja2yOUr+kr4aaz438R/Dj4f+IfiZ4Ks/hr8SNd8E+FNZ+IHw607xVa+OtP8A+N9U0KwvvFfgqx8bWWnaPZeMbPwrr09xYxarDaWsWoJALhIYlkCL21FFfHH/BRWBrr/AIJ9/t1Wy/euf2OP2nYF/wB6b4J+N4x+rV6X+ynH5X7Lv7N0X/PP4B/B6P8A74+Hnh1f6V75RRRRRX5d/wDBZ77Qf+Ca/wC0attcWlrK8/wWj8++GbSOOT9oL4UpOZ/9kwMw+pFfqJRRRRRRX5w/8FN7WzvPhL+z3HfX11Y26/t8/sKO0ltL5Syn/hpbwBGYbp/+fUpIWb0ZFr9HqKKKKKK+HfFTv/w8m+A6bfk/4Ye/a1YsMY3/APC+v2KgFJ+90zjt1r7ioooooor8/v8Agp3M9v8Asha1Ojcx/Hj9jNjH/wA/C/8ADZfwCDW3PI88Hbx61+gNFFFFFFfmH/wU4aQW37AKxrI5f/gp5+yEr+Xn5Y11XxdLI79f3aqnNfp5RRU46D6Cvhz9tL/gnN+yj/wUCi+HEH7UHhDxr4sj+FEni2XwSvg/4y/GH4Rmxk8br4bTxEdRb4T+OPBba/8AaV8J2PlC/NyLXY/k+X5su/8ALH4+f8EJP2EPgH8E/iv8a/2dfiN+03+xv8VvhL8O/GvxJ8JfHjwp+2P+0TfW/grWvBXhvUvEOnat4t0n4j/Ebxl4X1vwZa3WnqdVsp7VRd2Bmh8xN+4frV/wT5+M/j79oz9hP9jn4+fFSyFj8SvjL+zL8EfiZ47RLO0023vPFXjP4deHtf1vWLHTbELaaZpuvX9899a2yBfs9vcJGVUqQPsCvnHwN+13+zJ8S/jP44/Z58BfHD4eeLPjX8OIDc+Mvh1ouv295r+kRR/Zvt/lov8Aouqy6JLeRRalHZS3Eml3Eiw3YglISoP2h/2wP2av2Tl8ISftE/Frw98LI/Hj6+nhKTX7fWp01k+F49Lm8QNC2kaXqS28elRa5atK83lKBOuCeceKftV/Fr4bfHf/AIJjftZ/Fr4P+LtC+Ivw38bfsbftM6j4U8W6FPLdaH4gsoPhP4+0yeW1n8uKZ4o9QsZoJBtV1eNhgEV9D/ssJ5X7MX7OMR6x/Ab4QIf+A/D7w8v9K95ri/iH8RPAvwk8E+JfiT8TfFugeBPAPg7TJta8UeLfFGpW2kaHommwMqPc31/dyRwx+ZNIkcSZMk00iRxq0jqp/n38D/8ABz1/wT68RfEnxV4O8UaJ8ZfBXg3TvEF3pvhP4s/8Ieninwz4m0Wyh/e+ItV8OaBdTfEDwzHdXI221sulX80qEM/lMTGPtvQ/+C5n/BJ7xDE0th+2j8Orbbw0GvaB8R/C18p99P8AE/grSL4f9+6/THwB8QfAvxW8G+HviJ8MvGPhn4geAvFunpqvhjxl4N1vTvEfhrXtOeSSL7ZpWtaVcXVhfQrPC8bmN22SoyNhlYDsK/ML/gspF5//AATi/aAiyo36p8CuXh+0qMftGfCRsmH+Pp+B5r9PaKK/Jj9u3/gtH+w5+wXp+r6Z4y+IUXxW+Lun28rxfBL4O3ek+KvGcFwn2iML4rv/AO0Lfwz4FhhuIlE0epXcWo+VIskFncDivmn/AIJw/wDBwj+y7+3x8SR8FfEvhDVP2aPi3rV0bb4e6F448Y6H4j8NfEe64aPw/wCHfF9pY6BHD41uYD5kGk3NnFLdKCtvJLLtjb9/aK/Nr/gp3Z3N78Lv2aUt7Yagsf8AwUF/YVmudLYhV1a3X9orwYjae7MQirNJIpJPHy1+ktFcv4q8beDvAthHq3jfxZ4Y8HaVLOttFqfirX9K8PWEly4JW3S81e7tLd52UEhQ24gdK/Br/gpx/wAHAP7P37ErwfDr4GReFf2l/jjex+dqVroniyOT4ZfDqzkXdBceLfFOgx6hDruu3D/LHoenXMd5jLTSwEIknyr+xJ/wdF/Az4laleeEv23fBtl+z1f+VJe6P8UPAVt4u8b/AA0uoord5J9L8Q+H7XTtb8ceGdQgli2R3FuNWs5zKPMa2WMu/wDS58FPjz8GP2kPAWn/ABQ+A3xN8G/FjwDqcsltbeJ/BOuWmtafHfQxwy3Ol6h9mkNxpOs2cVzG09ldJDdQiRfMjXcM+t18I+IXMv8AwU3+EEeJsWH7CP7Rz7v+XfOrftA/srggf9Nj/YvP+zivu6iq11d2tjbz3l7cwWdnawvcXN3dTJb21tBEpeWe4nmZIoYokBLsxAA618vfED9uj9jD4WaTrOs/ED9qv9n3wzaaBYzajqkN58W/A8+rRWsKlmNvoFjrV1ruoTuB8kNvbyzSHhVJ4r+Tf9vb/g6W8YXniXUfAv7APhnQPDPhLTLua0l+OPxc8P8A9s+JvFcsK27C68E/Dq6u7fSfDuiJOsqmXWPt97dwOha106QMD7H/AMEnf+DkLxX8X/if4X/Z5/b2i8C6de+OtYtvDvgX9oXwxpv/AAh2jW3iW/S4/s3Qfij4ft577Q7KHxBq0kdhaapYrYWmmSIo1KNVd7qL+w+vz/8A+CmccMn7LVus8BuIj+0z+xAXiFz9k3Z/bR+AgOZfx6f4V+gFFfnj/wAFBf8Agpj+zd/wTn8BW3iL4w6tfa/4+8S2GoT/AA4+DvhL7LP428b3Fn+6NyzXc0Gm+F/C8F2wW51W/kSJVSRbaO7uUFs/8gnxu/4Ol/26PGOqzj4P+Bfgb8DfDgnaTTbS48Nav8TfFzRNtC2+s+KfE3iDR/CcqRlSQyaLpj5bkkYAyfgz/wAHRH/BQLwXqUV18U/DXwK+OPhaa7jbUrK98F6l8O/E0UCBg9toXiDwjrtrpGjiXeCW1DRNUcbeMc5/r7/4Juf8FP8A9n//AIKXfDLVPF/wsTVPBvxC8FnTbb4pfBrxZc6fceLPA93qkLy6bfxXmmu9h4k8Kax5En2HUoBEZAm2eC2m/cj9JK/Lf/gp9NcRXH/BOYQLu87/AIKkfsrQz/7MB0v4nu7fgUFfqRRRU6/dH0H8q/Fn/gtTpH7Gen/CL4VfE/8AbA+NH7Yfw4j8L+LtY8CfBv4c/sVfFb4j+CfjB8fPiT8TotC+x/Dvw54E+G9zbah8T/E7x+D1fTo7l4rXTEkuXeeFJ3Y/gDo/gH/gmT4f8ReENV/4KI/sZf8ABej4Z/s+XPizQpLDx/8A8FF/Hnxg8efsj6LqtxqFvB4Yb472Pgb4veJT4Kt9R1x4YjBrtjPo7Di/kFqZQf7h9BttDs9D0a08MQaTa+GrXSdOtvD1roMVnBodtocFnDFpMGjQ6eq2EOkw2CxrbrABCsIUINuK1CCQQCVJBAYYJGe43BlyD6giv54/2XvhV8YbMf8ABOn9lGX9mX4pfCn4j/sI/F3xf46+OP7R194KstJ+CnibwjZfDj4reAvEifCr4qTSw3fxLb9pfVviLo97dadZot1bRwTSapGsmnAH9J/2zNQ+Jvxg/Ze8Q+A/g74V+KGj6p8XPjV4G/Zs8cajH4Uk0rxl4M+Dnin9oTRvg78evito+i6/AGu/DS/CU61qmlaxD8g0q7t9XiZY4zjU/bC8D+Dvhd/wTe/as+Hvw/8ADek+D/Avgb9iv49eGfCnhfQbSPT9H0Dw/onwR8W2OmaVp1pFhILWys4VRec8ZYkkmvd/2acf8M5fADacj/hSfwqwfUf8ILoOD+Ir2yv86j/g4++MP7Zdh+3D4s+Bvx5+I2oat+z/AKfp+l/Ev9nnwXodkfDPgQ+APEIltrXVNT0jT8XPjfxl4b8QaXf6fd6hqRvUW7tQYUhgZI1/nSn1ebTZZL61jmtfON6v1s+e9a2la19tlsZbpPJi88ieeG3G2ADJHAPQ1/SP/wAER/8AgtDo37B1/qfwA/aB/tS4/ZT8b+Kn1vTfFFlY/bdV+CXjPVjFY6pr0llawPqet+D/ABI8Vt/adjEXe1e3bUNPikLX0M3+hFFLFPFHNDJHLDLGksUsTq8UkTqGSSORSVeN1OQRwRX5i/8ABZcRH/gm9+0KLiSSKE3/AMEBJJE22RUP7Q/wmB2t2J6fQ1+n9Ff5ov8AwWB+Pn/BTHwr+0n8Xf2bv2vv2i/iJ4m0Xw9rZ1Lwr4c8K31l4D+FHi34bavfz3nw+8R2fgf4ZjRdEvkvtMkSOV9WafWI5zLZSOX35/BrWL+SUxuieXF5/wDqeM59+TjrS6Zfpp80V5HPNa6lBc2d/Y3Fufst1b3mnci/z2xiv7uf+CGX/Bdn4tftNfFT4c/sPftTaRb+MvH+v6B4jsfhv8e7KWGw8Q+KLzwF4U17xrfaR8UtHRV0vU9dn8JaFOYtUsfs108lqv260e4nkuR/XLX58f8ABRyOGb4a/s5QTGULP+3/APsERhYjguf+GpPhsxVj2QIpJ+lfoPRX+Tf/AMFM/wBuT4hftx/tWfFf4ueNNf1m/wDC8XivxL4d+Engi+1O9k0L4e/DLR9S1PTPDOhaHp6LFDp13daZYQajrMiIn9q6vI8m0bsD8/bTzI7eOec+X/z7wzY6EcYNZEd9I7mOF/3f2gmaGE+/bGexr+jz/g2j+MnxB+HP/BS7wR8JvCurXMvgT49eAvifoHxK8OvPeNYtZeBPAWv/ABD8K+MLjTpj5cfiKx8TeGorFL0ruMWrXsY/1pr/AEba+BNevY5f+Co/wp07DedYfsC/tAXoP8Pl6t+0R+zRAwz/AHg+i1990V/l7f8ABZ79rb9oL4uft2/tXeGviT498Val4G+GHx8+Inwu8DfDVdZ8Q2XgLw94U+Huval4P8MyWnhj+110V9Q1PRtIXU55QqtqWqXlwcAfKPxkfXEd5JEtIbfJ/wCeHI6c/wAOTms43CO/n3aSy/8APHk/6RxzkjJHWtvTbg2p/cSHzMd8f8feeOpJGK/2DP2OtR8eax+yN+yxq3xTbV2+J2qfs5fBDUviM2vo0WvN47vvhn4ZufFza3G3zJq58QS3BuQeRPurxT/gpWkM37LkNtKtwHuf2lv2ILe2ltV3XFvdSftpfAHyp4RkfMoBH41980V/kwf8FJ/2rfiD+1L+2R8efiz45u7t7/VfiL4l8L+EdGvVjB8G+AvCniPVvDfhjwjbGNIBZW2k6Vp1u88xXfPcO8r5d2J/OnUL6QSkpdzeXEc85P06k9KswajfSW/lQXBcczgZ7/rxxX9HX/Br/ffES4/4KZaSnhCK7Xwyfgj8Vf8AhbbWk121oPBIs9FOgHUhKNj5+JA0PbnpJtr/AEcq/Lr/AIKdrG1z/wAE6PMiWXH/AAVE/ZYZA38Ei6V8Tykq/wC1GTkV+otFFTr90fQfyr8H/wDgqL478F/sz/t/f8Eqv21/2hrV7L9k74RP+1v8HfHXxXvtHbUPCH7O3xe/aH8JfCzR/g98TPHF5bJdXOheH/FA8Iav4aOsSwrZ6RNqKGWaIXWT7B+3x/wUz/4JveEP2NvjldeMf2lf2d/jNpHxG+EfjvwT4U+EPw8+J/gT4reMPjnrXjXwzqvhzQ/h74I8E+Ddb17WfFN34u1O/jsg8MDWtsJTLcSwwo8i/Tn/AATR+GnxP+DX/BPH9h/4TfGmG7svit8N/wBlT4EeC/Hmk34i/tDw94i8O/Dbw7pd94V1F4JriG41Hwm1sNNuJlkdZ5rVpAzbsn7door5D/4KBnb+wV+24w42/siftJn7hk6fBrxof9X/AMtP93v0r0n9mElv2av2eGYYY/A34Skj0J8A+HyR+Br3Ov5cf+Do/wDY0u/jH+yr8P8A9qzwlppu/FP7L/iG50zxvDZWaPe6j8HfidfaNpOq3s88TpeXaeDvFtjp00cAWRY7XUb6b5Arbv4ALawF3d/Z/wB9c+cfs9vFDuItx3z9P6V1umeA/G9y8slt4M8UzWvnm4nnt/CWvt9ns/X5dJ4r+pn/AIN9/wDgkhpH7SHi64/bH/aM8IDUPgL8NfE02nfBnwF4hsdlt8S/iR4a1TT7qXW/FOlara3Ca/4B8ByxNDLC4WDV9eAEjMtjeW8/949fmN/wWPtnvf8AgnP8fLSNtklxrHwHiRvRn/aO+EYHpX6c0V+an7c3/BJr9jP/AIKGa7pHjL9oLwn4xX4h+HvCKeBdD+IXgHx94g8JeINP8L2+r6pr9jpcmnLNqHg3WIdO1nXLy5g/tDSrzZJcvnI2gf54f/BU3/gmL8Tv+CZHx0XwH451S58c/C7x1YXGtfBz4v2uj3Wl6b4y0azeKLX/AA/eabBPcw6N438NPMjalYjUWDwTQzqWiniY/K/7Kn7HX7QX7anxc0L4N/s5eAdQ+IHinVIP7Q1O5jt2tdC8FeHv7Sj0i+8Y+OPEIDaX4b8L6XczxxR3V0pkmmkSKJHkeONv9FD/AIJb/wDBDf8AZ2/4Jyz6T8VdU1jUvjf+1IuhT6XdfFbxBG2neHPBKatYy6ZrmnfCbwXDLJaeHo7/AEaRdPudUu3u9VvLWN1WS0gup7Rv2/r8/wD/AIKIf8iJ+zRzGP8AjYB+wsf3i7s/8ZIeBeE/uyHse1foBRX+YV+1F/wQo/4KU/AXxb4jsbT9mjxn8Z/DkGpSJoPxC+A8UHxL0/xVpKuQb/8A4RLw/InjXRiRjFlqekxOTyu5ea/I7x54B8Y/Dq4g8P8AxC8IePfBGu2NxeW1zovjvwr4h8IajZuGOdPaPxRpKuhQjnODX0F+w5+xv8UP21f2jvhp8CPgzahtd8Xa0txrniaezXUdC+HfhLTXRPEfj/xJuK+bYaTp8iTxRKfMuJ5EijDSOoP+jT/wTD/4I3fs9f8ABM2DXvFvhrX9e+L/AMe/Geijw54u+MfivT7XRVh8Opq8+qp4b8BeDrW81WHwVoV3ILZ75ZL/AFO9vrm2EslzsKwp+vdfn3rc3m/8FVPhpBj/AI8f+Cffxvm3fau+rftG/s9pt+x9R/yBv9b3+7X6CUV/En/wW5/4INftMfE/9or4hftU/sfeFLP4ueGvjDfw+LPH3wxsdZ0/R/Hvg7x60EFrrupaPb+JPEOk2HjHwz4v1BRfmO2uY7yxup7gC38mNJG/Cnw5/wAEI/8Agq94stLk6T+xX8R7e7DCa4fxX4n+FfhMZGMjTh4x+IXh0yPn0ya42X/gi5/wVI0PX49D1T9hf47Xd3CPJLaZpGj65oRBHDr4q0TxJLoDLkdnNf0Pf8Elf+Da/wAZeHfiB4c/aN/4KG2Gj6ZZ+FdT0/xD4H/Zjt7/AETxhPrWu2GoPqFjrHxh1vTTd+Fk8PWhRHTw/Yid9TMudSe28mSzuP7Vq/P7/gpXp39r/s8eCtL+1Q2f9oftgfsG2puLg4RRL+2d8DAeecNxn6Cv0Bor+Hf/AILm/wDBB34rr8Q/H37ZP7GHhO7+IXg7xjqmqeOPi/8ABbwzb+f458I+J9XuFuNf8XeA/DtoqS+MvCF7cSSXlxpNgBqun3EjPbRyWwbyP5F/Cvwd+Jvj/wAaab8OPBXw+8YeM/iJq+pWeh6V4H8MeFPEOteK77XWOE02LwzDpL648m44Chc1/WT+wz/waj+MfHXhnT/H/wC3X8V9b+Dk+sfZb2z+Cnwfk8M61470rTnt5F2eMfH+p2Gs+CNC15SVJstN0zWYEQ/PcLJuiX+pL/gn5/wS5/ZR/wCCbHhnxVo/7PegeIr7xN46lsT40+J3xD1XT/EPxC8Q2Ol+Y+maE+o6Xovh7R9I8N6fdXE08djp9jaQNPKXkEjLGU/RWvy3/wCCncxi1L/gm2BB53nf8FSf2ZYf9b5Xk58J/GJzP/002BPu/wAWa/Uiiipx0H0H8q/Kv/gpf+0z8ZPAOq/s3/se/s1fAv4IfHv9oH9uLWfi14c8P+Hf2nNW1PTf2dPDnwy+DvgW38Y/FnxZ8V7PQNH8QeIvE2kQWutaZYQaVaWjNdy3/LZRYZvwz/ZX0/8AaT1j9ov4w/D/APYs/wCCH3/BIv8AZi/a4/Y01fw54Y/aW+OHxCuF0DwBp/xN+Iekan4o8FR/sy3Xwk+EVz8TLLwh4h+EtzpOuzzyHTnT+1jZyRkxJcXH9a3wij+KkPwn+GEXx0n8DXXxui+HngqP4x3XwxXW1+Gtz8VE8N6avxCn+Hq+JYbfxGvgabxaLxtJGoRpeiwMXnqJdwHodFFfIH/BQptn7Av7cL5Ubf2QP2l2y/3Rj4L+NTlv9kd69J/Zbl+0fsyfs6T9fO+BPwil/wC/nw/8PP35/ir3aqt5Z2mo2l1p+oWttfaffW1xZ31jewRXNneWdzG8Fza3VtOjw3FtcQuySI6lGViCCDXlvw9+AHwI+Ed1cX3wo+Cnwk+GN7doYrq7+Hvw38HeC7q5jb70dxceG9G02WZG7hmINeuUUV+Y/wDwWMj83/gnd8b4/tH2Td4o/Z2/0j/nnj9pz4NNnqOuMfjX6cUUV5t8T/g58IfjboUHhb4z/Cv4b/Fzwxa3yapbeHfif4H8MePdCt9Siikgj1CDSPFel6tp8V9HDM6LMsYkCswBwTXLfBv9mT9nP9neXxLN8A/gN8H/AIKyeMpNOl8WN8Kvhz4S8ADxE+kJcx6UdYXwtpOlpfjTkvJvJEgYRmV9uNzZ9yor86v+CkBm/wCEQ/ZKESFt3/BRX9h3zWHl/JCvxx0B3c+ZyR8u35fm54r9FaKKxfEPhvw94t0i98P+K9B0XxPoOpR+TqOh+IdLsda0i/h3BvKvdN1KC5srqPIztkRhmuT+HPwf+Evwd0/UNJ+Efwu+HXws0rVbtb/VNM+HPgjw14I0/Ur5EKLe39l4Z0zTLa8u1RiolkVnAOM16NRX5uXX2o/8FfNByP8AQl/4JueLNrf9PT/tO+C/MH/fqNa/SOiiiiiivgX/AIKKp5nwl+CMeM+Z+3V+wIn5/te/B/8Awr76ooqCO3t4ZJ5ooIYpbqRZLmWOJEkuJI0WFJJ3QBpnSKNUBbJCgDoKnoor8vf+CmXmHXP+Ca/lrv8A+NoX7O/mf7MY+HfxyLP/AMBIBr9QqKKnHQfQfyr8aP8AgpR8Ev25tY/aw/4J/ftY/sRfCT4P/GrXf2X9A/bA8OeNvBXxg+LNz8JdLntv2gvDHwb8OaDeadrFn4d8SXV5LZL4L1CSRFhAUrGCcPkfC/wV8Pf8F6/gx+0V+2b+0LZfsIfsU67qP7Y3jr4OeONa8M3X7Zeu2Vl4HuPg/wDAjwP8C7PT9L1GH4XTz6zDrdh4Jj1CV5Y4WhmuGiUMqhj/AEf/AAv1H4gav8NPh5q3xZ8N6F4N+KmqeBvCWo/Evwh4X1yXxP4Z8K/EC90DT7nxl4b8O+JZ7TT5vEOhaH4jlubW0vnt4Gu7eJJTGhfaO5oor42/4KLb/wDh3z+3ZsG5/wDhjf8Aae2A9C3/AApPxvtB+pr079lSJ4f2Xv2boZVKSQ/AT4PRSI3VXj+Hvh1WU+6sK97ooooor80v+Cvm7/h398YNgRj/AMJt+zVuEieYhi/4ak+C3nbk/iHk7q/S2iiiiiivzk/4KRxWsvhn9jv7UzBo/wDgo3+xRNa7bf7Rm6T4u2Plgj/lipQtmT+DrX6N0UUUUUV+alzdFv8AgsNo9j5s22L/AIJq+I7ryMf6Pm4/ai8Kw+aG4/fYtcEf3a/Suiiiiiivgj/gofMIfhl8AQf+W/7e/wDwT/hH1P7Xfwjk/wDadfe9FFFFFFfmH/wUoeVPE/8AwTQEXWT/AIKefAtH/wCuX/Cof2hGk6/7K1+nlFFTjoPoK+Uv2oPid+1b8K4vBPiD9m79l7w3+1PoTy+IIfij4OX466D8FvippMKf2HJ4W1L4YQeO/DFx8M/HDzL/AGmmoWWs+I/C5iZbRoLiYPMsfy74e/4K+fsjadrem+DP2mU+Ln7A/wAQtWu/sOn+Fv24/hjrPwK0DVbsEoE8N/HK4m179mzxkss6tHGdI8Z3zSOuAuSM/pvpGr6T4h0nS9f0DVNO1zQtc06y1fRda0i9ttS0nV9J1K2ivdO1TS9RspZ7PUNO1CznSWCeJ3ilidWVipBrRoor48/4KHc/sAfty9ef2PP2menX/ki3jbp718+/s3fsp/Ge6/Z/+BWoy/8ABQ39tPSmvPhB8Mb5tEtND/YrbT9Le58HeH7ltPtW8R/sfeKPEhtrTb5QW81S+l27t8rklj7XP+yd8ZLjfv8A+Cjn7bSeZ/zw8PfsHW+3v8nk/sPoFzULfsifFp4/Lb/gov8AtwH/AKaLY/sSxyf99RfsYR1Wl/Y4+KU33/8Agov+3UP+uUn7HMH6wfsexGpof2P/AIrwvvX/AIKL/tyOfSaD9i2df++Zf2NXFFz+yD8Wblg3/Dxf9uODHa2tf2KYVP1A/YyOart+xx8WWZGH/BR/9uxNn8Kw/sUbW/3w37GDbq+Bf+CoP7MPxJ8LfsV/ErxBrf7ff7Yfi3TNL8c/AKf/AIR/xXp/7I40Ka41D9on4WaXaXFz/wAIR+yr4G8VXU+iXeopfafHFqsIXULW3ZiyqUb78b9j/wCMDxXEf/Dx/wDbgH2xNlzImmfsRiQZHLWL/wDDGJOnNnkGPkVDD+xr8WoRCP8Ah5H+3fJ5Mfl5mi/YnkM3/TScn9i4GST34qT/AIY5+LQ37f8Ago/+3UN8XlH9x+xO+3Of3ieb+xfJtl560sn7HXxbljCN/wAFIP26gR/y0jtv2JIZD9TD+xagP5U9P2PfiygOf+Cjf7c7krtJe3/YpJ/3hj9jFQG98Utx+yD8Xp5Z5j/wUg/blgE3WK3039hqKKL/AK4IP2Jz5f4Gp4/2Rfi7HOk5/wCCjP7b8nlpsEMmm/sQNA3Odzxj9iwb39ya+If29P2dfifoPh39mCXV/wBvL9r/AMTJq/7d37IuiWlpqvh79kKCystW1n4o6XYQeIVPgn9lDwftv9HlAvrKPUjqGiR6hDFv0+TKgfdq/ssfGYSyS/8ADxb9tBvM/wCWbeFv2CTFHn/nmg/YeG38Saq3n7KPxtvJ/Mb/AIKQ/tp20QGBBZeEf2A7UfUun7DZLdfSnWv7KPxmtN+3/go5+2xNv/5+vDX7BFzt/wBzzf2HG20P+yp8bVeae1/4KO/toRXMi7U+0eEP2B7y0jJI5FlJ+w+kRH+6VPvUkP7Lnx2+U3n/AAUb/a/nIi8s/Z/Af7BNluY9ZML+xTLhvzqrJ+yz+0Az3Jj/AOCkn7XkYmh8uIH4c/sDuLdunmKf+GLlDP8AgD71HJ+y3+0cJA1v/wAFKf2rQm3HlXnwp/YKuBu/vhrX9jrTmH0JIr4Ck/Z4+N4/4KjxeHx+3/8AtF/8JnL/AME+9Q1OPxofhn+yMdVt9Cj/AGktLhGhLo4/Z6HgSTTW1CQuZW0T+1sfJ9v8sbK+5Lj9lb9rN/IFt/wU7/aLgVP9eZfgR+w9cST/AEdf2Z7dIvwQ0+f9lX9qxopBb/8ABTr9pWK5JzHJP8Dv2FriBB6Nbp+y1bSvz/01FQWv7Kf7WiL/AKZ/wVB/aRuJOeYPgP8AsMWsXP8A0zf9mC6k4/36sH9lb9qv91t/4Kc/tLDZ/rM/A/8AYWPnfXH7LI2fhio7j9lb9rFxi2/4Kd/tHQn+9L8Cf2HLj9F/ZkgqhL+yj+2GYwsH/BUj9oNJAuDJP+zz+xBNub+9sj/ZztgD7dKyl/ZP/blC4b/gqp8XXfZt3f8ADLP7HYG//nptHwm/TOK+OP24v2ZP2xdJ+G/wXu9Z/wCCknxD8Xqn7Z/7F1jYWfiP9nL9mHT7W11vWv2nPhzoPhLxIkHhT4f+GbnVdS8H+I9Ys9W+ySziyvv7O8t4lVmNfXC/snft+tdxyXH/AAVb+IQs0fdJa2P7If7J9vPMn/PM3d74K1JE+qxA1Lc/smft5uuLT/gq58VIG87fuuP2TP2QLrMP/PHCfDO1G7/a6+1B/ZP/AG+Aqgf8FWPiMzDzNzS/sjfsmsH3f6v5Y/AsO3y/rz3rPl/ZN/4KHmSZ4v8AgrF4wjD/AOpif9jH9lyWOH6ldEilk/FxUo/ZW/4KKYtM/wDBVXVyYf8Aj7P/AAxV+zoPtn0G8/Zv+A5rPP7KP/BSoRssf/BWKYyH7ss/7C37P0hX6pDrdsrfpVF/2U/+Cohbcn/BWbQF+bOxv+CfnwbePb/d2j4mJJ+O/NfH37THwR/bR+H3xk/4Jsaj+0P+3ToX7Rfg2+/4KK/Cmx0/wXZfst/D74IXFn4jh+Efxx1uy10eKPD/AIw8UX92LfSdE1CzNmI4knGoZLbo0r+hOiip1+6PoP5V8gftL/t3/sw/sd+PPgR4G/aW+Itv8Iof2jdV8YeHPhx8QfF1heaf8KI/F/g+PwvcN4Q8afEd4z4Z8Aa54ns/E/naOdYms7O/Gn3iLOssSRy/Tet6H4Z8beHr7QfEej6F4u8KeIrA22paNren6fr/AIe13S7tFc299p1/Dd6bqdhcxkHZIjxuMHBFL4Z8M+G/BXhvw94N8G+HtD8JeEPCWh6T4Z8KeFPDOk2Gg+G/DPhvQbC30rQ/D3h7Q9Kt7TS9F0PRdLtIra0tLaKOC2gjSONFRQBuUUV8gf8ABQobv2Bf24R6/sgftLD8/gx41Feqfs14/wCGdPgFjp/wpT4V4+n/AAguhYr2uiiiiivzO/4K/wAUdx/wT/8AizFMnmW8nxB/ZdFwmM7rf/hrH4HGf/yCGr9MaKKKKKK/Oz/govaLeaJ+xuGTf9n/AOCif7Ht2B6Nb/EKRw//AAA8/hX6J0UUUUUV+aH2pG/4LGGy/wCWkH/BNAXX/ALr9qRov52dfpfRRRRRRX56/wDBR25S3+Hn7M6v5f8ApP8AwUK/YDt18zu7ftR/DmQbP+mmIuPev0Koooooor8yv+CjdrNd+N/+CZMUSFlj/wCCmPwou5pMgeWln8Av2mLjGMc7ymK/TWiipx0H0H8q/M7/AIKbftWfs0/A74ZeGvgd8dP2f/Fv7Y/jP9qv/hKfB3wf/Y18DfCdfi34h+PmoeE4dCvPE0V1pur2M3gnwv4R8F/8JFpt5q2uazc20Gk28q3MQlljCV+Uf7GH/BLP/grn8M/hxrcPgn/go9ff8E1vhj4m8VXniP4XfsIeDPhj8Of+Cgvhz9mTwTe21t/Z3w4svj5+042p67NPp8oka503QAvhy1uWd7OSVZSV/pN+GHh/xl4T+Gvw88K/EX4gT/Fj4g+GvA3hLw/47+Kdz4Z0PwVc/Evxlo2gafp3if4gXHg3wwkfhvwlP4y1u2n1F9M09VsbBrkwQARRpXc0UV8e/wDBQ1in7AX7cjgkFf2Pf2mWBHUEfBbxqQR7jFe1/AaD7L8Dfgxak7vs3wo+HcG7+95PhDR48/jtr1eiiiiivzX/AOCugLfsFfFFQ0ql/iR+yvHuim+zyAv+1r8DU+Wb/lnndz6jiv0ooooooor8+/8AgoTGZNE/ZCXsP+Cgn7IsjfSL4gtIP/H1FfoJRRRRRRX5o2qb/wDgsXrknln9x/wTR8LR+bsOB9r/AGo/GLeWX6Zb7HkD2NfpdRRRRRRX51/8FI3KeBf2VUEvkmX/AIKKfsGRj/poF/aM8GTNF/20SIiv0Uoooooor81f+CiEZk8c/wDBM0AyDZ/wUq+GMh8vuE/Z2/ahYh/+mZ71+lVFFTjoPoK/BL/gql+2+P2IP+CgX/BKrxf4pu/jZqHwf8XeDv8AgoJpHxM8AfA/wR4t+Juv+N5tP8E/s9zeCG1L4f8Ag23vNV8QWHhnxDe/bA5ieOybM5xtyNT/AIiDP2M/+iHf8FDv/EC/2gv/AJmK/Zj4XfELRvi38M/h38VvDlj4i0zw98TfAvhH4haDpvi/QdQ8K+LNO0bxpoGn+JNLsfFHhfV4oNV8N+IrSx1JI72wukS4s7lXhkUOjCu6r5b+F/7an7Lvxp+LXjP4G/DD4xeHvF3xQ8BNqi+IPDllZeILWCY6DqC6V4i/4RjxHqWj2Phbxz/wjepsLfUjod7qP9nTMEufKYgHsP2hP2lPgp+yt4EtfiX8e/GyeA/Bd94n0XwbY6t/YHinxPNeeJ/EX2r+xtGtNH8HaH4g1y5ur77FKV2WzIAhLFRXzX+1z8V/A3xr/wCCX37YPxW+GWsXeveCPFn7GX7Ud/4d1e60DxF4ZvLuCy+EvxC0u4aTQPFmjaL4i02aK/sJY/Lu7OGQlNwUqQT9kfCaFLb4V/DS3j+5b/D/AMGwp/uReHNNRf0WvQK8T/aD/aK+C/7K/wALte+M3x78faN8Ovh54e8qK71nVmuJri/1K5WVrDQPDujWEN3rPiXxJqpgcWunWFvcXlwUbZGQrEfwhf8ABRb/AIOIP2l/2n/FU3hD9lLX/Gv7LXwD0TVJ7nS9e8KavdaF8cPHqWckcmna34u8XaLJMPAlnmEsuhafLNbfvXW/uL+PYsfEfsvf8HJv7f3wZ1vwzZfFDxX4U/aZ+G2iDTPD+reHPiToen+GvG+oacmxDqVn8VPB2iHWm1941Ia+1iz1uNyd0kTsd1f0V/s5f8HOP/BPr4rWn2b44R/EP9l7xFG5SQ+KfDmsfEnwPPI5UW8Nh4t+G+i6rrMc0uSW/tLRNMjTH3z2/ob8O+IdC8X6Bofizwvq9hr/AIZ8T6NpniHw7r2k3UV9pWt6FrVlBqWk6tpl7A0kF5p2pWFzHNDKhKSRuGBINfnR/wAFfnt0/YF+Jv2tbx7d/in+yVHKtgM3ZWX9sH4DRjyf9rc4z7Zr9MaK/K7/AIKtf8FRfAH/AATJ+DGg+ML7w3a/Ez4u/EbWp/Dvwu+Ff/CSW3hs6g1rYXN3rHjTxHqUltfz2HgvwwVgjuGihae7uruC2iKeY88P8pPiL/g6z/b8mmjg0r4Q/sj+HxcTFIR/whvxb1+4ReMCS+k+M1nYOee0C181ax/wcq/8FXtQ1GaC3+K3wj8OwTMIz/ZPwU8FyfYXJP8Ax6NrkfiI3BA/56GQV/Q5/wAEJP8AgtN8Qf20PEfjH9nf9sHxl4MuPja0Fh4j+CviHS/DFj4RuPiXpXl+JLzxn4X1CDw63/CGx+JfCWn6bZ3dnBBHZ3l9YyXUphYWrNX9QFfBH7fSq2lfsj7v4f29/wBlZl/3h4yu8V970V+cn/BQb/gqN+y//wAE3NC8F3vx1u/F3iDxb8RptVXwV8NPhrpWk67431TT9EtxJqviO9h1vXvDei6H4WsLueC2a7u72Npbi4CW8c5jn8r8mNE/4Osf2Hp5ZR4n+Bf7UGh220m3m0rRvhX4inmYdpLGT4n6JPCvv859q+ffjp/wdqfCzTnksv2bP2VvG3isRXcUMviL44+LdH+H1okQGZ2tfDfgmL4hzaid/CN/acHAyV54/UP/AII8/wDBZXw7/wAFP7D4m+FfFfw90X4M/Gj4ZDR9W/4RTTvGSeJNJ8feC9bl1W0Hinwkt9YaZq8TeH9S0k2+rWjJdLZG7tG+0P55SP8AcGvzYsY3k/4LBeKZvMHl2X/BNjwFG0X8W/U/2oPiQyyfTbpJFfpPRX4+/wDBUX/gst+z3/wTGs/D/hvxLoOsfGD45eMNKn17w58IPCms2GhvY6DDcraQ+IPHvim6tNXXwXo2r3gkt9OYWF9c39xC6RQFUd1/Hz4A/wDB2v8ACnxF40j0P9pX9lvxH8LPB2oThLbx18LfHy/FS40PfEg2+JPB+reEvAt09nBOGMtzYX11chD8tkWHzf1TfBb4/fBL9ozwdp/xB+BPxU8DfFfwhqVhpeox6x4J8Q6fra2cWs2S39hba3ZWszaj4d1Z7Zv3tjqENte27q0c0SSI6j1+ivze/wCClsNnN4N/ZCN4rsYf+CkP7CU1ps7Xg+PPh1Imb/YCyNmv0hor8xf2+/8Agrj+xv8A8E6YrPRvjZ4v1jxB8UtZ03+1tB+Dfw40lfEfjy802QSrDrGrCe607QPCuiPOgHnajewTTIxe2huAjgfz/wDi3/g778M2ty8Xgj9hDXdXtDKEg1DxX+0Ppnh24fpkvo2lfCDxNKMj/p44r1X4N/8AB3D+zh4k1W20746fsrfFn4U2lw0iHWfAPjfwt8XbSEbf3d1dWer6V8J7+Gzkk+8VWV0HQOeD/SP+yx+2V+zR+2p4Ek+Iv7NPxZ8OfEzQbOSC112009rrTvE/hPULhJHj03xb4S1m3sPEXhy7kMMgi+1W0cVz5TtA8qLur6er81P+Chdz5PxG/wCCY8P2cT/aP+Cj3gP5ym/ydn7N37T48wH+AjzM59q/Suiipx0H0H8q/F3/AIKX/EP4Tfs6ftbf8E9P2ufiV4V/az8fah8END/bF8OeFvA/7MP7Mfi79oO11ofGHwl8HvDGuXfxC1nwleJcfDsaNFplvNpIlt7gaw5ulUxi1Ytxf/D/AG/Zx/6M7/4Kn/8Aivn42f8AyJX7K/C/x/pvxX+Gfw7+KWjaP4q8O6R8SvAvhLx/pXh/x14evvCXjfQtN8Y6Bp/iKx0fxj4U1MLqXhnxVplrqKwahp9wBPZXaSQyfMhrsrqJrm2ubdJ5rVp4JoVubcqLi3aWNo1ngLq6iaFjuUkEZHQ1+DP7Lfwp+OVhqf8AwT2/ZQ1P9mPx98LJP+CdnijxL4j+Lf7S2sWfha1+EnxK06H4MfFz4R6ZY/ArX9P8QTeKPHX/AA0Hq3xJg8Ta/HLbWDaG1i0OqQz3RiWvq39t3w9oP7RP7P3wq8U/Ez9kb9pz4l+Fvh7+0pYeNdf+C/gDxRofgn446Xa+BJviZ8P9I+IFpovhXxsk/jrwtf3eo22pR6RpXiLTtUudG1KG9Zh9nuLF/Jdc0b4u+H/+COn7cEXxj07xz4b+3fs+ft5698K/BvxZ8Q/8JB8Wvh/8ANa8JfFTVPgp4H+KviFvEHjSebxv4a8AXNnBdLPrGrXdlGIra5uGuYZUT9aPhz/yT3wJxt/4ozwv8vp/xJLHj8K7OvyQ/wCCmf8AwR/+Bn/BS4eGvE/jDx/8SfhX8XPAnhvUPC/gzxv4V1JfEPheLS9Rvvt72nif4YeIpW8P6sltNNcPHc6ZPomqM84E15LDFFCv8Bn/AAUK/wCCZ/7UX/BPD4o2Hg74w6DZ614N8Z/arb4dfFTwb9vuPAXxEh04BnstMt3xNoHiSHzQdR0S/H210KyQvLC8UsnhP7J37Cv7VX7cHje88LfszfB/xP8AEWfT7qys/EXiZEj0L4f+DLloLqKF/F/jbV3XRdBiube0laC2ud88/lssSO3Ff2HfsX/8Gtnwe8G6Xo3iL9uj4oaj8avFlu9pNc/Cv4SaprXgz4TW627yebpmseNLm0034p+Nba4UofOt38NtEdyhH4ev6lPht8OfBfwh+H3gr4V/DnQ4vDPgH4deGNF8GeDPD0F1f30OieGvDmnwaVo2lx3uq3d/ql2tlYWyR+bcTzTybdzuzEk/DH/BWKE3P7D/AI3t1uYrRpvjN+x3EtzOf3UJk/bK+AKB3znA5r9HaK/je/4L3f8ABGz9tT9pn9ozxN+2H8ALmH46eELr4e+FdHvPg0msHTPiJ8P7XwNoVzZX+mfDvR9TaTS/FWka9M9zqghsZk1aTUdUu4YbCdmSSb+KDU9L1bS9VurW9sJ7S7024aC50y+gNneW99p2f7R0/UNOK5Gqiq2mWGqarqMVvb2k15qd3cf6Npdlb/ar25JYHbp2n4zkkV/Sz/wRz/4Inft2eOf2kfgj8f8A4yeCvHH7LnwU+DXxC8LfFeHWvGsF14R+I/jTWPBOqWmqaX4d8GeAdVkj8YaTb67cQxCbUdZsoNOOllxCJpQsMn+hhXwX+3q/+gfsh2/meX9r/b0/ZhTH/PX7L4j1XU/L/H7Bn8K+9KK/gE/4Oj/ht8TtI/bx+GXxM8cXVi/wq+InwK0zwn8G9cZby3GhX3gHU9c1Hx14MdoVVbrWT4g8Uwal5g3AR6zaJuyhA/lx1KW6tdWMdjcY+y3HF753c4+Y5z3qtbxxyWtz5l1EZcXn7mb/AJdwcj3Iyfyr+i//AIIN/wDBJz9pj9oL9oz4C/tg65o+q/CL9nT4JfEfwx8XdF+IOvafcW938W9Y8CeJrK+tfCXw0017i2fWtC8Sajob22ra6P8AiXW9iskMXnXAEJ/0Y6/NnTRN/wAPgfGreX/o/wDw7Z+GA83/AKbf8NP/ABdPl/8AfHNfpNRX+Tl/wVB8X/FX4o/8FAP2uPEvxyj1Gw8d2/xt8faDNpOoH9/omgeDta1HQPh54f0yMs5h0XRfA+kWa6dKXYzAb2ZiSx/NK4vY/P8Akklj8r/nvF16f411Ph7WLjSdWsbm0vry12ndfT6ZfnSry/z3/tDjGa/0U/8Ag2mh/bNh/ZE+Ih/aZsfiLYfCm68eaJqH7Ltv8Vm18+LD4L1Hw79u8Uz6IvionXE+G11fXFjJo+0DT2uHvns8xMTX9H1fnD/wUnmij8KfsfiXy/3v/BR39h6KPzJvK/et8aNIdQn/AD1kOzhO9fo9RX+PR+2D4z+JnxF/aU+Ovi74yzanL8Vdc+LHxB/4T+DUPtg1LSfEtv4h1PSNQ0O6tL9ES1g8NhBpyRKgSJbQKAAMV8nut1atG7xnyges3HHXj61Kj3V7L5x/1sf7jyec9jzjNf0gf8Gyfw3/AGlfEv8AwUU8HfEv4aaf4psfgz4E8M+ONM/aO8YW0N6ng288L+IfBHi1fBvgbX7r7VBpmo6zqvxETTLvTrBMyxrpzaj5bLAxP+kFX5vft8hJvix/wTBs2eRN/wDwUQ0K9UR9HbT/ANkz9rSYLJ/sktX6Q0UVOOg+gr8hv+Cjnx//AGsH+Pn7Hf7Av7FHjvwn8EfjB+1nB8bPiF48/aW8W+BNO+Ka/AX4Hfs+aV4Nu/E+p+E/htr1xa+G/Ffj/wAd+JfHumaVpn9pGewtUFwZolaSG5t/Yv8AgovoH7c1r8MNP+NH7Cnx48GfDzx18A9M8c/EPxT8Dvif8LtG8efDT9qjQtM0ey1aL4YeJvFJu9O8b/Cu+RdDuBp2saFcRym4vTFcYhPmxfR37J/7Qnh39rP9mL9n/wDac8KaZc6H4f8Aj58Hvh78WdP8P311Fe3/AIbXx14X03xBceGtQvYI4be71Dw5d30ljcSxqsck0DMo2kV9BUUV8S/8FLLv7F/wTs/bvnxKf+MP/wBpCACD/Wg3fwg8X2qlfoZ8n2r6k+GsnmfDrwA5BVn8E+FXKt94btCsCQfcE121FeT/ABq+Bfwe/aN+H2sfCr46/Djwl8VPh5rxjfUfCvjLSLfVtNa6txItrqNn5yi50vWLESube9tZIbu3Zi0Uik5rd+G3wx+HXwc8FaF8OPhP4G8K/DjwD4Ys1svD/g/wXoWneHfD2k2yksyWel6VBbWsckz5eR9pklkYu5ZmJPd0V+cf/BVy5S0/Yu8RTvDcXgT48fsZEabbD97qh/4bK+AhOnDg4FwoP5V+jlFFflr+2p/wRu/YI/bv1ebxh8YPhNJ4Z+J11LbNqXxZ+Eepp8PvH2uQ25kPk+Jru0sr7QfFcsisqfa9T0+71CKKJI4biONdh97/AGPf+Cfn7JX7Cfg6Dwl+zh8IfD3hO7MLRa34+1G2i174n+LJJltxcy+JvHuowvr19BO9ssgsopINNgfJgt4gzA/Z9Ffnh/wUKufs6/sP/f8A3/8AwUP/AGZrf5Bn/WyeMvv+ifLya/Q+ivDvj3+zV8Af2o/B8fgH9of4Q+A/jB4Str5NUsNJ8c+H7LWho+qxqFTVtBvZkGpaDqqx5T7TZSwTmMlC20kH8/vjJ/wQx/4Jf/GD4fa14FX9lbwB8M72/sbiLR/HnwntbnwR428Lam8Rjtda0vU9LnS3vrm2kw7QX8N3aTkYlifOa+X/ANjf/g29/Ya/Zm8UQ+P/AIo3HiX9q3xppWoPc+GofinY6VpPw30KGKe1m06Sb4Z6Du0jxRqlssUsU39sXF/pc0cvy6fE6h6/oJtLS0sbS2sLC2t7OxsreG0s7O0hjtrW0tbeNYbe1tbeFUhgt4IUCIiAKqgADFWa/NvSpY/+Hvvj2Hzf3v8Aw7b+Ecvk7ukf/DT3xsUy7fduM+1fpJRX5T/8FAv+COX7Gn/BRO4bxf8AFPw1rvgX40W+l2+maZ8bfhbqkXh7xlJBplvcR6HaeLLG6tNQ8NeOdM0qaSMouo2b3sUEKw211bp0/IjwD/waQ/spaVr8uofE39qH49+PdEeJguj+GfD/AMOvh7fNcEfu5brWrjR/GwmhQ9US1iZv7461+7X7GX/BMz9iz9gfSpLb9m/4LaJ4d8T3sXla38T/ABJPdeM/iproeFYJ47zx14ikvdX07TrmNRv0/TTY6XuywtlYkn71or8zv+CmQiOl/sLGWe3hx/wUx/Y3MYuOfPkPjDVwsEP/AE8NnKe61+mNFfip+39/wQf/AGI/2+PFmofFbWtP8UfBD446s4l1/wCJ3wiudN08eObhFt47ef4i+D9XsNQ8PeJLy0hifbf2q6bq8jODLeSqiIv5Bap/waB/D26YjTv26vGdpEbiaQx3/wABND1HMEucQFrb4qaTyufvY59K+vf2av8Ag1i/YJ+EOoaV4h+Nni74q/tLa7p1xdTzaNrV7p/wz+GuoJNC0VvDd+GPBKf8JncLaSN5gSbxPNC7AB4ym5W/or+GPwq+GfwV8GaP8OvhF4B8IfDLwH4fh8jR/CHgbw9pfhjw9YKQokkh0zSba1tTc3DKGlmZTNNJlpGZiSfQK/ND9veMt8d/+CWT75ECft6XGRGDhi37Jn7T5CuccKdvNfpfRRT95Hp+v+Nfkr/wVE+B1l431D9m34/fDb9rr4TfsZftkfsz+IPiPrn7NvxD+NeqaC3wv8e6B488P6L4f+M3wY+JPg7W/EWgT+Jvh5490+x0X7feaeZdU0G6tbS8tR5oCSfm147+MH/BUD9qrwV4g/Z/+N37ff8AwRE/ZT+EvxD03UvBPxZ+Mf7LHxr+IXxQ+O954E1zTrmx8S2vwh8P/FzVvDXgvwdqfiLTLt7BdR1G7ub7S/Ne5tlMsURb+in9nf4b/C74OfAL4J/Cf4IS2dx8G/hr8KPh94H+FV3Yavb+IbS++HnhfwppWj+D9Rg8Q2ryW/iBdR0CzgnN+jOt6ZDNubfk+x0UV+fH/BWK9Fh/wTO/bsnL7BJ+y98YLEv/AHf7T8H6npuf/JuvtD4ZyNL8OPh/K5DPL4I8KSOw6Mz6DYMzfiTXb0UUUUV+dH/BU8wf8MkpFcQTSpc/tM/sPW4mg/1liz/to/ARvtwz2iVSv/A6/Reiiiiiivzi/wCCit4tpcfsFKdmbz/go7+zbZpv/vPpnxGlOz/b2QnFfo7RRRRRRX5u6Lz/AMFfPiSfK+7/AME3fgj++/veZ+07+0D+6/4D5efxr9IqKKKKKK/NT/gpO4Mf7BVmYxJ9v/4KV/spoCf+WZ0+58aa55n1H9k4/Gv0roooooor80P26k+0/tJf8EqbPzpUz+214svxFHjbL/Zv7Hn7Tlx+86HamfyY1+l9FFFfi7/wWA1n9lw/8Mt/Dv4sf8E8NA/4KX/tI/Frxp8R/Cf7Kn7P2q6Z8Obe4hfSvCWneM/jJ4rufiF8UwPCPw38IaR4e8NaW2qXsxdpZfsi+WUV5IvxT8I6n+wJ8aPE/i34Q/s2f8GvMXxY/aJ+CaWel/tU/C3xN4S/ZZ+Dnh39nzx1rN9rsfhvwNqHxW8c6mnhXx7qPiXQdFTXLO50iOS1k0W/tpt3mmaCH+uP9mrSdS0H9nP4A6FrXwd0L9nfWNF+Cnwr0nVv2f8AwvrOieI/DPwM1LTvAuhWd98HfDviHwyqeHNd0L4ZXUL6JaXmngWN1b2KSwARMle10UV+Zf8AwWZuTaf8Esv255RJ5Rb4AeLbbfnGPtptLMj/AIEJ8fjX3x8LI/J+GPw5h27fK8B+EI9v93Z4f09dv4YrvKKKKKK/PT/gp7bPe/suaVZxxee11+1N+w1B5f8Ae3/tn/Acf+hYr9C6KKKKKK/Mn/gpUrtqP/BOzYu8j/gpt+zYzcbtqL4X+LZZ/baO9fptRRRRRRX5zeH7R3/4K3/Fu+BXy7f/AIJz/s72jj+IPeftMftQzIf90ixb8q/RmiiiiiivzP8A+CkkRudR/wCCeFoBCftH/BTD9nCU+b98DS/C3xb1smD/AKbf8Sz/AL43V+mFFFFFFFfmZ+2xPexftg/8EmYrcZt5/wBqr41C7/cedwn7Fn7RXl/N1iwrP83br2r9M6KKK/Df/gp9L+0V8MP24P8AgmZ+1f8ABT9kX40ftd+F/gL4Y/bn8OfEjwp8EpfB0PiPQLj40eDfgX4d8FXtzJ408Q+HtNW3vLjQNQb5ZWYpayDGcA/Dv7Ov7X/7c/wX/aw/4KD/ALQerf8ABFn9vPVNA/bE+JXwB8ceEtD0q8+BEes+FrP4Rfs0/Dv4JaxZ+JpLn4lQ2U13q2v+Dri+tzbSSoLa4XeVk3KP6cfhf4t1jx98NPh5478ReBvEfww8QeNfA3hLxbrvw18YNpz+Lfh5rHiPQNP1jU/A3il9Iur7SX8R+Er28ewvjazzW5ubd/Ld02se5oor8u/+C1a7/wDglR+3Mp7/AAM1w/iNR0oj9RT/AId/sIeLtQ8EeBtYf/goV/wUOhkvfCHhi++yJ8VPgwbS3a40DT28hIX/AGegJI4jyPMDMzEls5IrvD+wn458vy1/4KKf8FBk/wCmg8dfs7PJ3/il/ZlkHOfSrsX7EXj6OOCP/h4d+3xJ5K7fMk8Wfs1PJN/tTt/wy+PMb34pZ/2I/Htw25/+Chn7e6ZbdiDxX+zXbrz22w/swJ8vtR/wxH49CbE/4KGft7p/t/8ACV/s1yv+c/7MEooX9iX4hKu0/wDBRH9vd/dvEf7L+781/ZXWoB+w54+Vdq/8FE/2/gfKMW8+Lv2anfBOd/739l2RfN/2sZr4d/b/AP2O/HXhj4CeFZG/b+/bq1251f8Aar/Y18PWMPiDxd8CZrO01Hxb+1X8GfDljrqw+HPgH4RvZ77wNNeLrmkx/bI4/wC19Pt3k38ivuCb9hrx9Oct/wAFFP8AgoAntD4w/ZrhH/kL9l5apX37BXjS/aVn/wCCiv8AwUNtjLL5pFj8RfgDZqhxjZGsH7NihIv9npS2v7Bvjq0zs/4KNf8ABQyXP/P149/Z0uvy8/8AZjerVt+wv48tgg/4eK/8FArjZ3uvGf7N0rN/vn/hmBS/41Tt/wBgvx1bTPOv/BRz/godIXxmK48e/s53MC4/uRT/ALMLqtJdfsFeOLqeW5P/AAUb/wCCh0EsilQLb4gfs8QwRZx80Vqv7MptFYY6mM1Vtv8Agn74xtpnn/4eN/8ABRm4LrtMd18UfgVLCvukQ/ZwVVb6V8Bft7/sj+NvBOvfsB2s37dv7cnjE+N/+CinwQ8J2o8V+L/gTdr4Yvrn4ffF7Xm8X6ENE/Z80EDXtOTwuRbx6h/aGkILiUNYvvXb+jEv7DnjyVt3/DxH9v6I+kXjH9m5V/75P7MLLWbc/sDeM7q5F0//AAUX/wCCiCELt+z2/wASfgJa2p92ht/2bYwW9+tUpv8Agnx4vmbc3/BR3/go6v7zzMRfFf4HQL3+TbD+zmg8vnpUo/4J/eMVjdE/4KNf8FGFLtnzD8UfgVJIg/up537OMiKp/wB3PvT/APh3/wCK/NSX/h4l/wAFFTs6xf8AC2Pgv5Unf50/4Z7HX2Ipbf8AYC8YQSTSf8PFv+Cik3nf8s7j4n/AyWOL/rip/ZyBT8Sa27f9iT4i2jRNB/wUW/b5xGu0pceIP2WL9Jf9qX+0P2Ubt2b3yK+HNA/ZP8fyf8FOfi14bT9vP9suG+tv2FP2f9en8RrqP7N7+Kr3T9T+PH7QOi2Xhx2k/ZvPhaLw3ojeE5bmFxpB1aS91a8ka/KyNGfvN/2Ovis+P+NjX7dK4/ur+xgM/X/jDY5qvN+xh8UJizf8PGv28IpG/jhu/wBj9Mf7sJ/ZANsP++Klsv2NfilZQiL/AIeL/t2XmM5kvZv2Oppmz6uP2Poxx7YoX9jX4pJJ5if8FGP27B/0zab9jqaP/vmf9j2X09auTfsifFub/nIx+3BF6eTYfsSR4/L9i85/HNSD9k340RktD/wUf/bZUmXzcS+Gv2C7lP8ArnsuP2H3Ij9s1PH+yv8AG7yPKn/4KPftnyyd5k8G/sAwN/wEJ+w0cfma/Pn9vH9nb4v+HNW/YJ0+T9vv9qvWZvEf/BQL4T6HpOoeKfBn7H+p3HhvW5vhv8bdZtvEekP4Z/Zl8Cq2pxW2mTWaQ6r/AGxozJc4fTZflx94/wDDIvxxuP3+pf8ABSn9taa8MXlM+m+Ff2EtHsSe0i6dD+xZNEkuf4t+aZH+yF8cLQj+z/8AgpR+2sozkjVPDP7DWsE+o3XH7GUZ5qgn7IX7RsFwk0H/AAU3/a+2J1hu/h5+w/eI/wDveZ+yUoP5Gn/8Mp/tUi38sf8ABT39p4zbs+e/wS/YQf5O6GIfsoKCffOaqf8ADKf7XYYEf8FQv2ido35V/gB+w8xO4/Jlk/Zrj/1fsBnvXL3H7Jf7eralNcWn/BV74pQac8W2HTrj9kv9kC6khm/56m+j+Gto8qY/hKD61Suv2S/+ChckOy0/4Ky+Orab/nvP+xt+yldZ/wC2S+FrZRXyT4/+C37W3wr/AG2f+CYuqfH/APbuu/2l/D+rftN/G210nwfqP7O3wn+Dq6HfT/sdfH57bVLDVfhyo1LUfslhHdaeUu3kVzqYYngCv3vooor8lP2tPiX8RfDX/BWL/gkJ8NvDvj3xloPw7+J/g3/go9dfEnwJo3ibWdM8HfEC68CfCL4M6n4IuPGnhqzvYdG8UT+D9S1W6uNLe9hnawmuZXgKNIxPmXwl+P3xm8c/Cz/guy/if4i+JL+f9nv9pH9pT4efBO8guk0rUfhr4P8ADv7F3wT8c6FpHhfUtIisL+y/sjxf4ov9Qt7gyNdRXNyzLIMLt+1/+CY/ijxN44/4Js/8E9vGnjTxFrvi/wAY+L/2Hv2TvFHizxZ4o1fUPEHibxR4m8QfAXwDq2veIvEWvatcXeqa3rut6pdy3N3d3Mstxc3ErySOzsSfuGiivzY/4LE6emqf8Euv26LV/up+zt48vj9dKsV1RR+LWYr7p+FPHwu+G3/Yg+Dv/Ud06u+ooooor8y/+CrF49h8BvgncAyPAP26/wBhUXttbW32u+vrQftO/Dp5LPT4e93Kyrj/AGQR3r9NKKKKKKK/Lz/gpg0A8Sf8EzvOm8st/wAFQ/gIsIAyZpz8KP2gdsPtuXcT7A1+odFFFFFFfnd4Qnupf+Crv7QELxMLW1/YB/ZRSKb+Fnm/aC/a/lYf72S3/fNfojRRRRRRX5j/APBR54k8a/8ABMYTbvLf/gpv8JI/lh87Ep/Z/wD2nTbZH/LNTchAX/hr9OKKKKKKK/NX9sYhv21P+CTUfl7yf2jv2ipQfL37BF+xD+0KC2f4MeZnPtX6VUUUV+Df/BazwH+yz4h1v9k3xp+0H8GP+CnPxM8X+BX+OsHwm8V/8Ez/AAz8VtW8XfD6LxXY/C6x8fxfEXXfhRfadqvh2z8X2On6fDpizTKt6lpfooKpIK/n68E6Z/wSZ8Ur8YNN+HHwE/4Ob/EaX3j/AMQ+Gvj9p/gnS/2pNXW8+Kf/AAjegWXivw/8YLXQvGNwLjx//wAIfc6VDqFprSnUf7NktEmTyGhB/tg/ZB8O+A/CH7Jn7L3hP4WeFfiF4F+GPhf9nb4KeHfhz4I+Lml32h/Ffwb4D0T4a+GdN8IeFfidomqf8TPR/iF4e8P21vaa1a3H7+31KGaOT51NfRNFFfB3/BUfT49T/wCCbf7eVtLE0yJ+yP8AH+/aNfvN/ZXwy8SaoMd/lNnmvrL4UAj4W/DUN94eAPBwP1Hh3Td36139FFFFFfmP/wAFW4Dc/Aj4HQR/a2uZf28v2DltIbL/AF1zdH9qH4ceVB1HyuAf+BAV+nFFFFFFFfl9/wAFKlum8Xf8ExhaHDj/AIKgfBNpDu2/6MvwS/aRa7H+1m2D8V+oNFFFFFFfnn4I8h/+CqX7SDJ/rof2Ev2Q4pvq/wAdP2xpV/8AHCtfoZRRRRRRX5nf8FDbb7f8TP8Agl/YxuzXZ/4KReCdRSzVsC4ttG/Zi/ap1G/uGTgkWEEHm5zwcdc1+mNFFFFFFfm1+1xNN/w3H/wSftkBML/HH9qK6nPZTb/sT/HKCI59S16R+NfpLRRRXw5+2h/wUS/Zk/YEi+HM/wC0fq/xI0qP4qyeLIvB5+HvwV+LPxhM7+Cl8Nvrw1Zfhb4P8WN4eEa+K7PyDfeQLvdJ5O/yZdn89/8AwTt/4LSfsNfs/wDiT/gozqPxTv8A9ojQbT9oT/gpL8a/2ifhXJbfsiftN60de+E/jL4UfAXwt4f8QXSaP8LL19FurzWvA2pRNY3ggvYlhDvEqSRlv6ufhf8AEbwv8Yfhp8O/i34Im1K48F/FLwL4S+I3hC41nRNX8NaxP4X8b6Bp/ibQJtW8Oa/Z6dr2galLpWpxNPZXtvBd2spaKaNJFZR3VeOeDP2h/gL8RfiB4w+FHgD40fC3xt8Tvh8boeOvh/4U8e+GPEHjDwi1hfppeor4g8PaVqd1qmlyaXqksdtdrNEjWlxKkUwSR1U3/it8cvgt8CdGsfEXxt+Lfw1+EOganeTafpesfEzxv4a8DaZqeo21nPqM2n6beeJdS0yC/wBQjsLWSXyIWeXYhO3ivkL9v34heBPiT/wS2/bU+Inw68XeG/iD4F8U/sX/ALRmo+F/GHgvXbTxJ4a160uPhP4wtbbUNH1zw9dXVnf28d4NpkhkYK8bKx+VhX2v8L2R/hp8O2i3+U3gXwk0fmjbJsOgaeU8wdn29R613VflV+2t/wAFl/2Ef2HI73RvHvxTg+JHxRtb2/0n/hTnwYk0vx144sta0+5ewvNM8WzxapZ+Fvh/dWWpBYZoNb1GxvcsTFbzFGA/jf8A20P+Dir9ub9obx/pWtfBTxvN+yR8LPB3iOx1bwt4O+H9xb6p4r1i+00zSW138TPGV7DPB42sbg3LI+kppyeGriJYxc2c8iecfsH9mf8A4Om/2kfANlpekftRfBzwP8f9GhW0s7rx34LuZPg/8RJbrUJTJHealplrpXiH4a62YYj5axWdtoKkAEtnOf6If2Zf+C8v/BOL9p/XPAPgfw/8T/Fvw7+KfxL8XaB4C8J/C74l/DnxVY+ItT8Y+Jr9dK0bQ7fxD4PsvGfw4nkvdSbyVmj1x4VPLug5r9k6/NT/AIKh2a33wq/Zoge+n04N/wAFEf2Aj9pt/wDWKT+098Po/wAcb9w/2lFfpXRXzr+0f+1r+zX+yJ4Obx5+0p8Z/A3wi8PSQ3cumjxRqw/4SDxIdPa1+32vgzwdpkeoeMPG+o2QvYmmtdHsL65jjcO0YTLV/If+3j/wdKeMfEVxrnw9/YK8GQfD/wANyJdaVJ8fPippVnq/xCkuClhL/aPgf4aTy3fhHwz9lniuYQ+vSatPdW8qs1np8ygj4q/Zt/4Oev8AgoJ8HbeHQvjJZfDH9qHQYpbGxi1DxvosPgH4gW0Om27288MXir4d2+j6FfyXnyvJPf6LqFy8igmXLPu/qd/4J+f8F3v2Mv27tV8PfDae/wBS/Z+/aB8QypZaX8I/ibeWcth4t1V2mKad8NviPYRweGvGN68KxlLC4TStamkl2xWLqpev2yr8x/8Ago+bceNv+CYv2iV4v+NnHwkEOxN/m3B/Z+/adEUTf3UYnJPbFfpxRX5W/t7/APBY/wDYn/4J2a/pPgb40+JvF3i/4p6rZ22qn4T/AAf0HTfF/jnR9FvxJ/Z2seJl1fxB4W8M+F7bVWj22UN/qVve3wYPbwyx5cfmVN/wdh/sA/2TPNa/Bb9rO41yE7P7Nfwl8JbbSGkEhGf7fb4xyOsJi53fZN27jbj5q8n+F3/B2t+z14l+IVtoXxR/ZX+JHw3+G9zcXETfEHw58Q9F+Jus2ECxSNbahqPgK38IeFJWspLgKsjWWp3skaEsqyEbD/RT+y1+3Z+yH+2ro13rX7MHx78C/Ff+zoGutV0DTLu90Txzolks1vANQ1/4d+KrLQvHmiaZLcXUccV3d6dFbzSNtjdmBFfW1fnv8O4w3/BUT9qqcMCY/wBi/wDYygK79zLv+Ln7ZMv3P4QcfnX6EUUV8AftOf8ABUn9gH9j+S/0346/tPfDfQfFmm3Vzp958O/DOo3HxF+JlpqdtF5v9nan4A+H1t4l8U6BPOzBEk1K2s7cucGQYJH49+N/+DrT9hTQ71bbwb8Fv2mfG9qJ/Jk1e50n4V+ErFwD/rrO31P4oXWtT2+McyWkJz27196fsT/8F1v2Af24vFmm/DTwd428S/CT4t67qI0nwz8NfjlpGl+EdW8Z6j5UTGy8FeIdF17xR4I13U5Z5PLh0ttTg1yVhkWAFfsfX5tft8Sqnxj/AOCXUe/a8n/BQixwn99F/ZF/a03n/gJYfnX6S0V+bHxq/wCCwP8AwTS/Z68b6h8OPit+138M9F8baQbhNY0HQYfFXxBn0i7tOLvS9Wufh14d8V2Gl63bNxJYzzR3iNw0QNfn38ef+Dnr/gml8J7m1034c6j8Wf2i9RuFjmml+HPgG68K6HpsTFvNF/qPxdufh9qLzxjBAtbG6jbP3xXz34I/4OzP2Odc122s/G37PH7Qfgzw3PPJbS+IdPl8A+L57SUIWhln0aHxJozSWbvgPJFcSsg5CP0r+l74LfGj4X/tE/C3wb8aPgx4x0vx78MviBpS614U8VaR9pW11KzFxNZ3CSW17Baahp2oaff2s1tdWtzDDdWl1DJFLGkiMo9Rr81/2sGB/b1/4JUReVvI+Jf7WtwJd2PK8r9kn4gQH5f4/M+149q/Siiiivx9/wCClvxe/avvfj9+wf8AsNfsj/GvTP2YPF/7Yer/ALRnifx5+0hL8NvCnxd8S/Dz4a/s0fD7wv4q1jQvBPgLx1u8IXviL4ga1470+x+23iSHTreOSWJfMwRpf8EvfjV+1Xq3jP8Abe/Y+/bI+KXh79oD4xfsQ/Gr4f8Ag2x/aH8P/D/QfhTL8Xvhp8Z/gz4P+M/gDVfFHw88Kz3Hhfw94y0aw8RzWd8NP8u2dUiG15EkuJ/1wrA8UaXqGu+GfEWi6VrN14d1PWNC1fS9N8Q2SCS80K/1CwuLSz1m0QvGHutLuJlnjG8ZeMcivwj/AGOPhz8TLfV/+CYvwWvv2V/iV8JvGn7B3h/4r6P+0L8WPFPgyLQPhzqFvefBnxp8JtV0z4b/ABCZRYfGJPjz8V9f0/xjPNpM8oD6YL7UFExRW/QP9q7x3pXxK/ZhbS7/AOHn7Ytl4A+ON1q/wy8c3/wP+EkN38dfhr4O1SLxBo2teINW+GHjDwt4q8dW+ga/HpzWDT6R4b1TVobLU47y3SD5Lhfk/WPBnxA8E/8ABCj9ob4feP8AQtb8IXfgn9hv9rPwb4J0Lxbo+h6F460j4K+FPh/8UdF+B1t8QdD0B7jRNH+Iq/BnTtE/tqGBi0OpmXzCtwJFX9dfhw/mfDzwG/8Af8GeF3/760OxP9a7Ovwe/bN/4N6P2E/2q9X8YfELwbYeJf2afjP4tu9V1y/8a/C68e88Hat4o1rUI9Q1PXfFXwq1q5PhrU2vXV/Oj0ybRHlZ97SFhz/Kp44/4Nxv+Clnh/42a/8AC3wn8LvDnxD8LJFDf6H8b7Txr4G8M/CzU9El1CKNrye28S67J4s8O61bSPiTR4dMmvGVGkRnhxK3234G/wCDUb9rLWkRfip+1L+z14Qil2iWHwZ4Y+IXxKvEgIBMN5d6vY/Cm21SVMn70ajPRiK/TL9hv/g2y0b9j/8Aak+Dv7Sev/tcz/Fm3+EGv3Piqw+Hdv8AAaDwNZanr58L6/ommXbeI7n4weOJLGPTNY1qPUdsdiWle1VN0e7zF/p/r84/+CmEl1H8Ov2XGtoZZlP/AAUZ/wCCfYuxDy62o/aj+HjM5GRlRKEBHvX6OUV/Px/wUh/4ID/Dz/goP8cvEv7Q9z+098Yfhv471/w1omgweHtR0/SviP8ADjQG8O6TY6NYR+FvD17f+G9W8MaJdxWTXt3ZWmoKk2sTyXwYSMUP8tnj3/g26/4Kp6B411nwnoHwa8E/EfwZot2sWlePvCfxn+F3h/R/F1sygtfadpnjXxb4b8XaQEPyta3+lKmejEc1t/Bv/g2S/wCCnfxB8Wy2fj3wv8NP2evC8Y3nxD8Rfin4V8dvGoIBh0bRPhJqfxAu765bgqLiXS4+OXWv2S/Zb/4NZoPhX8Xfhj8UfjT+12njPTfhx4y8L+Pf+EI+GvwmvvBeoa3rfhS9OsaXayfEfVfiRrN9pujw6yqu0cOkCWWIFVkiY71/rtr8zf8AgotDNN4//wCCYHk7P3X/AAUy+Gs0vmf88U/Zo/atMm3n7/PHvX6ZUV/mWf8ABf79j79o/wDZ+/bu+KHxW+Ll/deN/A37Tfj3xn8S/g58QYptXudMuNCa8igg+Ess135jaH4q+GXh02GnJah3R9JFq1mRyq/iSvhrxZ/bOn+D9O8Na03iPxRcWcGn6CugX/8AbOtPqd6P7O0/TtNUHV9UZnPAAO7NfYE3/BMX/gox4c0y513V/wBhr9rSx02G2LXdzP8AAP4iQra2Sffa5hTwg8kKr/tqBX6yf8EI/wBi39vnwd/wUb/Z4+Mqfs3/ABl+H/w48Iax4uHxI+IfxM+HXiXwD4Mt/Aet/Dnxhoev6S1/4r0PRf7e17V5dcszp9ppyy+RqDrNIoiV5E/0Za/P34SvJJ/wUu/baBVPLg/ZZ/YLiDDO/c/jb9tuYBvr5jfhiv0Cqtd3dvYWl1fXkyW9pZ2813dXEpxHBb28TTTzSHskUaFj7Cv88D/gqF/wX3/aV/aq8SeLvhz+zr4o8W/s9/sxW19e6Ja2/g+7Gi/FX4p6E0D2z6x468Z6U8ut+F9O12GWVBoGny2tq8MhiuJtS2h6/m2vpfMfeksss2ozie4nmuP9L7HPoetdV4M+FvxB+KniHRtE+G/w/wDHXxH1bUCLW00XwH4V8Q+L9Su75iAqx6To+ka7JIxz0ANft5+xJ/wb7/8ABQb9oL4k+EL74hfCfxX+yn8Iodc0+/134kfFJdL8P+MtGsLK+bUrk+Gvhhc62PHt94iVojHpzXFhpNqknkGSaJVdh/pSwReTDFD5ks3lRJF5s7eZNLsUL5ksmF3yPjLHHJr8xf2+o45f2lP+CTgVtt6v7c3iWSH3tY/2Qv2l3v164+ZFSv1Ar46/4KDn40j9h/8Aar/4Z4gv7r4zH4HeP18CW+k8a9LqD6FcrfL4Z/jbxZ/YzXJ0oL87al5AX5sV/ko65ew2y3sVzbyx3MbXnnwzZ+1wE4GNQHGORXPjwx4ru7ObxOPDOuy+HoocXGuDR7/+x7c5HB1MD+ydxx61ilnfNxGIraz9fPBAzj3H51/phf8ABuB8B/i/8Cv+CbOhr8YtI1Lw5f8Axb+LHjH4zeB/DeswzWeq6P8ADfxR4b8D6H4bnvdMuP3ulnxTceFrvXIoujQaqko4kr966/Mb9q25gH/BRb/glHZSY8+bXf21ry3y2Dmz/Ztkt5cL/H+71H8K/TmiiivxC/4LIfCX4Mi8/Zd/a9+Kf/BRrxH/AME2tb/Zo1D4w+Bvh98UfDOheBvE+q+MtS/aM0rwJpeu+EdM8O+MtG8R3Ou6u2kfDYtHbabp93c/ZnuJmVI4TIv4+fHvwF8H/wDgnL+1B8Q9D+NH/Bx1+1N8Hv2lv2l7L4ffEj4yNp37M3w98eXl9o3hHw7YfCn4d+NPipqPgj4NeMvCfwt0Ww8PaTaabbXGpyaRDJAEuJA0Za4r+vf4NaXe6J8IPhVouo/E27+NeoaR8NvA2l3/AMZb9dDW++Ld7YeGNLtLr4m3q+GEi8NLd+PZ4m1WQacq2Ie7PkARbBXpNFFfF/8AwUgyP+CeH7ee1Wc/8MYftRAIuMsf+FH+OcKM926V9I/Cjb/wq34a7fu/8IB4O2/7v/CO6dj9K7+iiiiivzk/4KWSMvw//ZWhVrhftX/BRr/gn/GRBL5SusH7TfgK/ZbnP+tt8WmSndwtfo3RRRRRRX5u/wDBQKGOf4j/APBMaJ8Z/wCHkPguaMHPL2n7K37W91xx1Cwmv0iorPvtL03VEt01TTrHUUtLuG+tY7+0t7xLa+tiTb3lutykiwXcBc7JFw6Z4NOm03Trm9stRuNPsp9Q0wXI02+ntYJbywF7GIbwWV1JG09oLuFQsvlld6jByKvUUV8AfCAs3/BSb9uT5g0a/sy/sCIV3klJf+Eu/bdcjYeF3Rupr7/or4d1r/gmb/wTw8R+Kb/xpr37EP7LOseJdVnjutT1DUvgd8PbtL+7iQRrd3VhNoD6bNdsq5eVoS8jfMxLHNdvon7CP7D/AIauIrvw5+xr+yp4fuoP9Tc6J+zz8ItKuIf+uU1j4Qt5I/wIr6b0vSdL0Owt9L0XTbDSdMtVMdrp2mWdvYWNtGSWZILS0iit4ULMThVAya0KK/Mn9uwF/wBp3/glBH9i+1L/AMNqeOJjPni0MH7HP7TBD/Us27/gFfptRXzl4l/Y+/ZK8Z+Mb74ieMf2W/2dPFnxA1O6W+1Pxz4k+CXw013xhqN6ihEvL3xNqnhm61u7ulRQBJJMzgDANfQkVnaQWsdhBbW8FjDbpZw2UUMUdpFaxRCGO1it0UQx2yQqEEYUKFGAMV8zWn7D/wCxZp/imLxxYfsg/svWXjWHU/7ah8X2nwA+FFt4oi1kS/aBq0fiCHwkmrR6kJ/n+0CXzd/O7NfUVFflz+1ZtP8AwUu/4JPq0W/Cft3TJLkjypF+AvhmIcZ+bfHO49q/Uaiiivwj/wCCput6J8A/22P+CbX7cHxx+Gfj/wCI/wCyZ+zvpX7WPhT4h694C+GPiL4wN8Aviz8X/D3wotvhN8afEvgTwfpXiHxRPoS6f4Q13RV1e10+4k0W7v4dhV7xQ3wj8Ef2kfhB8Iv2d/8Agpl/wUw/bJ+GXxK0a9/4KwftGeMPg9+yz8DfFXwl8aaz8aP2gfgP8PPgnf8Awr/ZX+Dtt8LdN0bWdY0fWfivpek+Ip0tb5bbTDb3yXElyYJY5n/eT/gmR8Kfid8C/wDgnd+xJ8HPjPHd23xU+Gn7L/wX8G+OtLv5luL3w9r+heA9Fsrvwpd3CSTRT3HhIRLpjujvGzWhKsVwT9zUUV8af8FGSi/8E9v272kUtGP2NP2n2kVG2OyD4JeOC6q/OxmXoex5r6B+DM/2r4P/AAouc5+0/DXwLP6Z83wvpcmcf8Cr0qiiiiivzT/4KcX1hZ+Df2OY76fyZL//AIKWfsGWOnp/z9X7fHjw/cpB+FvbSyf9s6/SyiiiiiivzX/b6tXuPjP/AMEtWDbUt/8AgoVaTOP7zL+x/wDtcsn/AKCfzr9KKKKKKKK/Oz4GiI/8FLf+CgrJ/rP+FB/sARzf7y3P7XboD/wCUV+idFFFFFFfm7+2Rhv2wf8AglEhl2H/AIah+O0gj/56+X+w3+0x/wCgbv1r9IqKKKKKK/L39p+SNv8Agp1/wSvtpHcN/wAIn+37fQov3Xltvhb8JbQl/ZYdSfHua/UKiiiikKqSGKgsuSpIBKkjBIJ5GQaWiiivi7/gpC1sv/BPH9vI3kxgtm/Yz/afS4mVPNaKKX4JeN42dY+TIyhjhe54r6A+B8Rh+C3wgiLFzH8L/AEbO33nZPCmkqWJ55OK9Rooooor8q/+CrMRm0H9gPbIY3i/4KqfsKTKO0oT4k3nmRtweDCWI9wK/VSiiiiiivzC/b8u4/8Ahon/AIJQ6YIfMubr9vDV9Rjl/wCeUGl/sh/tPR3XP+0NRQ/hX6e0UUUUUV+aP7PAb/h55/wUqczGQf8ACnv+CfMYXtBt0b9pxzBz3Jk8zj/npX6XUUUUUUV+X37Zc1xJ+3p/wSG0yOEvby/Hf9qnVbiUHAhbSv2J/jVFBuGCD5jamw/Cv1Boooooor8s/wBpt1f/AIKof8Es4YxG8sPw0/4KF3M43Rebb2svgb4CwpLtb96FmuQFyvcc1+plFFFFFFFFFeZfGf4SeDfj18JPiX8E/iFBqVz4F+LHgjxL8PvF0Gj6pd6Jqsvh7xXpN1o2qpp+q2Lx3Vhdm0vH2SLkBvvKy5U/Bdp/wSn+E+nQ6XZaT+1B/wAFFNG07QbKx0vQNM0r9vL9oey0zQ9K0vTv7M0vS9Ms4/F3lQWOm2oHkoQwUjkkcVXT/glP4Cjubq8X9s3/AIKdC6vVVLmc/t+fHnfIkYCorY10DaqjAGMAV3dv/wAE3/h/b24gX9p//gofKw/5eLj9vn9pue4PPcyeP2j7/wB2vnTwx8Hf2LvGPxq8Rfs4+Fv+Chf7aevfGvw3Pqn9r+BbH/goD+0dNqVvc+HRG/inw1p9+fFw0PV9d8NrMDq+nWlzcatpADG4S22Nt9TH7Bfwi8QeJfE3w80X9ub/AIKAp438IWOg+I/FnhzRP2+fjFeeJ/DOj+O5fEEXhG81rTLvxDqU+k6ZrsnhXUf7OM0SCf7DNsLhHrO1X/glNo+olDB/wUE/4Kr6N5bFsaX+3H46If2kGp6RqW9fauKvf+CPzzQulp/wVB/4K9W08oCSTzftqX10rRDdlfJXwLbBH+fhlZSPeuU1j/giF4D8XaT4a0n4ift8f8FNfiXH4Q8SaL498NSeOP2sdY1v/hH/AIleGNx8LePdCSfw9/xK9c8Ms7GxdCXgY5Dk813Sf8EiI0l85f8Agpz/AMFfSdmxo3/baupIW6/OY5Ph4wV+eq7a3dM/4JQ2GnyPJc/8FFf+CsGsl4zGE1P9tvxKqKD/ABIml+F9MUS8fe616Fp3/BNzw/p0bx/8Nm/8FHr7e+8Saj+2f8TLqRP9lGLINnsQa+efhl8F/wBkf4/ePfF3w1+Dv/BTL9tzxp44+Htp9q8SeFvCX7c3xNurqOxtNVXSbvxFpk92ZU8WaJp+sOun3t9pk97ptreSx29wyTyRqeG8L+DP2I/Hfx0H7NHg7/gqb+3drPxw1K88c6fafDnTf2zvi7Pq8974AGpnxpDpt/d6PJpsreGV0i6aTyrtl2wsV3jFfUmof8ExdIvo5Ui/br/4KcaaZo44/N0/9tLxyJI/LB+eL7Zpl6iSSfxHbzWN/wAOrNLClR/wUF/4KpDP8X/DbXi9m/AvoLYrM1j/AIJAfBnxPq3gDxF4w/ah/wCChPjLxP8ACzUdW1r4feKPE/7ZnxR1XXvCOu65pk2iahrugX8kqHSdYm0a4ktPPtxE3kSOpzubPtVt+wT9m8JWXgX/AIbQ/b6m8M2WpyauVuP2ihN4tvL2VQCt98VG8Hf8Lcn0tQPlsP7eXTwefIzXO/8ADtjwr5yzf8Nhf8FHsr/yz/4bg+M3kt/vR/2rg1R1L/gnt4Q8N6Vq+tav+3D/AMFCtD0LS7K61jWdW1r9tLxzDpej6VpiSajqOoX2qa0skWmabZ2kLvPNJKiRQISWVQTXjfwb/Zw/Zz/aT0vxNqHwD/4KZ/t2fELTvD+oW2n+IH8Gftq+KNXudCur5Z5bC8dNT0m71ODS/EVvbSy6ddfNp19BG0tizIpauW+Dfwr/AGRf2h/GPiL4bfBn/gp7+3V8U/FvgzT9bude0zwr+2Z8SZvsFnomuW3hnVr461p2j2GmaoNN169itzJFc3A8xxjcnNfRv/Ds/wAP7t3/AA2x/wAFKv8AW+bj/htb4mbemNmMf6v2qWL/AIJraRE0zf8ADbn/AAUmlMvTzf2x/HLrF/1xQWCov4g1xXh3/gkn8N/Cfjbxf8SNA/a7/wCChth48+ISeF4PHnitP2rvED6x4xs/BOnajpfhCw8QzyaA8d/Z+GrTVrkWalQYvPcA4Yg9tN/wTasZ02P+3L/wUlX975u6H9rjxFA/f93vh8Pxt5XPSkl/4JuW0olA/br/AOCk0Jl/ii/ax1UGP/rkH8LOqflWh/w7xmC2o/4bt/4KLf6K27f/AMNJWRa49rrPw/xMvtgV8/eEfhR+zb8SYfHcvgL/AIK8/tW+KIfgrZ6z4g+J1x4Z/bM+GGrJ4F0azS9mv/EHjeaLwLdJp3hjS4tMuSt3d405Rby5dvLbbnfCL4d/szftLa1rPgb4Ff8ABYL9rX4reJdE8ORa7rXhn4bftd/DDXPFOj+H5b06aniK9tdH+G0viC0szqN7FE1w5EaytCrY3qH9svP+CaetXKIlr/wUi/4KeaVs6vZ/tD+ALl5P99tY+CWq/piseL/gl/4jiYMv/BT/AP4KqN7S/tDfCqVfxWT9n1hWHrX/AAST0zxH4s+G3jvXf+Cg3/BSPV/GXwh1fxDrvw58R6p8avhBqWo+F9U8V+G7jwh4iubE6h+z7dW5/tfw3eTWkyPG0bQzOu0B2z2dx/wTY8V3c0c91/wU0/4KeyGNNvlW/wAbfgvpkL/7UkWj/s56eHf3q7b/APBOXxXBn/jZT/wUwmz/AM/Hxm+C8pH0z+zsCKZrH7C2seE9M1/xX4m/4Kdf8FCdC0DSbC71vXtd8RfGH9nnTPD+gaVptvJd6hqd7d6h+zjbaTo2l2FpC0s8shjhjjUsxABNc34B/Zv8N/EHwm3xS+Gv/BXD9sb4hfDzRxr8V1458JfHf9k3xt4Ctm0lGl1s6n4i039nbVfDzSeHY1ZpvPmP2VQTKABWt4E/ZlsfjR4ah8WfCj/gq/8AtmfEvwlFe3+lp4q+GXxg/ZK8W6BJq1rsXULOXXfCv7Nmo6bPqFh5677dpSIS6kxg4p7/ALA/7SE2sfaZ/wDgqn+2p/Yf/QMttC/ZptNS7/8AMXT4HNaYP/XjXUzfsEfEO6Xbcf8ABSn/AIKJbu7WvjT9mOwz9Psf7K0TD8zWj8MP+CfHh7wH8evh1+0X42/aa/av/aG8efCbwl4/8H/Dq0+O/jX4Za74b8K2fxNg0Gz8X6lp+neB/hD8P7w63qdl4ctoGnluZA0YIdXIQr+glFFFFFFFFFFFFFfzx/Bv4OftE+D/AAJ+yH+yNrvwG+MOnav+wF+0d4v/AGjviN+0Jpen6E3w7+N/gvw3Z/HjUtFg+COvDxZLrnjT4pftBn4tWVtq+kajHp0lmL3WI9RuN6Bbj9PfhX4G8cR/tx/tDfGi88JanoHw7+KH7K37Hui6JqGrQ2ltf3PjHwZ41/aq1PxL4ev7aO5nnsta8LaT4500XsOXjX7XFh2Nfb9FFFFcF8U/B998Qvhj8R/AOl6/ceFdT8ceA/F/hDTvFFnF5934bvvEvh/UdGtNetoPMh8650e4vVuEXeu5owMjrX4l/s3eBvjdqGof8E+NI179k/4w/CL/AIdZ/CL4ieGPixrmvad4OjsvjX4kb9nWX4EWPgb9mm50Px5L/wALP0Tx5rUZ8WXGoX1tp+mo2mafG8x1CUJF6NpnhL42eL/+Cn3hH43eGvgp+058NvE2gy/EP4QfH/xR8V9d8AeJv2U739lHTdB1e78AW/7P19p9w+vT+L/ip8StJ8J+J5YbA2uo6NeLqVtrMBiXY37V0UUUV8Y/8FDvgZ45/aU/Yr/aI+CHw0GnTeOvHngN7Xw1pmr3q6bpXiPUNJ1jS/ER8H6jqUjxw6dZ+NLfSH0mS4kZYoFvS8jKik18v/AnWPiF8Sf2t/E/7Xw/Zg+PfwS+Hnhz9mT4V/smR/Dv4heE/DvhTx94x8d3/wAdLvxNe6t4f8I2viybTrn4T/A/R/EZ3eILmW0tpoNS1F9NS4it5i/zr+wX8JviZ8Jf2hfg38N/gb8Ef2xv2b/gX8NvBvx+tf2vfCn7Snj/AFn4kfCHxn4613VPDsnwin+AviTXPFWueGvEOsweI49TvZtY8F6VommSaIBHqA+0SwxN+/FFFFFeKftJeCPFvxM/Z2+Pfw38A6kmjeO/iD8Fvin4I8FaxLdPYx6V4t8V+Btd0Hw3qMl7EryWcdjrF/DKZVBaMJuAJFfkP+ypfeO/E3xv/Yt8RWn7KX7QPwC0n9j79ku+/ZY+OuqeMvg/eeE7Xxd4r+JWqfAXwd4I+HPw+uI3fUvib8Jfhn4k+G2r+Jb3xDapLoOj6dJBemRBeSGvsf8AYJ+GGu3+tftJfta/GD4daj4T+Pf7QXx7+K2gwTeMPDbaD408M/s6fBvxvqnwo+AngK1tbpjcaf4YvPCPgmLxIzRrCus3+tyai4kWW3Zf0coooor85P8AgqD4C8feOv2ePB1x4O8D+Iviz4Z+HX7SX7O/xc+NvwT8JaPF4o8S/Gn4FfDf4maP4k+Ivw90TwZcL9j8d6jcWlrDqUehzvGmqtp32cFndEf4C8T/AA40n4/6J/wUz1KX4S/ta/C/4D/tk+FPhLo3wy0Dwj+zn470j4keIfiP8Afh14y8X/En4pH4P+JvD2laZomkfEpYPC+gQR+Kf7Efx6+ltpgJ8yN69m/4JW+CvjLP8fP24/jv4qsbzR/g98WH/Z38OfD6K9/Ze8S/sff8JV43+F/gfxHpPxO8VW3wL8ceKPEXjTw41hcavYaTd6lebE1i8s5Eh2pp4jX9sqKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK/FH/gsv8AtnfEz9lSx/Y18EeFvjVpv7Ifwz/ae/aLn+E/xu/be1zwZ4f8b6X+zp4SsPBGteKdJtLKy8bafqvw+0HxV8VNc09NLstb1+1u9K0a2hu7qaEiPzoNz9jnTP2xPDHx68E3/gn9vDwx/wAFNP8Agn38SPAnjceKPin4pu/2en+KnwI+K3h9dG1LwU/hrx1+z9ovhLw38WfBHxIs7+5tbmzn02W80eaCKcXCwP5cnxz/AMEc/wDgo7+1V8fP24f28/2af2s/iFb+PfDcXxN/aC8VfsZahH4L8BeDR4e+Gv7On7T3j/4A/FX4V3V74N8L+GJPFur+ErfUfBOoxT35v9Saz1J5JriQ79vjth/wWK+P+i/8FXf2473xj46e4/4Jr/szfst/tmeKvBnwq0zwZ8P4td8Y+Pv2E4vglo3xl8dab4+fwwfHGonVfin4m8UeGtKtf7aOmSS6b81sHCtX1x+zx8Df+Crn7Z3wE8B/tc/ED/gp34w/ZJ8f/HfwTpPxe+FH7O/wI/Z0/Z08ZfBL4FeEfiHodnr/AIF8I+PJvi14H8U/EL4367a+Hb20n1i4utZ0oRX808Fn5KRpM/iEP/BYr9qy2/Y6vvhze+AfhTcf8FPrX/gpSv8AwSHsTGmqwfAXUvjvciLXLH9piXQvt8niS2+FB+FMja5JYecJG1WI24WO3dY19v8A2hPgJ/wVj/ZA+AHj/wDaw+Hf/BUXxb+1T8UPgd4H1n4vfEn9n746fs2fs6eFfgL8c/Dnw80C71/xj4H8GRfCnwZ4V+JPwX1PUdEsLmTSbiz8Qag0t9HBDds6yyXCebft2/Gn9tLxV/wTo8c/8Fbf2QP+CiHxQ+AXw31T9jvwJ+038Ov2XT+zn+yR8RfDuh3Gr/Dnw54hutA1r4g/EL4S+K/H1ze3Wo6hKb/N7LHb3fmR24SFUQUf2mPFP/BR79jT/glH+0n+3LqX/BTH4jfHP4jR/so+AvHPw20PxZ+zD+x74P0T4X+PvGHiX4Z6jc+K9NPgj4Nab/wlU1nomqX+mR2esxXuntDeNM0P2iOGRPqv9sv/AIKy/sweCP8Agnt+0t8R/gV+31+yVJ+0z4X/AGTfid4w+FMXh746/s/eMvFj/GbSvhbq2reEk0b4fXmu69p/iXXW8XwwiHSZNMu4rqfEDW0it5Z/NX9ov/gqf+318CP2gv8Agnn468Pa63xS/Zs0n/glT8OP28f+Cgfwn0/4f+AZfFvjjwRqfjb4ffDL4w/GHwLf6X4Og8UaZ4n+FsHxTg8YHRtL1DStEm0/QrqN4VjJU/rh8Y/2u/Hcn/BRn/gkx8Mvgx8UNL1T9mX9r/4J/tvfE7xnbaHpPhHXtF+KenfDn4a/Azxj8FPFGkeL73Rb/wASaXplgnjq7vIDpF/ZQahFeL9qWdEiCfkH+xx+0j+1T+1P8MPGHxV+LX/Bwf4W/ZV8VR/H/wDaL8A23wP1T4Mf8E2Ybzwj4V+Gnxn8Z+BvB/mH4j/Diy8Y3C6l4b0O2uFmvFaSZZA++QHefrj9vT9r345/sgfs7/8ABPH4fz/t/aY/hH9qT4+a58Ofjp/wVt1n4VfA+90XwH4CudK8ZeP/AAZNpHhfwd4dX9nHwzrvjxVs/C+l6/d6dcaNYWGl3Oo3UUk/m3UX01+xzpn7Ynhj49eCb/wT+3h4Y/4Kaf8ABPv4keBPG48UfFPxTd/s9P8AFT4EfFbw+ujal4Kfw146/Z+0Xwl4b+LPgj4kWd/c2tzZz6bLeaPNBFOLhYH8uT89vh7+1X8Yf22/2ov22/hp4x/4K7v/AME2/jJ+z1+1n8Uv2f8A4Afsc6F4L/Zfs9X1H4feAH0+w8A/G/xppH7RHg7X/HPx9h+MMd2dUNtoWo22iW9s0KW0iiZJT+/H7GF9+1de/s3fDpP23dE8BaN+07pcfiHQfiXN8MLyG78C+KJNB8U61pHhjx3oMUFzcrpP/CeeDrPT9Wu9PJQaffXk1uqJHGiD6jooooooooooooooooooooooooooooooooor8tP+Ck/xf+OfwU1P9njxRZfsy3H7W37Dmua58RvBX7d3wp8FfB9vjt8ZdK8N694c0ub4N/EjwT8L/tLDxl4K8IeNtPuh4usorHVL7+zrqGe3tmMEhH5Gfss/Bz4BePf+CtH7LX7Rn/BKr9jv48/sk/Abwl4J/aEg/b58f6p8AfiZ+yJ+zl8aNA8TfDxNG+CPwv8ADPwh+IWm+CLDxl4/8JfFaVNYnk0PQLK106FFuLiW8ZYlg8Fuvgd+2X+zJ+zTq37dvwP/AGW/jT41/ak/ZG/4K/8A/BRzxv4V+CMXwt8ez+PPjL+y3+2D418VfDTxJd+FPAVroD+KfG3hbVdUv/CXi2xvrG2ubJ7DQZL5GMURnj9U0f8A4JNfEhPHf7MX7H3iLw142i8OeLP+CFf7aX7P/wAcf2gz4S1zxB4B8P8A7V/7TPxO8AeOPHuq654ztrP/AIR2TxbqfxT8T654hsdOnvFvL22tWZA6Ru4+uf2W/wDgqR8Q/wBlf9m/4X/st/thfsCf8FArX9rr9n74a+G/g7d+G/gX+yx8RPjx8Of2gr34YeG9N8M6Z48+CPxh8CWuo/DvWdC8Z6Vp9tfXTajqGmrpF3PPbzHEG9vlq9/4Jr/tv6f+x9Zftlp8L9E1H/godb/8FaT/AMFh9W/ZatfFOiSlvDN3A/w9uv2PNO+JDXUGhS69YfAsxTfbd0kT69E9nE8qskr/AFT+0r/wU6+I37Y37NfxQ/Za/Yy/YT/b2tv2sP2hPh74i+Co0/8AaH/Zb+JHwA+GX7Ns/wATvDV54d8QeP8A44/F/wAb6fB4D0ix8C+H9TvL+yg0u71e41i+tYbe3ik89TXuH7c37IviT4Tf8EDPjL+xJ8F/D3i74weLPhj+wXpPwF8F6L4I8J6r4h8bfErXPBPgbQfCwvNG8HeHLXU9X1LXvFV5pkl41paRTSGWZgoOKyv+Cpvwk+K3xC/4IIfGr4O+Afhl8QvHHxd1X9kf4QeG9L+FnhDwX4k8S/EfUvEWmTfDQ6loFh4H0XTb3xNea1p4sJzPax2rTxeTJuUbGwz9vT/gmv8AsqXX/BMD9rSy+EX7An7Plz8erj9ir4vW3w7tPhx+yv8ADib4vT/FKT4Oa3H4bt/Bdv4Y8CN4zk8fyeJjEtiliDqJvygiHm4ryr9lT4CfFSy/4KG/sH+IvHXwX+IFp8NdB/4N7tJ+BHxD13xZ8OfEdv4G0X4oT/GH9ny71P4NeM9T1jR00DTvHs2h6Tfyz+Hb501FrS2uGe38uOQj5O/Z/wD2Mv2rf2T/APgtL+xf+z3Z/DH4ieNf+Ce37MfhP9vL4g/srfH228O+KvEHhP4T/DX9rDwb4LuY/wBlb4g+OE0680TQb/4PfEj4c6lbeGhqeotfajoGtWCJxEsMfgH7CVp+yH8BPg/43+Hv7af/AAQs/az+Onx5j/aV/aj8S6t8UJ/+COHiH48J4o8KeL/jz488SeAbqy+KHiL4Z3mo+I9MTwfqNmts3mPDFCFSImMKa/Xb41ftJ+L/AAn8Cf2MfiR8Gf8AgnX8StZ/4Jq6zffFz4Wftf8A7Gut/sXX+iftKfCrwLaodB+Cvjbwj+yNfW+mSxfDPw5458N3l1rmm2+kX1xceHdSsr+ytWQMw+H/ANln4OfALx7/AMFaP2Wv2i/+CVX7Hfx5/ZJ+A3hLwT+0JB+3z4/1T4A/Ez9kT9nL40aB4m+HiaN8Efhh4Z+EPxC03wRYeMvH/hL4rSprE8mh6BZWunQotxcS3jLEtv13xd+MHwB8fH4u/Ab/AILc/wDBMDxv8afjt4P+IfxG8MfCL4wfBf8AYF+Jfxz8JfH74I3fiXUr/wCD+qfBH4y/CrQ/GPijwJ45fwbdQ2urabNrWhXVlewM7tE8rxRfor/wQ6+DPx3+BH7AfhbwJ8c9B+JXgaE/FT40eIvgX8KPjPrn/CRfF74Nfsz+IvH+q6j8Dvhf8StUOo6sw8WeHPB8kbz2jXDvpsdwlkywtbmCP9e6KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK/9k="
     * /><br/>
     * Image modified from "Krumbein, W.C., and Sloss, L.L., Stratigraphy and Sedimentation, 1956,
     * Freeman and Company, San Francisco CA" <br/>
     * <br/>
     * NOTE: this descriptor is very sensitive to noise, as a single erroneous surface point may
     * dramatically decrease (resp. increase) the radius of the inscribed (resp. bounding) sphere,
     * and may result in both cases to a much lower roundness than otherwise expected. For a more
     * robust measure, consider using the sphericity (see {@link #computeSphericity()}).
     * 
     * @return
     */
    private double[][] computeRoundness()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Mesh3D contour = (Mesh3D) det;
                Point3d center = contour.getMassCenter(false);
                double minRadius = contour.mesh.getMinDistanceTo(center, null);
                double maxRadius = contour.boundingSphere.getRadius();
                result[trackIndex][det.getT()] = 100.0 * minRadius / maxRadius;
            }
        }
        
        return result;
    }
    
    /**
     * <b>Description: </b>The sphericity (as defined by H. Wadell in 1935) is based on the ratio
     * between the volume and surface area of the object. Given a sphere with a volume equal to the
     * volume of the object, the sphericity is the ratio between this sphere's surface area and the
     * real surface area of the object.<br/>
     * <b>Calculation: </b>Assume the object has a volume <code>V</code> and surface area
     * <code>S</code>. A perfect sphere of volume <code>V</code> would have a surface area of
     * <code>"S' = cubeRoot(36 * Pi * V)"</code>. The sphericity is thus calculated as
     * <code>sph = S'/S</code>, and expressed as a percentage.<br/>
     * <b>Interpretation: </b>sphericity, roundness and smoothness are particularly challenging to
     * interpret when their values differ from 1 (the theoretical value for a perfect sphere). The
     * figure below summarises the differences one may expect:<br/>
     * <br/>
     * <img alt="" src=
     * "data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/4gfMSUNDX1BST0ZJTEUAAQEAAAe8YXBwbAIgAABtbnRyR1JBWVhZWiAH0AACAA4ADAAAAABhY3NwQVBQTAAAAABub25lAAAAAAAAAAAAAAAAAAAAAAAA9tYAAQAAAADTLWFwcGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAVkZXNjAAAAwAAAAG9kc2NtAAABMAAABi5jcHJ0AAAHYAAAADh3dHB0AAAHmAAAABRrVFJDAAAHrAAAAA5kZXNjAAAAAAAAABVHZW5lcmljIEdyYXkgUHJvZmlsZQAAAAAAAAAAAAAAFUdlbmVyaWMgR3JheSBQcm9maWxlAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAbWx1YwAAAAAAAAAeAAAADHNrU0sAAAAqAAABeGRhREsAAAA0AAABomNhRVMAAAAsAAAB1nB0QlIAAAAqAAACAnVrVUEAAAAsAAACLGZyRlUAAAAqAAACWGh1SFUAAAAuAAACgnpoVFcAAAAQAAACsG5iTk8AAAAsAAACwGNzQ1oAAAAkAAAC7GhlSUwAAAAgAAADEGl0SVQAAAAuAAADMHJvUk8AAAAkAAADXmRlREUAAAA6AAADgmtvS1IAAAAYAAADvHN2U0UAAAAuAAAD1HpoQ04AAAAQAAAEAmphSlAAAAAWAAAEEmVsR1IAAAAkAAAEKHB0UE8AAAA4AAAETG5sTkwAAAAqAAAEhGVzRVMAAAAoAAAErnRoVEgAAAAkAAAE1nRyVFIAAAAiAAAE+mZpRkkAAAAsAAAFHGhySFIAAAA6AAAFSHBsUEwAAAA2AAAFgnJ1UlUAAAAmAAAFuGFyRUcAAAAoAAAF3mVuVVMAAAAoAAAGBgBWAWEAZQBvAGIAZQBjAG4A/QAgAHMAaQB2AP0AIABwAHIAbwBmAGkAbABHAGUAbgBlAHIAZQBsACAAZwByAOUAdABvAG4AZQBiAGUAcwBrAHIAaQB2AGUAbABzAGUAUABlAHIAZgBpAGwAIABkAGUAIABnAHIAaQBzACAAZwBlAG4A6AByAGkAYwBQAGUAcgBmAGkAbAAgAEMAaQBuAHoAYQAgAEcAZQBuAOkAcgBpAGMAbwQXBDAEMwQwBDsETAQ9BDgEOQAgBD8EQAQ+BEQEMAQ5BDsAIABHAHIAYQB5AFAAcgBvAGYAaQBsACAAZwDpAG4A6QByAGkAcQB1AGUAIABnAHIAaQBzAMEAbAB0AGEAbADhAG4AbwBzACAAcwB6APwAcgBrAGUAIABwAHIAbwBmAGkAbJAadShwcJaOgnJfaWPPj/AARwBlAG4AZQByAGkAcwBrACAAZwByAOUAdABvAG4AZQBwAHIAbwBmAGkAbABPAGIAZQBjAG4A/QAgAWEAZQBkAP0AIABwAHIAbwBmAGkAbAXkBegF1QXkBdkF3AAgAEcAcgBhAHkAIAXbBdwF3AXZAFAAcgBvAGYAaQBsAG8AIABnAHIAaQBnAGkAbwAgAGcAZQBuAGUAcgBpAGMAbwBQAHIAbwBmAGkAbAAgAGcAcgBpACAAZwBlAG4AZQByAGkAYwBBAGwAbABnAGUAbQBlAGkAbgBlAHMAIABHAHIAYQB1AHMAdAB1AGYAZQBuAC0AUAByAG8AZgBpAGzHfLwYACAARwByAGEAeQAg1QS4XNMMx3wARwBlAG4AZQByAGkAcwBrACAAZwByAOUAcwBrAGEAbABlAHAAcgBvAGYAaQBsZm6QGnBwXqZjz4/wZYdO9k4AgiwwsDDsMKQw1zDtMNUwoTCkMOsDkwO1A70DuQO6A8wAIAPAA8EDvwPGA68DuwAgA7MDugPBA7kAUABlAHIAZgBpAGwAIABnAGUAbgDpAHIAaQBjAG8AIABkAGUAIABjAGkAbgB6AGUAbgB0AG8AcwBBAGwAZwBlAG0AZQBlAG4AIABnAHIAaQBqAHMAcAByAG8AZgBpAGUAbABQAGUAcgBmAGkAbAAgAGcAcgBpAHMAIABnAGUAbgDpAHIAaQBjAG8OQg4bDiMORA4fDiUOTA4qDjUOQA4XDjIOFw4xDkgOJw5EDhsARwBlAG4AZQBsACAARwByAGkAIABQAHIAbwBmAGkAbABpAFkAbABlAGkAbgBlAG4AIABoAGEAcgBtAGEAYQBwAHIAbwBmAGkAaQBsAGkARwBlAG4AZQByAGkBDQBrAGkAIABwAHIAbwBmAGkAbAAgAHMAaQB2AGkAaAAgAHQAbwBuAG8AdgBhAFUAbgBpAHcAZQByAHMAYQBsAG4AeQAgAHAAcgBvAGYAaQBsACAAcwB6AGEAcgBv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"
     * /><br/>
     * Image modified from "Krumbein, W.C., and Sloss, L.L., Stratigraphy and Sedimentation, 1956,
     * Freeman and Company, San Francisco CA"
     * 
     * @return
     */
    private double[][] computeSphericity()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Mesh3D contour = (Mesh3D) det;
                double surfaceArea = contour.getDimension(1);
                double volume = contour.getDimension(2);
                
                result[trackIndex][det.getT()] = 100.0 * Math.cbrt(36 * Math.PI * volume * volume) / surfaceArea;
            }
        }
        
        return result;
    }
    
    // private void computeSpeedOverCurvature()
    // {
    // JFileChooser jfc = new JFileChooser();
    // if (jfc.showSaveDialog(getPanel()) == JFileChooser.APPROVE_OPTION)
    // {
    // XlsManager xls;
    // try
    // {
    // xls = new XlsManager(jfc.getSelectedFile());
    //
    // for (TrackSegment ts : trackPool.getTrackSegmentList())
    // {
    // int trackIndex = trackPool.getTrackIndex(ts);
    // xls.createNewPage("Cell " + trackIndex);
    // // System.out.println("Cell " + trackIndex);
    //
    // ArrayList<Detection> dets = ts.getDetectionList();
    //
    // int col = 0;
    //
    // for (int d = 1; d < dets.size(); d++, col += 2)
    // {
    // int row = 0;
    //
    // xls.setLabel(col, 0, "dK [" + (d - 1) + "-" + d + "]");
    // xls.setLabel(col + 1, 0, "dV [" + (d - 1) + "-" + d + "]");
    // // System.out.println("Interval " + (d - 1) + "-" + d + ":");
    // Mesh3D prevCt = (Mesh3D) dets.get(d - 1);
    // Mesh3D currCt = (Mesh3D) dets.get(d);
    //
    // Vector3d globalDisplacement = new Vector3d(currCt.getMassCenter());
    // globalDisplacement.sub(prevCt.getMassCenter());
    //
    // // register prevCt towards currCt
    // for (Point3d prevPt : prevCt)
    // if (prevPt != null) prevPt.add(globalDisplacement);
    //
    // Point3d currPt = new Point3d();
    //
    // for (Point3d prevPt : prevCt)
    // {
    // if (prevPt == null) continue;
    //
    // row++;
    //
    // currCt.getMinDistanceTo(prevPt, currPt);
    //
    // double dV = currPt.distance(prevPt) / trackPool.getDisplaySequence().getTimeInterval();
    // double dK = (prevCt.getCurvature(prevPt) - currCt.getCurvature(currPt)) /
    // trackPool.getDisplaySequence().getTimeInterval();
    //
    // xls.setNumber(col, row, dK);
    // xls.setNumber(col + 1, row, dV);
    // // System.out.println(dK + "\t" + dV);
    // }
    // // System.out.println();
    // // System.out.println();
    // // System.out.println();
    // }
    //
    // // System.out.println();
    // }
    //
    // xls.SaveAndClose();
    // }
    // catch (IOException e)
    // {
    // e.printStackTrace();
    // }
    //
    // }
    //
    // }
    
    @Override
    public void displaySequenceChanged()
    {
        
    }
    
    @Override
    public String getMainPluginClassName()
    {
        return ActiveContours.class.getName();
    }
    
    /**
     * Compute the best fitting ellipsoid for the given component.<br>
     * This method is adapted from Yury Petrov's Matlab code and ported to Java by the BoneJ project
     * 
     * @param in_points
     *            the points to fit
     * @param out_center
     *            (set to null if not wanted) the calculated ellipsoid center
     * @param out_radii
     *            (set to null if not wanted) the calculated ellipsoid radius in each
     *            eigen-direction
     * @param out_eigenVectors
     *            (set to null if not wanted) the calculated ellipsoid eigen-vectors
     * @param out_equation
     *            (set to null if not wanted) an array of size 9 containing the calculated ellipsoid
     *            equation
     * @throws IllegalArgumentException
     *             if the number of points in the component is too low (minimum is 9)
     * @throws SingularMatrixException
     *             if the component is flat (i.e. lies in a 2D plane)
     */
    public static void computeEllipse(Iterable<Point3d> in_points, Point3d out_center, Point3d out_radii, Vector3d[] out_eigenVectors, double[] out_equation)
            throws IllegalArgumentException
    {
        ArrayList<double[]> rawMatrix = new ArrayList<double[]>();
        for (Point3d p : in_points)
        {
            rawMatrix.add(new double[] { p.x * p.x, p.y * p.y, p.z * p.z, p.x * p.y * 2.0, p.x * p.z * 2.0, p.y * p.z * 2.0, p.x * 2.0, p.y * 2.0, p.z * 2.0 });
        }
        
        double[][] _rawMatrix = rawMatrix.toArray(new double[rawMatrix.size()][9]);
        Matrix D = new Matrix(_rawMatrix);
        Matrix ones = ones(_rawMatrix.length, 1);
        
        Matrix V;
        try
        {
            V = D.transpose().times(D).inverse().times(D.transpose().times(ones));
        }
        catch (RuntimeException e)
        {
            throw new SingularMatrixException("The component is most probably flat (i.e. lies in a 2D plane)");
        }
        
        double[] v = V.getColumnPackedCopy();
        
        double[][] a = { { v[0], v[3], v[4], v[6] }, { v[3], v[1], v[5], v[7] }, { v[4], v[5], v[2], v[8] }, { v[6], v[7], v[8], -1.0D } };
        Matrix A = new Matrix(a);
        Matrix C = A.getMatrix(0, 2, 0, 2).times(-1.0D).inverse().times(V.getMatrix(6, 8, 0, 0));
        Matrix T = Matrix.identity(4, 4);
        T.setMatrix(3, 3, 0, 2, C.transpose());
        Matrix R = T.times(A.times(T.transpose()));
        double r33 = R.get(3, 3);
        Matrix R02 = R.getMatrix(0, 2, 0, 2);
        EigenvalueDecomposition E = new EigenvalueDecomposition(R02.times(-1.0D / r33));
        Matrix eVal = E.getD();
        Matrix eVec = E.getV();
        Matrix diagonal = diag(eVal);
        
        if (out_radii != null) out_radii.set(Math.sqrt(1.0D / diagonal.get(0, 0)), Math.sqrt(1.0D / diagonal.get(1, 0)), Math.sqrt(1.0D / diagonal.get(2, 0)));
        
        if (out_center != null) out_center.set(C.get(0, 0), C.get(1, 0), C.get(2, 0));
        
        if (out_eigenVectors != null && out_eigenVectors.length == 3)
        {
            out_eigenVectors[0] = new Vector3d(eVec.get(0, 0), eVec.get(0, 1), eVec.get(0, 2));
            out_eigenVectors[1] = new Vector3d(eVec.get(1, 0), eVec.get(1, 1), eVec.get(1, 2));
            out_eigenVectors[2] = new Vector3d(eVec.get(2, 0), eVec.get(2, 1), eVec.get(2, 2));
        }
        if (out_equation != null && out_equation.length == 9) System.arraycopy(v, 0, out_equation, 0, v.length);
    }
    
    private static Matrix diag(Matrix matrix)
    {
        int min = Math.min(matrix.getRowDimension(), matrix.getColumnDimension());
        double[][] diag = new double[min][1];
        for (int i = 0; i < min; i++)
        {
            diag[i][0] = matrix.get(i, i);
        }
        return new Matrix(diag);
    }
    
    private static Matrix ones(int m, int n)
    {
        double[][] array = new double[m][n];
        for (double[] row : array)
            Arrays.fill(row, 1.0);
        return new Matrix(array, m, n);
    }
}
