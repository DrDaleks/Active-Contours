package plugins.adufour.activecontours;

import icy.gui.util.GuiUtil;
import icy.math.ArrayMath;
import icy.plugin.interface_.PluginBundled;
import icy.util.XLSUtil;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.vecmath.Point3d;
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

import plugins.adufour.activecontours.Mesh3D.Vertex;
import plugins.adufour.quickhull.QuickHull3D;
import plugins.fab.trackmanager.PluginTrackManagerProcessor;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;

/**
 * WARNING: this is an experimental class and may change in the future, use at your own risks!
 * 
 * @author Alexandre Dufour
 */
public class DeformationProfiler extends PluginTrackManagerProcessor implements PluginBundled
{
    private enum Descriptors
    {
        None, Perimeter, Volume, Roundness, RadiiVar, Convexity, Roughness
    }
    
    private JComboBox jComboDescriptors     = new JComboBox(Descriptors.values());
    
    private JButton   jButtonSaveToVTK      = new JButton("Export meshes to VTK files");
    
    private JButton   jButtonSaveShapeToXLS = new JButton("Export shape measures to XLS");
    
    private JPanel    chartPanel            = new JPanel();
    
    private double    tScale;
    
    public DeformationProfiler()
    {
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
                        exportMeshToVTK();
                    }
                }.run();
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
                        saveShapeToXLS();
                    }
                }.run();
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
        
        int dim = 3; // FIXME
        
        ChartPanel chart = null;
        
        Descriptors descriptor = (Descriptors) jComboDescriptors.getSelectedItem();
        
        switch (descriptor)
        {
        case None:
            return;
            
        case Perimeter:
            double[][] dim1 = computeDimension(1);
            chart = createChartPanel(dim1, (dim == 2) ? "Perimeter" : "Surface", "Time (sec.)", "\u03BC" + (dim == 2 ? "m" : "m\u00B2"));
            break;
        
        case Volume:
            double[][] dim2 = computeDimension(2);
            chart = createChartPanel(dim2, (dim == 2) ? "Area" : "Volume", "Time (sec.)", "\u03BC" + (dim == 2 ? "m\u00B2" : "m\u00B3"));
            break;
        
        case Roundness:
            double[][] roundness = computeRoundness();
            chart = createChartPanel(roundness, "Roundness", "Time (sec.)", "%");
            break;
        
        case RadiiVar:
            double[][] radiivar = computeRadiiVar();
            chart = createChartPanel(radiivar, "Radii Variance", "Time (sec.)", "Var.");
            break;
        
        case Convexity:
            double[][] convexHullDiff = computeConvexity();
            chart = createChartPanel(convexHullDiff, "Convexity", "Time (sec.)", "%");
            break;
        
        case Roughness:
            double[][] roughness = computeRoughness();
            chart = createChartPanel(roughness, "Roughness", "Time (sec.)", "A.U.");
            break;
        
        default:
            System.out.println("Measure " + descriptor.toString() + " not implemented yet");
            return;
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
        // FIXME super.setPanelHeight(355);
        super.panel.updateUI();
    }
    
    private void saveShapeToXLS()
    {
        JFileChooser jfc = new JFileChooser();
        
        if (jfc.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) return;
        
        try
        {
            WritableWorkbook wb = XLSUtil.createWorkbook(jfc.getSelectedFile());
            
            int dim = 3; // FIXME
            
            saveToXls(wb, (dim == 1) ? "Perimeter" : "Surface", computeDimension(1));
            saveToXls(wb, (dim == 1) ? "Area" : "Volume", computeDimension(2));
            saveToXls(wb, "Roundness", computeRoundness());
            saveToXls(wb, "Convexity", computeConvexity());
            
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
    
    // private void exportMeshToXLS()
    // {
    // JFileChooser jfc = new JFileChooser();
    //
    // if (jfc.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) return;
    //
    // try
    // {
    // WritableWorkbook book = XLSUtil.createWorkbook(jfc.getSelectedFile());
    //
    // int cpt = 1;
    // for (Detection det : trackPool.getAllDetection())
    // {
    // if (!(det instanceof Mesh3D)) continue;
    //
    // Mesh3D mesh = (Mesh3D) det;
    //
    // WritableSheet sheet = XLSUtil.createNewPage(book, "Object " + (cpt++));
    // XLSUtil.setCellString(sheet, 0, 0, "X");
    // XLSUtil.setCellString(sheet, 1, 0, "Y");
    // XLSUtil.setCellString(sheet, 2, 0, "Z");
    // XLSUtil.setCellString(sheet, 3, 0, "NX");
    // XLSUtil.setCellString(sheet, 4, 0, "NY");
    // XLSUtil.setCellString(sheet, 5, 0, "NZ");
    //
    // int row = 1;
    //
    // for (Vertex v : mesh.vertices)
    // {
    // if (v == null) continue;
    //
    // XLSUtil.setCellNumber(sheet, 0, row, v.position.x);
    // XLSUtil.setCellNumber(sheet, 1, row, v.position.y);
    // XLSUtil.setCellNumber(sheet, 2, row, v.position.z);
    //
    // XLSUtil.setCellNumber(sheet, 3, row, v.normal.x);
    // XLSUtil.setCellNumber(sheet, 4, row, v.normal.y);
    // XLSUtil.setCellNumber(sheet, 5, row, v.normal.z);
    //
    // row++;
    // }
    // }
    //
    // book.close();
    // }
    // catch (IOException e)
    // {
    // e.printStackTrace();
    // }
    // catch (WriteException e)
    // {
    // e.printStackTrace();
    // }
    // }
    
    private void exportMeshToVTK()
    {
        JFileChooser jfc = new JFileChooser();
        
        if (jfc.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) return;
        
        File f = jfc.getSelectedFile();
        
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
                        BufferedWriter writer = new BufferedWriter(new FileWriter(f.getAbsolutePath() + fileName));
                        ((Mesh3D) det).saveToVTK(writer);
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
    
    //
    // private void exportMeshToOFF()
    // {
    // JFileChooser jfc = new JFileChooser();
    //
    // if (jfc.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) return;
    //
    // File f = jfc.getSelectedFile();
    //
    // int cpt = 1;
    // for (TrackSegment ts : trackPool.getTrackSegmentList())
    // {
    // String filePrefix = "_#" + (cpt < 10 ? "0" : "") + cpt;
    //
    // for (Detection det : ts.getDetectionList())
    // {
    // String fileName = filePrefix + "_T" + (det.getT() < 10 ? "0" : "") + det.getT() + ".off";
    // PrintStream ps;
    // try
    // {
    // ps = new PrintStream(f.getAbsolutePath() + fileName);
    // Mesh c3d = (Mesh) det;
    // c3d.exportToOFF(ps);
    // }
    // catch (FileNotFoundException e)
    // {
    // e.printStackTrace();
    // }
    // }
    //
    // cpt++;
    // }
    // }
    //
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
    
    /**
     * Computes the coefficient of variation of shape for all objects over time. This coefficient is
     * given as the variance of all distances between the mass center and the vertices, normalized
     * over the average mesh radius
     * 
     * @return
     */
    private double[][] computeRadiiVar()
    {
        double[][] result = new double[trackPool.getTrackSegmentList().size()][trackPool.getDisplaySequence().getSizeT()];
        
        for (TrackSegment ts : trackPool.getTrackSegmentList())
        {
            int trackIndex = trackPool.getTrackIndex(ts);
            
            for (Detection det : ts.getDetectionList())
            {
                if (!(det instanceof Mesh3D)) continue;
                
                Mesh3D contour = (Mesh3D) det;
                int n = (int) contour.getDimension(0);
                double[] distances = new double[n];
                int cpt = 0;
                
                Point3d center = contour.getMassCenter();
                
                for (Point3d p : contour)
                {
                    if (p == null) continue;
                    distances[cpt++] = p.distance(center);
                }
                
                result[trackIndex][det.getT()] = ArrayMath.var(distances, true) / ArrayMath.mean(distances);
                
//                // compute the surface to area ratio instead (aka sphericity)
//                
//                double surfaceArea = contour.getDimension(1);
//                double volume = contour.getDimension(2);
//                // Sphericity:
//                // (36 * Pi * Vol^2)^(1/3) / Area 
//                result[trackIndex][det.getT()] = Math.pow(36 * Math.PI * volume * volume, 1.0/3.0) / surfaceArea;
//                
//                // try this: cubic root of (the volume divided by the volume of the circumscribing sphere)
//                // from here: http://homepage.smc.edu/grippo_alessandro/roundness.html
//                
//                double r = contour.getMaxDistanceTo(contour.getMassCenter());
//                double sphVol = 4*Math.PI*r*r*r/3;
//                result[trackIndex][det.getT()] = Math.cbrt(volume / sphVol);
            }
        }
        
        return result;
    }
    
    /**
     * <b>Description: </b>The roughness provides a dimension-less measure of the local
     * (high-frequency) deviations of the surface from a hypothetical "smoothed" version. The term
     * "local" implies that a sampling distance is defined, such that a the roughness is first
     * computed locally for each surface point, while the global roughness measure is given as the
     * average roughness over the entire surface.<br/>
     * <b>Calculation: </b>.<br/>
     * <b>Interpretation: </b>Values close to 0 indicate a smooth surface, while increasing values
     * indicate a rough surface. This measure is dimension-less, however it is locally normalized
     * with respect to the local surface, therefore it can be reliably compared across multiple
     * objects, regardless of their size, as long as the sampling and sampling length are identical.
     * 
     * @param samplingDistance
     *            the distance over which local roughness measures should be computed (typical
     *            values could be 3 to 5 times the surface {@link ActiveContour#sampling})
     * @return
     */
    private double[][] computeRoughness()
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
                
                final Mesh3D contour = (Mesh3D) det;
                final double samplingDistance = contour.sampling.getValue() * 3;
                
                results.add(service.submit(new Runnable()
                {
                    @Override
                    public void run()
                    {
                        Point3d center = contour.getMassCenter();
                        
                        double[] localRoughness = new double[(int) contour.getDimension(0)];
                        
                        int n = 0;
                        
                        // compute local roughness for each mesh point
                        for (Vertex v1 : contour.vertices)
                        {
                            if (v1 == null) continue;
                            
                            // 1) Build a neighborhood for each vertex w.r.t. the sampling distance
                            ArrayList<Vertex> neighborhood = new ArrayList<Vertex>();
                            neighborhood.add(v1);
                            
                            // FIXME sweep neighbors instead of entire mesh to speed up
                            // for (Vertex v2 : contour.vertices)
                            // {
                            // if (v2 == null || v1 == v2 || neighborhood.contains(v2)) continue;
                            //
                            // if (v1.distanceTo(v2) <= samplingDistance) neighborhood.add(v2);
                            // }
                            
                            for (Integer n1 : v1.neighbors)
                            {
                                Vertex vn1 = contour.vertices.get(n1);
                                if (vn1 == null || vn1 == v1 || neighborhood.contains(vn1)) continue;
//                                 if (v1.distanceTo(vn1) <= samplingDistance)
                                {
                                    neighborhood.add(vn1);
                                    for (Integer n2 : vn1.neighbors)
                                    {
                                        Vertex vn2 = contour.vertices.get(n2);
                                        if (vn2 == null || vn2 == v1 || vn2 == vn1 || neighborhood.contains(vn2)) continue;
//                                         if (v1.distanceTo(vn2) <= samplingDistance)
                                        {
                                            neighborhood.add(vn2);
                                            for (Integer n3 : vn2.neighbors)
                                            {
                                                Vertex vn3 = contour.vertices.get(n3);
                                                if (vn3 == null || vn3 == v1 || vn3 == vn1 || vn3 == vn2 || neighborhood.contains(vn3)) continue;
//                                                if (v1.distanceTo(vn3) <= samplingDistance)
                                                {
                                                    double d = samplingDistance;
                                                    neighborhood.add(vn3);
                                                }
                                            }
                                            
                                        }
                                    }
                                    
                                }
                            }
                            
                            // 2) Compute local roughness for the current neighborhood
                            
                            // 2.a) compute the distance from each point to the center
                            double[] distancesToCenter = new double[neighborhood.size()];
                            for (int i = 0; i < distancesToCenter.length; i++)
                                distancesToCenter[i] = neighborhood.get(i).distanceTo(center);
                            // 2.b) subtract their average
                            // ArrayMath.subtract(distancesToCenter,
                            // ArrayMath.mean(distancesToCenter), distancesToCenter);
                            
                            localRoughness[n++] = ArrayMath.std(distancesToCenter, true) / ArrayMath.mean(distancesToCenter);
                        }
                        
                        result[trackIndex][contour.getT()] = ArrayMath.mean(localRoughness);
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
     * <b>Description: </b>The roundness (following the ISO 1101 standard) is based on the ratio
     * between the smallest inscribed and largest circumscribed spheres needed to best fit the
     * interior and exterior of the object.<br/>
     * <b>Calculation: </b>The roundness is approximated here by the ratio between the smallest and
     * largest distance from the mass center to any point on the surface.<br/>
     * <b>Interpretation: </b>A value of 1 corresponds to a perfect sphere, while a value close to 0
     * describes flat objects (e.g. star-shaped objects with flat branches)
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
                Point3d center = contour.getMassCenter();
                double minRadius = contour.getMinDistanceTo(center, null);
                double maxRadius = contour.boundingSphere.getRadius();
                result[trackIndex][det.getT()] = minRadius / maxRadius;
            }
        }
        
        return result;
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
                    hullVolume += surf * cross.x * (v1.x + v2.x + v3.x);
                }
                
                result[trackIndex][det.getT()] = 100 * contour.getDimension(2) / hullVolume;
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
}
