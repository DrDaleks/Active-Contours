package plugins.adufour.activecontours;

import icy.image.IcyBufferedImage;
import icy.main.Icy;
import icy.painter.Overlay;
import icy.painter.Overlay.OverlayPriority;
import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROIUtil;
import icy.sequence.DimensionId;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.swimmingPool.SwimmingObject;
import icy.system.IcyHandledException;
import icy.system.SystemUtil;
import icy.system.thread.ThreadUtil;
import icy.type.DataType;
import icy.util.ShapeUtil.BooleanOperator;
import icy.util.StringUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.vecmath.Point3d;

import plugins.adufour.activecontours.SlidingWindow.Operation;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzException;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDimensionPicker;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.filtering.Convolution1D;
import plugins.adufour.filtering.ConvolutionException;
import plugins.adufour.filtering.Kernels1D;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarROIArray;
import plugins.adufour.vars.util.VarException;
import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackManager;
import plugins.fab.trackmanager.TrackSegment;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.nchenouard.spot.Detection;

public class ActiveContours extends EzPlug implements EzStoppable, Block
{
    private final double                   EPSILON               = 0.0000001;
    
    private final EzVarBoolean             showAdvancedOptions   = new EzVarBoolean("Show advanced options", false);
    
    public final EzVarSequence             input                 = new EzVarSequence("Input");
    private Sequence                       inputData;
    private Sequence                       currentFrame_float;
    
    public final EzVarDouble               init_isovalue         = new EzVarDouble("Isovalue", 1, 0, 1000000, 0.01);
    
    public final EzVarDouble               regul_weight          = new EzVarDouble("Contour smoothness", 0.05, 0, 1.0, 0.01);
    
    public final EzGroup                   edge                  = new EzGroup("Find bright/dark edges");
    public final EzVarDimensionPicker      edge_c                = new EzVarDimensionPicker("Find edges in channel", DimensionId.C, input);
    public final EzVarDouble               edge_weight           = new EzVarDouble("Edge weight", 0, -1, 1, 0.1);
    
    public final EzGroup                   region                = new EzGroup("Find homogeneous intensity areas");
    public final EzVarDimensionPicker      region_c              = new EzVarDimensionPicker("Find regions in channel", DimensionId.C, input);
    public final EzVarDouble               region_weight         = new EzVarDouble("Region weight", 1.0, 0.0, 1.0, 0.1);
    
    public final EzVarDouble               balloon_weight        = new EzVarDouble("Contour inflation", 0, -0.5, 0.5, 0.001);
    
    public final EzVarDouble               axis_weight           = new EzVarDouble("Axis constraint", 0, 0.0, 1, 0.1);
    
    public final EzVarBoolean              coupling_flag         = new EzVarBoolean("Multi-contour coupling", true);
    
    public final EzGroup                   evolution             = new EzGroup("Evolution parameters");
    public final EzVarSequence             evolution_bounds      = new EzVarSequence("Bound field to ROI of");
    public final EzVarDouble               contour_resolution    = new EzVarDouble("Contour resolution", 2, 0.1, 1000.0, 0.1);
    public final EzVarInteger              contour_minArea       = new EzVarInteger("Contour min. area", 10, 1, 100000000, 1);
    public final EzVarDouble               contour_timeStep      = new EzVarDouble("Evolution time step", 0.1, 0.1, 10, 0.01);
    public final EzVarInteger              convergence_winSize   = new EzVarInteger("Convergence window size", 50, 10, 10000, 10);
    public final EzVarEnum<Operation>      convergence_operation = new EzVarEnum<SlidingWindow.Operation>("Convergence operation", Operation.values(), Operation.VAR_COEFF);
    public final EzVarDouble               convergence_criterion = new EzVarDouble("Convergence criterion", 0.001, 0, 0.1, 0.001);
    
    public final EzVarBoolean              output_rois           = new EzVarBoolean("Regions of interest (ROI)", true);
    
    public final EzVarBoolean              tracking              = new EzVarBoolean("Track objects over time", false);
    
    private final Sequence                 edgeData              = new Sequence("Edge information");
    private final Sequence                 contourMask_buffer    = new Sequence("Mask");
    
    private Sequence                       region_data;
    private HashMap<ActiveContour, Double> region_cin            = new HashMap<ActiveContour, Double>(0);
    private double                         region_cout;
    
    private VarROIArray                    roiInput              = new VarROIArray("input ROI");
    private VarROIArray                    roiOutput             = new VarROIArray("Regions of interest");
    
    private boolean                        globalStop;
    
    private TrackGroup                     trackGroup;
    
    private ActiveContoursOverlay          overlay;
    
    private ExecutorService                multiThreadService    = Executors.newFixedThreadPool(SystemUtil.getAvailableProcessors());
    
    public TrackGroup getTrackGroup()
    {
        return trackGroup;
    }
    
    @Override
    public void initialize()
    {
        addEzComponent(showAdvancedOptions);
        
        addEzComponent(input);
        
        // regul
        regul_weight.setToolTipText("Higher values result in a smoother contour, but may also slow its growth");
        addEzComponent(regul_weight);
        
        // edge
        edge.setToolTipText("Sets the contour(s) to follow image intensity gradients");
        edge_weight.setToolTipText("Negative (resp. positive) weight pushes contours toward decreasing (resp. increasing) intensities");
        edge.addEzComponent(edge_c, edge_weight);
        addEzComponent(edge);
        
        // region
        region.setToolTipText("Sets the contour(s) to isolate homogeneous intensity regions");
        region_weight.setToolTipText("Set to 0 to deactivate this parameter");
        region.addEzComponent(region_c, region_weight);
        addEzComponent(region);
        
        // balloon force
        balloon_weight.setToolTipText("Positive (resp. negative) values will inflate (resp. deflate) the contour");
        addEzComponent(balloon_weight);
        
        // axis contraint
        axis_weight.setToolTipText("Higher values restrict the evolution along the principal axis");
        addEzComponent(axis_weight);
        
        // coupling
        coupling_flag.setToolTipText("Prevents multiple contours from overlapping");
        addEzComponent(coupling_flag);
        
        // contour
        contour_resolution.setToolTipText("Sets the contour(s) precision as the distance (in pixels) between control points");
        
        contour_minArea.setToolTipText("Contours with a surface (in pixels) below this value will be removed");
        showAdvancedOptions.addVisibilityTriggerTo(contour_minArea, true);
        
        contour_timeStep.setToolTipText("Defines the evolution speed (warning: keep a low value to avoid vibration effects)");
        
        convergence_winSize.setToolTipText("Defines over how many iterations the algorithm should check for convergence");
        showAdvancedOptions.addVisibilityTriggerTo(convergence_winSize, true);
        
        convergence_operation.setToolTipText("Defines the operation used to detect convergence");
        showAdvancedOptions.addVisibilityTriggerTo(convergence_operation, true);
        
        convergence_criterion.setToolTipText("Defines the value of the criterion used to detect convergence");
        
        evolution_bounds.setNoSequenceSelection();
        evolution_bounds.setToolTipText("Bounds the evolution of the contour to all ROI of the given sequence (select \"No sequence\" to deactivate)");
        showAdvancedOptions.addVisibilityTriggerTo(evolution_bounds, true);
        
        evolution.addEzComponent(evolution_bounds, contour_resolution, contour_minArea, contour_timeStep, convergence_winSize, convergence_operation, convergence_criterion);
        addEzComponent(evolution);
        
        contour_resolution.addVarChangeListener(new EzVarListener<Double>()
        {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue)
            {
                convergence_winSize.setValue((int) (100.0 / newValue));
            }
        });
        
        // output
        output_rois.setToolTipText("Clone the original sequence and with results overlayed as ROIs");
        addEzComponent(output_rois);
        
        tracking.setToolTipText("Track objects over time and export results to the track manager");
        addEzComponent(tracking);
        
        setTimeDisplay(true);
    }
    
    @Override
    public void execute()
    {
        roiOutput.setValue(null);
        inputData = input.getValue(true);
        
        globalStop = false;
        
        int startT = inputData.getFirstViewer() == null ? 0 : inputData.getFirstViewer().getPositionT();
        int endT = tracking.getValue() ? inputData.getSizeT() - 1 : startT;
        
        ThreadUtil.invokeNow(new Runnable()
        {
            @Override
            public void run()
            {
                trackGroup = new TrackGroup(inputData);
                trackGroup.setDescription("Active contours (" + new Date().toString() + ")");
                if (tracking.getValue())
                {
                    SwimmingObject object = new SwimmingObject(trackGroup);
                    Icy.getMainInterface().getSwimmingPool().add(object);
                }
            }
        });
        
        // replace any ActiveContours Painter object on the sequence by ours
        for (Overlay overlay : inputData.getOverlays())
            if (overlay instanceof ActiveContoursOverlay) inputData.removeOverlay(overlay);
        
        overlay = new ActiveContoursOverlay(trackGroup);
        overlay.setPriority(OverlayPriority.TOPMOST);
        inputData.addOverlay(overlay);
        
        if (getUI() != null)
        {
            roiInput.setValue(new ROI[0]);
            
            if (inputData.getFirstViewer() != null)
            {
                startT = inputData.getFirstViewer().getPositionT();
            }
        }
        
        for (int t = startT; t <= endT; t++)
        {
            if (inputData.getFirstViewer() != null) inputData.getFirstViewer().setPositionT(t);
            
            // initialize contours
            initContours(t, t == startT);
            
            if (Thread.currentThread().isInterrupted()) break;
            
            // if (firstRun)
            // {
            // // the thread pool now is warmed up
            // // and the JIT did its business
            // // => restart at full speed
            // firstRun = false;
            // execute();
            // return;
            // }
            
            // evolve contours on the current image
            evolveContours(t);
            
            if (Thread.currentThread().isInterrupted()) break;
            
            // store detections and results
            storeResult(t);
            
            if (Thread.currentThread().isInterrupted()) break;
        }
        
        if (getUI() != null)
        {
            getUI().setProgressBarValue(0.0);
            
            if (output_rois.getValue())
            {
                for (ROI roi : roiOutput.getValue())
                    inputData.addROI(roi, false);
            }
            
            if (tracking.getValue() && !isHeadLess())
            {
                ThreadUtil.invokeLater(new Runnable()
                {
                    public void run()
                    {
                        TrackManager tm = new TrackManager();
                        tm.reOrganize();
                        tm.setDisplaySequence(inputData);
                    }
                });
            }
        }
        else
        {
            // possibly block mode, remove the painter after processing
            // if (inputData != null) inputData.removePainter(painter);
        }
    }
    
    private void initContours(final int t, boolean isFirstImage)
    {
        // Viewer v = inputData.getFirstViewer();
        //
        // int z = (v == null) ? 0 : inputData.getFirstViewer().getPositionZ();
        region_cin.clear();
        region_cout = 0.0;
        
        currentFrame_float = SequenceUtil.convertToType(SequenceUtil.extractFrame(inputData, t), DataType.FLOAT, true, true);
        
        // 1) Initialize the edge data
        {
            if (edge_c.getValue() >= currentFrame_float.getSizeC())
            {
                throw new IcyHandledException("The selected edge channel is invalid.");
            }
            
            if (region_c.getValue() >= currentFrame_float.getSizeC())
            {
                throw new IcyHandledException("The selected region channel is valid.");
            }
            
            Sequence gradX = SequenceUtil.getCopy(inputData.getSizeC() == 1 ? currentFrame_float : SequenceUtil.extractChannel(currentFrame_float, edge_c.getValue()));
            
            Sequence gradient = Kernels1D.GRADIENT.toSequence();
            Sequence gaussian = Kernels1D.CUSTOM_GAUSSIAN.createGaussianKernel1D(0.5).toSequence();
            
            // smooth the signal first
            try
            {
                Convolution1D.convolve(gradX, gaussian, gaussian, null);
            }
            catch (ConvolutionException e)
            {
                throw new EzException("Cannot smooth the signal: " + e.getMessage(), true);
            }
            
            // clone into gradY
            Sequence gradY = SequenceUtil.getCopy(gradX);
            Sequence gradZ = (inputData.getSizeZ() == 1) ? null : SequenceUtil.getCopy(gradX);
            
            // compute the gradient in each direction
            try
            {
                Convolution1D.convolve(gradX, gradient, null, null);
                Convolution1D.convolve(gradY, null, gradient, null);
                if (gradZ != null) Convolution1D.convolve(gradZ, null, null, gradient);
            }
            catch (ConvolutionException e)
            {
                throw new EzException("Cannot compute the gradient information: " + e.getMessage(), true);
            }
            
            // combine the edge data as multi-channel data
            if (gradZ == null)
            {
                float[][] edgeSlice = new float[][] { gradX.getDataXYAsFloat(0, 0, 0), gradY.getDataXYAsFloat(0, 0, 0) };
                edgeData.setImage(0, 0, new IcyBufferedImage(inputData.getWidth(), inputData.getHeight(), edgeSlice));
            }
            else for (int z = 0; z < edgeData.getSizeZ(); z++)
            {
                float[][] edgeSlice = new float[][] { gradX.getDataXYAsFloat(0, z, 0), gradY.getDataXYAsFloat(0, z, 0), gradZ.getDataXYAsFloat(0, z, 0) };
                edgeData.setImage(0, z, new IcyBufferedImage(inputData.getWidth(), inputData.getHeight(), edgeSlice));
            }
        }
        
        // 2) initialize the region data
        {
            region_data = (inputData.getSizeC() == 1) ? currentFrame_float : SequenceUtil.extractChannel(currentFrame_float, region_c.getValue());
            
            // initialize the mask buffer (used to calculate average intensities inside/outside
            if (isFirstImage)
            {
                for (int z = 0; z < inputData.getSizeZ(); z++)
                    contourMask_buffer.setImage(0, z, new IcyBufferedImage(inputData.getWidth(), inputData.getHeight(), 1, DataType.UBYTE));
            }
        }
        
        // 3) Initialize the contours
        if (isFirstImage)
        {
            if (roiInput.getValue().length == 0)
            {
                if (isHeadLess()) throw new VarException("Active contours: no input ROI");
                
                ArrayList<ROI2D> roiFromSequence = inputData.getROI2Ds();
                
                if (roiFromSequence.isEmpty()) throw new EzException("Please draw or select a ROI", true);
                
                roiInput.setValue(roiFromSequence.toArray(new ROI2D[roiFromSequence.size()]));
            }
            
            ArrayList<Future<?>> tasks = new ArrayList<Future<?>>(roiInput.getValue().length);
            
            for (ROI roi : roiInput.getValue())
            {
                if (!(roi instanceof ROI2D))
                {
                    System.err.println("Warning: skipped non-2D ROI");
                    continue;
                }
                
                final ROI2D roi2d = (ROI2D) roi;
                
                Runnable initializer = new Runnable()
                {
                    public void run()
                    {
                        if (roi2d instanceof ROI2DArea)
                        {
                            // special case: check if the area has multiple components => split them
                            ROI2DArea area = (ROI2DArea) roi2d;
                            
                            BooleanMask2D[] components = area.getBooleanMask(true).getComponents();
                            
                            for (BooleanMask2D comp : components)
                            {
                                ROI2DArea roi = new ROI2DArea(comp);
                                final ActiveContour contour = new Polygon2D(ActiveContours.this, contour_resolution, contour_minArea, new SlidingWindow(
                                        convergence_winSize.getValue()), roi);
                                contour.setX(roi.getBounds2D().getCenterX());
                                contour.setY(roi.getBounds2D().getCenterY());
                                contour.setT(t);
                                TrackSegment segment = new TrackSegment();
                                segment.addDetection(contour);
                                synchronized (trackGroup)
                                {
                                    trackGroup.addTrackSegment(segment);
                                }
                                synchronized (region_cin)
                                {
                                    region_cin.put(contour, 0.0);
                                }
                            }
                        }
                        else
                        {
                            final ActiveContour contour = new Polygon2D(ActiveContours.this, contour_resolution, contour_minArea,
                                    new SlidingWindow(convergence_winSize.getValue()), roi2d);
                            contour.setX(roi2d.getBounds2D().getCenterX());
                            contour.setY(roi2d.getBounds2D().getCenterY());
                            contour.setT(t);
                            TrackSegment segment = new TrackSegment();
                            segment.addDetection(contour);
                            synchronized (trackGroup)
                            {
                                trackGroup.addTrackSegment(segment);
                            }
                            synchronized (region_cin)
                            {
                                region_cin.put(contour, 0.0);
                            }
                        }
                    }
                };
                
                tasks.add(multiThreadService.submit(initializer));
            }
            
            try
            {
                for (Future<?> future : tasks)
                    future.get();
            }
            catch (InterruptedException e)
            {
                // restore the interrupted flag
                Thread.currentThread().interrupt();
                return;
            }
            catch (Exception e)
            {
                // stopExecution();
                if (e.getCause() instanceof EzException) throw (EzException) e.getCause();
                e.printStackTrace();
                throw new RuntimeException(e.getCause());
            }
        }
        else
        {
            for (TrackSegment segment : trackGroup.getTrackSegmentList())
            {
                Detection previous = segment.getDetectionAtTime(t - 1);
                
                if (previous == null) continue;
                
                ActiveContour clone = ((ActiveContour) previous).clone();
                clone.setT(t);
                segment.addDetection(clone);
            }
        }
    }
    
    private void evolveContours(final int t)
    {
        // retrieve the contours on the current frame and store them in currentContours
        final HashSet<ActiveContour> allContours = new HashSet<ActiveContour>(trackGroup.getTrackSegmentList().size());
        
        for (TrackSegment segment : trackGroup.getTrackSegmentList())
        {
            Detection det = segment.getDetectionAtTime(t);
            if (det != null) allContours.add((ActiveContour) det);
        }
        
        if (allContours.size() == 0) return;
        
        // get the bounded field of evolution
        ROI field;
        
        Sequence boundSource = evolution_bounds.getValue();
        
        if (boundSource == null || !boundSource.getDimension2D().equals(inputData.getDimension2D()) || boundSource.getROI2Ds().size() == 0)
        {
            field = new ROI2DRectangle(0, 0, inputData.getWidth(), inputData.getHeight());
        }
        else
        {
            field = ROIUtil.merge(boundSource.getROIs(), BooleanOperator.OR);
        }
        
        int nbConvergedContours = 0;
        
        long iter = 0;
        
        final HashSet<ActiveContour> evolvingContours = new HashSet<ActiveContour>(allContours.size());
        
        while (!globalStop && nbConvergedContours < allContours.size())
        {
            nbConvergedContours = 0;
            
            // update region information every 10 iterations
            if (region_weight.getValue() > EPSILON && iter % 10 == 0) updateRegionInformation(allContours, t);
            
            // take a snapshot of the current list of evolving (i.e. non-converged) contours
            evolvingContours.clear();
            for (ActiveContour contour : allContours)
            {
                Double criterion = contour.convergence.computeCriterion(convergence_operation.getValue());
                
                if (criterion != null && criterion < convergence_criterion.getValue() / 10)
                {
                    nbConvergedContours++;
                    continue;
                }
                
                // if the contour hasn't converged yet, store it for the main loop
                evolvingContours.add(contour);
            }
            
            if (getUI() != null) getUI().setProgressBarValue((double) nbConvergedContours / allContours.size());
            
            if (evolvingContours.size() == 0) break;
            
            // re-sample the contours to ensure homogeneous resolution
            resampleContours(evolvingContours, allContours, t);
            
            // compute deformations issued from the energy minimization
            deformContours(evolvingContours, allContours, field);
            
            // compute energy
            // computeEnergy(mainService, allContours);
            
            iter++;
            
            if (Thread.currentThread().isInterrupted()) globalStop = true;
            
            overlay.painterChanged();
        }
    }
    
    /**
     * Deform contours together (coupling involved)
     * 
     * @param service
     * @param evolvingContours
     * @param allContours
     */
    private void deformContours(final HashSet<ActiveContour> evolvingContours, final HashSet<ActiveContour> allContours, final ROI field)
    {
        
        if (evolvingContours.size() == 1 && allContours.size() == 1)
        {
            // no multi-threading needed
            
            ActiveContour contour = evolvingContours.iterator().next();
            
            if (Math.abs(edge_weight.getValue()) > EPSILON) contour.computeEdgeForces(edge_weight.getValue(), edgeData);
            
            if (regul_weight.getValue() > EPSILON) contour.computeInternalForces(regul_weight.getValue());
            
            if (region_weight.getValue() > EPSILON) contour.computeRegionForces(region_data, region_weight.getValue(), region_cin.get(contour), region_cout);
            
            if (axis_weight.getValue() > EPSILON) contour.computeAxisForces(axis_weight.getValue());
            
            if (Math.abs(balloon_weight.getValue()) > EPSILON) contour.computeBalloonForces(balloon_weight.getValue());
            
            contour.move(field, contour_timeStep.getValue());
        }
        else
        {
            final double wRegul = regul_weight.getValue();
            final double wEdge = edge_weight.getValue();
            final double wRegion = region_weight.getValue();
            final double wAxis = axis_weight.getValue();
            final double wBall = balloon_weight.getValue();
            
            ArrayList<Future<ActiveContour>> tasks = new ArrayList<Future<ActiveContour>>(evolvingContours.size());
            
            for (final ActiveContour contour : evolvingContours)
            {
                tasks.add(multiThreadService.submit(new Callable<ActiveContour>()
                {
                    public ActiveContour call()
                    {
                        if (wRegul > EPSILON) contour.computeInternalForces(wRegul);
                        
                        if (Math.abs(wEdge) > EPSILON) contour.computeEdgeForces(wEdge, edgeData);
                        
                        if (wRegion > EPSILON) contour.computeRegionForces(region_data, wRegion, region_cin.get(contour), region_cout);
                        
                        if (wAxis > EPSILON) contour.computeAxisForces(wAxis);
                        
                        if (Math.abs(wBall) > EPSILON) contour.computeBalloonForces(wBall);
                        
                        if (coupling_flag.getValue())
                        {
                            // Don't move the contours just now: coupling feedback must be computed
                            // against ALL contours (including those which have already converged)
                            for (ActiveContour otherContour : allContours)
                            {
                                if (otherContour == null || otherContour == contour) continue;
                                
                                contour.computeFeedbackForces(otherContour);
                            }
                        }
                        else
                        {
                            // move contours asynchronously
                            contour.move(field, contour_timeStep.getValue());
                        }
                        
                        return contour;
                    }
                }));
            }
            
            try
            {
                for (Future<ActiveContour> future : tasks)
                    future.get();
            }
            catch (InterruptedException e)
            {
                // reset the interrupted flag
                Thread.currentThread().interrupt();
                return;
            }
            catch (ExecutionException e)
            {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
            
            if (coupling_flag.getValue())
            {
                // motion is synchronous, and can be done now
                for (ActiveContour contour : evolvingContours)
                    contour.move(field, contour_timeStep.getValue());
            }
        }
    }
    
    private void resampleContours(final HashSet<ActiveContour> evolvingContours, final HashSet<ActiveContour> allContours, final int t)
    {
        final VarBoolean loop = new VarBoolean("loop", true);
        
        final VarBoolean change = new VarBoolean("change", false);
        
        while (loop.getValue())
        {
            loop.setValue(false);
            
            if (evolvingContours.size() == 1)
            {
                // no multi-threading needed
                
                ActiveContour contour = evolvingContours.iterator().next();
                ReSampler reSampler = new ReSampler(trackGroup, contour, evolvingContours, allContours);
                if (reSampler.call())
                {
                    change.setValue(true);
                    loop.setValue(true);
                }
            }
            else
            {
                ArrayList<Future<Boolean>> tasks = new ArrayList<Future<Boolean>>(evolvingContours.size());
                
                for (final ActiveContour contour : evolvingContours)
                    tasks.add(multiThreadService.submit(new ReSampler(trackGroup, contour, evolvingContours, allContours)));
                
                try
                {
                    for (Future<Boolean> resampled : tasks)
                    {
                        if (resampled.get())
                        {
                            change.setValue(true);
                            loop.setValue(true);
                        }
                    }
                }
                catch (InterruptedException e)
                {
                    // reset the interrupted flag
                    Thread.currentThread().interrupt();
                    return;
                }
                catch (ExecutionException e)
                {
                    throw new RuntimeException(e);
                }
                catch (RuntimeException e)
                {
                    throw e;
                }
                finally
                {
                    tasks.clear();
                }
            }
        }
    }
    
    private void updateRegionInformation(Collection<ActiveContour> contours, int t)
    {
        int sizeZ = inputData.getSizeZ();
        
        ArrayList<Future<?>> tasks = new ArrayList<Future<?>>(contours.size());
        
        for (final ActiveContour contour : contours)
            tasks.add(multiThreadService.submit(new Runnable()
            {
                public void run()
                {
                    region_cin.put(contour, contour.computeAverageIntensity(region_data, contourMask_buffer));
                }
            }));
        
        try
        {
            for (Future<?> task : tasks)
                task.get();
        }
        catch (InterruptedException e)
        {
            // reset the interrupted flag
            Thread.currentThread().interrupt();
            return;
        }
        catch (ExecutionException e)
        {
            throw new RuntimeException(e);
        }
        
        double outSum = 0, outCpt = 0;
        for (int z = 0; z < sizeZ; z++)
        {
            byte[] _mask = contourMask_buffer.getDataXYAsByte(0, z, 0);
            float[] _data = region_data.getDataXYAsFloat(0, z, 0);
            
            for (int i = 0; i < _mask.length; i++)
            {
                if (_mask[i] == 0)
                {
                    double value = _data[i];
                    outSum += value;
                    outCpt++;
                }
                else _mask[i] = 0;
            }
        }
        region_cout = outSum / outCpt;
    }
    
    private void storeResult(int t)
    {
        ArrayList<TrackSegment> segments = trackGroup.getTrackSegmentList();
        
        ArrayList<ROI> rois = null;
        if (output_rois.getValue()) rois = new ArrayList<ROI>(Arrays.asList(roiOutput.getValue()));
        
        for (int i = 1; i <= segments.size(); i++)
        {
            TrackSegment segment = segments.get(i - 1);
            
            ActiveContour contour = (ActiveContour) segment.getDetectionAtTime(t);
            if (contour == null) continue;
            
            // store detection parameters
            Point3d center = new Point3d();
            for (Point3d p : contour)
                center.add(p);
            center.scale(1.0 / contour.getDimension(0));
            contour.setX(center.x);
            contour.setY(center.y);
            
            // output as ROIs
            if (output_rois.getValue())
            {
                ROI roi = contour.toROI();
                roi.setColor(contour.getColor());
                roi.setName("[T=" + StringUtil.toString(t, 1 + (int) Math.round(Math.log10(inputData.getSizeT()))) + "] Object #" + i);
                rois.add(roi);
            }
        }
        
        if (output_rois.getValue() && rois.size() > 0) roiOutput.setValue(rois.toArray(new ROI2D[rois.size()]));
    }
    
    @Override
    public void clean()
    {
        if (inputData != null) inputData.removeOverlay(overlay);
        
        // contoursMap.clear();
        // contours.clear();
        // trackGroup.clearTracks();
        if (region_weight.getValue() > EPSILON) region_cin.clear();
        
        // meanUpdateService.shutdownNow();
        multiThreadService.shutdownNow();
    }
    
    @Override
    public void stopExecution()
    {
        globalStop = true;
    }
    
    @Override
    public void declareInput(VarList inputMap)
    {
        inputMap.add("input sequence", input.getVariable());
        inputMap.add("Input ROI", roiInput);
        inputMap.add("regularization: weight", regul_weight.getVariable());
        inputMap.add("edge: weight", edge_weight.getVariable());
        edge_c.setActive(false);
        edge_c.setValues(0, 0, 16, 1);
        inputMap.add("edge: channel", edge_c.getVariable());
        inputMap.add("region: weight", region_weight.getVariable());
        region_c.setActive(false);
        region_c.setValues(0, 0, 16, 1);
        inputMap.add("region: channel", region_c.getVariable());
        
        inputMap.add("balloon: weight", balloon_weight.getVariable());
        
        coupling_flag.setValue(true);
        inputMap.add("contour resolution", contour_resolution.getVariable());
        contour_resolution.addVarChangeListener(new EzVarListener<Double>()
        {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue)
            {
                convergence_winSize.setValue((int) (100.0 / newValue));
            }
        });
        
        // inputMap.add("minimum object size", contour_minArea.getVariable());
        inputMap.add("region bounds", evolution_bounds.getVariable());
        inputMap.add("time step", contour_timeStep.getVariable());
        // inputMap.add("convergence window size", convergence_winSize.getVariable());
        inputMap.add("convergence value", convergence_criterion.getVariable());
        output_rois.setValue(true);
        inputMap.add("tracking", tracking.getVariable());
    }
    
    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add(roiOutput);
    }
    
}
