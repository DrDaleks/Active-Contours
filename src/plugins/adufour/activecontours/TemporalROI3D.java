package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.canvas.IcyCanvas2D;
import icy.canvas.IcyCanvas3D;
import icy.painter.OverlayEvent;
import icy.painter.OverlayListener;
import icy.roi.ROI;
import icy.roi.ROI3D;
import icy.roi.ROI4D;
import icy.roi.ROIEvent;
import icy.roi.ROIListener;
import icy.sequence.Sequence;
import icy.system.IcyExceptionHandler;
import icy.type.point.Point5D;
import icy.type.rectangle.Rectangle3D;
import icy.type.rectangle.Rectangle4D;
import icy.util.XMLUtil;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.Semaphore;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * A temporal sequence of {@link ROI3D} objects.
 * 
 * @author Alexandre Dufour
 */
class TemporalROI3D extends ROI4D implements ROIListener, OverlayListener, Iterable<ROI3D>
{
    public static final String              PROPERTY_USECHILDCOLOR = "useChildColor";
    
    protected final TreeMap<Integer, ROI3D> roiSequence            = new TreeMap<Integer, ROI3D>();
    
    protected final Class<? extends ROI3D>  roiClass;
    protected boolean                       useChildColor;
    protected Semaphore                     modifyingSlice;
    
    /**
     * Creates a new 3D ROI based on the given 2D ROI type.
     */
    public TemporalROI3D(Class<? extends ROI3D> roiClass)
    {
        super();
        
        this.roiClass = roiClass;
        useChildColor = false;
        modifyingSlice = new Semaphore(1);
    }
    
    @Override
    protected ROIPainter createPainter()
    {
        return new ROI3DSequencePainter();
    }
    
    /**
     * Create a new empty 2D ROI slice.
     */
    protected ROI3D createSlice()
    {
        try
        {
            return roiClass.newInstance();
        }
        catch (Exception e)
        {
            IcyExceptionHandler.showErrorMessage(e, true, true);
            return null;
        }
    }
    
    /**
     * Returns <code>true</code> if the ROI directly uses the 2D slice color draw property and
     * <code>false</code> if it uses the global 3D ROI color draw property.
     */
    public boolean getUseChildColor()
    {
        return useChildColor;
    }
    
    /**
     * Set to <code>true</code> if you want to directly use the 2D slice color draw property and
     * <code>false</code> to keep the global 3D ROI color draw property.
     * 
     * @see #setColor(int, Color)
     */
    public void setUseChildColor(boolean value)
    {
        if (useChildColor != value)
        {
            useChildColor = value;
            propertyChanged(PROPERTY_USECHILDCOLOR);
            // need to redraw it
            getOverlay().painterChanged();
        }
    }
    
    /**
     * Set the painter color for the specified ROI slice.
     * 
     * @see #setUseChildColor(boolean)
     */
    public void setColor(int t, Color value)
    {
        final ROI3D slice = getSlice(t);
        
        modifyingSlice.acquireUninterruptibly();
        try
        {
            if (slice != null) slice.setColor(value);
        }
        finally
        {
            modifyingSlice.release();
        }
    }
    
    @Override
    public void setColor(Color value)
    {
        beginUpdate();
        try
        {
            super.setColor(value);
            
            if (!getUseChildColor())
            {
                modifyingSlice.acquireUninterruptibly();
                try
                {
                    for (ROI3D slice : roiSequence.values())
                        slice.setColor(value);
                }
                finally
                {
                    modifyingSlice.release();
                }
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setOpacity(float value)
    {
        beginUpdate();
        try
        {
            super.setOpacity(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setOpacity(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setStroke(double value)
    {
        beginUpdate();
        try
        {
            super.setStroke(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setStroke(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setCreating(boolean value)
    {
        beginUpdate();
        try
        {
            super.setCreating(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setCreating(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setReadOnly(boolean value)
    {
        beginUpdate();
        try
        {
            super.setReadOnly(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setReadOnly(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setFocused(boolean value)
    {
        beginUpdate();
        try
        {
            super.setFocused(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setFocused(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setSelected(boolean value)
    {
        beginUpdate();
        try
        {
            super.setSelected(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setSelected(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    @Override
    public void setC(int value)
    {
        beginUpdate();
        try
        {
            super.setC(value);
            
            modifyingSlice.acquireUninterruptibly();
            try
            {
                for (ROI3D slice : roiSequence.values())
                    slice.setC(value);
            }
            finally
            {
                modifyingSlice.release();
            }
        }
        finally
        {
            endUpdate();
        }
    }
    
    /**
     * Returns <code>true</code> if the ROI stack is empty.
     */
    public boolean isEmpty()
    {
        return roiSequence.isEmpty();
    }
    
    /**
     * @return The size of this ROI stack along Z.<br>
     *         Note that the returned value indicates the difference between upper and lower bounds
     *         of this ROI, but doesn't guarantee that all slices in-between exist (
     *         {@link #getSlice(int)} may still return <code>null</code>.<br>
     */
    public int getSizeT()
    {
        if (roiSequence.isEmpty()) return 0;
        
        return (roiSequence.lastKey().intValue() - roiSequence.firstKey().intValue()) + 1;
    }
    
    /**
     * Returns the ROI slice at given T position.
     */
    public ROI3D getSlice(int t)
    {
        return roiSequence.get(Integer.valueOf(t));
    }
    
    /**
     * Returns the ROI slice at given T position.
     */
    public ROI3D getSlice(int t, boolean createIfNull)
    {
        ROI3D result = getSlice(t);
        
        if ((result == null) && createIfNull)
        {
            result = createSlice();
            if (result != null) setTimePoint(t, result);
        }
        
        return result;
    }
    
    /**
     * Sets the slice for the given Z position.
     */
    public void setTimePoint(int t, ROI3D roi)
    {
        // set Z, T and C position
        roi.setT(t);
        roi.setC(getC());
        // listen events from this ROI and its overlay
        roi.addListener(this);
        roi.getOverlay().addOverlayListener(this);
        
        roiSequence.put(Integer.valueOf(t), roi);
        
        // notify ROI changed
        roiChanged();
    }
    
    /**
     * Removes slice at the given Z position and returns it.
     */
    public ROI3D removeSlice(int z)
    {
        // remove the current slice (if any)
        final ROI3D result = roiSequence.remove(Integer.valueOf(z));
        
        // remove listeners
        if (result != null)
        {
            result.removeListener(this);
            result.getOverlay().removeOverlayListener(this);
        }
        
        // notify ROI changed
        roiChanged();
        
        return result;
    }
    
    /**
     * Removes all slices.
     */
    public void clear()
    {
        for (ROI3D slice : roiSequence.values())
        {
            slice.removeListener(this);
            slice.getOverlay().removeOverlayListener(this);
        }
        
        roiSequence.clear();
    }
    
    /**
     * Called when a ROI slice has changed.
     */
    protected void sliceChanged(ROIEvent event)
    {
        if (modifyingSlice.availablePermits() <= 0) return;
        
        final ROI source = event.getSource();
        
        switch (event.getType())
        {
        case ROI_CHANGED:
            roiChanged();
            break;
        
        case FOCUS_CHANGED:
            setFocused(source.isFocused());
            break;
        
        case SELECTION_CHANGED:
            setSelected(source.isSelected());
            break;
        
        case PROPERTY_CHANGED:
            final String propertyName = event.getPropertyName();
            
            if ((propertyName == null) || propertyName.equals(PROPERTY_READONLY)) setReadOnly(source.isReadOnly());
            if ((propertyName == null) || propertyName.equals(PROPERTY_CREATING)) setCreating(source.isCreating());
            break;
        default:
            break;
        }
    }
    
    /**
     * Called when a ROI slice overlay has changed.
     */
    protected void sliceOverlayChanged(OverlayEvent event)
    {
        switch (event.getType())
        {
        case PAINTER_CHANGED:
            // forward the event to ROI stack overlay
            getOverlay().painterChanged();
            break;
        
        case PROPERTY_CHANGED:
            // forward the event to ROI stack overlay
            getOverlay().propertyChanged(event.getPropertyName());
            break;
        }
    }
    
    @Override
    public Rectangle4D computeBounds4D()
    {
        Rectangle3D r3 = null;
        
        for (ROI3D roi : roiSequence.values())
        {
            final Rectangle3D bnd3d = roi.getBounds3D();
            
            // only add non empty bounds
            if (!bnd3d.isEmpty())
            {
                if (r3 == null) r3 = (Rectangle3D) bnd3d.clone();
                else r3.add(bnd3d);
            }
        }
        
        // create empty 2D bounds
        if (r3 == null) r3 = new Rectangle3D.Double();
        
        final int minT;
        final int sizeT;
        
        if (!roiSequence.isEmpty())
        {
            minT = roiSequence.firstKey().intValue();
            sizeT = getSizeT();
        }
        else
        {
            minT = 0;
            sizeT = 0;
        }
        
        return new Rectangle4D.Double(r3.getX(), r3.getY(), r3.getZ(), minT, r3.getSizeX(), r3.getSizeY(), r3.getSizeZ(), sizeT);
    }
    
    @Override
    public boolean contains(double x, double y, double z, double t)
    {
        final ROI3D roi = getSlice((int) t);
        
        if (roi != null) return roi.contains(x, y, z);
        
        return false;
    }
    
    @Override
    public boolean contains(double x, double y, double z, double t, double sizeX, double sizeY, double sizeZ, double sizeT)
    {
        final Rectangle4D bounds = getBounds4D();
        
        // easy discard
        if (!bounds.contains(x, y, z, t, sizeX, sizeY, sizeZ, sizeT)) return false;
        
        for (int localT = (int) t; localT < (int) (t + sizeZ); localT++)
        {
            final ROI3D roi = getSlice(localT);
            if (roi == null || !roi.contains(x, y, z, sizeX, sizeY, sizeZ)) return false;
        }
        
        return true;
    }
    
    @Override
    public boolean intersects(double x, double y, double z, double t, double sizeX, double sizeY, double sizeZ, double sizeT)
    {
        final Rectangle4D bounds = getBounds4D();
        
        // easy discard
        if (!bounds.intersects(x, y, z, t, sizeX, sizeY, sizeZ, sizeT)) return false;
        
        for (int localT = (int) t; localT < (int) (t + sizeT); localT++)
        {
            final ROI3D roi = getSlice(localT);
            if (roi != null && roi.intersects(x, y, z, sizeX, sizeY, sizeZ)) return true;
        }
        
        return false;
    }
    
    @Override
    public boolean hasSelectedPoint()
    {
        // default
        return false;
    }
    
    @Override
    public double computeNumberOfContourPoints()
    {
        // 3D edge points = first slice points + inter slices edge points + last slice points
        // TODO: only approximation, fix this to use real 3D edge point
        double perimeter = 0;
        
        if (roiSequence.size() <= 2)
        {
            for (ROI3D slice : roiSequence.values())
                perimeter += slice.getNumberOfPoints();
        }
        else
        {
            final Entry<Integer, ROI3D> firstEntry = roiSequence.firstEntry();
            final Entry<Integer, ROI3D> lastEntry = roiSequence.lastEntry();
            final Integer firstKey = firstEntry.getKey();
            final Integer lastKey = lastEntry.getKey();
            
            perimeter = firstEntry.getValue().getNumberOfPoints();
            
            for (ROI3D slice : roiSequence.subMap(firstKey, false, lastKey, false).values())
                perimeter += slice.getNumberOfContourPoints();
            
            perimeter += lastEntry.getValue().getNumberOfPoints();
        }
        
        return perimeter;
    }
    
    @Override
    public double computeNumberOfPoints()
    {
        double volume = 0;
        
        for (ROI3D slice : roiSequence.values())
            volume += slice.getNumberOfPoints();
        
        return volume;
    }
    
    @Override
    public boolean canTranslate()
    {
        // only need to test the first entry
        if (!roiSequence.isEmpty()) return roiSequence.firstEntry().getValue().canTranslate();
        
        return false;
    }
    
    // called when one of the slice ROI changed
    @Override
    public void roiChanged(ROIEvent event)
    {
        // propagate children change event
        sliceChanged(event);
    }
    
    // called when one of the slice ROI overlay changed
    @Override
    public void overlayChanged(OverlayEvent event)
    {
        // propagate children overlay change event
        sliceOverlayChanged(event);
    }
    
    @Override
    public Iterator<ROI3D> iterator()
    {
        return roiSequence.values().iterator();
    }
    
    @Override
    public boolean loadFromXML(Node node)
    {
        beginUpdate();
        try
        {
            if (!super.loadFromXML(node)) return false;
            
            // we don't need to save the 2D ROI class as the parent class already do it
            clear();
            
            for (Element e : XMLUtil.getElements(node, "slice"))
            {
                // faster than using complete XML serialization
                final ROI3D slice = createSlice();
                
                // error while reloading the ROI from XML
                if ((slice == null) || !slice.loadFromXML(e)) return false;
                
                setTimePoint(slice.getT(), slice);
            }
        }
        finally
        {
            endUpdate();
        }
        
        return true;
    }
    
    @Override
    public boolean saveToXML(Node node)
    {
        if (!super.saveToXML(node)) return false;
        
        for (ROI3D slice : roiSequence.values())
        {
            Element sliceNode = XMLUtil.addElement(node, "slice");
            
            if (!slice.saveToXML(sliceNode)) return false;
        }
        
        return true;
    }
    
    public class ROI3DSequencePainter extends ROIPainter
    {
        ROI3D getSliceForCanvas(IcyCanvas canvas)
        {
            final int t = canvas.getPositionT();
            
            return t >= 0 ? getSlice(t) : null;
        }
        
        @Override
        public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
        {
            if (isActiveFor(canvas))
            {
                if (canvas instanceof IcyCanvas3D)
                {
                    // TODO
                    
                }
                else if (canvas instanceof IcyCanvas2D)
                {
                    // forward event to current slice
                    final ROI3D slice = getSliceForCanvas(canvas);
                    
                    if (slice != null && slice.getT() == canvas.getPositionT()) slice.getOverlay().paint(g, sequence, canvas);
                }
            }
        }
        
        @Override
        public void keyPressed(KeyEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.keyPressed(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().keyPressed(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void keyReleased(KeyEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.keyReleased(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().keyReleased(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseEntered(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseEntered(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseEntered(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseExited(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseExited(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseExited(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseMove(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseMove(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseMove(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseDrag(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseDrag(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseDrag(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mousePressed(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mousePressed(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mousePressed(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseReleased(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseReleased(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseReleased(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseClick(MouseEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseClick(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseClick(e, imagePoint, canvas);
            }
        }
        
        @Override
        public void mouseWheelMoved(MouseWheelEvent e, Point5D.Double imagePoint, IcyCanvas canvas)
        {
            // send event to parent first
            super.mouseWheelMoved(e, imagePoint, canvas);
            
            // then send it to active slice
            if (isActiveFor(canvas))
            {
                // forward event to current slice
                final ROI3D slice = getSliceForCanvas(canvas);
                
                if (slice != null) slice.getOverlay().mouseWheelMoved(e, imagePoint, canvas);
            }
        }
    }
    
}
