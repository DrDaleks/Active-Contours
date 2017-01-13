package plugins.adufour.roi;

import icy.canvas.IcyCanvas;
import icy.canvas.IcyCanvas2D;
import icy.painter.OverlayEvent;
import icy.painter.OverlayListener;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROI3D;
import icy.roi.ROIEvent;
import icy.roi.ROIListener;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import icy.type.rectangle.Rectangle5D;
import icy.util.XMLUtil;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.lang.reflect.Method;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.concurrent.Semaphore;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

public class TemporalROI<R extends ROI> extends ROI implements ROIListener, OverlayListener, Iterable<R>
{
    public static final String          PROPERTY_USECHILDCOLOR = "useChildColor";
    
    protected final TreeMap<Integer, R> slices                 = new TreeMap<Integer, R>();
    
    private int                         channel;
    
    protected boolean                   useChildColor;
    protected Semaphore                 modifyingSlice;
    
    public TemporalROI()
    {
        super();
        useChildColor = false;
        modifyingSlice = new Semaphore(1);
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
        final R slice = getSlice(t);
        
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
                    for (R slice : slices.values())
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
                for (R slice : slices.values())
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
                for (R slice : slices.values())
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
                for (R slice : slices.values())
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
                for (R slice : slices.values())
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
                for (R slice : slices.values())
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
                for (R slice : slices.values())
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
    
    public int getC()
    {
        return channel;
    }
    
    public void setC(int value)
    {
        beginUpdate();
        try
        {
            this.channel = value;
            
            modifyingSlice.acquireUninterruptibly();
            
            try
            {
                for (R slice : slices.values())
                    setC(slice, value);
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
        return slices.isEmpty();
    }
    
    public int getFirstT()
    {
        return slices.firstKey();
    }
    
    public int getLastT()
    {
        return slices.lastKey();
    }
    
    /**
     * @return The size of this ROI stack along Z.<br>
     *         Note that the returned value indicates the difference between upper and lower bounds
     *         of this ROI, but doesn't guarantee that all slices in-between exist (
     *         {@link #getSlice(int)} may still return <code>null</code>.<br>
     */
    public int getSizeT()
    {
        if (slices.isEmpty()) return 0;
        
        return (slices.lastKey().intValue() - slices.firstKey().intValue()) + 1;
    }
    
    /**
     * Returns the ROI slice at given T position.
     */
    public R getSlice(int t)
    {
        return slices.get(Integer.valueOf(t));
    }
    
    /**
     * Removes all slices.
     */
    public void clear()
    {
        for (R slice : slices.values())
        {
            slice.removeListener(this);
            slice.getOverlay().removeOverlayListener(this);
        }
        
        slices.clear();
    }
    
    /**
     * Sets the slice for the given Z position.
     */
    public void setSlice(int t, R slice)
    {
        setT(slice, t);
        setC(slice, getC());
        
        // listen events from this ROI and its overlay
        slice.addListener(this);
        slice.getOverlay().addOverlayListener(this);
        
        slices.put(Integer.valueOf(t), slice);
        
        // notify ROI changed
        roiChanged(true);
    }
    
    /**
     * calls {@link ROI2D#setC(int)} or {@link ROI3D#setC(int)} via reflection (because {@link ROI}
     * has not such method)
     * 
     * @param roi
     * @return
     */
    private void setC(R roi, int c)
    {
        try
        {
            Method setC = roi.getClass().getMethod("setC", Integer.TYPE);
            setC.invoke(roi, c);
        }
        catch (Exception e)
        {
            throw new UnsupportedOperationException("Cannot assign channel to " + roi, e);
        }
    }
    
    /**
     * calls {@link ROI2D#setT(int)} or {@link ROI3D#setT(int)} via reflection (because {@link ROI}
     * has not such method)
     * 
     * @param roi
     * @return
     */
    private void setT(R roi, int t)
    {
        try
        {
            Method setT = roi.getClass().getMethod("setT", Integer.TYPE);
            setT.invoke(roi, t);
        }
        catch (Exception e)
        {
            throw new UnsupportedOperationException("Cannot assign time point to " + roi, e);
        }
    }
    
    /**
     * calls {@link ROI2D#getT()} or {@link ROI3D#getT()} via reflection (because {@link ROI} has
     * not such method)
     * 
     * @param roi
     * @return
     */
    private int getT(R roi)
    {
        try
        {
            Method getT = roi.getClass().getMethod("getT");
            return (Integer) getT.invoke(roi);
        }
        catch (Exception e)
        {
            throw new UnsupportedOperationException("Cannot get time point from " + roi, e);
        }
    }
    
    /**
     * Removes slice at the given Z position and returns it.
     */
    public R removeSlice(int t)
    {
        // remove the current slice (if any)
        final R result = slices.remove(t);
        
        // remove listeners
        if (result != null)
        {
            result.removeListener(this);
            result.getOverlay().removeOverlayListener(this);
        }
        
        // notify ROI changed
        roiChanged(true);
        
        return result;
    }
    
    @Override
    public Iterator<R> iterator()
    {
        return slices.values().iterator();
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
                String roiClassName = XMLUtil.getAttributeValue(e, "className", null);
                
                // faster than using complete XML serialization
                @SuppressWarnings("unchecked")
                final R slice = (R) Class.forName(roiClassName).newInstance();
                
                // error while reloading the ROI from XML
                if (slice == null || !slice.loadFromXML(e)) return false;
                
                setSlice(getT(slice), slice);
            }
        }
        catch (Exception e1)
        {
            throw new RuntimeException("Couldn't load temporal ROI from XML", e1);
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
        
        for (R slice : slices.values())
        {
            Element sliceNode = XMLUtil.addElement(node, "slice");
            XMLUtil.setAttributeValue(sliceNode, "className", slice.getClassName());
            
            if (!slice.saveToXML(sliceNode)) return false;
        }
        
        return true;
    }
    
    @Override
    public double computeNumberOfPoints()
    {
        double interior = 0;
        
        for (R slice : slices.values())
            interior += slice.getNumberOfPoints();
        
        return interior;
    }
    
    @Override
    public double computeNumberOfContourPoints()
    {
        double contour = 0;
        
        for (R slice : slices.values())
            contour += slice.getNumberOfContourPoints();
        
        return contour;
    }
    
    @Override
    public boolean hasSelectedPoint()
    {
        // default
        return false;
    }
    
    // LISTENERS //
    
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
            roiChanged(true);
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
    protected ROIPainter createPainter()
    {
        return new TemporalROIPainter();
    }
    
    protected class TemporalROIPainter extends ROIPainter
    {
        @Override
        public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
        {
            if (canvas instanceof IcyCanvas2D && isActiveFor(canvas))
            {
                // forward event to current slice
                final R slice = getCurrentSlice(canvas);
                
                if (slice != null && getT(slice) == canvas.getPositionT()) slice.getOverlay().paint(g, sequence, canvas);
            }
        }
        
        protected R getCurrentSlice(IcyCanvas canvas)
        {
            final int t = canvas.getPositionT();
            
            return t >= 0 ? getSlice(t) : null;
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
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
                final R slice = getCurrentSlice(canvas);
                
                if (slice != null) slice.getOverlay().mouseWheelMoved(e, imagePoint, canvas);
            }
        }
    }
    
    @Override
    public int getDimension()
    {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public boolean isActiveFor(IcyCanvas canvas)
    {
        int currentT = canvas.getPositionT();
        return currentT == -1 || (slices.containsKey(currentT));
    }
    
    @Override
    public Rectangle5D computeBounds5D()
    {
        Rectangle5D bounds = null;
        
        for(R roi : slices.values())
        {
            if (bounds == null)
            {
                bounds = roi.getBounds5D();                
            }
            else
            {
                bounds = bounds.createUnion(roi.getBounds5D());
            }
        }
        
        return bounds;
    }
    
    @Override
    public boolean canSetBounds()
    {
        return false;
    }
    
    @Override
    public boolean canSetPosition()
    {
        return false;
    }
    
    @Override
    public void setBounds5D(Rectangle5D bounds)
    {
        System.out.println("setting bounds 5D");
    }
    
    @Override
    public void setPosition5D(Point5D position)
    {
        System.out.println("setting position 5D");
    }
    
    @Override
    public boolean contains(double x, double y, double z, double t, double c)
    {
        return slices.containsKey(t) && slices.get(t).contains(x, y, z, t, c);
    }
    
    @Override
    public boolean contains(double x, double y, double z, double t, double c, double sizeX, double sizeY, double sizeZ, double sizeT, double sizeC)
    {
        if (t < getFirstT() || t + sizeT > getLastT()) return false;
        
        for (R slice : slices.values())
            if (!slice.contains(x, y, z, t, c, sizeX, sizeY, sizeZ, sizeT, sizeC)) return false;
        
        return true;
    }
    
    @Override
    public boolean intersects(double x, double y, double z, double t, double c, double sizeX, double sizeY, double sizeZ, double sizeT, double sizeC)
    {
        if (t > getLastT() || t + sizeT < getFirstT()) return false;
        
        for (R slice : slices.values())
            if (!slice.intersects(x, y, z, t, c, sizeX, sizeY, sizeZ, sizeT, sizeC)) return false;
        
        return true;
    }
}
