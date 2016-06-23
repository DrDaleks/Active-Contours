package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.Overlay;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceEvent.SequenceEventSourceType;
import icy.sequence.SequenceEvent.SequenceEventType;
import icy.sequence.SequenceListener;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.ConcurrentModificationException;

import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;

public class ActiveContoursOverlay extends Overlay implements SequenceListener
{
    private TrackGroup trackGroup;
    
    public ActiveContoursOverlay(TrackGroup trackGroup)
    {
        super("Active contours");
        this.trackGroup = trackGroup;
    }
    
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        if (!Arrays.asList(sequence.getListeners()).contains(this)) sequence.addListener(this);
        
        int t = canvas.getPositionT();
        
        ArrayList<TrackSegment> segments = new ArrayList<TrackSegment>(trackGroup.getTrackSegmentList());
        
        for (int i = 1; i <= segments.size(); i++)
        {
            TrackSegment segment = segments.get(i - 1);
            
            if (segment == null) continue;
            
            for (int d = 0; d < segment.getDetectionList().size(); d++)
            {
                ActiveContour contour = (ActiveContour) segment.getDetectionList().get(d);
                
                if (contour == null) continue;
                
                contour.paint(g, sequence, canvas);
                
                if (g != null && contour.getT() == t)
                {
                    // in 2D, draw the contour number in its center (and mind the zoom factor)
                    float f = (float) canvas.canvasToImageLogDeltaX(18);
                    
                    float x = (float) contour.getX();
                    float y = (float) contour.getY();
                    
                    // adjust the text positioning
                    x -= (i < 10 ? f / 2 : f);
                    y += f / 2;
                    
                    g.drawString("" + i, x, y);
                }
            }
        }
        
    }
    
    @Override
    public void remove()
    {
        ArrayList<TrackSegment> segments = new ArrayList<TrackSegment>(trackGroup.getTrackSegmentList());
        
        for (int i = 1; i <= segments.size(); i++)
        {
            TrackSegment segment = segments.get(i - 1);
            
            if (segment == null) continue;
            
            try
            {
                ArrayList<Detection> detections = segment.getDetectionList();
                
                if (detections != null) for (Detection detection : detections)
                {
                    ((ActiveContour) detection).clean();
                }
            }
            catch (ConcurrentModificationException e)
            {
                // segment has probably changed while looking for a detection to paint
                // => ignore and wait for the next repaint
            }
        }
        
        super.remove();
        trackGroup.getSequence().removeOverlay(this);
    }
    
    @Override
    public void sequenceChanged(SequenceEvent sequenceEvent)
    {
        if (sequenceEvent.getSourceType() != SequenceEventSourceType.SEQUENCE_OVERLAY) return;
        
        if (sequenceEvent.getSource() == this && sequenceEvent.getType() == SequenceEventType.REMOVED)
        {
            sequenceEvent.getSequence().removeListener(this);
            remove();
        }
    }
    
    @Override
    public void sequenceClosed(Sequence sequence)
    {
        sequence.removeListener(this);
        remove();
    }
}
