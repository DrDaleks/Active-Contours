package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.Overlay;
import icy.sequence.Sequence;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.ConcurrentModificationException;

import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;

public class ActiveContoursOverlay extends Overlay
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
            
            try
            {
                for (Detection detection : segment.getDetectionList())
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
    }
}
