package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.Overlay;
import icy.sequence.Sequence;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.HashMap;

import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;

public class ActiveContoursOverlay extends Overlay
{
    private HashMap<Integer, ArrayList<ActiveContour>> contoursMap;
    
    private TrackGroup                                 trackGroup;
    
    public ActiveContoursOverlay(HashMap<Integer, ArrayList<ActiveContour>> contours)
    {
        super("Active contours");
        contoursMap = contours;
    }
    
    public ActiveContoursOverlay(TrackGroup trackGroup)
    {
        super("Active contours");
        this.trackGroup = trackGroup;
    }
    
    @Override
    public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
    {
        if (g == null) return;
        
        int t = canvas.getPositionT();
        
        if (trackGroup == null)
        {
            if (contoursMap.containsKey(t)) for (ActiveContour contour : contoursMap.get(t))
                contour.paint(g, sequence, canvas);
        }
        else
        {
            ArrayList<TrackSegment> segments = new ArrayList<TrackSegment>(trackGroup.getTrackSegmentList());
            
            for (int i = 1; i <= segments.size(); i++)
            {
                TrackSegment segment = segments.get(i - 1);
                
                ActiveContour contour = (ActiveContour) segment.getDetectionAtTime(t);
                if (contour == null) continue;
                contour.paint(g, sequence, canvas);
                
                // draw the contour number in its center (and mind the zoom factor)
                float f = (float) canvas.canvasToImageLogDeltaX(18);
                g.drawString("" + i, (float) contour.getX() - (i < 10 ? f / 2 : f), (float) contour.getY() + f / 2);
            }
        }
    }
}
