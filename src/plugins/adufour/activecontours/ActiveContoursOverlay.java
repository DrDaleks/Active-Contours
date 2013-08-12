package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.Overlay;
import icy.sequence.Sequence;
import icy.util.GraphicsUtil;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.HashMap;

import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;

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
        int t = canvas.getPositionT();
        
        if (trackGroup == null)
        {
            if (contoursMap.containsKey(t)) for (ActiveContour contour : contoursMap.get(t))
                contour.paint(g, sequence, canvas);
        }
        else
        {
            ArrayList<TrackSegment> segments = new ArrayList<TrackSegment>(trackGroup.getTrackSegmentList());
            
            int id = 1;
            for (TrackSegment segment : segments)
            {
                ArrayList<Detection> detections = new ArrayList<Detection>(segment.getDetectionList());
                
                for (Detection det : detections)
                {
                    if (det.getT() == canvas.getPositionT())
                    {
                        ((ActiveContour) det).paint(g, sequence, canvas);
                        GraphicsUtil.drawCenteredString(g, "" + id, (int) Math.round(det.getX()), (int) Math.round(det.getY()), false);
                    }
                }
                
                id++;
            }
        }
    }
}
