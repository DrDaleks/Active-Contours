package plugins.adufour.activecontours;

import icy.canvas.IcyCanvas;
import icy.painter.PainterAdapter;
import icy.sequence.Sequence;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.HashMap;

import plugins.fab.trackmanager.Detection;
import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;

public class ActiveContoursPainter extends PainterAdapter
{
	private HashMap<Integer, ArrayList<ActiveContour>>	contoursMap;
	
	private TrackGroup									trackPool;
	
	public ActiveContoursPainter(HashMap<Integer, ArrayList<ActiveContour>> contours)
	{
		contoursMap = contours;
	}
	
	public ActiveContoursPainter(TrackGroup trackPool)
	{
		this.trackPool = trackPool;
	}
	
	@Override
	public void paint(Graphics2D g, Sequence sequence, IcyCanvas canvas)
	{
		super.paint(g, sequence, canvas);
		
		int t = canvas.getT();
		
		if (trackPool == null)
		{
			if (contoursMap.containsKey(t))
				for (ActiveContour contour : contoursMap.get(t))
					contour.paint(g, sequence, canvas);
		}
		else
		{
			ArrayList<TrackSegment> segments = trackPool.getTrackSegmentList();
			
			for (int i = 1; i <= segments.size(); i++)
			{
				TrackSegment segment = segments.get(i - 1);
				ArrayList<Detection> detections = segment.getDetectionList();
				for (int d = 0; d < detections.size(); d++)
				{
					Detection det = detections.get(d);
					
					if (det.getT() == canvas.getT())
					{
						((ActiveContour) det).paint(g, sequence, canvas);
						g.drawString("#" + i, (float) det.getX(), (float) det.getY());
					}
				}
			}
		}
	}
}
