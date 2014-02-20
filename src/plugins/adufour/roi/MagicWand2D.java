package plugins.adufour.roi;

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginROI;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import plugins.adufour.activecontours.ActiveContours;
import plugins.kernel.roi.roi2d.ROI2DEllipse;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

public class MagicWand2D extends Plugin implements PluginROI
{
    @Override
    public String getROIClassName()
    {
        return ROI2DPolygon.class.getName();
    }
    
    @Override
    public ROI createROI(Point5D pt)
    {
        ActiveContours ac = new ActiveContours();
        Sequence s = getActiveSequence();
        
        double radius = 10.0;
        
        ac.input.setValue(s);
        ac.fastInit = true;
        ac.contour_timeStep.setValue(2.0);
        ac.convergence_criterion.setValue(0.1);
        ac.evolution_bounds.setNoSequenceSelection();
        ac.roiInput.add(new ROI2DEllipse(pt.getX() - radius, pt.getY() - radius, pt.getX() + radius, pt.getY() + radius));
        ac.region_weight.setValue(1.0);
        ac.edge_weight.setValue(0.0);
        ac.contour_resolution.setValue(radius / 2);
        ac.output_rois.setValue(true);
        
        try
        {
            ac.run();
            ac.clean();
            
            ROI[] roi = ac.roiOutput.getValue();
            if (roi.length > 0) return roi[0];
            
            throw new NullPointerException();
        }
        catch (Exception e)
        {
            new AnnounceFrame("Unable to find an object", 5000);
            ac.clean();
            System.err.println(e.getMessage());
            e.printStackTrace();
            return null;
        }
    }
    
    @Override
    public ROI createROI()
    {
        // TODO Auto-generated method stub
        return null;
    }
}
