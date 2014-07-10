package plugins.adufour.activecontours;

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginROI;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import plugins.adufour.activecontours.ActiveContours.ExportROI;
import plugins.kernel.roi.roi2d.ROI2DEllipse;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

public class MagicWand2D extends Plugin implements PluginROI
{
    public MagicWand2D()
    {
        
    }
    
    @Override
    public String getROIClassName()
    {
        return ROI2DPolygon.class.getName();
    }
    
    @Override
    public ROI createROI(Point5D pt)
    {
        double radius = 10.0;
        ROI2DEllipse defaultROI = new ROI2DEllipse(pt.getX() - radius, pt.getY() - radius, pt.getX() + radius, pt.getY() + radius); 
        
        ActiveContours ac = new ActiveContours();
        Sequence s = getActiveSequence();
        
        
        ac.input.setValue(s);
        ac.contour_timeStep.setValue(2.0);
        ac.convergence_criterion.setValue(0.1);
        ac.evolution_bounds.setNoSequenceSelection();
        ac.roiInput.add(defaultROI);
        ac.region_weight.setValue(1.0);
        ac.edge_weight.setValue(0.0);
        ac.contour_resolution.setValue(radius / 10);
        ac.output_rois.setValue(ExportROI.ON_INPUT);
        
        try
        {
            ac.run();
            ac.clean();
            
            ROI[] roi = ac.roiOutput.getValue();
            if (roi.length > 0) return roi[0];
        }
        catch (Exception e)
        {
            new AnnounceFrame("Unable to find an object", 5000);
            ac.clean();
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
        
        return defaultROI;
    }
    
    @Override
    public ROI createROI()
    {
        return null;
    }
}
