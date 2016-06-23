package plugins.adufour.activecontours;

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginROI;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import plugins.adufour.activecontours.ActiveContours.ExportROI;
import plugins.adufour.activecontours.ActiveContours.ROIType;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DEllipse;
import plugins.kernel.roi.roi3d.ROI3DArea;

public class MagicWand extends Plugin implements PluginROI
{
    public MagicWand()
    {
    
    }
    
    @Override
    public String getROIClassName()
    {
        return ROI2DArea.class.getName();
    }
    
    @Override
    public ROI createROI(Point5D pt)
    {
        Sequence s = getActiveSequence();
        if (s == null)
        {
            new AnnounceFrame("No active sequence");
            return null;
        }
        
        int radius = 5;
        
        ROI roi = null;
        
        if (s.getSizeZ() > 1)
        {
            // 3D
            roi = new ROI3DArea();
            
            int rx = radius, rrx = rx * rx;
            int ry = radius, rry = ry * ry;
            int rz = radius, rrz = rz * rz;
            
            for (int z = -rz; z <= rz; z++)
                for (int y = -ry; y <= ry; y++)
                    for (int x = -rx; x <= rx; x++)
                    {
                        double xr2 = rrx == 0 ? 0 : x * x / rrx;
                        double yr2 = rry == 0 ? 0 : y * y / rry;
                        double zr2 = rrz == 0 ? 0 : z * z / rrz;
                        
                        if (xr2 + yr2 + zr2 <= 1.0)
                        {
                            int px = (int) Math.round(pt.getX() + x);
                            int py = (int) Math.round(pt.getY() + y);
                            int pz = (int) Math.round(pt.getZ() + z);
                            
                            ((ROI3DArea) roi).addPoint(px, py, pz);
                        }
                    }
        }
        else
        {
            roi = new ROI2DEllipse(pt.getX() - radius, pt.getY() - radius, pt.getX() + radius, pt.getY() + radius);
        }
        
        ActiveContours ac = new ActiveContours();
        
        ac.input.setValue(s);
        ac.contour_timeStep.setValue(2.0);
        ac.convergence_criterion.setValue(0.1);
        ac.evolution_bounds.setNoSequenceSelection();
        ac.roiInput.add(roi);
        ac.region_weight.setValue(1.0);
        ac.edge_weight.setValue(0.0);
        if (s.getSizeZ() > 1)
        {
            ac.contour_resolution.setValue(s.getPixelSizeZ() * 0.75);
        }
        else
        {
            ac.contour_resolution.setValue(radius / 3.0);
        }
        ac.output_rois.setValue(ExportROI.ON_INPUT);
        ac.output_roiType.setValue(ROIType.POLYGON);
        
        try
        {
            ac.run();
            ac.clean();
            
            ROI[] rois = ac.roiOutput.getValue();
            if (rois.length > 0) return rois[0];
        }
        catch (Exception e)
        {
            new AnnounceFrame("Unable to find an object", 5000);
            ac.clean();
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
        
        return roi;
    }
    
    @Override
    public ROI createROI()
    {
        return null;
    }
}
