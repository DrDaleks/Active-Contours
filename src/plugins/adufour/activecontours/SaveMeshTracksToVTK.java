package plugins.adufour.activecontours;

import java.io.File;
import java.util.ArrayList;

import icy.plugin.abstract_.Plugin;
import icy.util.StringUtil;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.Var;
import plugins.adufour.vars.lang.VarFile;
import plugins.adufour.vars.lang.VarString;
import plugins.fab.trackmanager.TrackGroup;
import plugins.fab.trackmanager.TrackSegment;
import plugins.nchenouard.spot.Detection;

public class SaveMeshTracksToVTK extends Plugin implements Block
{
    Var<TrackGroup> meshTracks = new Var<TrackGroup>("Mesh tracks", TrackGroup.class);
    VarFile folder = new VarFile("Folder", null);
    VarString prefix = new VarString("Prefix", "Mesh");
    
    @Override
    public void declareInput(VarList inputMap)
    {
        inputMap.add("ROI", meshTracks);
        inputMap.add("folder", folder);
        inputMap.add("prefix", prefix);
    }
    
    @Override
    public void declareOutput(VarList outputMap)
    {
    }
    
    @Override
    public void run()
    {
        ArrayList<TrackSegment> segments = new ArrayList<TrackSegment>(meshTracks.getValue(true).getTrackSegmentList());
        
        for (int i = 1; i <= segments.size(); i++)
        {
            TrackSegment segment = segments.get(i - 1);
            
            if (segment == null) continue;
            
            for(Detection detection : segment.getDetectionList())
            
            if (detection instanceof Mesh3D)
            {
                String meshID = "_#" + StringUtil.toString(i, 3);
                String timeID = "_T" + StringUtil.toString(detection.getT(), 3);
                String fileName = prefix.getValue() + meshID + timeID + ".vtk";
                File vtkFile = new File(folder.getValue(true).getPath() + File.separator + fileName);
                ((Mesh3D) detection).mesh.saveToVTK(vtkFile);
            }
        }
    }
    
}
