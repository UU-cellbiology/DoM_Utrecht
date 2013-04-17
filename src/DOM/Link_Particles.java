package DOM;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.PlugIn;

public class Link_Particles implements PlugIn {
	
	ImagePlus imp;
	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLLinker smlLink;
	Overlay ovTracks;
	

	public void run(String arg) {
		
		double [] frames;
		IJ.register(Link_Particles.class);
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.columnExists(13))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for particle linking, please load one.");
			return;
		}
		
		
		//show 'Link Particles' dialog window
		if (!dlg.LinkParticles()) return;
		
		//sort table by frame number		
		Sort_Results.sorting_external_silent(sml, 13, true);
		
		imp = null;
		ovTracks = null;
		//whether or not show tracks
		if(dlg.bShowTracks)
		{
			imp = IJ.getImage();		

			if (imp==null )
			{
			    IJ.noImage();
			    return;
			}
			else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 ) 
			{
			    IJ.error("8 or 16 bit greyscale image required");
			    return;
			}
			frames = sml.ptable.getColumnAsDoubles(13);
			if(frames[frames.length-1]>imp.getStackSize())
			{
			    IJ.error("Results table is probably different from images stack.\n Maximum frame number in table is bigger than the stack size.");
			    return;
				
			}
			ovTracks = imp.getOverlay();
			if (ovTracks==null)
				ovTracks = new Overlay();
			else
			{
				ovTracks.clear();
			}
		}
		
		
		
		smlLink = new SMLLinker(sml, dlg, ovTracks, imp);
		
		IJ.showStatus("Linking trajectories...");
		
		smlLink.Link_Closest();
		
		//getting length (number of detections per track)		
		Sort_Results.sorting_external_silent(sml, 18, true);
		
		smlLink.calculate_Tracks_Lengths();
		smlLink.add_Tracks_Lengths();	
		
		//marking last particles
		smlLink.mark_Last_Detections(dlg.nLastDetections);
		
		//adding tracks to overlay
		if(dlg.bShowTracks)
		{
			smlLink.addTracksToOverlay();
		}
		
		IJ.showStatus("Linking trajectories... Done.");
		
		if(dlg.bShowTracks)
		{
			smlLink.addTracksToOverlay();
			imp.setOverlay(smlLink.ovTracks);
			imp.updateAndRepaintWindow();
			imp.show();
		}
		sml.showTable();
		
	}
	

}
