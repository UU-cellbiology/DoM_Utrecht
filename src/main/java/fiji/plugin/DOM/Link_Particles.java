package fiji.plugin.DOM;

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
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for particle linking, please load one.");
			return;
		}
		
		
		//show 'Link Particles' dialog window
		if (!dlg.LinkParticles()) return;
		
		LogLinkParameters(dlg);
			
		//sort table by frame number		
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_FrameN, true);
		
		imp = null;
		ovTracks = null;
		//whether or not show tracks
		if(dlg.bShowTracks || dlg.bShowParticlesLink)
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
			frames = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
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
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_TrackID, true);
		
		smlLink.calculate_Tracks_Lengths();
		smlLink.add_Tracks_Lengths();	
		
		//marking particles in each track in reverse order
		//in the Results table
		smlLink.mark_Particles_in_Reverse();
		
		//adding tracks to overlay
		if(dlg.bShowTracks)
		{
			smlLink.addTracksToOverlay(true);
		}
		//adding particles to overlay
		if(dlg.bShowParticlesLink)
		{
			smlLink.addParticlesToOverlay();		
		}
		
		IJ.showStatus("Linking trajectories... Done.");
		
		if(dlg.bShowTracks || dlg.bShowParticlesLink)
		{
					
			imp.setOverlay(smlLink.ovTracks);
			imp.updateAndRepaintWindow();
			imp.show();
		}
		
		sml.showTable();
		
	}
	void LogLinkParameters(SMLDialog dlg)
	{
		//let's log stuff
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Linking parameters");
		if(dlg.nLinkFP==0)
		{
			IJ.log("For linking use only true positives");
		}
		else
		{
			IJ.log("For linking use all particles");
		}
		IJ.log("Max distance to search over one frame: "+String.format("%d",(int)dlg.dLinkDistance)+" pixels");
		if(dlg.nLinkTrace==0)
		{
			IJ.log("Measure distance from initial position");
		}
		else
		{
			IJ.log("Measure distance from next detected position");
		}
		IJ.log("Maximum linking gap "+String.format("%d",(int)dlg.nLinkFrameGap)+" frames");
	}
	

}
