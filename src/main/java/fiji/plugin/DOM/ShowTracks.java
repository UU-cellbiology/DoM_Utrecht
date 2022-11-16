package fiji.plugin.DOM;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.PlugIn;

public class ShowTracks implements PlugIn{
	
	
	ImagePlus imp;
	SMLAnalysis sml = new SMLAnalysis();
	int indX, indY, indFrame, indTrack, indTrLength;
	@Override
	
	public void run(String arg) {
		imp = IJ.getImage();		
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		
		if (sml.ptable.getCounter()==0)
		{
			IJ.error("Error. No Results table.");
			return;
		}
		indX = getResultsColIndex(sml,"X_(px)");
		if(indX<0)
		{
			IJ.error("Error. Results table does not contain tracks information in DoM format.");
			return;
		}
		indY = getResultsColIndex(sml,"Y_(px)");
		if(indY<0)
		{
			IJ.error("Error. Results table does not contain tracks information in DoM format.");
			return;
		}
		indFrame = getResultsColIndex(sml,"Frame");
		if(indFrame<0)
			indFrame = getResultsColIndex(sml,"Frame_Number");
		if(indFrame<0)
		{
			IJ.error("Error. Results table does not contain tracks information in DoM format.");
			return;
		}
		indTrack = getResultsColIndex(sml,"Track_ID");
		if(indTrack<0)
		{
			IJ.error("Error. Results table does not contain tracks information in DoM format.");
			return;
		}
		indTrLength = getResultsColIndex(sml,"Track_Length");
		if(indTrLength<0)
		{
			IJ.error("Error. Results table does not contain tracks information in DoM format.");
			return;
		}
		//sort table by frame number		
		Sort_Results.sorting_external_silent(sml, indFrame, true);
		//and by Track ID now
		Sort_Results.sorting_external_silent(sml, indTrack, true);
		Overlay ovTracks = imp.getOverlay();
		if (ovTracks==null)
			ovTracks = new Overlay();
		else
		{
			ovTracks.clear();
		}
		SMLLinker smlLink = new SMLLinker(sml,ovTracks,  indX, indY, indFrame, indTrack, indTrLength);
		smlLink.addTracksToOverlay(false);
		imp.setOverlay(smlLink.ovTracks);
		imp.updateAndRepaintWindow();
		imp.show();
	}
	
	public static int getResultsColIndex(final SMLAnalysis sml, final String colTitle)
	{
		
		for(int i=0;i<sml.ptable.getHeadings().length;i++)
			if(sml.ptable.getHeadings()[i].equals(colTitle))
				return i;
		return -1;
	}

}
