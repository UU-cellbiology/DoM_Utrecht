package DOM;



import ij.IJ;
import ij.plugin.PlugIn;


public class Reconstruct_Image implements PlugIn{

	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLReconstruct smlViewer;
	
	
	public void run(String arg) 
	{
		double [] xloc;
		double [] yloc;
		double xlocavg, ylocavg, pxsize;
		int i, sz;
		IJ.register(Reconstruct_Image.class);

		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.columnExists(13))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
	
		//calculate average localization precision 
		
		xloc = sml.ptable.getColumnAsDoubles(7);
		yloc = sml.ptable.getColumnAsDoubles(8);
		sz = xloc.length; 
		xlocavg=0; ylocavg = 0;
		for (i=0; i<sz; i++)
		{
			xlocavg+=xloc[i];
			ylocavg+=yloc[i];
		}
		pxsize =  sml.ptable.getValueAsDouble(3, 0)/sml.ptable.getValueAsDouble(1, 0);
		xlocavg = pxsize*xlocavg/(double)sz;
		ylocavg = pxsize*ylocavg/(double)sz;
		//show dialog with options
		if (!dlg.ReconstructImage(xlocavg,ylocavg)) return;
		
		smlViewer = new SMLReconstruct("Reconstruction", sml, dlg);
		
		
		//smlViewer.clear();
		if(dlg.bDrift)
		{			
			smlViewer.sortbyframe();
			smlViewer.DriftCorrection();
			smlViewer.imp.setTitle("Reconstructed Image (DriftCorrection frames="+dlg.nDriftFrames+" pixels="+dlg.nDriftPixels);
			//smlViewer.correctDriftCOM();
		}	
		//for(int i = 1; i<48000; i+=5000)
		//{
			smlViewer.draw_unsorted(1, smlViewer.nframes);
			smlViewer.imp.show();
			//smlViewer.draw(i, i+5000);
		//}
		//new ImagePlus("Results", smlViewer.imp).show();
		
		
	}
}
