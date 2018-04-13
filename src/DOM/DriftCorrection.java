package DOM;

import java.util.Arrays;

import ij.IJ;
import ij.plugin.PlugIn;

public class DriftCorrection implements PlugIn{


	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLReconstruct smlViewer;
	
	@Override
	public void run(String arg) 
	{

		String imagename; 
		double [] xloc;	
		double [] yloc;
		double [] xnm;
		double [] ynm;
		double [] falsepos;
		double [] frames;
		double xmax, ymax;
		double xlocavg, ylocavg;
		double fminframe, fmaxframe;
		int i, sz;
		double dPatCount;
		
		IJ.register(DriftCorrection.class);
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for drift correction, please load one.");
			return;
		}
		
		//calculate average localization precision 	
		falsepos = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Fp);
		xloc = sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
		yloc = sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);		
		sz = xloc.length; 
		xlocavg=0; ylocavg = 0;
		dPatCount=0;
		double pxsize =  sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0);
		for (i=0; i<sz; i++)
		{
			if(xloc[i]*pxsize<1 && yloc[i]*pxsize<1)
			{
				xlocavg+=xloc[i];
				ylocavg+=yloc[i];
				dPatCount++;
			}
		}
		
		xlocavg = xlocavg/dPatCount;
		ylocavg = ylocavg/dPatCount;
		
		//calculate min and max frame number		
		frames = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		Arrays.sort(frames);
		fminframe = frames[0];
		fmaxframe = frames[frames.length-1];
		
		//calculate max x and y coordinates
		xnm = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
		ynm = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);		
		xmax = 0;
		ymax = 0;
		for (i=0; i<sz; i++)
		{
			//check if localization precision is less than 1 pixels
			if(xloc[i]*pxsize<1 && yloc[i]*pxsize<1)
			{
				if(xnm[i]>xmax)
					xmax=xnm[i];
				if(ynm[i]>ymax)
					ymax=ynm[i];
			}
		}
		
		//show dialog with options
		if (!dlg.DriftCorrection(xlocavg,ylocavg,fminframe,fmaxframe, xmax, ymax)) return;
		
		
		imagename = "Reconstructed Image (Drift Corrected. Window size: " + String.format("%d",dlg.nDriftFrames)+" frames. Max shift: "+ String.format("%.2f",dlg.nDriftMaxDistnm)+ " nm)";
		//create reconstruction object
		smlViewer = new SMLReconstruct(imagename, sml, dlg);
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");		
		IJ.log("Calculating drift correction ");
		IJ.log("Batch size: "+ String.format("%d",dlg.nDriftFrames)+" frames");
		IJ.log("Pixel size of intermediate reconstructions: "+ String.format("%.2f",dlg.nDriftScale)+ " nm");
		if(dlg.nDriftMethod==0)
			IJ.log("Drift correction method: Direct CC first batch and all others");
		if(dlg.nDriftMethod==1)
			IJ.log("Drift correction method: Direct CC consecutive batches");
		
		//IJ.log("Maximum shift between windows: "+ String.format("%.2f",dlg.nDriftMaxDistnm)+ " nm");
		if(dlg.bDriftUpdateTable)
			IJ.log("Results Table update: On");
		else
			IJ.log("Results Table update: Off");
		
		
		//make drift correction
		// x and y coordinates arrays are updated inside
		smlViewer.DriftCorrection();
		// make table update if needed
		if(dlg.bDriftUpdateTable)
			smlViewer.DriftUpdateResultTable(sml, dlg);
		
		if(dlg.bDriftMakeReconstruction)
		{
			if(dlg.bFramesInterval)
			{	smlViewer.draw_unsorted((int)dlg.nFrameMin, (int)dlg.nFrameMax);}
			else
			{	smlViewer.draw_unsorted(1, smlViewer.nframes);}
			
			if(dlg.bCutoff)
			{
				imagename += " (Cutoff localization=" + dlg.dcutoff +" nm)";
			}
			smlViewer.imp.setTitle(imagename);
			smlViewer.imp.show();
			
		}
		
	}
}
