package DOM;
//in detect molecules on line 329
//add width and height of original image to the results table
//sml.ptable.addValue("Original_image_size",imp.getWidth());
//sml.ptable.addValue("Original_image_size",imp.getHeight());


import java.util.Arrays;

import ij.IJ;
import ij.plugin.PlugIn;


public class Reconstruct_Image implements PlugIn{

	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLReconstruct smlViewer;
	
	
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
		IJ.register(Reconstruct_Image.class);

		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}

		
		imagename = "Reconstructed Image";
		
		//calculate average localization precision 	
		falsepos = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Fp);
		xloc = sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
		yloc = sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);		
		sz = xloc.length; 
		xlocavg=0; ylocavg = 0;
		dPatCount=0;
		for (i=0; i<sz; i++)
		{
			if(falsepos[i]<0.2)
			{
				xlocavg+=xloc[i];
				ylocavg+=yloc[i];
				dPatCount++;
			}
		}
		//pxsize =  sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0);
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
			//use only true positives to estimate image borders
			if(falsepos[i]<0.2)
			{
				if(xnm[i]>xmax)
					xmax=xnm[i];
				if(ynm[i]>ymax)
					ymax=ynm[i];
			}
		}
		
		//show dialog with options
		if (!dlg.ReconstructImage(xlocavg,ylocavg,fminframe,fmaxframe, xmax, ymax)) return;
		
		
		//print parameters values to Log window
		LogParameters(dlg);
		
		//add frame interval to image name
		if(dlg.bFramesInterval)
		{
			imagename += " (frames " +Math.round(dlg.nFrameMin)+" till "+Math.round(dlg.nFrameMax)+")";

		}
		
		//create reconstruction object
		smlViewer = new SMLReconstruct(imagename, sml, dlg);
		
				
		if(dlg.bAveragePositions)
		{
			smlViewer.averagelocalizations(sml);
			
		}
		
		
		//if(!dlg.b3D)//if you don't want a z-stack
		//{
			if(dlg.bFramesInterval)
			{	smlViewer.draw_unsorted((int)dlg.nFrameMin, (int)dlg.nFrameMax);}
			else
			{	smlViewer.draw_unsorted(1, smlViewer.nframes);}
			
			if(dlg.bCutoff)
			{
				imagename += " (Cutoff localization=" + dlg.dcutoff +" nm)";
			}
		//}
		//else
		//{//create z-stack
			//int zStep = (int)dlg.dDistBetweenZSlices;
			
//			if(dlg.bFramesInterval)
	//		{	smlViewer.draw_zstack((int)dlg.nFrameMin, (int)dlg.nFrameMax,zStep);}
		//	else
			//{	smlViewer.draw_zstack(1, smlViewer.nframes,zStep);}
//		}
		
		smlViewer.imp.setTitle(imagename);
		smlViewer.imp.show();
		
		
	}
	void LogParameters(SMLDialog dlg)
	{
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Reconstruct image parameters");
		if(dlg.nRecParticles==0)
			IJ.log("Use only true positives");
		else
			IJ.log("Use all particles");
		IJ.log("Pixel size: "+ String.format("%.1f",dlg.dRecPixelSize)+ " nm");
		
		if(dlg.nSDIndex==0)
			IJ.log("Plot spots with width of localization precision");
		else
			IJ.log("Plot spots with fixed width, SD = "+ String.format("%.1f",dlg.dFixedSD)+" nm");
		if(dlg.bCutoff)
			IJ.log("Cut-off localization precision: " + String.format("%.1f",dlg.dcutoff)+" nm");
		else
			IJ.log("No localization cut-off");
		if(dlg.bTranslation)
		{
			IJ.log("Translate image in "+ String.format("%.1f",dlg.dTranslationX)+" nm in X and "+String.format("%.1f",dlg.dTranslationY)+" nm in Y");
		}
		if(dlg.bFramesInterval)
		{
			IJ.log("Use frame interval: "+String.format("%d",(int)dlg.nFrameMin)+"-"+String.format("%d",(int)dlg.nFrameMax));
		}
		if(dlg.bAveragePositions)
		{
			IJ.log("Average localizations in consecutive frames within 1 pixel: on");
			if(dlg.bUpdateAveragePositions)
				IJ.log("Update Results table with average localizations: on");
			
		}
	}
}
