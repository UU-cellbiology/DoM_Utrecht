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
		double [] falsepos;	
		double [] yloc;
		double [] frames;
		double xlocavg, ylocavg, pxsize;
		double fminframe, fmaxframe;
		int i, sz;
		double dPatCount;
		IJ.register(Reconstruct_Image.class);

		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.columnExists(13))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}

		imagename = "Reconstructed Image";
		
		//calculate average localization precision 	
		falsepos = sml.ptable.getColumnAsDoubles(6);
		xloc = sml.ptable.getColumnAsDoubles(7);
		yloc = sml.ptable.getColumnAsDoubles(8);		
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
		pxsize =  sml.ptable.getValueAsDouble(3, 0)/sml.ptable.getValueAsDouble(1, 0);
		xlocavg = pxsize*xlocavg/dPatCount;
		ylocavg = pxsize*ylocavg/dPatCount;
		
		//calculate min and max frame number		
		frames = sml.ptable.getColumnAsDoubles(13);
		Arrays.sort(frames);
		fminframe = frames[0];
		fmaxframe = frames[frames.length-1];
		
		//show dialog with options
		if (!dlg.ReconstructImage(xlocavg,ylocavg,fminframe,fmaxframe)) return;
		
		//check some parameters for consistensy 
		if(dlg.bFramesInterval)
		{
			if(dlg.nFrameMax>fmaxframe ||  dlg.nFrameMax<fminframe)
			{
				IJ.error("Maximum frame number is out of range!");
				return;
			}
			if(dlg.nFrameMin>fmaxframe ||  dlg.nFrameMin<fminframe)
			{
				IJ.error("Minimum frame number is out of range!");
				return;
			}
			if(dlg.nFrameMin>dlg.nFrameMax)
			{
				IJ.error("Minimum frame should be less then maximum frame number!");
				return;
			}
			imagename += " (frames " +dlg.nFrameMin+" till "+dlg.nFrameMax+")";

		}
		
		//create reconstruction object
		smlViewer = new SMLReconstruct(imagename, sml, dlg);
		
		
		//smlViewer.clear();
		if(dlg.bDrift)
		{			
			smlViewer.sortbyframe();
			smlViewer.DriftCorrection();
			imagename += " (DriftCorrection frames="+dlg.nDriftFrames+" pixels="+dlg.nDriftPixels+")";
			//
			//smlViewer.correctDriftCOM();
		}	

		if(dlg.bAveragePositions)
		{
			smlViewer.averagelocalizations();
			
		}
		if(!dlg.b3D)//if you don't want a z-stack
		{
			if(dlg.bFramesInterval)
			{	smlViewer.draw_unsorted((int)dlg.nFrameMin, (int)dlg.nFrameMax);}
			else
			{	smlViewer.draw_unsorted(1, smlViewer.nframes);}
			
			if(dlg.bCutoff)
			{
				imagename += " (Cutoff localization=" + dlg.dcutoff +"px)";
			}
		}
		else
		{//create z-stack
			int zStep = (int)dlg.dDistBetweenZSlices;
			
			if(dlg.bFramesInterval)
			{	smlViewer.draw_zstack((int)dlg.nFrameMin, (int)dlg.nFrameMax,zStep);}
			else
			{	smlViewer.draw_zstack(1, smlViewer.nframes,zStep);}
		}
		
		smlViewer.imp.setTitle(imagename);
		smlViewer.imp.show();
		
		
	}
}
