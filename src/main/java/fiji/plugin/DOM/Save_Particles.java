package fiji.plugin.DOM;

import ij.IJ;

import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
/** 
 * old DOM format (less than v.1) reverse compartibility
 * 
 * */
public class Save_Particles implements PlugIn {

	SMLAnalysis sml = new SMLAnalysis();
	ImageProcessor ip;
	ImagePlus imp;
	
	public void run(String arg) {
		
		IJ.register(Save_Particles.class);
		GenericDialog dgSaveParticles = new GenericDialog("Save Particles");
		
		String [] SaveFileFormatOptions = new String [] {
		"DoM format (all columns) .tif","QuickPALM compartible .tif (truncated)"};
		int FFindex;
		
		dgSaveParticles.addChoice("File format:", SaveFileFormatOptions, Prefs.get("SiMoLoc.FFStore", "DoM format (all columns) .tif"));
		dgSaveParticles.showDialog();
		if (dgSaveParticles.wasCanceled())
            return;
		FFindex = dgSaveParticles.getNextChoiceIndex();
		Prefs.set("SiMoLoc.FPStore", SaveFileFormatOptions[FFindex]);
		
		SaveDialog sd = new SaveDialog("File to save particles into", "Particles Table", ".tif");
        String path = sd.getDirectory();
        if (path==null)
        	return;
        
        String filename = path+sd.getFileName();
        
		
		double [] s = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_AmplFit);
		double [] x = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		double [] y = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		double [] x_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
		double [] y_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
		double [] z_ = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Znm);
		double [] left = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_Fp);
		double [] right = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);//transform to px
		double [] up = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);//transform to px
		double [] down = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_BGfit);
		double [] xsym = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_IntegrInt);
		double [] ysym = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		double [] wmh = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_chi);
		double [] frame = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		double [] iter =    sml.ptable.getColumnAsDoubles(DOMConstants.Col_IterN);
		double [] sdx =     sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);//transform to px
		double [] sdy =     sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);//transform to px
		double [] amploc =  sml.ptable.getColumnAsDoubles(DOMConstants.Col_Amp_error);
		
		double pxscale =  sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0);
		
		
		
		if (FFindex ==0)
			ip = new FloatProcessor(18, s.length);
		else
			ip = new FloatProcessor(14, s.length);
		imp = new ImagePlus("Particles Table", ip);
		
		IJ.showStatus("Generating Particles Table Image...");
		for (int n=0; n<s.length; n++)
		{
			IJ.showProgress(n, s.length);
			ip.setf(0, n, (float) s[n]);
			ip.setf(1, n, (float) x[n]);
			ip.setf(2, n, (float) y[n]);
			ip.setf(3, n, (float) x_[n]/1000);
			ip.setf(4, n, (float) y_[n]/1000);
			ip.setf(5, n, (float) z_[n]/1000);
			ip.setf(6, n, (float) left[n]);
			ip.setf(7, n, (float)(right[n]*pxscale));
			ip.setf(8, n, (float)(up[n]*pxscale));
			ip.setf(9, n, (float) down[n]);
			ip.setf(10, n, (float) xsym[n]);
			ip.setf(11, n, (float) ysym[n]);
			ip.setf(12, n, (float) wmh[n]);
			ip.setf(13, n, (float) frame[n]/1000000);
			if(FFindex ==0)
			{
				ip.setf(14, n, (float) iter[n]);
				ip.setf(15, n, (float)(sdx[n]*pxscale));
				ip.setf(16, n, (float)(sdy[n]*pxscale));
				ip.setf(17, n, (float) amploc[n]);
			}
		}
		IJ.showStatus("Saving Particles Table Image...");
		IJ.save(imp, filename);
		imp.close();
    }
}
