package fiji.plugin.DOM;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
/** 
 * old DOM format (less than v.1) reverse compartibility
 * 
 * */
public class Load_Particles implements PlugIn{

	ImagePlus imp;
	ImageProcessor ip;
	SMLAnalysis sml = new SMLAnalysis();
	
	public void run(String arg) {
		
		IJ.register(Load_Particles.class);
		GenericDialog dgLoadParticles = new GenericDialog("Load Particles");
		
		String [] LoadFileFormatOptions = new String [] {
		"DoM old format (all columns) .tif","QuickPALM compartible .tif (truncated)"};
		int FFindex;		
		String [] LoadFPOptions = new String [] {
		"Load all particles","Load only true positives","Load true and half positives"};
		int FPindex;
		double dFPThreshold =2.0;

		dgLoadParticles.addChoice("File format:", LoadFileFormatOptions, Prefs.get("SiMoLoc.FFStore", "DoM format (all columns) .tif"));
		dgLoadParticles.addChoice("What particles to load?", LoadFPOptions, Prefs.get("SiMoLoc.FPStore", "Load all particles"));
				
		dgLoadParticles.showDialog();
		if (dgLoadParticles.wasCanceled())
            return;
		FFindex = dgLoadParticles.getNextChoiceIndex();
		FPindex = dgLoadParticles.getNextChoiceIndex();
		
		Prefs.set("SiMoLoc.FFStore", LoadFileFormatOptions[FFindex]);
		Prefs.set("SiMoLoc.FPStore", LoadFPOptions[FPindex]);
		imp = IJ.openImage();
		if(imp == null)
			return;
		ip = imp.getProcessor();
		sml.ptable.reset();
		switch(FPindex){
			case 0: dFPThreshold = 2.0; break;
			case 1: dFPThreshold = 0.4; break;
			case 2: dFPThreshold = 0.6; break;
		}
		
				
		IJ.showStatus("Loading Particles Table...");
		int nParticles = ip.getHeight();
		int nCheck = ip.getWidth();
		
		if(nCheck<18 && FFindex ==0)
		{
			IJ.error("The format of loading file is wrong.");
			FFindex =1;
			return;
		}
			
		float pxscale = (float) ip.getf(3,1)*1000/(float) ip.getf(1,1);
		boolean bQuality;
		for (int n=0;n<nParticles;n++)
		{
			IJ.showProgress(n, nParticles);
			bQuality=true;
			if((float) ip.getf(6,n)>dFPThreshold)
				bQuality=false;
			//localizations make sense
			if((float) ip.getf(7,n)*pxscale>1000)
				bQuality=false;
			if((float) ip.getf(8,n)*pxscale>1000)
				bQuality=false;
			
			if(bQuality)
			{
				sml.ptable.incrementCounter();
				
				sml.ptable.addValue("X_(px)",               (float) ip.getf(1,n));
				sml.ptable.addValue("Y_(px)",               (float) ip.getf(2,n));
				sml.ptable.addValue("Frame_Number",         (float) ip.getf(13,n)*1000000);
				sml.ptable.addValue("X_(nm)",               (float) ip.getf(3,n)*1000);
				sml.ptable.addValue("X_loc_error_(nm)",     (float) ip.getf(7,n)*pxscale);
				sml.ptable.addValue("Y_(nm)",               (float) ip.getf(4,n)*1000);
				sml.ptable.addValue("Y_loc_error_(nm)",     (float) ip.getf(8,n)*pxscale);
				sml.ptable.addValue("Z_(nm)",               (float) ip.getf(5,n)*1000);
				sml.ptable.addValue("Z_loc_error_(nm)",     0.0);

				sml.ptable.addValue("Amplitude_fit",        (float) ip.getf(0,n));
				sml.ptable.addValue("Amp_error",            (float) ip.getf(17,n));
				
				sml.ptable.addValue("BGfit",                (float) ip.getf(9,n));
				sml.ptable.addValue("BGfit_error",          0.0);
				sml.ptable.addValue("SD_X_(nm)",            (float) ip.getf(15,n)*pxscale);
				sml.ptable.addValue("SD_X_error(nm)",       0.0);
				sml.ptable.addValue("SD_Y_(nm)",            (float) ip.getf(16,n)*pxscale);
				sml.ptable.addValue("SD_Y_error(nm)",       0.0);
				sml.ptable.addValue("False_positive",       (float) ip.getf(6,n));
				sml.ptable.addValue("IntegratedInt",        (float) ip.getf(10,n));
				sml.ptable.addValue("SNR",       			(float) ip.getf(11,n));
				sml.ptable.addValue("chi2_fit",             (float) ip.getf(12,n));
				sml.ptable.addValue("Iterations_fit",       (float) ip.getf(14,n));
			
			}
		}
		imp.close();
		sml.showTable();
	}
}
