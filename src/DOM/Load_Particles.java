package DOM;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class Load_Particles implements PlugIn{

	ImagePlus imp;
	ImageProcessor ip;
	SMLAnalysis sml = new SMLAnalysis();
	
	public void run(String arg) {
		
		
		GenericDialog dgLoadParticles = new GenericDialog("Load Particles");
		
		String [] LoadFileFormatOptions = new String [] {
		"DoM format (all columns) .tif","QuickPALM compartible .tif (truncated)"};
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
			IJ.error("The format of loading file and loading options are different.\nIt could be due to plug-in version mismatch.\nPlug-in will load truncated QuickPALM format version.");
			FFindex =1;
		}
			
				
		for (int n=0;n<nParticles;n++)
		{
			IJ.showProgress(n, nParticles);
			if((float) ip.getf(6,n)<dFPThreshold)
			{
				sml.ptable.incrementCounter();
				sml.ptable.addValue("Amplitude_fit",        (float) ip.getf(0,n));
				sml.ptable.addValue("X_(px)",               (float) ip.getf(1,n));
				sml.ptable.addValue("Y_(px)",               (float) ip.getf(2,n));
				sml.ptable.addValue("X_(nm)",               (float) ip.getf(3,n)*1000);
				sml.ptable.addValue("Y_(nm)",               (float) ip.getf(4,n)*1000);
				sml.ptable.addValue("Z_(nm)",               (float) ip.getf(5,n)*1000);
				sml.ptable.addValue("False_positive",       (float) ip.getf(6,n));
				sml.ptable.addValue("X_loc_error_(px)",     (float) ip.getf(7,n));
				sml.ptable.addValue("Y_loc_error_(px)",     (float) ip.getf(8,n));
				sml.ptable.addValue("BGfit",                (float) ip.getf(9,n));			
				sml.ptable.addValue("IntegratedInt",        (float) ip.getf(10,n));
				sml.ptable.addValue("SNR",       			(float) ip.getf(11,n));
				sml.ptable.addValue("chi2_fit",             (float) ip.getf(12,n));
				sml.ptable.addValue("Frame_Number",         (float) ip.getf(13,n)*1000000);
				if(FFindex ==0)
				{
					sml.ptable.addValue("Iterations_fit",            (float) ip.getf(14,n));
					sml.ptable.addValue("SD_X_fit, pix",             (float) ip.getf(15,n));
					sml.ptable.addValue("SD_Y_fit, pix",             (float) ip.getf(16,n));
					sml.ptable.addValue("Amp_loc_error",             (float) ip.getf(17,n));
					
				}
			}
		}
		imp.close();
		sml.showTable();
	}
}
