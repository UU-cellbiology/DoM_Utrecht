package DOM;


import java.text.DecimalFormat;

import ij.Prefs;
import ij.gui.GenericDialog;

public class SMLDialog {

	
	int nWidth;  //image width of original image
	int nHeight; //image height of original image
	int nSlices; //number of slices
	
	//finding particles
	double dPSFsigma; //standard deviation of PSF approximated by Gaussian
	int nKernelSize; //size of Gaussian kernel for particles enhancement
	double dPixelSize; //size of pixel in nm of original image
	int nThreads; //number of threads for calculation
	int nAreaCut; //threshold of particles area
	boolean bShowParticles; //whether or not add overlay to detected particles
	boolean bIgnoreFP; //whether or not to ignore false positives
	
	
	//reconstructing image
	double dRecPixelSize;  //pixel size in nm of reconstructed image
	int nRecWidth;         //width of original image for reconstruction	
	int nRecHeight;        //height of original image for reconstruction
	int nIntIndex;         //parameter of intensity reconstruction
	int nSDIndex;          //parameter of intensity reconstruction
	double dFixedSD; 	   //value of SD in case of fixed value
	double dMagnification; //magnification coefficient
	boolean bCutoff; 	   //whether or not apply cut-off in localization
	double  dcutoff; 	   //cut-off value in original pixels
	
	
	//drift correction parameters
	boolean bDrift;   //whether or not make drift correction
	int nDriftFrames; //number of frames for averaging per one time interval during drift correction
	int nDriftPixels; //maximal shift in pixels per one time period defined by previous parameter 
	boolean bShowIntermediate;     //show intermediate reconstructions
	boolean bShowCrossCorrelation; //show cross-correlation 
	
	//dialog showing options for particle search algorithm		
	public boolean findParticles() {
		GenericDialog fpDial = new GenericDialog("Find Particles");
		fpDial.addNumericField("PSF standard devation, pix", Prefs.get("SiMoLoc.dPSFsigma", 2), 3);
		fpDial.addNumericField("Gaussial kernel size, \nodd number from 7(fast) till 13 (slow)  ", Prefs.get("SiMoLoc.nKernelSize", 7), 0);
		fpDial.addNumericField("Pixel size, nm", Prefs.get("SiMoLoc.dPixelSize", 66), 2);
		fpDial.addNumericField("Number of parallel threads", Prefs.get("SiMoLoc.nThreads", 50), 0);
		fpDial.addCheckbox("Mark detected particles? (better not use this feature on big datasets)", Prefs.get("SiMoLoc.bShowParticles", false));
		fpDial.addCheckbox("Ignore false positives?", Prefs.get("SiMoLoc.bIgnoreFP", false));
		
		fpDial.showDialog();
		if (fpDial.wasCanceled())
            return false;
		
		dPSFsigma = fpDial.getNextNumber();
		Prefs.set("SiMoLoc.dPSFsigma", dPSFsigma);
		nKernelSize = (int) fpDial.getNextNumber();
		Prefs.set("SiMoLoc.nKernelSize", nKernelSize);
		dPixelSize = fpDial.getNextNumber();
		Prefs.set("SiMoLoc.dPixelSize", dPixelSize);
		nThreads = (int) fpDial.getNextNumber();
		Prefs.set("SiMoLoc.nThreads", nThreads);
		bShowParticles = fpDial.getNextBoolean();
		Prefs.set("SiMoLoc.bShowParticles", bShowParticles);
		bIgnoreFP = fpDial.getNextBoolean();
		Prefs.set("SiMoLoc.bIgnoreFP", bIgnoreFP);		
		return true;
	}
	
	
	public boolean ReconstructImage(double xlocavg_, double ylocavg_) //dialog showing options for reconstruction image		
	{
		GenericDialog dgReconstruct = new GenericDialog("Reconstruct Dataset");
		String [] RecIntOptions = new String [] {
				"Integrated spot intensity","Amplitude of Gaussian fitting"};
		String [] RecSDOptions = new String [] {
				"Localization precision","Constant value"};
		
		dgReconstruct.addNumericField("Pixel size of reconstructed image, nm", Prefs.get("SiMoLoc.Rec_PixSize", 30), 2);
		dgReconstruct.addMessage("Average localization precision in X: " + new DecimalFormat("#.##").format(xlocavg_) + " nm, in Y: " +  new DecimalFormat("#.##").format(ylocavg_) +" nm.");
		dgReconstruct.addNumericField("Original image width, px", Prefs.get("SiMoLoc.Rec_ImWidth", 512), 0);
		dgReconstruct.addNumericField("Original image height, px", Prefs.get("SiMoLoc.Rec_ImHeight", 512), 0);
		dgReconstruct.addChoice("As spot's intensity use:", RecIntOptions, Prefs.get("SiMoLoc.Rec_Int", "Integrated spot intensity"));
		dgReconstruct.addChoice("As spot's SD use:", RecSDOptions, Prefs.get("SiMoLoc.Rec_SD", "Localization precision"));
		dgReconstruct.addNumericField("Value of SD in case of constant (in original pixels):", Prefs.get("SiMoLoc.Rec_SDFixed", 2), 2);
		dgReconstruct.addCheckbox("Apply cut-off for localization precision", Prefs.get("SiMoLoc.applycutoff", false));
		dgReconstruct.addNumericField("Cut-off particles with localization less than (in original pixels): ", Prefs.get("SiMoLoc.cutoff", 0.3), 2);
		dgReconstruct.addMessage("\n\n");
		dgReconstruct.addCheckbox("Apply drift correction", Prefs.get("SiMoLoc.drift", false));		
		dgReconstruct.addNumericField("Number of frames for averaging:", Prefs.get("SiMoLoc.drift_frames", 1000), 0);
		dgReconstruct.addNumericField("Max shift in pixels:", Prefs.get("SiMoLoc.drift_pixels", 10), 0);
		dgReconstruct.addCheckbox("Show intermediate reconstructions:", Prefs.get("SiMoLoc.drift_intermediate_reconstr", false));		
		dgReconstruct.addCheckbox("Show cross-correlation images:", Prefs.get("SiMoLoc.drift_cross_correlation", false));
		dgReconstruct.addMessage("\n\n");		
		dgReconstruct.showDialog();
		if (dgReconstruct.wasCanceled())
            return false;
		
        dRecPixelSize = dgReconstruct.getNextNumber();
        Prefs.set("SiMoLoc.Rec_PixSize", dRecPixelSize);
		nRecWidth = (int) dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.Rec_ImWidth", nRecWidth);
		nRecHeight = (int) dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.Rec_ImHeight", nRecHeight);
		nIntIndex = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.Rec_Int", RecIntOptions[nIntIndex]);
		nSDIndex = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.Rec_SD", RecSDOptions[nSDIndex]);
		dFixedSD = dgReconstruct.getNextNumber();
	    Prefs.set("SiMoLoc.Rec_SDFixed", dFixedSD);
		bCutoff = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.applycutoff", bCutoff);
		dcutoff = dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.cutoff", dcutoff);
		bDrift = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.drift", bDrift);
		nDriftFrames = (int) dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.drift_frames", nDriftFrames);
		nDriftPixels = (int) dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.drift_pixels", nDriftPixels);
		bShowIntermediate = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.drift_intermediate_reconstr", bShowIntermediate);
		bShowCrossCorrelation = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.drift_cross_correlation", bShowCrossCorrelation);
		
		return true;
	}
	

}
