package DOM;


import java.awt.MenuItem;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;

import ij.IJ;
import ij.Menus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.util.StringSorter;
import ij.util.Tools;

public class SMLDialog {

	/** image width of original image */
	int nWidth;  
	/** image height of original image */ 
	int nHeight; 
	/** number of slices */
	int nSlices; 
	
	//finding particles
	/** task to perform: finding particles, finding and fitting or fitting */
	int nDetectionType;
	/** standard deviation of PSF approximated by Gaussian */
	double dPSFsigma;
	/** size of Gaussian kernel for particles detection */
	int nKernelSize; 
	/** size of pixel in nm of original image */
	double dPixelSize;
	/** number of threads for detection/fitting */
	int nThreads; 
	/** threshold of particle's area */
	int nAreaCut; 
	/** whether or not add overlay of detected particles */
	boolean bShowParticles;
	/** whether or not to ignore false positives */
	boolean bIgnoreFP; 
	/**  wheter or not to use GPU acceleration using the OpenCL framework */
	boolean bUseGPUAcceleration;
	/**  number of spots per batch for GPU processing */
	int nBatchSize; 
	/** size of local workgroups on GPU */
	int nGroupSize; 
	/** fixed number of iterations on GPU */
	int nIterations; 
	/** fixed number of iterations on GPU */
	boolean bUseMLE; 

	//reconstructing image
	/** what type of particles (all/TP/FP) use for reconstruction */
	int nRecParticles;       
	/** pixel size in nm of reconstructed image */
	double dRecPixelSize;    
	/** width of original image for reconstruction in nm */
	double nRecWidth;
	/** height of original image for reconstruction in nm */
	double nRecHeight;          
	/** parameter of intensity reconstruction */
	//int nIntIndex;           
	/** parameter of intensity reconstruction */
	int nSDIndex;            
	/** value of SD in case of fixed value */
	double dFixedSD; 	     
	/** magnification coefficient */
	double dMagnification;   
	/** whether or not apply cut-off in localization */
	boolean bCutoff; 	     
	/** cut-off localization value in original pixels */
	double  dcutoff; 	     
	/** apply translation during reconstruction */
	boolean bTranslation;    
	/** translation in X direction */
	double dTranslationX;    
	/** translation in Y direction */
	double dTranslationY;    
	/** use all frames (false) or some interval defined by nFrameMin and nFrameMax */
	boolean bFramesInterval; 
	/** minimal frame to use in reconstruction */
	double nFrameMin;
	/** maximal frame to use in reconstruction */
	double nFrameMax;
	/** do average localizations in consecutive frames within 1 pixel for single photoactivation events **/
	boolean bAveragePositions; //
	/** update Results table after averaging consecutive localizations**/
	boolean bUpdateAveragePositions; 
	
	//drift correction parameters
	//
	/** pixel size of image for drift reconstruction, nm */
	double nDriftScale;
	/** what kind of particles use for drift correction */
	//int nDriftParticles; 
	
	/** whether object was created from Drift Correction or Reconstruct Image menus*/
	boolean bDrift=false;
	/** number of frames for averaging per one time interval during drift correction */
	int nDriftFrames; 
	/** Maximal shift in pixels per one time period defined by previous parameter for drift correction */
	int nDriftPixels;
	/** Maximal shift in nm per one time period defined by previous parameter for drift correction */
	double nDriftMaxDistnm;
	
	/** update results table with results of drift correction*/
	boolean bDriftUpdateTable;     
	/** make reconstruction after Drift Correction*/
	boolean bDriftMakeReconstruction;     	
	/** show intermediate reconstructions during drift correction */
	boolean bShowIntermediate;     
	/** show cross-correlation map of drift correction*/
	boolean bShowCrossCorrelation; 
	
	//3D reconstruction parameters
	/** whether or not to make 3D reconstruction */
	boolean b3D;		
	/** 3D rendering type 0 - stack, 1 - color coded */
	int n3DRenderType; 
	/** z-distance between slices in the z rended stack in nm */
	double dDistBetweenZSlices;	
	/** (re)calculate the z-values based on calibration-file */
	boolean bCalculateZValues;
	
	//z-calibration curve parameters
	/** what kind of input data  use for z-calibration */ 
	int sZcUse;
	/** degree of polynomial used for astigmatism z fitting*/
	//int fitPolynomialDegree;
	/** distance between frames in calibration z stack, in nm*/
	int zCalDistBetweenPlanes;
	/** Rsquared threshold */
	double zCalRsquareThreshold;
	/** adjust LUT for particle number in Z plane */
	boolean bDynamicZscale;
	/** name of LUT for colorcoded reconstruction of Z*/
	String sZLUTName;
	boolean bZCalWobbling;
	
	//color calibration parameters
	/** maximum distance between particle image in both channel
	 * */
	double dCCDist;
	/** SNR threshold for particles used for color calibration
	 * */
	double dCCSNR;
	/**  Show distortion map */
	boolean bCCShowMap;
	/**  Show points used for calibration */
	boolean bCCShowPoints;
	/**  Put coordinates of points used for calibration to table */
	boolean bCCShowResults;
	
		
	//particle linking parameters
	/** what kind of particles use for linking */
	int nLinkFP; 
	/** maximum distance between particles to link in one frame */
	double dLinkDistance; 
	/** whether to measure distance from initial spot position or 'moving particle' mode */
	int nLinkTrace; 
	/** maximum linking gap in frames */
	int nLinkFrameGap;
	/** whether to show linked tracks in overlay */
	boolean bShowTracks;
	/** add detected particles to overlay after linking*/
	boolean bShowParticlesLink;
	
	/** Dialog showing options for z axis calibration 
	 * 		
	 * @return
	 */
		public boolean zCalibration() {
		
			String [] zcUseOptions = new String [] {"Image stack","Particle table"};//,"Polynomial coefficients"};
			
			GenericDialog zcDial = new GenericDialog("Z Calibration");
			zcDial.addChoice("Make z-calibration from: ",zcUseOptions, Prefs.get("SiMoLoc.ZC_Use", "Image stack"));
			//zcDial.addChoice("Degree of polynomial fit", zcPolDegreeOptions, Prefs.get("SiMoLoc.ZC_fitPolDegree", "1"));
			zcDial.addNumericField("Spacing between Z planes: ", Prefs.get("SiMoLoc.ZC_distBetweenPlanes", 20), 0, 5, " nm");
			zcDial.addNumericField("Upper R^2 threshold (0 - no filter): ", Prefs.get("SiMoLoc.ZC_zCalRsquareThreshold", 0), 2, 5, " ");
			zcDial.addCheckbox("Account for \"wobbling\" in X and Y", Prefs.get("SiMoLoc.bZCalWobbling", true));
			zcDial.setResizable(false);
			zcDial.showDialog();
			if (zcDial.wasCanceled())
	            return false;
			
			sZcUse = zcDial.getNextChoiceIndex();
			Prefs.set("SiMoLoc.ZC_Use", zcUseOptions[sZcUse]);
			zCalDistBetweenPlanes = (int)zcDial.getNextNumber();
			Prefs.set("SiMoLoc.ZC_distBetweenPlanes", zCalDistBetweenPlanes);
			zCalRsquareThreshold = zcDial.getNextNumber();
			Prefs.set("SiMoLoc.ZC_zCalRsquareThreshold", zCalRsquareThreshold);
			bZCalWobbling = zcDial.getNextBoolean();
			Prefs.set("SiMoLoc.bZCalWobbling", bZCalWobbling);
			
			return true;
		}
		
	/** Dialog showing options for particle search algorithm
	 * 		
	 * @return
	 */
	public boolean findParticles() {
		
		String [] DetectionType = new String [] {
				"Detect molecules (no fitting)","Detect molecules and fit", "Fit detected molecules"};
		GenericDialog fpDial = new GenericDialog("Find Particles");
		fpDial.addChoice("Task:", DetectionType, Prefs.get("SiMoLoc.DetectionType", "Detect molecules and fit"));
		fpDial.addNumericField("PSF standard devation", Prefs.get("SiMoLoc.dPSFsigma", 1.7), 2,4," pixels");
		//fpDial.addNumericField("Gaussial kernel size, \nodd number from 7(fast) till 13 (slow)  ", Prefs.get("SiMoLoc.nKernelSize", 7), 0);		
		fpDial.addNumericField("Pixel size", Prefs.get("SiMoLoc.dPixelSize", 66), 2, 4, "nm");
		fpDial.addNumericField("Parallel threads number", Prefs.get("SiMoLoc.nThreads", 500), 0);
		fpDial.addNumericField("Fitting iterations number", Prefs.get("SiMoLoc.nIterations", 3), 0);
		fpDial.addCheckbox("Mark detected particles in overlay?", Prefs.get("SiMoLoc.bShowParticles", false));
		fpDial.addCheckbox("Ignore false positives?", Prefs.get("SiMoLoc.bIgnoreFP", false));
		
		fpDial.setInsets(15, 20, 0); // extra space on top
		
		//*GPU part -> rewrite to new format of table
		//fpDial.addCheckbox("Accelerate using GPU", Prefs.get("SiMoLoc.bUseGPUAcceleration", false));
		//fpDial.addNumericField("Batch size", Prefs.get("SiMoLoc.nBatchSize", 4096), 0);//, 6, "Max : "); //TODO: get maximum value of the GPU
		//fpDial.addNumericField("Group size", Prefs.get("SiMoLoc.nGroupSize", 128), 0);//, 6, "Max : ");  //TODO: get maximum value of the GPU
		
		
		//!!! TODO fpDial.addCheckbox("Use log MLE instead of Chi^2", Prefs.get("SiMoLoc.bUseMLE", false));
		
		fpDial.setResizable(false);
		fpDial.showDialog();
		if (fpDial.wasCanceled())
            return false;
		
		nDetectionType = fpDial.getNextChoiceIndex();
		Prefs.set("SiMoLoc.DetectionType", DetectionType[nDetectionType]);
		dPSFsigma = fpDial.getNextNumber();
		Prefs.set("SiMoLoc.dPSFsigma", dPSFsigma);
		
		//calculate Gaussian kernel dimensions in pixels for spot detection
		nKernelSize = (int) Math.ceil(3.0*dPSFsigma);
		if(nKernelSize%2 == 0)
			 nKernelSize++;

		dPixelSize = fpDial.getNextNumber();
		Prefs.set("SiMoLoc.dPixelSize", dPixelSize);
		nThreads = (int) fpDial.getNextNumber();
		Prefs.set("SiMoLoc.nThreads", nThreads);
		nIterations = (int)fpDial.getNextNumber();
		Prefs.set("SiMoLoc.nIterations", nIterations);
		bShowParticles = fpDial.getNextBoolean();
		Prefs.set("SiMoLoc.bShowParticles", bShowParticles);
		bIgnoreFP = fpDial.getNextBoolean();
		Prefs.set("SiMoLoc.bIgnoreFP", bIgnoreFP);
		
		//*GPU part -> rewrite to new format of table
		//bUseGPUAcceleration = fpDial.getNextBoolean();
		//Prefs.set("SiMoLoc.bUseGPUAcceleration", bUseGPUAcceleration);
		//nBatchSize = (int)fpDial.getNextNumber();
		//Prefs.set("SiMoLoc.nBatchSize", nBatchSize);
		//nGroupSize = (int)fpDial.getNextNumber();
		//Prefs.set("SiMoLoc.nGroupSize", nGroupSize);
		//bUseMLE = false;// TODO: set to fpDial.getNextBoolean();
		//Prefs.set("SiMoLoc.bUseMLE", bUseMLE);
		
		return true;
	}
	/** Dialog showing options for reconstruction image
	 * 
	 * @param xlocavg_ average localization precision x
	 * @param ylocavg_ average localization precision y
	 * @param fminframe minimal frame
	 * @param fmaxframe maximal frame 
	 * @param xmax largest x
	 * @param ymax largest y
	 * @return user choice
	 */
	
	public boolean ReconstructImage(double xlocavg_, double ylocavg_, double fminframe, double fmaxframe, double xmax, double ymax) 		
	{
		int LutChoice;
		GenericDialog dgReconstruct = new GenericDialog("Reconstruct Dataset");
		//String [] RecIntOptions = new String [] {
				//"Normalized probability", "Integrated spot intensity","Amplitude of Gaussian fitting"};
		String [] RecSDOptions = new String [] {
				"Localization precision","Constant value"};
		String [] RecFPOptions = new String [] {
				"Only true positives", "All particles"};
		String [] Rec3DOptions = new String [] {
				"Z-stack", "Colorcoded (added colors)", "Colorcoded (max proj)"};
		String [] sLabelsCheckbox = new String [] {
				"Calculate z-values based on calibration","Scaled colorcode for depth"};
		boolean [] sLabelsDefault = new boolean [] {
				Prefs.get("SiMoLoc.recalZvalues", false), Prefs.get("SiMoLoc.bDynamicZscale", false)};
		String [] luts = IJ.getLuts();
		
		dgReconstruct.addChoice("For reconstruction use:", RecFPOptions, Prefs.get("SiMoLoc.Rec_FP", "Only true positives"));
		dgReconstruct.addMessage("Average localization precision in X: " + new DecimalFormat("#.##").format(xlocavg_) + " nm, in Y: " +  new DecimalFormat("#.##").format(ylocavg_) +" nm.");
		dgReconstruct.addNumericField("Pixel size of reconstructed image", Prefs.get("SiMoLoc.Rec_PixSize", 30), 2,6,"nm");
		dgReconstruct.addNumericField("Width of original image", xmax,0, 5, " pixels");
		dgReconstruct.addNumericField("Height of original image", ymax, 0, 5, " pixels");
		dgReconstruct.addChoice("SD of spots:", RecSDOptions, Prefs.get("SiMoLoc.Rec_SD", "Localization precision"));
		dgReconstruct.addNumericField("Value of SD in case of constant:", Prefs.get("SiMoLoc.Rec_SDFixed", 60), 2,6," nm");
		dgReconstruct.addCheckbox("Cut-off for localization precision:", Prefs.get("SiMoLoc.applycutoff", false));
		dgReconstruct.addNumericField("Cut particles off with localization less than: ", Prefs.get("SiMoLoc.cutoff", 20), 2,4," nm ");
		dgReconstruct.addMessage("\n");
		
		
		//dgReconstruct.addMessage("\n");
		dgReconstruct.addCheckbox("Reconstruct with translation (in nm):", Prefs.get("SiMoLoc.bTranslate", false));
		dgReconstruct.addNumericField("X_Offset:", Prefs.get("SiMoLoc.dTransX", 0), 2,6, "nm");
		dgReconstruct.addNumericField("Y_Offset:", Prefs.get("SiMoLoc.dTransY", 0), 2,6, "nm");
		//dgReconstruct.addMessage("\n");
		dgReconstruct.addCheckbox("Use frame interval:", Prefs.get("SiMoLoc.bFramesInterval", false));
		if(Prefs.get("SiMoLoc.bFramesInterval", false))
			dgReconstruct.addStringField("Range:", Prefs.get("SiMoLoc.sFrameRange", "1 - 2"));		
		else
			dgReconstruct.addStringField("Range:", new DecimalFormat("#").format(fminframe) + "-" +  new DecimalFormat("#").format(fmaxframe));		
		dgReconstruct.addMessage("\n");
		dgReconstruct.addCheckbox("Average localizations in consecutive frames within 1 pixel?", Prefs.get("SiMoLoc.bAveragePositions", false));
		dgReconstruct.addCheckbox("Update Results table with average localizations?", Prefs.get("SiMoLoc.bUpdateAveragePositions", false));
		
		dgReconstruct.addMessage("~~~~~~~~");
		dgReconstruct.addCheckbox("3D-reconstruction", Prefs.get("SiMoLoc.Render3D", false));
		dgReconstruct.addChoice("Render as:", Rec3DOptions, Prefs.get("SiMoLoc.n3DRenderType", "Z-stack"));
		dgReconstruct.addNumericField("Z-distance between slices (Z-stack):", Prefs.get("SiMoLoc.distZSlices", 100), 0,4," (nm)");
		dgReconstruct.addCheckboxGroup(1,2,sLabelsCheckbox,sLabelsDefault);
		dgReconstruct.addChoice("LUT for colorcode render:",luts,Prefs.get("SiMoLoc.zLutChoice","Fire"));
		
		//dgReconstruct.addCheckbox("Calculate z-values based on calibration", Prefs.get("SiMoLoc.recalZvalues", false));
		//dgReconstruct.addMessage("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		
		dgReconstruct.showDialog();
		if (dgReconstruct.wasCanceled())
            return false;
		
		nRecParticles = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.Rec_FP", RecFPOptions[nRecParticles]);
        dRecPixelSize = dgReconstruct.getNextNumber();
        Prefs.set("SiMoLoc.Rec_PixSize", dRecPixelSize);
        nRecWidth = (int) dgReconstruct.getNextNumber();
        Prefs.set("SiMoLoc.Rec_ImWidth", nRecWidth);
        nRecHeight = (int) dgReconstruct.getNextNumber();
        Prefs.set("SiMoLoc.Rec_ImHeight", nRecHeight);
		
		nSDIndex = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.Rec_SD", RecSDOptions[nSDIndex]);
		dFixedSD = dgReconstruct.getNextNumber();
	    Prefs.set("SiMoLoc.Rec_SDFixed", dFixedSD);
		bCutoff = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.applycutoff", bCutoff);
		dcutoff = dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.cutoff", dcutoff);
		
		
		bTranslation = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.bTranslate", bTranslation);
		dTranslationX =  dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.dTransX", dTranslationX);
		dTranslationY =  dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.dTransY", dTranslationY);
		
		bFramesInterval = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.bFramesInterval", bFramesInterval);			
			
		//range of frames		
		String sFrameRange = dgReconstruct.getNextString();
		Prefs.set("SiMoLoc.sFrameRange", sFrameRange);	
		String[] range = Tools.split(sFrameRange, " -");
		double c1 = dgReconstruct.parseDouble(range[0]);
		double c2 = range.length==2?dgReconstruct.parseDouble(range[1]):Double.NaN;
		nFrameMin = Double.isNaN(c1)?fminframe:(int)c1;
		nFrameMax = Double.isNaN(c2)?nFrameMin:(int)c2;
		if (nFrameMin<fminframe) nFrameMin = fminframe;
		if (nFrameMax>fmaxframe) nFrameMax = fmaxframe;
		if (nFrameMin>nFrameMax) {nFrameMin=fminframe; nFrameMax=fmaxframe;}				
		
		bAveragePositions = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.bAveragePositions", bAveragePositions);
		
		
		bUpdateAveragePositions= dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.bUpdateAveragePositions", bUpdateAveragePositions);
		
				
		//values for 3D reconstruction
		b3D = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.Render3D", b3D);
		n3DRenderType = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.n3DRenderType", Rec3DOptions[n3DRenderType]);
		dDistBetweenZSlices =  dgReconstruct.getNextNumber();
		Prefs.set("SiMoLoc.distZSlices", dDistBetweenZSlices);
		bCalculateZValues = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.recalZvalues", bCalculateZValues);
		bDynamicZscale = dgReconstruct.getNextBoolean();
		Prefs.set("SiMoLoc.bDynamicZscale", bDynamicZscale);	
		LutChoice = dgReconstruct.getNextChoiceIndex();
		Prefs.set("SiMoLoc.zLutChoice", luts[LutChoice]);
		sZLUTName = luts[LutChoice];
		
		return true;
	}
	/** Dialog showing parameters for color correction
	 * 
	 * 
	 */
	
	public boolean dglColorCalibration() 
	{
		GenericDialog ccDial = new GenericDialog("Color Calibration parameters");
		ccDial.addNumericField("Max distance between particle images in both channels:",Prefs.get("SiMoLoc.dCCDist", 7),0,3," original pixels");
		ccDial.addNumericField("SNR threshold filter:",Prefs.get("SiMoLoc.dCCSNR", 10),1,4," ");
		ccDial.addCheckbox("Show points used for calibration", Prefs.get("SiMoLoc.bCCShowPoints", true));
		ccDial.addCheckbox("Show final distortion map", Prefs.get("SiMoLoc.bCCShowMap", true));
		ccDial.addCheckbox("Provide reference points coordinates in Results", Prefs.get("SiMoLoc.bCCShowResults", true));
		
		ccDial.showDialog();
		if (ccDial.wasCanceled())
            return false;
		
		dCCDist =  ccDial.getNextNumber();
		Prefs.set("SiMoLoc.dCCDist", dCCDist);
		dCCSNR =  ccDial.getNextNumber();
		Prefs.set("SiMoLoc.dCCSNR", dCCSNR);
		bCCShowPoints = ccDial.getNextBoolean();
		Prefs.set("SiMoLoc.bCCShowMap", bCCShowPoints);
		bCCShowMap = ccDial.getNextBoolean();
		Prefs.set("SiMoLoc.bCCShowMap", bCCShowMap);
		bCCShowResults = ccDial.getNextBoolean();
		Prefs.set("SiMoLoc.bCCShowMap", bCCShowResults);
		return true;
	
	}
	/** Dialog showing parameters for drift correction
	 * 
	 * 
	 */
	public boolean DriftCorrection(double xlocavg_, double ylocavg_, double fminframe, double fmaxframe, double xmax, double ymax) 		
	{
		
		
		GenericDialog dgDriftCorrection = new GenericDialog("Calculate drift correction (correlation based)");
		dgDriftCorrection.addMessage("Average localization precision in X: " + new DecimalFormat("#.##").format(xlocavg_) + " nm, in Y: " +  new DecimalFormat("#.##").format(ylocavg_) +" nm.");
		dgDriftCorrection.addNumericField("Pixel size for intermediate reconstruction: ", Prefs.get("SiMoLoc.nDriftScale", 10), 2,6,"nm");
		dgDriftCorrection.addNumericField("Window size:", Prefs.get("SiMoLoc.drift_frames", 1000), 0,6," frames");
		//dgDriftCorrection.addNumericField("Maximum shift between windows:", Prefs.get("SiMoLoc.drift_nm", 50), 0, 3, " nm" );
		dgDriftCorrection.addCheckbox("Apply correction to Results Table", Prefs.get("SiMoLoc.drift_apply", true));
		dgDriftCorrection.addCheckbox("Show corrected image (using settings from Reconstruct)", Prefs.get("SiMoLoc.drift_reconstruct", true));
		dgDriftCorrection.addMessage("~~~~~~~~~~~~~~~~~~~~~");
		dgDriftCorrection.addCheckbox("Show intermediate reconstructions", Prefs.get("SiMoLoc.drift_intermediate_reconstr", false));		
		//dgDriftCorrection.addCheckbox("Show cross-correlation maps", Prefs.get("SiMoLoc.drift_cross_correlation", false));
		dgDriftCorrection.addMessage("~~~~~~~~~~~~~~~~~~~~~~~");
		//Get and show values from "Reconstruct window"
		String sRenderParameters="Rendering parameters (from Reconstruct menu)\n";

        //Take values of cut-off and false/true positives and SD rendering from ImageJ REGISTRY
        //also FRAME INTERVAL!!!
        String sParticles = Prefs.get("SiMoLoc.Rec_FP", "Only true positives");
        if(sParticles.equals("Only true positives"))
        {
        	nRecParticles=0;
        	sRenderParameters=sRenderParameters+"Use: only true positives\n";
        }
        else
        {
        	nRecParticles=1;
        	sRenderParameters=sRenderParameters+"Use: all particles\n";
        	
        }
        
        bCutoff=Prefs.get("SiMoLoc.applycutoff", false);
        if(!bCutoff)
        	sRenderParameters=sRenderParameters+"Cut-off localizations: off\n";        	
        else
        {
        	dcutoff=Prefs.get("SiMoLoc.cutoff", 20);
        	sRenderParameters=sRenderParameters+"Cut-off localizations: on, by "+new DecimalFormat("#.##").format(dcutoff)+" nm\n";
        }
        
        String sSDoption = Prefs.get("SiMoLoc.Rec_SD", "Localization precision");        
        if(sSDoption.equals("Localization precision"))
        {
        	sRenderParameters=sRenderParameters+ "SD of spots: " + sSDoption + "\n";
        	nSDIndex=0;
        }
        else
        {        	
        	nSDIndex=1;
        	dFixedSD=Prefs.get("SiMoLoc.Rec_SDFixed", 60);
        	sRenderParameters=sRenderParameters+ "SD of spots: " + sSDoption +", "+new DecimalFormat("#").format(dFixedSD)+ " nm \n";
        }
        
        
        
        bFramesInterval=Prefs.get("SiMoLoc.bFramesInterval", false);
        if(bFramesInterval)
        {
        	String sFrameRange = Prefs.get("SiMoLoc.sFrameRange", "1 - 2");
        	
        	String[] range = Tools.split(sFrameRange, " -");
			double c1 = dgDriftCorrection.parseDouble(range[0]);
			double c2 = range.length==2?dgDriftCorrection.parseDouble(range[1]):Double.NaN;
			nFrameMin = Double.isNaN(c1)?fminframe:(int)c1;
			nFrameMax = Double.isNaN(c2)?nFrameMin:(int)c2;
			if (nFrameMin<fminframe) nFrameMin = fminframe;
			if (nFrameMax>fmaxframe) nFrameMax = fmaxframe;
			if (nFrameMin>nFrameMax) {nFrameMin=fminframe; nFrameMax=fmaxframe;}	
			sRenderParameters=sRenderParameters+ "Frame Interval: " + sFrameRange + " \n";
        }

		dgDriftCorrection.addMessage(sRenderParameters);		
		
		dgDriftCorrection.showDialog();
		if (dgDriftCorrection.wasCanceled())
            return false;
		
		nDriftScale = (double) dgDriftCorrection.getNextNumber();
		Prefs.set("SiMoLoc.nDriftScale", nDriftScale);
		//Prefs.set("SiMoLoc.drift_frames", nDriftFrames);
		//nDriftParticles = dgDriftCorrection.getNextChoiceIndex();
		//Prefs.set("SiMoLoc.drift_pattype", RecFPOptions[nDriftParticles]);
		
		nDriftFrames = (int) dgDriftCorrection.getNextNumber();
		Prefs.set("SiMoLoc.drift_frames", nDriftFrames);
		
		//nDriftMaxDistnm = dgDriftCorrection.getNextNumber();
		//nDriftPixels = (int) Math.ceil(nDriftMaxDistnm /nDriftScale);
		//Prefs.set("SiMoLoc.drift_nm", nDriftMaxDistnm);


		bDriftUpdateTable = dgDriftCorrection.getNextBoolean();
		Prefs.set("SiMoLoc.drift_apply", bDriftUpdateTable);
		bDriftMakeReconstruction = dgDriftCorrection.getNextBoolean();
		Prefs.set("SiMoLoc.drift_reconstruct", bDriftMakeReconstruction);
		
		
		bShowIntermediate = dgDriftCorrection.getNextBoolean();
		Prefs.set("SiMoLoc.drift_intermediate_reconstr", bShowIntermediate);
		bShowCrossCorrelation = false;
		//bShowCrossCorrelation = dgDriftCorrection.getNextBoolean();
		//Prefs.set("SiMoLoc.drift_cross_correlation", bShowCrossCorrelation);
        
		// width and height of image
		nRecWidth = xmax + xlocavg_*3.0;
        nRecHeight = ymax + ylocavg_*3.0;
        
        //make sure we are in "Drift correction" mode, i.e. table sorting will happen
        bDrift = true;
        
        //if needed final reconstruction
        dRecPixelSize=Prefs.get("SiMoLoc.Rec_PixSize", 30);
        
		return true;
	}
	
	
	
	/** Dialog showing options for linking particles after detection
	 * 
	 * 
	 */
	
	public boolean LinkParticles() 
	{
		GenericDialog dgLink = new GenericDialog("Link Particles");
		String [] Link_Dist = new String [] {
				"Initial position", "Next detected position"};
		String [] LinkFPOptions = new String [] {
				"Only true positives", "All particles"};
		
		dgLink.addChoice("For linking use:", LinkFPOptions, Prefs.get("SiMoLoc.Link_FP", "All particles"));
		dgLink.addNumericField("Max distance to search over one frame:", Prefs.get("SiMoLoc.LinkDist", 5), 2,4," original pixels");
		dgLink.addChoice("Measure distance from:", Link_Dist, Prefs.get("SiMoLoc.LinkTrace", "Initial position"));
		dgLink.addNumericField("Maximum linking gap in frames:", Prefs.get("SiMoLoc.LinkFrames", 0), 0);
		dgLink.addCheckbox("Display tracks in overlay?", Prefs.get("SiMoLoc.bShowTracks", true));
		dgLink.addCheckbox("Show detected particles?", Prefs.get("SiMoLoc.bShowParticlesLink", false));
		dgLink.showDialog();
		if (dgLink.wasCanceled())
            return false;
		
		nLinkFP = dgLink.getNextChoiceIndex();
		Prefs.set("SiMoLoc.Link_FP", LinkFPOptions[nLinkFP]);
		dLinkDistance= dgLink.getNextNumber();
		Prefs.set("SiMoLoc.LinkDist", dLinkDistance);
		nLinkTrace = dgLink.getNextChoiceIndex();
		Prefs.set("SiMoLoc.LinkTrace", Link_Dist[nLinkTrace]);
		nLinkFrameGap= (int)dgLink.getNextNumber();
		Prefs.set("SiMoLoc.LinkFrames", nLinkFrameGap);
		bShowTracks = dgLink.getNextBoolean();
		Prefs.set("SiMoLoc.bShowTracks", bShowTracks);
		bShowParticlesLink = dgLink.getNextBoolean();
		Prefs.set("SiMoLoc.bShowParticlesLink", bShowParticlesLink);
		return true;		
	}


}
