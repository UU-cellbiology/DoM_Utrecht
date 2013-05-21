package DOM;


import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
//import ij.process.ShortProcessor;
import ij.process.ShortProcessor;
import ij.text.TextWindow;

public class SMLReconstruct {

	ImagePlus imp;
	ImageStack driftstack;
	ImageStack crosscorrstack;
	ImageStack zstack;
	ImageProcessor ip;
	//ResultsTable table;
	SMLDialog settings;
		
	double [] amp;
	double [] x;
	double [] y;
	double [] z_nm;
	double [] sdx;
	double [] sdy;
	double [] f;
	double [] fp;
	double [] driftx;
	double [] drifty;
	double [] approxx;
	double [] approxy;

	double max=0;
	double min=9999999;
	double nrange=0;
	int nframes = 0;
	int nParticlesCount = 0;

	int new_width;
	int new_height;
	boolean bFrameSorted;
	boolean bNormalized;
	double dFPThreshold; //threshold that determines what particles are used for reconstruction
	
	
	//parameters of drift correction
	double dDriftMagn; //magnification of reconstructed images used for drift correction (average precision)
	int drift_width;  //width of reconstructed images used for drift correction
	int drift_height; //width of reconstructed images used for drift correction
	int nIntervalsCount; //number of intervals (and images) used for drift correction
	

	
	//constructor for image reconstruction 
	
	SMLReconstruct(java.lang.String title, SMLAnalysis sml_, SMLDialog dlg_)
	{
		
		settings = dlg_;
		bFrameSorted = false;
		bNormalized = false;
		double pixelsize = sml_.ptable.getValueAsDouble(3, 0)/sml_.ptable.getValueAsDouble(1, 0);
		settings.dMagnification = pixelsize/settings.dRecPixelSize;
		
		new_width  = (int) (settings.nRecWidth*settings.dMagnification);
		new_height = (int) (settings.nRecHeight*settings.dMagnification);

		imp = new ImagePlus();
		imp.setTitle("Reconstructed image");
		
		//threshold for particles reconstruction (include false positives or not)
		switch (dlg_.nRecParticles)
		{ //only true positives
			case 0:
				dFPThreshold = 0.3;
				break;
		 //true and half positives
			case 1:
				dFPThreshold = 0.7;
				break;
		 //all particles
			case 2:
				dFPThreshold = 1.1;
				break;
			default:
				dFPThreshold = 0.3;
				break;
		
		}
		
		
		// load data
		sml_.ptable_lock.lock();
		if(settings.nIntIndex==0 || settings.nIntIndex==1) //integrated spot intensity as amplitude
			amp = sml_.ptable.getColumnAsDoubles(10);
		else 					  //fitted Gaussian peak amplitude as amplitude
			amp = sml_.ptable.getColumnAsDoubles(0);
		//use normalized 2D gaussian to plot each particle 
		if(settings.nIntIndex==0)
			bNormalized = true;
		
		//frame number
		f   = sml_.ptable.getColumnAsDoubles(13);
		nParticlesCount = f.length;
		//coordinates
		x   = sml_.ptable.getColumnAsDoubles(1);		
		y   = sml_.ptable.getColumnAsDoubles(2);

		//calculate z-values if necessary
		if(settings.bCalculateZValues)
			calcZValues(sml_);
		
		z_nm = sml_.ptable.getColumnAsDoubles(5);
		
		if (settings.bTranslation)			
		{
			for (int i=0; i<nParticlesCount; i++)
			{
				x[i] += settings.dTranslationX;
				y[i] += settings.dTranslationY;
			}
		
		}
		//false positive mark
		fp  = sml_.ptable.getColumnAsDoubles(6);
		//localization precision
		sdx = sml_.ptable.getColumnAsDoubles(7);
		sdy = sml_.ptable.getColumnAsDoubles(8);
		
		sml_.ptable_lock.unlock();
		
		
		// load min & max values
		for (int n=0;n<f.length;n++)
		{
			if (f[n]>nframes) nframes=(int) f[n];
			if (amp[n]>max) max=amp[n];
			if (amp[n]<min) min=amp[n];
		}
		nrange = max-min;
	}
	
	/** Reconstruction drawing function used by the "Reconstruct Image" plugin.
	 *  Slow one since it assumes that Results table is not sorted. 
	 * @param fstart show only particle after this frame
	 * @param fstop show only particle before this frame
	*/
	void draw_unsorted(int fstart, int fstop)
	{	
		FloatProcessor ipf = new FloatProcessor(new_width, new_height);
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, sdxmag, sdymag;
		int xmag, ymag, xsq, ysq;
		int i,j;
		double preG=0.25*Math.PI;
		double dCutoff = 1.0;
		double dNorm = 1.0;
		//double dVerif = 0.0;
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;


		IJ.showStatus("Reconstructing Image...");
		
		for (int n=0;n<nParticlesCount;n++)
		{
			if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold)
			{
				IJ.showProgress(n, nParticlesCount);
				xmag=(int) Math.round(x[n]*settings.dMagnification);
				ymag=(int) Math.round(y[n]*settings.dMagnification);
				

				if(sdx[n]<dCutoff && sdy[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						sdx[n] = settings.dFixedSD;
						sdy[n] = settings.dFixedSD;
					}
					xsq = (int) Math.round(3*sdx[n]*settings.dMagnification);
					ysq = (int) Math.round(3*sdy[n]*settings.dMagnification);
					xpeak=x[n]*settings.dMagnification;
					ypeak=y[n]*settings.dMagnification;
					sdxmag=sdx[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					sdymag=sdy[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					if(bNormalized)
					{
						dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/sdxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/sdxmag);
						dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/sdymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/sdymag);
						dNorm = 1/(dErrx*dErry);
					}
					//dVerif=0;
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak)/sdxmag) - ErrorFunction.erf2((1+i-xpeak)/sdxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/sdymag) - ErrorFunction.erf2((1+j-ypeak)/sdymag);
								if(bNormalized)
									{new_i = old_i + dNorm*dErrx*dErry;}
									//dVerif+=dNorm*dErrx*dErry;}
								else
									{new_i = old_i + preG*amp[n]*sdxmag*sdymag*dErrx*dErry;}
								
								ipf.setf(i, j, (float)new_i);
							}
						}
				}		
			}
		}

		imp.setProcessor(ipf.convertToShort(true));
		IJ.run(imp, "Set Scale...", "distance=1 known="+settings.dRecPixelSize+" pixel=1 unit=nm");
		//imp.addSlice(null, ipf.convertToShort(true));
		IJ.run(imp, "Enhance Contrast", "saturated=0.35");
	}
	
	/** Reconstruction drawing function used by the "Reconstruct Image" plugin.
	 *  Slow one since it assumes that Results table is not sorted. 
	 * @param fstart show only particle after this frame
	 * @param fstop show only particle before this frame
	 * @param fstep distance between z-slices
	 * @param 
	*/
	void draw_zstack(int fstart, int fstop, double zstep)
	{	
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, sdxmag, sdymag;
		int xmag, ymag, xsq, ysq;
		int i,j;
		double preG=0.25*Math.PI;
		double dCutoff = 1.0;
		double dNorm = 1.0;
		//double dVerif = 0.0;

		
		int nSlices;		//total number of slices
		int sliceNumber;	//slice number of particle
		double[] zValues;	//array with cloned z values
		double zmin, zmax;
		
		//create new array to sort without changing original
		zValues = z_nm.clone();
		Arrays.sort(zValues);
		//zmin = zValues[0];
		zmin = 0;
		zmax = zValues[zValues.length-1];
		
		
		/*DEBUGGING!*/
		zmin = 0;
		zmax = 800;
		/*DEBUGGING!*/
		
		nSlices = (int) Math.ceil((zmax-zmin)/zstep);
		//in case all z-values are zero still one slice has to be drawn
		if(nSlices==0)
			nSlices++;
		
		FloatProcessor[] ipf = new FloatProcessor[nSlices];
		
		for(int k = 0; k<nSlices; k++){
			ipf[k] = new FloatProcessor(new_width, new_height);
		}
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;
		

		IJ.showStatus("Reconstructing Z-stack...");
		
		for (int n=0;n<nParticlesCount;n++)
		{
			if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold)
			{
				IJ.showProgress(n, nParticlesCount);
				xmag=(int) Math.round(x[n]*settings.dMagnification);
				ymag=(int) Math.round(y[n]*settings.dMagnification);
				
				//calculate sliceNumber for this particle
				sliceNumber = (int) Math.floor((z_nm[n]-zmin)/zstep);
				
				/*DEBUGGING*/
				if(sliceNumber>=nSlices)
					sliceNumber = nSlices-1;
				/*DEBUGGING*/
				
				if(sdx[n]<dCutoff && sdy[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						sdx[n] = settings.dFixedSD;
						sdy[n] = settings.dFixedSD;
					}
					xsq = (int) Math.round(3*sdx[n]*settings.dMagnification);
					ysq = (int) Math.round(3*sdy[n]*settings.dMagnification);
					xpeak=x[n]*settings.dMagnification;
					ypeak=y[n]*settings.dMagnification;
					sdxmag=sdx[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					sdymag=sdy[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					if(bNormalized)
					{
						dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/sdxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/sdxmag);
						dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/sdymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/sdymag);
						dNorm = 1/(dErrx*dErry);
					}
					//dVerif=0;
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
							{
								old_i=ipf[sliceNumber].getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak)/sdxmag) - ErrorFunction.erf2((1+i-xpeak)/sdxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/sdymag) - ErrorFunction.erf2((1+j-ypeak)/sdymag);
								if(bNormalized)
									{new_i = old_i + dNorm*dErrx*dErry;}
									//dVerif+=dNorm*dErrx*dErry;}
								else
									{new_i = old_i + preG*amp[n]*sdxmag*sdymag*dErrx*dErry;}
								
								ipf[sliceNumber].setf(i, j, (float)new_i);
							}
						}
				}		
			}
		}
		zstack = new ImageStack(new_width, new_height);
		
		for(int k=0; k<nSlices; k++){
			zstack.addSlice(null, ipf[k].convertToShort(true));
		}
		
		new ImagePlus("Z-stack (slices@"+zstep+"nm", zstack).show();

		//IJ.run(imp, "Set Scale...", "distance=1 known="+settings.dRecPixelSize+" pixel=1 unit=nm");
		//IJ.run(imp, "Enhance Contrast", "saturated=0.35");
	}
	
	/** Reconstruction drawing function used by the "Reconstruct Image" plugin during drift correction.
	 *  Faster one since it assumes that Results table is sorted by frame number. 
	 *  Puts resulting reconstructions to driftstack variable.
	 *  Uses different magnification determined for drift correction only
	 * @param fstart show only particle after this frame
	 * @param fstop show only particle before this frame
	*/
	void draw_sorted(int fstart, int fstop)
	{	
		FloatProcessor ipf = new FloatProcessor(drift_width, drift_height);
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, sdxmag, sdymag;
		int xmag, ymag, xsq, ysq;
		int i,j,n;
		double preG=0.25*Math.PI;
		boolean bContinue = true;
		double dCutoff = 1.0;
		double dNorm = 1.0;
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;
		
		n=0;
		while(f[n]<fstart)
			n++;
		
		while (bContinue) 
		{
			if (fp[n]<dFPThreshold)
			{
				
				xmag=(int) Math.round(x[n]*dDriftMagn);
				ymag=(int) Math.round(y[n]*dDriftMagn);
				
				if(sdx[n]<dCutoff && sdy[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						sdx[n] = settings.dFixedSD;
						sdy[n] = settings.dFixedSD;
					}
					xsq = (int) Math.round(3*sdx[n]*dDriftMagn);
					ysq = (int) Math.round(3*sdy[n]*dDriftMagn);
					xpeak=x[n]*dDriftMagn;
					ypeak=y[n]*dDriftMagn;
					sdxmag=sdx[n]*dDriftMagn*1.41421356; //number is just sqrt(2)
					sdymag=sdy[n]*dDriftMagn*1.41421356; //number is just sqrt(2)
					if(bNormalized)
					{
						dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/sdxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/sdxmag);
						dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/sdymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/sdymag);
						dNorm = 1/(dErrx*dErry);
					}					
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<drift_width) && (j<drift_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak)/sdxmag) - ErrorFunction.erf2((1+i-xpeak)/sdxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/sdymag) - ErrorFunction.erf2((1+j-ypeak)/sdymag);
								if(bNormalized)
									{new_i = old_i + dNorm*dErrx*dErry;}								
								else
									{new_i = old_i + preG*amp[n]*sdxmag*sdymag*dErrx*dErry;}
								
								
						
								ipf.setf(i, j, (float)new_i);
							}
						}
				}	
			}
			n++;
			if(n>=nParticlesCount)
				bContinue = false;
			else
			  if(f[n]>fstop)
				  bContinue = false;
		}

		
		driftstack.addSlice(null, ipf.convertToShort(true));
		
	}
	
	/** Main function performing drift correction  */
	void DriftCorrection()
	{
		double dXAver, dYAver; //average precision
		int i;
		int nCount=0;
		int nIntervalFrames;
		int [] xymax;
		String sDriftData ="";
		TextWindow DriftTable; 
		
		nIntervalFrames = settings.nDriftFrames;
		
		
		IJ.showStatus("Applying Drift correction...");
		
		//determine magnification of reconstructed images 
		dXAver=0; dYAver=0;
		for (i=0;i<nParticlesCount;i++)
		{
			if(fp[i]<0.3 && sdx[i]<1.0 && sdy[i]<1.0)
			{
				nCount++;
				dXAver+=sdx[i];
				dYAver+=sdy[i];
			}
		}
		dXAver /=nCount;
		dYAver /=nCount;
		dDriftMagn = 1.0/Math.max(dXAver, dYAver);
		drift_width  = (int) (settings.nRecWidth*dDriftMagn);
		drift_height = (int) (settings.nRecHeight*dDriftMagn);
		driftstack = new ImageStack(drift_width, drift_height);
		crosscorrstack = new ImageStack(2*settings.nDriftPixels+1, 2*settings.nDriftPixels+1);
		//calculating number of time intervals
		if(settings.bFramesInterval)
			nIntervalsCount = (int) Math.floor((settings.nFrameMax-settings.nFrameMin+1.0)/((double)nIntervalFrames));
		else
			nIntervalsCount = (int) Math.floor(((double)nframes)/((double)nIntervalFrames));
		
		if(nIntervalsCount <= 1)
		{
			IJ.error("Drift correction interval is too big. Drift correction is cancelled.");
			return;
		}
	
		//the rest of the table is truncated for now			
		driftx = new double [nIntervalsCount];
		drifty = new double [nIntervalsCount];
		
		//reconstructing images for cross-correlation calculation		
		IJ.showStatus("Applying Drift correction (images reconstruction)...");
		for (i=0; i<nIntervalsCount; i++)
		{
			IJ.showProgress(i, nIntervalsCount);
			if(settings.bFramesInterval)
				draw_sorted((int)settings.nFrameMin+i*nIntervalFrames, (int)settings.nFrameMin -1 +(i+1)*nIntervalFrames);
			else
				draw_sorted(i*nIntervalFrames+1, (i+1)*nIntervalFrames);			
		}
		if(settings.bShowIntermediate)
			new ImagePlus("Intermediate reconstructions (Drift frames="+settings.nDriftFrames+" pixels="+settings.nDriftPixels+")", driftstack).show();
		
		IJ.showStatus("Applying Drift correction (calculating cross-correlation)...");
		for (i=1; i<nIntervalsCount; i++)
		{
			IJ.showProgress(i-1, nIntervalsCount);			
			crosscorrstack.addSlice(null, crosscorrelation(driftstack.getProcessor(i),driftstack.getProcessor(i+1)));
		}
		if(settings.bShowCrossCorrelation)
			new ImagePlus("Cross Correlation (Drift frames="+settings.nDriftFrames+" pixels="+settings.nDriftPixels+")", crosscorrstack).show();
		
		driftx[0]=0;
		drifty[0]=0;
		 
		for (i=1; i<nIntervalsCount; i++)
		{
			xymax = getmaxpositions(crosscorrstack.getProcessor(i));
			driftx[i] = driftx[i-1]+((xymax[0]-settings.nDriftPixels)/dDriftMagn);
			drifty[i] = drifty[i-1]+((xymax[1]-settings.nDriftPixels)/dDriftMagn);
		}
		//applysimplecorrection();
		applylinearapproximationcorrection();		
		
		//cut out uncorrected frames
		nframes = nIntervalsCount*nIntervalFrames;
		
		//show results of drift
		for(i=0;i<nframes;i++)
		{
			if(settings.bFramesInterval)
			{
				sDriftData = sDriftData + Integer.toString(i+(int)settings.nFrameMin)+"\t"+Double.toString(approxx[i])+"\t"+Double.toString(approxy[i])+"\n";
			}
			else
			{
				sDriftData = sDriftData + Integer.toString(i+1)+"\t"+Double.toString(approxx[i])+"\t"+Double.toString(approxy[i])+"\n";				
			}
		}
		Frame frame = WindowManager.getFrame("Drift Correction (frames="+settings.nDriftFrames+" pixels="+settings.nDriftPixels+")");
		if (frame!=null && (frame instanceof TextWindow) )
		{
			DriftTable = (TextWindow)frame;
			DriftTable.getTextPanel().clear();
			DriftTable.getTextPanel().append(sDriftData);
			DriftTable.getTextPanel().updateDisplay();			
		}
			else
				DriftTable = new TextWindow("Drift Correction (frames="+settings.nDriftFrames+" pixels="+settings.nDriftPixels+")", "Frame_Number\tX_drift_(px)\tY_drift_(px)", sDriftData, 450, 300);			
		
		return;
		
	}

	void applysimplecorrection()
	{
		int i;
		int dIndex;
		int nIntervalFrames = settings.nDriftFrames;
		
		
		for(i=0;f[i]<nIntervalsCount*nIntervalFrames;i++)
		{
			dIndex = (int) Math.floor(((double)f[i])/((double)nIntervalFrames));
			x[i]+=driftx[dIndex];
			y[i]+=drifty[dIndex];
		}
			
		
	}
	
	void applylinearapproximationcorrection()
	{
		int i;
		int dIndex;
		int nIntervalFrames = settings.nDriftFrames;
		double coeffx, coeffy;
		
		
		approxx = new double [nIntervalsCount*nIntervalFrames];
		approxy = new double [nIntervalsCount*nIntervalFrames];
		
		
		//making linear approximation
		//first half interval
		for(i=0; i<nIntervalFrames*0.5; i++)
		{
			approxx[i] = 0.;
			approxy[i] = 0.;			
		}
		//middle part
		for(;i<nIntervalFrames*(nIntervalsCount-0.5); i++)
		{
			dIndex = (int) Math.floor((i-nIntervalFrames*0.5)/nIntervalFrames);
			coeffx = (driftx[dIndex+1]-driftx[dIndex])/nIntervalFrames;
			coeffy = (drifty[dIndex+1]-drifty[dIndex])/nIntervalFrames;
			approxx[i] = coeffx*(i-(dIndex+0.5)*nIntervalFrames) + driftx[dIndex];
			approxy[i] = coeffy*(i-(dIndex+0.5)*nIntervalFrames) + drifty[dIndex];			
			
		}
		dIndex = nIntervalsCount-1;
		//last part
		for(;i<nIntervalFrames*nIntervalsCount; i++)
		{
			approxx[i] = driftx[dIndex];
			approxy[i] = drifty[dIndex];		
		}
		
		if(settings.bFramesInterval)
		{
			i=0;
			while(f[i]<settings.nFrameMin)
				i++;
			for(;f[i]<(settings.nFrameMin+nIntervalsCount*nIntervalFrames-1);i++)
			{
				dIndex = (int) f[i]-(int)settings.nFrameMin;
				x[i]+=approxx[dIndex];
				y[i]+=approxy[dIndex];
			}						
		}
		else
		{
			for(i=0;f[i]<nIntervalsCount*nIntervalFrames;i++)
			{
				dIndex = (int) f[i];
				x[i]+=approxx[dIndex];
				y[i]+=approxy[dIndex];
			}
			
		}
		return;
			
		
	}
	
	ShortProcessor crosscorrelation (ImageProcessor ip1, ImageProcessor ip2)
	{
		int nMaxPix = settings.nDriftPixels;
		int i, j, m, n;
		int tot;
		int dMaskWidth, dMaskHeight;
		double dCC;
		
		
		FloatProcessor resultIP;
		ShortProcessor returnIP;
		ImageProcessor driftmask;
		
		tot = 2*nMaxPix+1;
		resultIP = new FloatProcessor(tot, tot);
		dMaskWidth = drift_width-tot+1;
		dMaskHeight = drift_height-tot+1;
		
		ip2.setRoi(nMaxPix,nMaxPix,dMaskWidth, dMaskHeight);
		driftmask = ip2.crop();
		ip2.resetRoi();
		
		for (i=0;i<tot;i++)
			for (j=0;j<tot;j++)
			{
				dCC=0;
				for(m=0; m<dMaskWidth; m++)
					for(n=0; n<dMaskHeight; n++)
						dCC+=ip1.get(m+i,n+j)*driftmask.get(m, n);
				resultIP.setf(i,j,(float)dCC);
			}
		
		returnIP = (ShortProcessor) resultIP.convertToShort(true);
		returnIP.smooth();
		return returnIP;		
	}
	

	
	/** Sorts particles data by frame number, ascending. */	
	void sortbyframe()
	{
		int i, nSize;
		nSize = f.length;
	
		//put all data together in one array
		double [][] data = new double[nSize][7];
		for (i=0; i<nSize; i++)
		{
			data[i][0] = f[i];
			data[i][1] = amp[i];
			data[i][2] = x[i];
			data[i][3] = y[i];
			data[i][4] = sdx[i];
			data[i][5] = sdy[i];
			data[i][6] = fp[i];
		}
		//sort
		Arrays.sort(data, new tableComparator(0, true));
		//put data back
		for (i=0; i<nSize; i++)
		{
			  f[i] = data[i][0] ;
			amp[i] = data[i][1];
			  x[i] = data[i][2];
		  	  y[i] = data[i][3];
			sdx[i] = data[i][4];
			sdy[i] = data[i][5];
			 fp[i] = data[i][6];
		}		
		bFrameSorted = true;
	}
	
	/** Cleans the reconstruction viewer image. */
	void clear()
	{
	
		for (int i=0;i<ip.getWidth();i++)
			for (int j=0;j<ip.getHeight();j++)
				ip.set(i, j, 0);
	
	}
	/** Smoothes correlation images with simple gaussian filter */
	/*void convolvecorrelation()
	{
		GaussianBlur lowpassGauss = new GaussianBlur();
		for (int i=1; i<nIntervalsCount;i++)
		lowpassGauss.blurGaussian(crosscorrstack.getProcessor(i), 1.0, 1.0, 0.03);
		
		
	}*/
	
	/** Gets x and y coordinates of maximum intensity pixel on image. */	
	int [] getmaxpositions(ImageProcessor ipp)
	{
		int [] results = new int [2];
		results[0]=0;
		results[1]=0;

		int s = 0, smax=-10;
		
		for (int i=0;i<ipp.getWidth();i++)
		{
			for (int j=0;j<ipp.getHeight();j++)
			{
				s=ipp.get(i, j);	
				if (s>smax)
				{
					smax=s;
					results[0]=i;
					results[1]=j;
				}
			}
		}
		return results;		
	}
	
	/** Function averages localizations in consecutive frames within 1 pixel. */	
	void averagelocalizations()
	{
		int nUniqFrames,i;
		int nPatNumber = f.length;
		int nCurrFrame, nCount, nFrameCount, nCurrentParticle, nCompareFrame, nCompareParticle,nFrameCompareCount;
		int nAveragesNumber;
		int [][] framesstat;
		double [][] amp_un;
		double [][] x_un;
		double [][] y_un;
		double [][] sdx_un;
		double [][] sdy_un;
		double [][] fp_un;
		boolean bContinue;
		double dDistance, dAverX, dAverY, dWeightX, dWeightY, dSumWeightX, dSumWeightY, dSumAmp, dWeightAmp, dSumWeightAmp ;
		ArrayList<double[]> Averaging = new ArrayList<double[]>();
		
		IJ.showStatus("Averaging single photoactivation events...");
		//sort table
		if(!bFrameSorted)
			this.sortbyframe();
		
		//calculate number of particles in each frame and
		//store it for each frame in framesstat array
		nUniqFrames = uniqFrames();
		framesstat = new int [nUniqFrames][2];
		nCurrFrame = (int) f[0];
		nCount = 0;
		nFrameCount = 0;
		framesstat[nFrameCount][0]=nCurrFrame;
		for (i=0;i<nPatNumber;i++)
		{
			if((int)f[i]==nCurrFrame)
				nCount++;
			else
			{	
				framesstat[nFrameCount][1]=nCount;
				nCurrFrame = (int)f[i];
				nFrameCount++;
				framesstat[nFrameCount][0]=nCurrFrame;
				nCount = 1;				
			}
			
		}
		//last element
		framesstat[nFrameCount][1]=nCount;
		framesstat[nFrameCount][0]=nCurrFrame;
		
		//resorting data to new arrays for convenience
				
		//1) allocating arrays with uneven rows 
		//  for particles from each frame separately
		
		amp_un = new double [nUniqFrames][];
		x_un   = new double [nUniqFrames][];
		y_un   = new double [nUniqFrames][];
		sdx_un = new double [nUniqFrames][];
		sdy_un = new double [nUniqFrames][];
		fp_un  = new double [nUniqFrames][];
		
		for(i=0;i<nUniqFrames;i++)
		{
			amp_un[i] = new double[framesstat[i][1]];
			x_un[i]   = new double[framesstat[i][1]];
			y_un[i]   = new double[framesstat[i][1]];
			sdx_un[i] = new double[framesstat[i][1]];
			sdy_un[i] = new double[framesstat[i][1]];
			fp_un[i]  = new double[framesstat[i][1]];
		}
		
		//2) putting all data to those arrays
		nFrameCount = 0;
		nCurrentParticle = -1;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			nCurrentParticle++;
			if(nCurrentParticle == framesstat[nFrameCount][1])
			{
				nCurrentParticle=0;
				nFrameCount++;
				IJ.showProgress(nFrameCount, nUniqFrames);
			}
			amp_un[nFrameCount][nCurrentParticle] = amp[nCount];
			x_un[nFrameCount][nCurrentParticle]   = x[nCount];
			y_un[nFrameCount][nCurrentParticle]   = y[nCount];
			sdx_un[nFrameCount][nCurrentParticle] = sdx[nCount];
			sdy_un[nFrameCount][nCurrentParticle] = sdy[nCount];
			fp_un[nFrameCount][nCurrentParticle]  = fp[nCount];					
		}
		
		//all right!
		//now let's check single photoactivation events 
		nFrameCount = 0;
		nCurrentParticle = -1;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			nCurrentParticle++;
			if(nCurrentParticle == framesstat[nFrameCount][1])
			{
				nCurrentParticle=0;
				nFrameCount++;
			}
			
			
			if(fp_un[nFrameCount][nCurrentParticle]<dFPThreshold && nFrameCount<(nUniqFrames-1)) 
			{ 	
				nCurrFrame = framesstat[nFrameCount][0];
				//accumulating information
				Averaging.clear();
				
				bContinue = true;
				nFrameCompareCount = nFrameCount+1;
				nCompareFrame = framesstat[nFrameCompareCount][0];
				nCompareParticle = -1;
				nAveragesNumber = 0;
				while(bContinue) 
				{
					nCompareParticle++;
					if(nCompareParticle == framesstat[nFrameCompareCount][1]) //end of particles list in current frame
					{
						if (nFrameCompareCount == nUniqFrames - 1)	//reached the end of total particle list, stop						
							{bContinue = false;}
						else
						{//jumping to next frame
						 // but let's check whether we found something while scanning previous one?
							//no, we didn't, stop scanning
							if(nCompareFrame-nCurrFrame>nAveragesNumber)
								{bContinue = false;}
							//yes, we did, let's jump then
							else
							{
								nFrameCompareCount++;
								nCompareFrame = framesstat[nFrameCompareCount][0];
								nCompareParticle = 0;
							}													
						}
							
					}
					//ok, let's see if all checks are passed
					if(bContinue && fp_un[nFrameCompareCount][nCompareParticle]<dFPThreshold)
					{//seems so. let's measure distances between molecules then
						dDistance = Math.sqrt(Math.pow(x_un[nFrameCount][nCurrentParticle]-x_un[nFrameCompareCount][nCompareParticle], 2) +Math.pow(y_un[nFrameCount][nCurrentParticle]-y_un[nFrameCompareCount][nCompareParticle], 2)); 
						//they are closer then one pixel
						if(dDistance<1.0)
						{
							//found same photoactivation event
							//let's store values
							Averaging.add(new double[] {amp_un[nFrameCompareCount][nCompareParticle], x_un[nFrameCompareCount][nCompareParticle], y_un[nFrameCompareCount][nCompareParticle], sdx_un[nFrameCompareCount][nCompareParticle], sdy_un[nFrameCompareCount][nCompareParticle]});
							nAveragesNumber++;
							//mark it as already used (fp)
							//debug stub
							fp_un[nFrameCompareCount][nCompareParticle]+=1.2;//mark it as an averaged one (>1)
							//exit from current frame (a bit lame one)
							nCompareParticle = framesstat[nFrameCompareCount][1] -1;
						}
					}
				}
				//was photoactivation event found?
				if(nAveragesNumber>0)
				{
					//adding particle itself to the list of averaged values
					Averaging.add(new double[] {amp_un[nFrameCount][nCurrentParticle], x_un[nFrameCount][nCurrentParticle], y_un[nFrameCount][nCurrentParticle], sdx_un[nFrameCount][nCurrentParticle], sdy_un[nFrameCount][nCurrentParticle]});
					nAveragesNumber++;
					dAverX=0; dAverY=0; dSumWeightX=0; dSumWeightY=0; dSumAmp=0; dSumWeightAmp=0; dWeightAmp=0; dWeightX=0; dWeightY=0;
					//let's calculate weighted average 
					for(double [] item:Averaging)
					{
						dWeightX = 1/(item[3]*item[3]);
						dWeightY = 1/(item[4]*item[4]);
						dWeightAmp  = 0.5*(dWeightX+dWeightY);
						dAverX  += dWeightX* item[1];
						dAverY  += dWeightY* item[2];
						dSumWeightX += dWeightX;
						dSumWeightY += dWeightY;
						dSumWeightAmp += dWeightAmp;
						dSumAmp += dWeightAmp*item[0];												
					}
					//weighted averages final
					dAverX = dAverX/dSumWeightX;
					dAverY = dAverY/dSumWeightY;
					dSumAmp = dSumAmp/dSumWeightAmp;
					//let's update values for current particle 
					amp_un[nFrameCount][nCurrentParticle] = dSumAmp;
					x_un[nFrameCount][nCurrentParticle] = dAverX;
					y_un[nFrameCount][nCurrentParticle] = dAverY;
					sdx_un[nFrameCount][nCurrentParticle] = 1/Math.sqrt(dSumWeightX);
					sdy_un[nFrameCount][nCurrentParticle] = 1/Math.sqrt(dSumWeightY);
					fp_un[nFrameCount][nCurrentParticle]-=1.2; //mark it as an averaged)

				}
			}//false positive check
		}
		//updating initial tables		
		nFrameCount = 0;
		nCurrentParticle = -1;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			nCurrentParticle++;
			if(nCurrentParticle == framesstat[nFrameCount][1])
			{
				nCurrentParticle=0;
				nFrameCount++;
			}
			if(fp_un[nFrameCount][nCurrentParticle]>1)
			{
				fp[nCount]+=1.2;
			}
			if(fp_un[nFrameCount][nCurrentParticle]<0)
			{
				amp[nCount] = amp_un[nFrameCount][nCurrentParticle];
				x[nCount]   = x_un[nFrameCount][nCurrentParticle] ;
				y[nCount]   = y_un[nFrameCount][nCurrentParticle];
				sdx[nCount] = sdx_un[nFrameCount][nCurrentParticle];
				sdy[nCount] = sdy_un[nFrameCount][nCurrentParticle];
				fp[nCount]  = fp_un[nFrameCount][nCurrentParticle];
			}
		}
		//show table with updated results
		//show results of drift
	/*	String sAverData="";
		TextWindow AverTable;
				for(i=0;i<nPatNumber;i++)
				{
						sAverData = sAverData + Double.toString(f[i]) + "\t" +Double.toString(amp[i])+"\t"+Double.toString(fp[i])+"\t"+Double.toString(x[i])+"\t"+Double.toString(y[i])+"\t"+Double.toString(sdx[i])+"\t"+Double.toString(sdy[i])+"\n";
					
				}
				Frame frame = WindowManager.getFrame("After Averaging");
				if (frame!=null && (frame instanceof TextWindow) )
				{
					AverTable = (TextWindow)frame;
					AverTable.getTextPanel().clear();
					AverTable.getTextPanel().append(sAverData);
					AverTable.getTextPanel().updateDisplay();			
				}
					else
						AverTable = new TextWindow("After Averaging", "FrameNumber\tAmplitude\tFalsePositive\tX_px\tY_px\tLocX_px\tLocY_px", sAverData, 450, 300);
		*/
		
		
	}
	
	/** Function caclulating number of unique frames */
	int uniqFrames()
	{
		int i, nCount=1, nLength = f.length;
		int nCurrFrame;
		boolean bContinue = true;
		
		i=0;
		nCurrFrame = (int)f[i];		
		while (bContinue)
		{
			i++;
			if(i==nLength)
				bContinue = false;
			else
			{
				if((int)f[i]!=nCurrFrame)
				{
					nCurrFrame=(int)f[i];
					nCount++;
					
				}
			}
			
		}
		
		return nCount;
		
	}


	/**function calculating z-values based on calibration curve*/
	
	void calcZValues(SMLAnalysis sml_)
	{
		double[] SDx = sml_.ptable.getColumnAsDoubles(15);
		double[] SDy = sml_.ptable.getColumnAsDoubles(16);
		
		//double SDx;
		//double SDy;
		
		double z;
		
		//calibration values, according to calibration curve: calValues[0]+calValues[1]*(SDx-SDy)+calValues[2]*(SDx-SDy)^2
		double[] calValues = {389.4,-101.3,-4.075};
		/*
		for(int i=0; i<x.length; i++){
			SDx = sml_.ptable.getValueAsDouble(15, i);
			SDy = sml_.ptable.getValueAsDouble(16, i);
			z = calValues[0]+calValues[1]*(SDx-SDy)+calValues[2]*(SDx-SDy)*(SDx-SDy);
			sml_.ptable.setValue(5,i, z);
		}
		*/
		
		for(int i=0; i<SDx.length; i++){
			z = calValues[0]+calValues[1]*(SDx[i]-SDy[i])+calValues[2]*(SDx[i]-SDy[i])*(SDx[i]-SDy[i]);
			if(z<0)
				z=0;
			sml_.ptable.setValue(5,i, z);
		}
		
		sml_.ptable.updateResults();
		sml_.ptable.show("Results");
	}

}