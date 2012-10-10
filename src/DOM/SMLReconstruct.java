package DOM;


import java.awt.Frame;
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
	ImageProcessor ip;
	//ResultsTable table;
	SMLDialog settings;
		
	double [] amp;
	double [] x;
	double [] y;
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
	
	
	//parameters of drift correction
	double dDriftMagn; //magnification of reconstructed images used for drift correction (average precision)
	int drift_width;  //width of reconstructed images used for drift correction
	int drift_height; //width of reconstructed images used for drift correction
	int nIntervalsCount; //number of intervals (and images) used for drift correction
	

	
	//constructor for offline image reconstruction 
	
	SMLReconstruct(java.lang.String title, SMLAnalysis sml_, SMLDialog dlg_)
	{
		
		settings = dlg_;
		double pixelsize = sml_.ptable.getValueAsDouble(3, 0)/sml_.ptable.getValueAsDouble(1, 0);
		settings.dMagnification = pixelsize/settings.dRecPixelSize;
		
		new_width  = (int) (settings.nRecWidth*settings.dMagnification);
		new_height = (int) (settings.nRecHeight*settings.dMagnification);

		imp = new ImagePlus();
		imp.setTitle("Reconstructed image");
		
		
		// load data
		sml_.ptable_lock.lock();
		if(settings.nIntIndex==0) //integrated spot intensity as amplitude
			amp = sml_.ptable.getColumnAsDoubles(10);
		else 					  //fitted Gaussian peak amplitude as amplitude
			amp = sml_.ptable.getColumnAsDoubles(0);
		
		x   = sml_.ptable.getColumnAsDoubles(1);		
		y   = sml_.ptable.getColumnAsDoubles(2);
		fp  = sml_.ptable.getColumnAsDoubles(6);
		sdx = sml_.ptable.getColumnAsDoubles(7);
		sdy = sml_.ptable.getColumnAsDoubles(8);
		f   = sml_.ptable.getColumnAsDoubles(13);
		sml_.ptable_lock.unlock();
		
		nParticlesCount = f.length;
		// load max & min values
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
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;


		IJ.showStatus("Reconstructing Image...");
		
		for (int n=0;n<nParticlesCount;n++)
		{
			if (f[n]>=fstart && f[n]<=fstop && fp[n]<0.3)
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
					
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak-0.5)/sdxmag) - ErrorFunction.erf2((0.5+i-xpeak)/sdxmag);
								dErry = ErrorFunction.erf2((j-ypeak-0.5)/sdymag) - ErrorFunction.erf2((0.5+j-ypeak)/sdymag);
								new_i = old_i + preG*amp[n]*sdxmag*sdymag*dErrx*dErry;
						
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
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;
		
		n=0;
		while(f[n]<fstart)
			n++;
		
		while (bContinue) 
		{
			if (fp[n]<0.3)
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
					
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<drift_width) && (j<drift_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak-0.5)/sdxmag) - ErrorFunction.erf2((0.5+i-xpeak)/sdxmag);
								dErry = ErrorFunction.erf2((j-ypeak-0.5)/sdymag) - ErrorFunction.erf2((0.5+j-ypeak)/sdymag);
								new_i = old_i + preG*amp[n]*sdxmag*sdymag*dErrx*dErry;
						
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
		nIntervalsCount = (int) Math.floor(((double)nframes)/((double)nIntervalFrames));
		//the rest of the table is truncated for now
			
		driftx = new double [nIntervalsCount];
		drifty = new double [nIntervalsCount];
		//reconstructing images for cross-correlation calculation
		
		IJ.showStatus("Applying Drift correction (images reconstruction)...");
		for (i=0; i<nIntervalsCount; i++)
		{
			IJ.showProgress(i, nIntervalsCount);
			draw_sorted(i*nIntervalFrames+1, (i+1)*nIntervalFrames);			
		}
		if(settings.bShowIntermediate)
			new ImagePlus("Intermediate reconstructions (Drift Correction)", driftstack).show();
		
		IJ.showStatus("Applying Drift correction (calculating cross-correlation)...");
		for (i=1; i<nIntervalsCount; i++)
		{
			IJ.showProgress(i-1, nIntervalsCount);			
			crosscorrstack.addSlice(null, crosscorrelation(driftstack.getProcessor(i),driftstack.getProcessor(i+1)));
		}
		if(settings.bShowCrossCorrelation)
			new ImagePlus("Cross Correlation (Drift Correction)", crosscorrstack).show();
		
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
			sDriftData = sDriftData + Integer.toString(i+1)+"\t"+Double.toString(approxx[i])+"\t"+Double.toString(approxy[i])+"\n";
		
		Frame frame = WindowManager.getFrame("Drift Correction");
		if (frame!=null && (frame instanceof TextWindow) )
		{
			DriftTable = (TextWindow)frame;
			DriftTable.getTextPanel().clear();
			DriftTable.getTextPanel().append(sDriftData);
			DriftTable.getTextPanel().updateDisplay();			
		}
			else
				DriftTable = new TextWindow("Drift Correction", "Frame_Number\tX_drift_(px)\tY_drift_(px)", sDriftData, 450, 300);			
		
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
		for(i=0; i<(int)(nIntervalFrames*0.5); i++)
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
		
		for(i=0;f[i]<nIntervalsCount*nIntervalFrames;i++)
		{
			dIndex = (int) f[i];
			x[i]+=approxx[dIndex];
			y[i]+=approxy[dIndex];
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
}





/** Makes drift correction based on center of mass calculation.
 * Using sliding average. Turns out to be not so good idea. */	
/*void correctDriftCOM()
{
	int nBlocks;
	int nDriftFr;
	int nCurrFrame;
	double dTotInt;
	double dCMx, dCMy;
	int nFrameTimer;
	int i,j;

	double max_x = -100;
	double max_y = -100;
	double min_x = 1000000000;
	double min_y = 1000000000;
	
	nDriftFr = settings.nDriftFrames;
	
	if (((double)nframes)%((double)nDriftFr)!=0)
		nBlocks = (int) Math.ceil(((double)nframes)/((double)nDriftFr));
	else
		nBlocks = nframes/nDriftFr;
	
	//center of masses 
	double [] cmx = new double [nBlocks];
	double [] cmy = new double [nBlocks];
	
	IJ.showStatus("Applying Drift Correction...");	
	
	nCurrFrame = (int)f[0];
	nFrameTimer = nDriftFr;
	dCMx = 0; dCMy = 0; dTotInt = 0;
	j=0;
	for(i=0;i<f.length; i++)
	{
		if(nCurrFrame != (int)f[i])
		{
			nFrameTimer--;
			nCurrFrame = (int)f[i];				
		}
		if(nFrameTimer == 0)
		{
			nFrameTimer = nDriftFr;
			cmx[j]=dCMx/dTotInt;
			cmy[j]=dCMy/dTotInt;
			j++;
			dCMx = 0; dCMy = 0; dTotInt = 0;				
		}
		if(fp[i]<0.3)
		{
			//dCMx += (amp[i]-min)*x[i]/nrange;				
			//dCMy += (amp[i]-min)*y[i]/nrange;
			//dTotInt += (amp[i]-min)/nrange;
			dCMx += x[i];				
			dCMy += y[i];
			dTotInt ++;

		}
	}

	cmx[nBlocks-1]=dCMx/dTotInt;
	cmy[nBlocks-1]=dCMy/dTotInt;

	for( i = 1;i<nBlocks; i++)
	{
		dCMx = cmx[i]-cmx[0];
		dCMy = cmy[i]-cmy[0];
		cmx[i]=dCMx;
		cmy[i]=dCMy;
		if(dCMx>max_x) max_x=dCMx;
		if(dCMy>max_y) max_y=dCMy;
		if(dCMx<min_x) min_x=dCMx;
		if(dCMy<min_y) min_y=dCMy;
	}
	cmx[0]=0;
	cmy[0]=0;
	
	//correcting
	nCurrFrame = (int)f[0];
	nFrameTimer = nDriftFr;
	j=0;
	for(i=0;i<f.length; i++)
	{
			if(nCurrFrame != (int)f[i])
			{
				nFrameTimer--;
				nCurrFrame = (int)f[i];				
			}
			if(nFrameTimer == 0)
			{
				nFrameTimer = nDriftFr;
				j++;
			
			}
			x[i]+=cmx[j];
			y[i]+=cmy[j];	
	}
		
}
*/
/** Makes drift correction based on center of mass calculation.
 * Using sliding average. Turns out to be not so good idea. */	
/*void correctDriftCOMSliding()
{

	int i,j;
	int nBlocks;
	int nDriftFr;
	double dTotInt;
	double dCMx, dCMy;
	double max_x = -100;
	double max_y = -100;
	double min_x = 1000000000;
	double min_y = 1000000000;
	boolean bFirstSwitch = true;
	int nFirstFrameIndex;
	int nCurrFrame;
	int nFrameTimer;
	 
	nDriftFr = settings.nDriftFrames;
	
	nBlocks = 1 + nframes - nDriftFr;
	double [] cmx = new double [nBlocks];
	double [] cmy = new double [nBlocks];
	IJ.showStatus("Applying Drift Correction...");
	
	nFirstFrameIndex = 0;						
	for(i = 0;i<nBlocks; i++)
	{
		IJ.showProgress(i, nBlocks);
		j = nFirstFrameIndex;
		nFrameTimer = nDriftFr;
		nCurrFrame = (int)f[j];
		bFirstSwitch = true;
		dCMx =0; dCMy = 0; dTotInt = 0;
		while(nFrameTimer>0)
		{
				if(fp[j]<0.3)
				{
					dCMx += (amp[j]-min)*x[j]/nrange;				
					dCMy += (amp[j]-min)*y[j]/nrange;
					dTotInt += (amp[j]-min)/nrange;
					//dCMx += x[j];				
					//dCMy += y[j];
					//dTotInt ++;

				}
				j++;
				if(j==f.length)
					nFrameTimer=0;
				else
				{
					if(nCurrFrame != (int)f[j])
					{
						nFrameTimer--;
						nCurrFrame = (int)f[j];
						if(bFirstSwitch)
						{
							bFirstSwitch = false;
							nFirstFrameIndex = j;
						}
					}					
				}
		}
		cmx[i]=dCMx/dTotInt;
		cmy[i]=dCMy/dTotInt;
		//if(cmx[i]>max_x) max_x=cmx[i];
		//if(cmy[i]>max_y) max_y=cmy[i];
		//if(cmx[i]<min_x) min_x=cmx[i];
		//if(cmy[i]<min_y) min_y=cmy[i];
	}
	
	for( i = 1;i<nBlocks; i++)
	{
		dCMx = cmx[i]-cmx[0];
		dCMy = cmy[i]-cmy[0];
		cmx[i]=dCMx;
		cmy[i]=dCMy;
		if(dCMx>max_x) max_x=dCMx;
		if(dCMy>max_y) max_y=dCMy;
		if(dCMx<min_x) min_x=dCMx;
		if(dCMy<min_y) min_y=dCMy;
	}
	cmx[0]=0;
	cmy[0]=0;
	
	//finding first frame to correct
	i=0;
	while((int)f[i]<=nDriftFr)
		i++;
	nCurrFrame= (int)f[i];
	j=1;
	for( ;i<f.length; i++)
	{
		if(nCurrFrame!= (int)f[i])
		{
			nCurrFrame = (int)f[i];
			j++;
		}
		
		x[i]+=cmx[j];
		y[i]+=cmy[j];			
	}
}
*/

