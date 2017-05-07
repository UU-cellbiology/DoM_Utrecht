package DOM;


import java.awt.Color;
import java.awt.Frame;
import java.text.DecimalFormat;
import java.util.ArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.text.TextWindow;

public class SMLReconstruct {

	ImagePlus imp;
	/** stack with intermediate reconstruction for drift correction */
	ImageStack driftstack;
	/** stack containing cross correlation maps of intermediate reconstructions */
	ImageStack crosscorrstack;
	ImageStack zstack;
	ImageProcessor ip;
	//ResultsTable table;
	SMLDialog settings;
		
	
	/** x coordinates of particles in nm*/
	double [] x;
	/** y coordinates of particles in nm*/
	double [] y;
	/** z coordinates of particles in nm*/
	double [] z;
	/** x localization error (nm) */
	double [] loc_errx;
	/** y localization error (nm) */
	double [] loc_erry;
	/** z localization error (nm) */
	double [] loc_errz;
	/** frame numbers of particle*/
	double [] f;
	/** false positive mark */
	double [] fp;
	/** drift x coordinates between reconstructions */
	double [] driftx;
	/** drift y coordinates between reconstructions */
	double [] drifty;
	/** linearly interpolated drift x coordinates for each frame */
	double [] approxx;
	/** linearly interpolated drift x coordinates for each frame */
	double [] approxy;
	
	/** max Z value*/
	double zmax;
	/** min Z value*/
	double zmin;

	double max=0;
	double min=9999999;
	/** total number of frames in particle table */
	int nframes = 0;
	/** total number of particles in table */
	int nParticlesCount = 0;

	int new_width;
	int new_height;
	boolean bFrameSorted;

	/** threshold that determines what particles are used for reconstruction */
	double dFPThreshold; 
	
	
	//parameters of drift correction
	/** magnification of reconstructed images used for drift correction (average precision) */
	double dDriftMagn;
	/** width of reconstructed images used for drift correction */
	int drift_width;  
	/** height of reconstructed images used for drift correction */
	int drift_height; 
	/** number of intervals (and images) used for drift correction */
	int nIntervalsCount; 
	

	
	/**
	 * constructor for image reconstruction 
	 * @param title title of image
	 * @param sml_ corresponding analysis object
	 * @param dlg_  dialog
	 */
	
	SMLReconstruct(java.lang.String title, SMLAnalysis sml_, SMLDialog dlg_)
	{
		
		settings = dlg_;
		bFrameSorted = false;
		
		settings.dMagnification = 1.0/settings.dRecPixelSize;
		
		new_width  = (int) Math.ceil(settings.nRecWidth*settings.dMagnification);
		new_height = (int) Math.ceil(settings.nRecHeight*settings.dMagnification);

		imp = new ImagePlus();
		imp.setTitle(title);
		
		//threshold for particles reconstruction (include false positives or not)
		switch (dlg_.nRecParticles)
		{ //only true positives
			case 0:
				dFPThreshold = 0.3;
				break;
		 //all particles
			case 1:
				dFPThreshold = 1.1;
				break;
			default:
				dFPThreshold = 0.3;
				break;
		
		}
		
		//sort table by frame if required
		if(settings.bAveragePositions || settings.bDrift)
		{
			IJ.showStatus("Sorting table by frame number (averaging preparation)...");
			Sort_Results.sorting_external_silent(sml_, DOMConstants.Col_FrameN, true);		
			bFrameSorted = true;
			IJ.showStatus("Sorting table by frame number (averaging preparation)...done");
			
		}
		
		
		// load data
		sml_.ptable_lock.lock();
	
		//frame number
		f   = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		nParticlesCount = f.length;
		//coordinates
		x   = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);		
		y   = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);	
		z   = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_Znm);
		
		//false positive mark
		fp  = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_Fp);
		//localization precision
		loc_errx = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
		loc_erry = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);
		loc_errz = sml_.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errZ);
		
		sml_.ptable_lock.unlock();
		
		
		// find maximum frame number
		for (int n=0;n<f.length;n++)
		{
			if (f[n]>nframes) nframes=(int) f[n];
			
		}
		if(dlg_.b3D)
		{
			getZmaxmin();
		}
	
	}
	
	/** Reconstruction drawing function used by the "Reconstruct Image" plugin.
	 *  Rather slow one, since it assumes that Results table is not sorted. 
	 * @param fstart show only particle after this frame
	 * @param fstop show only particle before this frame
	*/
	void draw_unsorted(int fstart, int fstop)
	{	
		FloatProcessor ipf = new FloatProcessor(new_width, new_height);
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, loc_errxmag, loc_errymag;
		int xmag, ymag, xsq, ysq;
		int i,j;
		double dCutoff = 1000.0; //default cut-off is 1000 nm
		double dNorm = 1.0;
		double loc_errxlocal,loc_errylocal;
		
		
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
				

				if(loc_errx[n]<dCutoff && loc_erry[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						loc_errxlocal  = settings.dFixedSD;
						loc_errylocal = settings.dFixedSD;
					}
					else
					{
						loc_errxlocal=loc_errx[n];
						loc_errylocal=loc_erry[n];
					}
					xsq = (int) Math.ceil(3*loc_errxlocal*settings.dMagnification);
					ysq = (int) Math.ceil(3*loc_errylocal*settings.dMagnification);
					xpeak=x[n]*settings.dMagnification;
					ypeak=y[n]*settings.dMagnification;
					loc_errxmag=loc_errxlocal*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					loc_errymag=loc_errylocal*settings.dMagnification*1.41421356; //last number is just sqrt(2)

					///calculate area under pixelated gaussian
					dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/loc_errxmag);
					dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/loc_errymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/loc_errymag);
					dNorm = 1/(dErrx*dErry);

					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+i-xpeak)/loc_errxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/loc_errymag) - ErrorFunction.erf2((1+j-ypeak)/loc_errymag);
								
								new_i = old_i + dNorm*dErrx*dErry;
								ipf.setf(i, j, (float)new_i);
							}
						}
				}		
			}
		}

		//imp.setProcessor(ipf.convertToShort(true));
		imp.setProcessor(ipf);
		IJ.run(imp, "Set Scale...", "distance=1 known="+settings.dRecPixelSize+" pixel=1 unit=nm");
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
		double old_i, new_i, dErrx, dErry, dErrz, xpeak, ypeak, zpeak, loc_errxmag, loc_errymag,loc_errzmag;
		int xmag, ymag, zmag, xsq, ysq, zsq;
		int i,j,k;
		
		double dCutoff = 1000.0; //default cut-off is 1000 nm
		double dNorm;
		

		
		int nSlices;		//total number of slices
		int sliceNumber;	//slice number of particle
		
		//double zmin, zmax;
		double zMagnification =1/zstep;
		
			
	
		//zmin=Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		//zmax=Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
	
		
		nSlices = (int) Math.ceil((zmax-zmin)/zstep);
		
		//in case all z-values are zero still one slice has to be drawn
		if(nSlices==0)
			nSlices++;
		
		FloatProcessor[] ipf = new FloatProcessor[nSlices];
		
		for(k = 0; k<nSlices; k++)
		{
			ipf[k] = new FloatProcessor(new_width, new_height);
		}
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;
		

		IJ.showStatus("Reconstructing Z-stack...");
		
		for (int n=0;n<nParticlesCount;n++)
		{
			//if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold)
			if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold && z[n]>zmin && z[n]<zmax)
			{
				IJ.showProgress(n, nParticlesCount);
				xmag=(int) Math.round(x[n]*settings.dMagnification);
				ymag=(int) Math.round(y[n]*settings.dMagnification);
				zmag=(int) Math.round(z[n]*zMagnification);
				
				//calculate sliceNumber for this particle
				//sliceNumber = (int) Math.floor((z[n]-zmin)/zstep);
				
				/*DEBUGGING*/
				//if(sliceNumber>=nSlices)
				//	sliceNumber = nSlices-1;
				/*DEBUGGING*/
				
				if(loc_errx[n]<dCutoff && loc_erry[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						loc_errx[n] = settings.dFixedSD;
						loc_erry[n] = settings.dFixedSD;
					}
					xsq = (int) Math.ceil(3*loc_errx[n]*settings.dMagnification);
					ysq = (int) Math.ceil(3*loc_erry[n]*settings.dMagnification);
					zsq = (int) Math.ceil(3*loc_errz[n]*zMagnification);
					xpeak=x[n]*settings.dMagnification;
					ypeak=y[n]*settings.dMagnification;
					zpeak= z[n]*zMagnification;
					loc_errxmag=loc_errx[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					loc_errymag=loc_erry[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					loc_errzmag=loc_errz[n]*zMagnification*1.41421356; //last number is just sqrt(2)
			
					dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/loc_errxmag);
					dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/loc_errymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/loc_errymag);
					dErrz = ErrorFunction.erf2((zmag-zsq-zpeak)/loc_errzmag) - ErrorFunction.erf2((1+zmag+zsq-zpeak)/loc_errzmag);
					dNorm = 1/(dErrx*dErry*dErrz);
			
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
							for(k=zmag-zsq;k<zmag+zsq+1;k++)
							{							
								if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
								{
									sliceNumber = (int)Math.round(k-(zmin/zstep));
									if(sliceNumber>=0 && sliceNumber<nSlices)
									{
										old_i=ipf[sliceNumber].getf(i, j);
				
										dErrx = ErrorFunction.erf2((i-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+i-xpeak)/loc_errxmag);
										dErry = ErrorFunction.erf2((j-ypeak)/loc_errymag) - ErrorFunction.erf2((1+j-ypeak)/loc_errymag);	
										dErrz = ErrorFunction.erf2((k-zpeak)/loc_errzmag) - ErrorFunction.erf2((1+k-zpeak)/loc_errzmag);
										new_i = old_i + dNorm*dErrx*dErry*dErrz;
										ipf[sliceNumber].setf(i, j, (float)new_i);
									}
								}
							}
				}		
			}
		}
		zstack = new ImageStack(new_width, new_height);
		
		for(k=0; k<nSlices; k++)
		{
			//zstack.addSlice(null, ipf[k].convertToShort(true));
			zstack.addSlice(null, ipf[k]);
		}
		imp = new ImagePlus("Z-stack" , zstack);
		//IJ.run(imp, "Set Scale...", "distance=1 known="+settings.dRecPixelSize+" pixel=1 unit=nm");
		String sProp = "channels=1 slices="+Integer.toString(nSlices) + " frames=1 unit=nm pixel_width=";
		sProp = sProp + Double.toString(settings.dRecPixelSize)+"  pixel_height="+ Double.toString(settings.dRecPixelSize);
		sProp = sProp + " voxel_depth=" + Double.toString(settings.dDistBetweenZSlices);
		IJ.run(imp, "Properties...", sProp);
	}
	
	/** Renders image in 3D color-coded mode
	 * @param fstart show only particle after this frame
	 * @param fstop show only particle before this frame
	 * @param 
	 * @param 
	*/
	void draw_colorcodedZ(int fstart, int fstop)
	{	

		Color c;

		int [] dZcolorScale = getZcumHistogram();
		
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, loc_errxmag, loc_errymag;
		int xmag, ymag, xsq, ysq;
		int i,j;
		
		double dCutoff = 1000.0; //default cut-off is 1000 nm
		double dNorm;		
		//double zmin, zmax;
		//zmin=zmin;
		//zmax=zmax;
		double zhue;

		int zOldInd;
		int zInd;
	
		float [][] zLUT;
	
		zLUT=getHSBLutTable();
	
		//zmin=Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		//zmax=Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		
		ColorProcessor imcol = new ColorProcessor(new_width, new_height);
		FloatProcessor[] ipf = new FloatProcessor[2];
		
		for(int k = 0; k<2; k++)
		{
			ipf[k] = new FloatProcessor(new_width, new_height);
		}
		
		if(settings.bCutoff)
			dCutoff = settings.dcutoff;
		

		IJ.showStatus("Reconstructing colorcoded Z projection...");
		
		for (int n=0;n<nParticlesCount;n++)
		{
			//if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold)
			if (f[n]>=fstart && f[n]<=fstop && fp[n]<dFPThreshold && z[n]>zmin && z[n]<zmax)
			{
				IJ.showProgress(n, nParticlesCount);
				xmag=(int) Math.round(x[n]*settings.dMagnification);
				ymag=(int) Math.round(y[n]*settings.dMagnification);
				

				
				if(loc_errx[n]<dCutoff && loc_erry[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						loc_errx[n] = settings.dFixedSD;
						loc_erry[n] = settings.dFixedSD;
					}
					xsq = (int) Math.ceil(3*loc_errx[n]*settings.dMagnification);
					ysq = (int) Math.ceil(3*loc_erry[n]*settings.dMagnification);

					xpeak=x[n]*settings.dMagnification;
					ypeak=y[n]*settings.dMagnification;
					loc_errxmag=loc_errx[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
					loc_errymag=loc_erry[n]*settings.dMagnification*1.41421356; //last number is just sqrt(2)
			
					dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/loc_errxmag);
					dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/loc_errymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/loc_errymag);
					dNorm = 1/(dErrx*dErry);
			
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<new_width) && (j<new_height)&&(i>0)&&(j>0))
							{
								
		
								dErrx = ErrorFunction.erf2((i-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+i-xpeak)/loc_errxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/loc_errymag) - ErrorFunction.erf2((1+j-ypeak)/loc_errymag);					
								
								zhue = (z[n]-zmin)/(zmax-zmin);
							
								
								zInd = dZcolorScale[(int)Math.round(zhue*255)];
	
								if(settings.n3DRenderType==2)
								{
									old_i=ipf[0].getf(i, j);
									new_i = dNorm*dErrx*dErry;
									if(new_i>old_i)
									{
										ipf[0].setf(i, j, (float)(new_i));
										ipf[1].setf(i, j, (float)(zInd));
									}
									
								}
								else
								{
							
									old_i=ipf[0].getf(i, j);
									new_i = dNorm*dErrx*dErry;
									ipf[0].setf(i, j, (float)(new_i+old_i));
									//calculate new color index
									//by weighted average
									zOldInd=(int)Math.round(ipf[1].getf(i, j));
									zhue = (zInd*new_i+zOldInd*old_i)/(new_i+old_i);
									ipf[1].setf(i, j, (int)Math.round(zhue));
								}
						
							}
						}
				}		
			}
		}
		
	
		int [] newrgb = new int[3];
		//generate final RGB image
		double intMax = ipf[0].getMax();
		for(i=0;i<new_width;i++)
			for(j=0;j<new_height;j++)
			{
				if(ipf[0].getf(i, j)>0)
				{
					
					zInd = (int)Math.round(ipf[1].getf(i, j));
					c = Color.getHSBColor(zLUT[zInd][0], zLUT[zInd][1], (float)(ipf[0].getf(i, j)/intMax));
					newrgb[0]=c.getRed();
					newrgb[1]=c.getGreen();
					newrgb[2]=c.getBlue();
					imcol.putPixel(i, j, newrgb);
				}
				
			}
		
		//show colorbar
		ColorProcessor imcolcode = new ColorProcessor(256, 40);
		for(i=0;i<256;i++)
			for(j=0;j<25;j++)
			{
				zInd = dZcolorScale[i];
				c = Color.getHSBColor(zLUT[zInd][0], zLUT[zInd][1], 1f);								
				newrgb[0]=c.getRed();
				newrgb[1]=c.getGreen();
				newrgb[2]=c.getBlue();
				imcolcode.putPixel(i, j, newrgb);
			}
		newrgb[0]=255;
		newrgb[1]=255;
		newrgb[2]=255;
		for(i=0;i<256;i++)
			for(j=25;j<40;j++)
				imcolcode.putPixel(i, j, newrgb);
		imcolcode.setJustification(ImageProcessor.LEFT_JUSTIFY);
		imcolcode.drawString(new DecimalFormat("#").format(zmin)+ " nm",0,40);
		imcolcode.setJustification(ImageProcessor.RIGHT_JUSTIFY);
		imcolcode.drawString(new DecimalFormat("#").format(zmax)+ " nm",250,40);		
		new ImagePlus("Z colorbar" , imcolcode).show();
		
	
		IJ.showProgress(nParticlesCount, nParticlesCount);
		
		imp = new ImagePlus("Colorcoded" , imcol);
		IJ.run(imp, "Set Scale...", "distance=1 known="+settings.dRecPixelSize+" pixel=1 unit=nm");
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
		double old_i, new_i, dErrx, dErry, xpeak, ypeak, loc_errxmag, loc_errymag;
		int xmag, ymag, xsq, ysq;
		int i,j,n;
		double loc_errxlocal,loc_errylocal;
	
		boolean bContinue = true;
		double dCutoff = 1000.0; //default cut-off is 1000 nm
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
				
				if(loc_errx[n]<dCutoff && loc_erry[n]<dCutoff)
				{
					if(settings.nSDIndex>0)
					{
						loc_errxlocal  = settings.dFixedSD;
						loc_errylocal = settings.dFixedSD;
					}
					else
					{
						loc_errxlocal=loc_errx[n];
						loc_errylocal=loc_erry[n];
					}
					xsq = (int) Math.ceil(3*loc_errxlocal*dDriftMagn);
					ysq = (int) Math.ceil(3*loc_errylocal*dDriftMagn);
					xpeak=x[n]*dDriftMagn;
					ypeak=y[n]*dDriftMagn;
					loc_errxmag=loc_errxlocal*dDriftMagn*1.41421356; //number is just sqrt(2)
					loc_errymag=loc_errylocal*dDriftMagn*1.41421356; //number is just sqrt(2)
					//if(bNormalized)
					//{
					dErrx = ErrorFunction.erf2((xmag-xsq-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+xmag+xsq-xpeak)/loc_errxmag);
					dErry = ErrorFunction.erf2((ymag-ysq-ypeak)/loc_errymag) - ErrorFunction.erf2((1+ymag+ysq-ypeak)/loc_errymag);
					dNorm = 1/(dErrx*dErry);
					//}					
					for(i=xmag-xsq;i<xmag+xsq+1;i++)
						for(j=ymag-ysq;j<ymag+ysq+1;j++)
						{
							if((i<drift_width) && (j<drift_height)&&(i>0)&&(j>0))
							{
								old_i=ipf.getf(i, j);
		
								dErrx = ErrorFunction.erf2((i-xpeak)/loc_errxmag) - ErrorFunction.erf2((1+i-xpeak)/loc_errxmag);
								dErry = ErrorFunction.erf2((j-ypeak)/loc_errymag) - ErrorFunction.erf2((1+j-ypeak)/loc_errymag);							
								new_i = old_i + dNorm*dErrx*dErry;
						
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

		
		//driftstack.addSlice(null, ipf.convertToShort(true));
		driftstack.addSlice(null, ipf);
		
	}
	
	/** Main function performing correlation based drift correction  */
	void DriftCorrection()
	{
				
		int i;

		int nIntervalFrames;
		//int [] xymax;
		String sDriftData ="";
		TextWindow DriftTable; 
		
		//some variables for timing
		long startTime;
		long reconstructionTime=0;
		long fullTime;
		double [] driftxfull;
		double [] driftyfull;
		
		nIntervalFrames = settings.nDriftFrames;
		
		
		IJ.showStatus("Applying Drift correction...");


		dDriftMagn = 1.0/settings.nDriftScale;

		drift_width  = (int) (settings.nRecWidth*dDriftMagn);
		drift_height = (int) (settings.nRecHeight*dDriftMagn);
		driftstack = new ImageStack(drift_width, drift_height);
		crosscorrstack = new ImageStack(drift_width, drift_height);
		
		//calculating number of time intervals
		if(settings.bFramesInterval)
			nIntervalsCount = (int) Math.floor((settings.nFrameMax-settings.nFrameMin+1)/((double)nIntervalFrames));
		else
			nIntervalsCount = (int) Math.floor(((double)nframes)/((double)nIntervalFrames));

		
		if(nIntervalsCount <= 1)
		{
			IJ.error("Drift correction window is larger than total number of frames. Drift correction is cancelled.");
			return;
		}
	
		
		driftx = new double [nIntervalsCount];
		drifty = new double [nIntervalsCount];
		driftxfull =  new double [nframes];
		driftyfull =  new double [nframes];
		//let's start measuring time		
		startTime = System.nanoTime();
		
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
		IJ.showProgress(nIntervalsCount, nIntervalsCount);
		
		reconstructionTime = System.nanoTime() - startTime;
		IJ.log("Intermediate reconstructions time: " + String.format("%.2f",((double)Math.abs(reconstructionTime))*0.000000001) + " s");
		//show them, if asked
		if(settings.bShowIntermediate)
			new ImagePlus("Intermediate reconstructions (Drift frames="+settings.nDriftFrames+" px size="+String.format("%d",((int)settings.nDriftScale))+" nm)", driftstack).show();
			//new ImagePlus("Intermediate reconstructions (Drift frames="+settings.nDriftFrames+" max shift="+String.format("%d",((int)settings.nDriftMaxDistnm))+" nm)", driftstack).show();
		
		
		IJ.showStatus("Applying Drift correction (calculating cross-correlation)...");
		/*
		for (i=1; i<nIntervalsCount; i++)
		{
			IJ.showProgress(i-1, nIntervalsCount);			
			crosscorrstack.addSlice(null, crosscorrelation(driftstack.getProcessor(i),driftstack.getProcessor(i+1)));
		}
		IJ.showProgress(nIntervalsCount, nIntervalsCount);
		
		fullTime = System.nanoTime() - startTime;
		IJ.log("Cross correlation calculation time: " + String.format("%.2f", ((double)Math.abs(fullTime-reconstructionTime))*0.000000001)+ " s");
		IJ.log("Total time: " + String.format("%.2f",((double)Math.abs(fullTime))*0.000000001) + " s");

		if(settings.bShowCrossCorrelation)
			new ImagePlus("Cross Correlation (Drift frames="+settings.nDriftFrames+" max shift="+String.format("%d",((int)settings.nDriftMaxDistnm))+" nm)", crosscorrstack).show();
		*/
		
		driftx[0]=0;
		drifty[0]=0;
		ImCrossCorrelation imCrCorr = new ImCrossCorrelation();
		double [] xymaxd;
		//define pixel with highest cross correlation value on correlation map
		for (i=1; i<nIntervalsCount; i++)
		{
			
			if(settings.bShowCrossCorrelation)
			{
				crosscorrstack.addSlice(imCrCorr.calcFFTCorrelationImage(driftstack.getProcessor(i),driftstack.getProcessor(i+1),1));
			}
			xymaxd=imCrCorr.calcShiftFFTCorrelationDouble(driftstack.getProcessor(i),driftstack.getProcessor(i+1),1);			
			driftx[i] = driftx[i-1]+(xymaxd[0]*settings.nDriftScale);
			drifty[i] = drifty[i-1]+(xymaxd[1]*settings.nDriftScale);
	
			IJ.showProgress(i-1, nIntervalsCount);			
		}
		IJ.showProgress(nIntervalsCount, nIntervalsCount);
		fullTime = System.nanoTime() - startTime;
		IJ.log("Cross correlation calculation time: " + String.format("%.2f", ((double)Math.abs(fullTime-reconstructionTime))*0.000000001)+ " s");
		IJ.log("Total time: " + String.format("%.2f",((double)Math.abs(fullTime))*0.000000001) + " s");

		
		if(settings.bShowCrossCorrelation)
			new ImagePlus("Cross Correlation (Drift frames="+settings.nDriftFrames+" px size="+String.format("%d",((int)settings.nDriftScale))+" nm)", crosscorrstack).show();
		
		//applysimplecorrection();
		
		//linear interpolation of drift
		applylinearapproximationcorrection();		
		
		int lastDriftFrame=nIntervalsCount*nIntervalFrames;
		
		if(!settings.bFramesInterval)
		{
			//get results		
		
			for(i=0;i<lastDriftFrame;i++)
			{
				driftxfull[i]=approxx[i];
				driftyfull[i]=approxy[i];
			}
			
			//fill in the rest of correction with last value				
			if(nframes>lastDriftFrame)
			{
				for(i=lastDriftFrame;i<nframes;i++)
				{
					driftxfull[i]=approxx[lastDriftFrame-1];
					driftyfull[i]=approxy[lastDriftFrame-1];
				}
			}

		}
		else
		{
			//only corrected segment
			for(i=(int) settings.nFrameMin;i<=(int)(lastDriftFrame+settings.nFrameMin-1);i++)
			{
				driftxfull[i-1]=approxx[(int) (i-settings.nFrameMin)];
				driftyfull[i-1]=approxy[(int) (i-settings.nFrameMin)];				
			}
			//filling the rest with last value
			if(settings.nFrameMax>lastDriftFrame+settings.nFrameMin-1)
			{
				i=(int)(lastDriftFrame+settings.nFrameMin);
				while (i<=settings.nFrameMax)
				{
					driftxfull[i-1]=approxx[lastDriftFrame-1];
					driftyfull[i-1]=approxy[lastDriftFrame-1];				
					i++;
				}
				
			}

		}
		//show results of drift
		for(i=0;i<nframes;i++)
		{
			sDriftData = sDriftData + Integer.toString(i+1)+"\t"+Double.toString(driftxfull[i])+"\t"+Double.toString(driftyfull[i])+"\n";						
		}
		//Frame frame = WindowManager.getFrame("Drift Correction (frames="+settings.nDriftFrames+" max shift="+String.format("%d",((int)settings.nDriftMaxDistnm))+" nm)");
		Frame frame = WindowManager.getFrame("Drift Correction (frames bin="+settings.nDriftFrames +" px size="+String.format("%d",((int)settings.nDriftScale))+" nm)");
		if (frame!=null && (frame instanceof TextWindow) )
		{
			DriftTable = (TextWindow)frame;
			DriftTable.getTextPanel().clear();
			DriftTable.getTextPanel().append(sDriftData);
			DriftTable.getTextPanel().updateDisplay();			
		}
			else
				DriftTable = new TextWindow("Drift Correction (frames="+settings.nDriftFrames+" px size="+String.format("%d",((int)settings.nDriftScale))+" nm)", "Frame_Number\tX_drift_(nm)\tY_drift_(nm)", sDriftData, 450, 300);
				//DriftTable = new TextWindow("Drift Correction (frames="+settings.nDriftFrames+" max shift="+String.format("%d",((int)settings.nDriftMaxDistnm))+" nm)", "Frame_Number\tX_drift_(nm)\tY_drift_(nm)", sDriftData, 450, 300);			
		
		return;
		
	}
	
	void DriftUpdateResultTable(SMLAnalysis sml_, SMLDialog dlg_)
	{
		int i;
		//let's restore scale
		double pxscale =  sml_.ptable.getValueAsDouble(DOMConstants.Col_X, 0)/sml_.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0);
		
		//lock table
		sml_.ptable_lock.lock();
		for(i=0;i<nParticlesCount;i++)
		{
			sml_.ptable.setValue(DOMConstants.Col_Xnm, i, x[i]);
			sml_.ptable.setValue(DOMConstants.Col_Ynm, i, y[i]);
			sml_.ptable.setValue(DOMConstants.Col_X, i, x[i]*pxscale);
			sml_.ptable.setValue(DOMConstants.Col_Y, i, y[i]*pxscale);
		}
		sml_.ptable_lock.unlock();
		//sml_.ptable.updateResults();
		sml_.showTable();
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
	/** function uses linear approximation for interpolation of intermediate values 
	 * for discrete set of drift correction values (corresponding to bins of frames)
	 * */
	void applylinearapproximationcorrection()
	{
		int i;
		int dIndex;
		int nIntervalFrames = settings.nDriftFrames;
		double coeffx, coeffy;
		boolean bStop;
		
		int lastDriftFrame=nIntervalsCount*nIntervalFrames;
		
		approxx = new double [lastDriftFrame];
		approxy = new double [lastDriftFrame];
		
		
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
		//last segment
		for(;i<lastDriftFrame; i++)
		{
			approxx[i] = driftx[dIndex];
			approxy[i] = drifty[dIndex];		
		}
		
		//int lastDriftFrame=nIntervalsCount*nIntervalFrames;
		//for(i=0;i<lastDriftFrame;i++)
		
		
		if(settings.bFramesInterval)
		{
			i=0;
			//skip beginning
			while(f[i]<settings.nFrameMin)
				i++;
			bStop = false;
			while (!bStop)
			{
				dIndex = (int) f[i]-(int)settings.nFrameMin;
				x[i]+=approxx[dIndex];
				y[i]+=approxy[dIndex];
				
				i++;
				if(i>=nParticlesCount)
					bStop=true;
				else
				{
					//or the end of drift correction "clean" interval
					if(f[i]>=(settings.nFrameMin+nIntervalsCount*nIntervalFrames))
						bStop=true;
				}
			}
			//filling the rest
			if(settings.nFrameMax>settings.nFrameMin+nIntervalsCount*nIntervalFrames-1)
			{
				dIndex = nIntervalsCount-1;
				bStop = false;
				while (!bStop)
				{
					x[i]+=driftx[dIndex];
					y[i]+=drifty[dIndex];
					i++;
					if(i>=nParticlesCount)
						bStop=true;
					else
					{
						//or the end of drift correction "clean" interval
						if(f[i]>(settings.nFrameMax))
							bStop=true;
					}

				}
			}
		}
		else
		{
			
			bStop = false; 
			i=0;
			while (!bStop)
			{
				dIndex = (int) (f[i]-1);
				x[i]+=approxx[dIndex];
				y[i]+=approxy[dIndex];
				i++;
				//let's check whether we reached the end of all particles list
				if(i>=nParticlesCount)
					bStop=true;
				else
				{
					//or the end of drift correction "clean" interval
					if(f[i]>=(lastDriftFrame+1))
						bStop=true;
				}
				
			}

			//fill in the rest of correction with last value				
			if(nframes>lastDriftFrame)
			{
				dIndex = nIntervalsCount-1;
				for(;i<nParticlesCount;i++)
				{
					x[i]+=driftx[dIndex];
					y[i]+=drifty[dIndex];
				}
			}
			
		}
		return;
			
		
	}
	
	/* Function calculates cross-correlation between two images //OBSOLETE	
	ShortProcessor crosscorrelation (ImageProcessor ip1, ImageProcessor ip2)
	{
		int nMaxPix = settings.nDriftPixels;
		int i, j, m, n;
		int tot;
		double dCC,dCC1,dCC2;
		double val1,val2;
		
		
		FloatProcessor resultIP;
		ShortProcessor returnIP;
		ImageProcessor extendedip1;
		
		
		tot = 2*nMaxPix+1;
		extendedip1 = new FloatProcessor(drift_width+tot-1, drift_height+tot-1);

		for (i=0;i<drift_width;i++)
			for (j=0;j<drift_height;j++)
				extendedip1.setf(i+nMaxPix, j+nMaxPix, ip1.get(i,j));

		resultIP = new FloatProcessor(tot, tot);
		
		
		for (i=0;i<tot;i++)
			for (j=0;j<tot;j++)
			{
				dCC=0;
				dCC1=0;
				dCC2=0;
				for(m=0; m<drift_width; m++)
					for(n=0; n<drift_height; n++)
					{
						val1=extendedip1.get(m+i,n+j);
						val2=ip2.get(m, n);
						dCC+=val1*val2;
						dCC1+=val1;
						dCC2+=val2;
					}
				if(dCC1!=0.0 && dCC2!=0.0)
					resultIP.setf(i,j,(float)(dCC/(dCC1*dCC2)));
				else
					resultIP.setf(i,j,0);
			}
		
		returnIP = (ShortProcessor) resultIP.convertToShort(true);
		returnIP.smooth();
		return returnIP;		
	}*/
	
	/** *
	 * builds cumulative distribution histogram 
	 * of z coordinate with 256 bins 
	 */
	int [] getZcumHistogram()
	{
		double [] cumHist = new double[256];
		int [] cumHistFinal = new int[256];
		//double zmin, zmax;
		double ztemp;
		int i;	
	
		//zmin=Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		//zmax=Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		
		if(settings.bDynamicZscale)
		{
			//build histogram
			for (i =0; i<nParticlesCount; i++)
			{
				if(z[i]>zmin && z[i]<zmax)
				{
					ztemp = (z[i]-zmin)/(zmax-zmin);
					cumHist[(int)Math.round(ztemp*255)]++;
				}
			}
			
			//build cumulative hist
			for(i=1;i<256;i++)
				cumHist[i]=(cumHist[i]+cumHist[i-1]);
			//build cumulative normalize
			for(i=1;i<256;i++)
				cumHistFinal[i]=(int)Math.round(255*cumHist[i]/cumHist[255]);
		}
		//simple linear stuff
		else
		{
			for(i=0;i<256;i++)
				cumHistFinal[i]=i;
			
		}
		return cumHistFinal;
		
	}
	
	/** function gets LUT specified by sZLUTName in settings
	 * and returns 256x3 table map in HSB format */
	float [][]  getHSBLutTable()
	//void  getHSBLutTable()
	{
		int i,j;
		//ColorProcessor ipLUT;
		float[] hsbvals = new float[4];
		//int [][] rgbtable = new int[256][4]; 
		int [] onepix; 
		float [][] HSBLutTable = new float[256][3];
		ByteProcessor ish = new ByteProcessor(256,10);
		for ( i=0; i<256; i++)
			for (j=0; j<10; j++)
				ish.putPixel(i, j, i);
		ImagePlus ccc = new ImagePlus("test",ish);
		ccc.show();
		IJ.run(settings.sZLUTName);
		IJ.run("RGB Color");
		//ipLUT= (ColorProcessor) ccc.getProcessor();
		ccc.setSlice(1);
		for(i=0;i<256;i++)
		{
			
			onepix= ccc.getPixel(i, 2);
			//rgbtable[i]=ccc.getPixel(i, 1);
			java.awt.Color.RGBtoHSB(onepix[0], onepix[1], onepix[2], hsbvals);
			HSBLutTable[i][0]=hsbvals[0];
			HSBLutTable[i][1]=hsbvals[1];
			HSBLutTable[i][2]=hsbvals[2];
		}
		/*
		int [] newrgb= new int[3];
		Color c;
		ColorProcessor imcolcode = new ColorProcessor(256, 40);
		for(i=0;i<256;i++)
			for(j=0;j<25;j++)
			{
				
					c = Color.getHSBColor(HSBLutTable[i][0], HSBLutTable[i][1], HSBLutTable[i][2]);
				
				
				newrgb[0]=c.getRed();
				newrgb[1]=c.getGreen();
				newrgb[2]=c.getBlue();
				imcolcode.putPixel(i, j, newrgb);
			}
		new ImagePlus("testafter",imcolcode).show();
		
		*/
		ccc.changes=false;
		ccc.close();
		return HSBLutTable;
		//return;
	}
	
	
	/** Cleans the reconstruction viewer image. */
	void clear()
	{
	
		for (int i=0;i<ip.getWidth();i++)
			for (int j=0;j<ip.getHeight();j++)
				ip.set(i, j, 0);
	
	}
	
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
	void averagelocalizations(SMLAnalysis sml_)
	{
		int nUniqFrames,i,j;
		int nPatNumber = f.length;
		int nCount;
		int nCurrFrame, nFrameCount, nCurrentParticle, nCompareFrame, nCompareParticle,nFrameCompareCount;
		int nAveragesNumber;
		int [][] framesstat;
		double [] temp;
		double [][][] res_table_unique;
		double [][] res_table;	
		
		boolean bContinue;
		double dDistance;
		/** array containing accumulating values of averaged numbers*/
		ArrayList<double[]> Averaging_tbl = new ArrayList<double[]>();
		
		IJ.showStatus("Averaging single photoactivation events...");
		
		//sort table by frame number
		//should be already sorted in constructor

		sml_.ptable_lock.lock();

		res_table = new double [nPatNumber][DOMConstants.Col_Total_No_Tracks];
		
		//read results table to double array with transpose operation (makes later code easier)
		for (i=0;i<DOMConstants.Col_Total_No_Tracks;i++)
		{
			temp = sml_.ptable.getColumnAsDoubles(i);
			for(j=0;j<nPatNumber;j++)
				res_table[j][i]=temp[j];
		}			
		sml_.ptable_lock.unlock();			

		
		//calculate number of particles in each frame and
		//store it for each frame in framesstat array
		nUniqFrames = uniqFrames();
		framesstat = new int [nUniqFrames][2];
		nCurrFrame = (int) f[0];
		nCount = 0;
		nFrameCount = 0;
		framesstat[nFrameCount][0] = nCurrFrame;
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
		
		res_table_unique = new double [nUniqFrames][][];
		
		for(i=0;i<nUniqFrames;i++)
		{
			nCount = framesstat[i][1];					
			res_table_unique[i] = new double [nCount][DOMConstants.Col_Total_No_Tracks];
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
			res_table_unique[nFrameCount][nCurrentParticle]= res_table[nCount];
									
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
			
			if(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]<dFPThreshold && nFrameCount<(nUniqFrames-1))
			//if(fp_un[nFrameCount][nCurrentParticle]<dFPThreshold && nFrameCount<(nUniqFrames-1)) 
			{ 	
				nCurrFrame = framesstat[nFrameCount][0];
				//accumulating information				
				Averaging_tbl.clear();
				
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
					//in the beginning look the next frame should be just after first
					if(nAveragesNumber==0 && (nCompareFrame-nCurrFrame>1.5))
						{bContinue = false;}
					//ok, let's see if all checks are passed
					
					if(bContinue && res_table_unique[nFrameCompareCount][nCompareParticle][DOMConstants.Col_Fp]<dFPThreshold)
					{//seems so. let's measure distances between molecules then
						dDistance = Math.sqrt(Math.pow(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_X] -
													   res_table_unique[nFrameCompareCount][nCompareParticle][DOMConstants.Col_X], 2) 
								             +Math.pow(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Y]-
								            		   res_table_unique[nFrameCompareCount][nCompareParticle][DOMConstants.Col_Y], 2)); 

						//they are closer then one pixel
						if(dDistance<1.0)
						{
							//found same photoactivation event
							//let's store values
							//mark it as already used (fp)
							res_table_unique[nFrameCompareCount][nCompareParticle][DOMConstants.Col_Fp]+=1.2;//mark it as an averaged one (>1), so skip in update
							Averaging_tbl.add(res_table_unique[nFrameCompareCount][nCompareParticle]);
							nAveragesNumber++;
							//exit from current frame (a bit lame one)
							nCompareParticle = framesstat[nFrameCompareCount][1] -1;
						}
					}
				}
				//was photoactivation event found?
				if(nAveragesNumber>0)
				{
					//adding particle itself to the list of averaged values
					//mark it as already used (fp)
					res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]+=1.2;//mark it as an averaged one (>1), so skip in update
					Averaging_tbl.add(res_table_unique[nFrameCount][nCurrentParticle]);
					nAveragesNumber++;
					
					//calculating weighted average for this particle
					res_table_unique[nFrameCount][nCurrentParticle] = weightedAverage(Averaging_tbl);
					//function will mark it as averaged inside (do subtraction -1.2);
				}
			}//false positive check
		}
		
		//updating initial arrays used for rendering picture		
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
			//this one was used in averaging
			if(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]>1)
			{
				fp[nCount]+=1.2;
			}
			//this one was is the cumulative averaged particle
			if(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]<0)
			{				
				x[nCount]   = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Xnm];
				y[nCount]   = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Ynm];
				z[nCount]   = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Znm];
				loc_errx[nCount] = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_loc_errX];
				loc_erry[nCount] = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_loc_errY];
				res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]+=1.2; //restore averaged value
				fp[nCount]  = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp];
				f[nCount] = res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_FrameN];
			}
		}
		
		///update results table
		if(settings.bUpdateAveragePositions)
		{
			sml_.ptable_lock.lock();
			//clear particle table
			sml_.ptable.reset();
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
				if(res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]<=1)
				{
					sml_.ptable.incrementCounter();
					
					sml_.ptable.addValue("X_(px)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_X]);							
					sml_.ptable.addValue("Y_(px)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Y]);
					sml_.ptable.addValue("Frame Number", res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_FrameN]);
					sml_.ptable.addValue("X_(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Xnm]);							
					sml_.ptable.addValue("X_loc_error(nm)", res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_loc_errX]);							
					sml_.ptable.addValue("Y_(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Ynm]);
					sml_.ptable.addValue("Y_loc_error(nm)", res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_loc_errY]);
					sml_.ptable.addValue("Z_(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Znm]);
					sml_.ptable.addValue("Z_loc_error(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_loc_errZ]);							
					sml_.ptable.addValue("Amplitude_fit",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_AmplFit]);
					sml_.ptable.addValue("Amp_error",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Amp_error]);
					sml_.ptable.addValue("BGfit",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_BGfit]);
					sml_.ptable.addValue("BGfit_error",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_BGfit_error]);
					sml_.ptable.addValue("SD_X_(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_SD_X]);
					sml_.ptable.addValue("SD_X_error(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_SD_X_err]);
					sml_.ptable.addValue("SD_Y_(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_SD_Y]);
					sml_.ptable.addValue("SD_Y_error(nm)",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_SD_Y_err]);
					sml_.ptable.addValue("False positive", res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_Fp]);							
					sml_.ptable.addValue("IntegratedInt",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_IntegrInt]);
					sml_.ptable.addValue("SNR", res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_SNR]);	
					//ptable.addValue("chi2_fit",chi2_fit);
					sml_.ptable.addValue("R2_fit",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_chi]);												
					//ptable.addValue("Iterations_fit",fdg.nIterations);
					sml_.ptable.addValue("Iterations_fit",res_table_unique[nFrameCount][nCurrentParticle][DOMConstants.Col_IterN]);
				}

			}
			
			sml_.ptable_lock.unlock();
			sml_.showTable();
		}
					
	}
	
	/** Function calculating weighted average of results
	 * */
	double [] weightedAverage(ArrayList<double[]> Averaging_tbl)
	{
		double [] dAverW;
		double dWeight;
		int nCount=0;
		dAverW = new double [DOMConstants.Col_Total_No_Tracks];
		
		//calculating averaged values
		for(double [] item:Averaging_tbl)
		{
			dWeight = Math.pow(item[DOMConstants.Col_loc_errX], -2);
			dAverW[DOMConstants.Col_X] += item[DOMConstants.Col_X]*dWeight; 				
			dAverW[DOMConstants.Col_Xnm] += item[DOMConstants.Col_Xnm]*dWeight;			
			dAverW[DOMConstants.Col_loc_errX] += dWeight;
			
			dWeight = Math.pow(item[DOMConstants.Col_loc_errY], -2);
			dAverW[DOMConstants.Col_Y] += item[DOMConstants.Col_Y]*dWeight;
			dAverW[DOMConstants.Col_Ynm] += item[DOMConstants.Col_Ynm]*dWeight;			
			dAverW[DOMConstants.Col_loc_errY] += dWeight;
			
			dWeight = Math.pow(item[DOMConstants.Col_Amp_error], -2);
			dAverW[DOMConstants.Col_AmplFit] += item[DOMConstants.Col_AmplFit]*dWeight;
			dAverW[DOMConstants.Col_Amp_error] += dWeight;
			
			dWeight = Math.pow(item[DOMConstants.Col_BGfit_error], -2);
			dAverW[DOMConstants.Col_BGfit] += item[DOMConstants.Col_BGfit]*dWeight;
			dAverW[DOMConstants.Col_BGfit_error] += dWeight;
			
			dWeight = Math.pow(item[DOMConstants.Col_SD_X_err], -2);
			dAverW[DOMConstants.Col_SD_X] += item[DOMConstants.Col_SD_X]*dWeight;
			dAverW[DOMConstants.Col_SD_X_err] += dWeight;
			
			dWeight = Math.pow(item[DOMConstants.Col_SD_Y_err], -2);
			dAverW[DOMConstants.Col_SD_Y] += item[DOMConstants.Col_SD_Y]*dWeight;
			dAverW[DOMConstants.Col_SD_Y_err] += dWeight;
			
			//TODO recalculate Z as weighted when Z error would be present
			/*if(Math.abs(item[DOMConstants.Col_Znm])>0)
			{
				dWeight = Math.pow(item[DOMConstants.Col_loc_errZ], -2);
				dAverW[DOMConstants.Col_Znm] += item[DOMConstants.Col_Znm]*dWeight;
				dAverW[DOMConstants.Col_loc_errZ] += dWeight;
			}
			*/
			//simple average part
			dAverW[DOMConstants.Col_Znm] += item[DOMConstants.Col_Znm];
			
			dAverW[DOMConstants.Col_IntegrInt] += item[DOMConstants.Col_IntegrInt];
			dAverW[DOMConstants.Col_SNR] += item[DOMConstants.Col_SNR];
			dAverW[DOMConstants.Col_chi] += item[DOMConstants.Col_chi];
			dAverW[DOMConstants.Col_IterN] += item[DOMConstants.Col_IterN];
			dAverW[DOMConstants.Col_FrameN] += item[DOMConstants.Col_FrameN];
			dAverW[DOMConstants.Col_Fp] += item[DOMConstants.Col_Fp]-1.2;//to restore real value
			nCount++;
		}

		//weighted average sum divide by sum of weights 
		dAverW[DOMConstants.Col_X]/=dAverW[DOMConstants.Col_loc_errX];
		dAverW[DOMConstants.Col_Xnm]/=dAverW[DOMConstants.Col_loc_errX];
		//for the SD (SEM) itself it is this formula
		dAverW[DOMConstants.Col_loc_errX] = 1/Math.sqrt(dAverW[DOMConstants.Col_loc_errX]);
		
		dAverW[DOMConstants.Col_Y]/=dAverW[DOMConstants.Col_loc_errY];
		dAverW[DOMConstants.Col_Ynm]/=dAverW[DOMConstants.Col_loc_errY];
		dAverW[DOMConstants.Col_loc_errY] = 1/Math.sqrt(dAverW[DOMConstants.Col_loc_errY]);
		
		dAverW[DOMConstants.Col_AmplFit]/=dAverW[DOMConstants.Col_Amp_error];		
		dAverW[DOMConstants.Col_Amp_error] = 1/Math.sqrt(dAverW[DOMConstants.Col_Amp_error]);

		dAverW[DOMConstants.Col_BGfit]/=dAverW[DOMConstants.Col_BGfit_error];		
		dAverW[DOMConstants.Col_BGfit_error] = 1/Math.sqrt(dAverW[DOMConstants.Col_BGfit_error]);
		
		dAverW[DOMConstants.Col_SD_X]/=dAverW[DOMConstants.Col_SD_X_err];		
		dAverW[DOMConstants.Col_SD_X_err] = 1/Math.sqrt(dAverW[DOMConstants.Col_SD_X_err]);
		
		dAverW[DOMConstants.Col_SD_Y]/=dAverW[DOMConstants.Col_SD_Y_err];		
		dAverW[DOMConstants.Col_SD_Y_err] = 1/Math.sqrt(dAverW[DOMConstants.Col_SD_Y_err]);
		
		//TODO recalculate Z as weighted when Z error would be present
		/*if(Math.abs(dAverW[DOMConstants.Col_Znm])>0)
		{
			dAverW[DOMConstants.Col_Znm]/=dAverW[DOMConstants.Col_loc_errZ];		
			dAverW[DOMConstants.Col_loc_errZ] = 1/Math.sqrt(dAverW[DOMConstants.Col_loc_errZ]);			
		}*/

		//just simple average for these values
		dAverW[DOMConstants.Col_Znm] /= (double)nCount;
		dAverW[DOMConstants.Col_IntegrInt] /= (double)nCount;
		dAverW[DOMConstants.Col_SNR] /= (double)nCount;
		dAverW[DOMConstants.Col_chi] /= (double)nCount;
		dAverW[DOMConstants.Col_IterN] /= (double)nCount;
		
		dAverW[DOMConstants.Col_Fp] /= (double)nCount;
		dAverW[DOMConstants.Col_Fp]= Math.round(dAverW[DOMConstants.Col_Fp]);
		dAverW[DOMConstants.Col_Fp]= dAverW[DOMConstants.Col_Fp] -1.2; //mark it as averaged
		//dAverW[DOMConstants.Col_Fp] = Averaging_tbl.get(nCount-1)[DOMConstants.Col_Fp]-1.2;
		dAverW[DOMConstants.Col_FrameN]/= (double)nCount;
		dAverW[DOMConstants.Col_FrameN] = Math.round(dAverW[DOMConstants.Col_FrameN]);
		//dAverW[DOMConstants.Col_FrameN] = Averaging_tbl.get(nCount-1)[DOMConstants.Col_FrameN];

		return dAverW;
	}
	
	/** Function calculating number of unique frames */
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
	/** function updating arrays of x and y with translation value */
	void addTranslation()
	{
		nParticlesCount = x.length;
		for (int i=0; i<nParticlesCount; i++)
		{
			x[i] += settings.dTranslationX;
			y[i] += settings.dTranslationY;
		}	
		IJ.log("Coordinates translated " + Double.toString(settings.dTranslationX) +" nm in X and "+ Double.toString(settings.dTranslationY) +" in Y");

	}
	/**function calculating zmin and zmax
	 * */
	void getZmaxmin()
	{
		zmax=(-1)*Double.MAX_VALUE;
		zmin=Double.MAX_VALUE;
		for (int i=0;i<z.length;i++)
		{
			if(z[i]>zmax)
				zmax=z[i];
			if(z[i]<zmin)
				zmin=z[i];
			
		}
	}

}