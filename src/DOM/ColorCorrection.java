package DOM;

import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.StackProcessor;


public class ColorCorrection implements PlugIn{
	
	ImagePlus imp;
	ImagePlus imp_refMax;
	ImagePlus imp_warpMax;
	
	ImageStack is_temp;
	//ImagePlus imp_tempFFT;
	//FloatProcessor [] is_refFFT;
	//FloatProcessor [] is_warpFFT;
	//FloatProcessor [] is_FFTMult;
	Roi refROI, warpedROI;
	SMLAnalysis sml;
	SMLDialog dlg = new SMLDialog();
	
	double [] xref, yref, SNRref, fref;
	double [] xwarp, ywarp, SNRwarp, fwarp;
	double [] dXYShift;
	/** list of colocalized particles coordinates 
	 * */
	double [] xy_col;
	/** Max SNR threshold*/
	double dCCSNR;
	/** Max distance between particles in pixels */
	double dCCDist;
	
	//final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Fit Z-calibration"); 
	@Override
	public void run(String arg) {
	
		int i;
		
		IJ.register(ColorCorrection.class);
	
		//show parameters dialog
		if(!dlg.dglColorCalibration())
			return;
		dCCSNR= dlg.dCCSNR;
		dCCDist = dlg.dCCDist;
		
		NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Choose image and ROI of REFERENCE channel");
		
		nb.addMessage("Choose calibration image and ROI region for REFERENCE channel");
		nb.showDialog();
		if(!nb.wasOKed())
			return;
		//let's start detection
		
		IJ.log("Creating color chromatic calibration from image stack");
		IJ.log("Step 1: Detecting molecules in reference channel");
		IJ.run("Detect Molecules");
		
		imp = IJ.getImage();
		
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 ) 
		{
		    IJ.error("8 or 16 bit greyscale image required");
		    return;
		}
		
		
		
		refROI=imp.getRoi();
		if(!(refROI==null))
		{
			imp=imp.duplicate();
		}
		ZProjector zProj = new ZProjector(imp);
		zProj.setMethod(ZProjector.MAX_METHOD);
	
		zProj.doProjection();
		imp_refMax = zProj.getProjection();
		
		//imp_refMax.show();
		
		sml = new SMLAnalysis();
		//get detection results
		xref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		yref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		SNRref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		fref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		
		nb = new NonBlockingGenericDialog("Choose image and ROI of WARPED channel");
		nb.addMessage("Choose calibration image and ROI region for WARPED channel");
		nb.showDialog();
		if(!nb.wasOKed())
			return;

		
		imp = IJ.getImage();
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 ) 
		{
		    IJ.error("8 or 16 bit greyscale image required");
		    return;
		}
		//let's start detection	
		IJ.log("Step 2: Detecting molecules in warped channel");
		IJ.run("Detect Molecules");
		
		warpedROI=imp.getRoi();
		if(!(warpedROI==null))
		{
			imp=imp.duplicate();
		}
		zProj = new ZProjector(imp);
		zProj.doProjection();
		imp_warpMax = zProj.getProjection();
		
		//imp_warpMax.show();
		sml = new SMLAnalysis();
		//get detection results
		xwarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		ywarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		SNRwarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		fwarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		

		IJ.log("Aligning maximum projection pictures...");
		
		//TODO make sizes of images equal
		dXYShift = calcShiftFFTCorrelation(imp_refMax.getProcessor(),imp_warpMax.getProcessor());
		
		//update warped particles coordinates
		for(i=0;i<xwarp.length;i++)
		{
			xwarp[i]+=dXYShift[0];
			ywarp[i]+=dXYShift[1];
		}
		
		IJ.showProgress(1.0);
		
		IJ.log("X shift: "+Double.toString(dXYShift[0])+ " px");
		IJ.log("Y shift: "+Double.toString(dXYShift[1]) + " px");
		IJ.log("Finding corresponding images of particles in both channels...");
		
		IJ.log("In total "+Integer.toString(CCfindParticles())+" particles found.");
		
		IJ.log("Calculating B-spline transform...");
	}
	
	
	//TODO add subpixel precision
	/** Function calculates shift in X and Y between images using
	 * cross correlation calculated through FFT transform
	 * */
	public double [] calcShiftFFTCorrelation(ImageProcessor ip1, ImageProcessor ip2)
	{
		double [] xyshift = new double [2];
		
		
		int i,j;
		float dRefRe,dRefIm, dWarpRe, dWarpIm;

		int imheight, imwidth;
		int imheightFFT, imwidthFFT;
	    FHT fht1,fht2;
	    
		imheight=ip1.getHeight();
		imwidth =ip1.getWidth();
		
	    fht1 = new FHT(pad(ip1));
	    fht2 = new FHT(pad(ip2));
	    fht1.originalHeight=imheight;
	    fht1.originalWidth=imwidth;
	    
	    fht2.originalHeight=imheight;
	    fht2.originalWidth=imwidth;
	    
	    fht1.transform();
	    fht2.transform();
	    ImageStack ct1 = fht1.getComplexTransform();
	    ImageStack ct2 = fht2.getComplexTransform();
	    ImageStack istemp;
		
	    FloatProcessor [] is_refFFT = new FloatProcessor[2];
		FloatProcessor [] is_warpFFT = new FloatProcessor[2];
		FloatProcessor [] is_FFTMult = new FloatProcessor[2];
		imheightFFT=ct1.getHeight();
		imwidthFFT =ct1.getWidth();
		for(i=0;i<2;i++)
		{
			is_refFFT[i]=(FloatProcessor) ct1.getProcessor(i+1);
		}
		for(i=0;i<2;i++)
		{
			is_warpFFT[i]=(FloatProcessor) ct2.getProcessor(i+1);
			is_FFTMult[i] = new FloatProcessor(imwidthFFT,imheightFFT);
		}
			
		//calculate convolution between two FFT
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dRefRe=is_refFFT[0].getf(i,j);
				dRefIm=is_refFFT[1].getf(i,j);
				
				dWarpRe=is_warpFFT[0].getf(i,j);
				//conjugate
				dWarpIm=(float) ((-1.0)*is_warpFFT[1].getf(i,j));
		
				
				is_FFTMult[0].setf(i,j,dRefRe*dWarpRe-dWarpIm*dRefIm);
				is_FFTMult[1].setf(i,j,dRefRe*dWarpIm+dWarpRe*dRefIm);
				
				
			}

		istemp = new ImageStack(imwidthFFT, imheightFFT);
		


		istemp.addSlice("Real", is_FFTMult[0]);
		istemp.addSlice("Imaginary", is_FFTMult[1]);
		
		
		//ct1 = doComplexInverseTransform(istemp,imwidth, imheight);
		ct1 = doComplexInverseTransform(istemp,0, 0);
		
		//ImagePlus imp1 = new ImagePlus("Complex of 1", ct1);
		//imp1.show();
		 
	    
		for(i=0;i<2;i++)
		{
			is_FFTMult[i]=(FloatProcessor) ct1.getProcessor(i+1);
		}
		
		int niMax=0,njMax=0;
		double nAbsMax=(-1.0)*Double.MAX_VALUE;
		double dVal;
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dVal=Math.sqrt(is_FFTMult[0].getf(i,j)*is_FFTMult[0].getf(i,j)+is_FFTMult[1].getf(i,j)*is_FFTMult[1].getf(i,j));
				
				if(dVal>0)
				{
					dVal++;
				}
				if(dVal>nAbsMax)
				{
					nAbsMax=dVal;
					niMax=i;
					njMax=j;
				}
			}
		
		if((niMax+1)>Math.floor(0.5*imwidthFFT))
		{
			niMax = niMax - imwidthFFT;
		}
		/*else
		{
			niMax = niMax;
		}*/
		if((njMax+1)>Math.floor(0.5*imheightFFT))
		{
			njMax = njMax - imheightFFT;
		}
		/*
		else
		{
			njMax = njMax;
		}*/
		
		xyshift[0]=(double)niMax;
		xyshift[1]=(double)njMax;
	    return xyshift;
		
	}
	
	
	/** function locating particles images in both channels
	 *  returns number of colocalized particles  
	 * */
	int CCfindParticles()
	{
		ArrayList<double[]> colocTable = new ArrayList<double[]>();
		double [] xyfound = new double [4];
		double dMinDist, dDist;
		int ind_warped;
		int i,j;
		int nRef, nWarp;
		nRef = xref.length;
		nWarp = xwarp.length;
		//looking for colocalization
		for (i=0;i<nRef;i++)
		{
			if(SNRref[i]>dCCSNR)
			{
				dMinDist = Double.MAX_VALUE;
				ind_warped = -1;
				for (j=0;j<nWarp;j++)
				{
					if(SNRwarp[j]>dCCSNR)
					{
						//not marked as used yet
						if(fwarp[j]<3)
						{
							//on the same frame
							if((int)fref[i]== (int)fwarp[j])
							{
								dDist = Math.sqrt(Math.pow(xref[i]-xwarp[j], 2)+Math.pow(yref[i]-ywarp[j], 2));
								if (dDist>dCCDist)
								{
									if(dDist<dMinDist)
									{
										dMinDist=dDist;
										ind_warped = j;
										
									}
								}
							}
						}
					}
				}
				//let's see if found anything
				if(ind_warped>=0)
				{

					//mark particle as used
					fwarp[ind_warped]=5;
					xyfound[0]=xref[i];
					xyfound[1]=yref[i];
					xyfound[2]=xwarp[ind_warped];
					xyfound[3]=ywarp[ind_warped];
					colocTable.add(xyfound);
				}
			}
		}
		return colocTable.size();
	}

	/** 
	 * some parts from ij.plugin.FFT code
	 *  
	 * */
    ImageProcessor pad(ImageProcessor ip) 
    {
        int originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) 
        {
        
            return ip;
        }
        maxN = i;
   
        ImageStatistics stats = ImageStatistics.getStatistics(ip, ImageStatistics.MEAN, null);
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(stats.mean);
        ip2.fill();
        ip2.insert(ip, 0, 0);
  
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
	/**	Swap quadrants 1 and 3 and 2 and 4 of the specified ImageProcessor 
	so the power spectrum origin is at the center of the image.
	<pre>
	    2 1
	    3 4
	</pre>
*/
    void swapQuadrants(ImageStack stack) {
        FHT fht = new FHT(new FloatProcessor(1, 1));
        for (int i=1; i<=stack.getSize(); i++)
            fht.swapQuadrants(stack.getProcessor(i));
    }
    
    /** Complex to Complex Inverse Fourier Transform
    *   @author Joachim Wesner
    */
    void c2c2DFFT(float[] rein, float[] imin, int maxN, float[] reout, float[] imout) 
    {
            FHT fht = new FHT(new FloatProcessor(maxN,maxN));
            float[] fhtpixels = (float[])fht.getPixels();
            // Real part of inverse transform
            for (int iy = 0; iy < maxN; iy++)
                  cplxFHT(iy, maxN, rein, imin, false, fhtpixels);
            fht.inverseTransform();
            // Save intermediate result, so we can do a "in-place" transform
            float[] hlp = new float[maxN*maxN];
            System.arraycopy(fhtpixels, 0, hlp, 0, maxN*maxN);
            // Imaginary part of inverse transform
            for (int iy = 0; iy < maxN; iy++)
                  cplxFHT(iy, maxN, rein, imin, true, fhtpixels);
            fht.inverseTransform();
            System.arraycopy(hlp, 0, reout, 0, maxN*maxN);
            System.arraycopy(fhtpixels, 0, imout, 0, maxN*maxN);
      }
    
    /** Build FHT input for equivalent inverse FFT
    *   @author Joachim Wesner
    */
    void cplxFHT(int row, int maxN, float[] re, float[] im, boolean reim, float[] fht) 
    {
            int base = row*maxN;
            int offs = ((maxN-row)%maxN) * maxN;
            if (!reim) {
                  for (int c=0; c<maxN; c++) {
                        int l =  offs + (maxN-c)%maxN;
                        fht[base+c] = ((re[base+c]+re[l]) - (im[base+c]-im[l]))*0.5f;
                  }
            } else {
                  for (int c=0; c<maxN; c++) {
                        int l = offs + (maxN-c)%maxN;
                        fht[base+c] = ((im[base+c]+im[l]) + (re[base+c]-re[l]))*0.5f;
                  }
            }
      }
    
    ImageStack doComplexInverseTransform(ImageStack stack, int width, int height ) 
    {
   
        int maxN = stack.getWidth();
        swapQuadrants(stack);
        float[] rein = (float[])stack.getPixels(1);
        float[] imin = (float[])stack.getPixels(2);
        float[] reout= new float[maxN*maxN];
        float[] imout = new float[maxN*maxN];
        c2c2DFFT(rein, imin, maxN, reout, imout);
        ImageStack stack2 = new ImageStack(maxN, maxN);
        swapQuadrants(stack);
        stack2.addSlice("Real", reout);
        stack2.addSlice("Imaginary", imout);
        stack2 = unpad(stack2,width,height);
        /*String name = WindowManager.getUniqueName(imp.getTitle().substring(10));
        ImagePlus imp2 = new ImagePlus(name, stack2);
        imp2.getProcessor().resetMinAndMax();
        imp2.show();*/
        return stack2;
    }
    
    ImageStack unpad(ImageStack stack, int width, int height) 
    {
      
 
        if (width==0 || height==0 || (width==stack.getWidth()&&height==stack.getHeight()))
            return stack;
        StackProcessor sp = new StackProcessor(stack, null);
        ImageStack stack2 = sp.crop(0, 0, width, height);
        return stack2;
    }
    
}
