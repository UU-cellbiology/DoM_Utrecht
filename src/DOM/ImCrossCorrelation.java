package DOM;


import java.awt.Rectangle;

import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.StackProcessor;

public class ImCrossCorrelation {

	
	public ImageProcessor calcFFTCorrelationImage(ImageProcessor ip1, ImageProcessor ip2, int nNorm)
	{
	
		
		int i,j;
		float dRefRe,dRefIm, dWarpRe, dWarpIm, dVal, fNorm;
		ImageProcessor paddedip1, paddedip2, normIP; 
		//ip1 and ip2 assumed to be same size
		int origW =ip1.getWidth();
		int origH =ip1.getHeight();
		int nPaddedS;
		
		int imheightFFT, imwidthFFT;
	    FHT fht1,fht2;
	   
	    paddedip1=padzeros(ip1);
	    paddedip2=padzeros(ip2);
	    //since now width and height is the same, just one measurement
	    nPaddedS=paddedip1.getWidth();
	    
	    fht1 = new FHT(paddedip1);
	    fht2 = new FHT(paddedip2);
	    fht1.setShowProgress(false);
	    fht2.setShowProgress(false);
	    
	    fht1.transform();
	    fht2.transform();
	    ImageStack ct1 = fht1.getComplexTransform();
	    ImageStack ct2 = fht2.getComplexTransform();
	    ImageStack istemp;
		
	    FloatProcessor [] is_refFFT = new FloatProcessor[2];
		FloatProcessor [] is_warpFFT = new FloatProcessor[2];
		FloatProcessor [] is_FFTMult = new FloatProcessor[2];
		
		ImageProcessor crossCorrIm;
		
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
		
		ct1 = doComplexInverseTransform(istemp,0, 0);
		
		crossCorrIm = ((FloatProcessor) ct1.getProcessor(1)).duplicate();

		fht1.swapQuadrants(crossCorrIm);
		
		
		// normalize
		if(nNorm==0)
		{
			normIP = calcNormCorrelationCoeff(paddedip1,paddedip2);
			for(i=0;i<imwidthFFT;i++)
				for(j=0;j<imheightFFT;j++)
				{
					dVal=crossCorrIm.getf(i,j)/normIP.getf(i,j);
					crossCorrIm.setf(i,j,dVal);
				}
		}
		else
		{
			fNorm = calcSimpleCorrelationCoeff(paddedip1,paddedip2);
			crossCorrIm.multiply(1.0/fNorm);
		}
		//new ImagePlus("xcorr", crossCorrIm.duplicate()).show();
		
		//crop to original size
		crossCorrIm.setRoi(new Rectangle((int)((nPaddedS-origW)*0.5),(int)((nPaddedS-origH)*0.5),origW,origH));
		
	    return crossCorrIm.crop();
		
	}

	/** Function calculates slow (direct) cross correlation between to images
	 * by shifting them and calculating result
	 * **/
	public ImageProcessor calcDirectCorrelationImage(ImageProcessor ipp1, ImageProcessor ipp2)
	{
		
		ImageProcessor ip1,ip2;

		//to have the same size as FFT cross correlation for comparison
		ip1=padzeros(ipp1);
		ip2=padzeros(ipp2);
		
		int origW = ipp1.getWidth();
	    int origH = ipp1.getHeight();
		int nCorrW = ip1.getWidth();
		
		
	    double dCC;//,dCC1,dCC2;
	    int i,j,m,n;
	    double val1,val2;
	    //double mean1,mean2;
	    
	    //mean1 = ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean;
	    //mean2 = ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean;
	    
	    FloatProcessor crosscorr= new FloatProcessor (nCorrW,nCorrW);
	    
	    int dx,dy;
	    int dMaxx=(int)Math.round(nCorrW*0.5-1);
	    int dMaxy=(int)Math.round(nCorrW*0.5-1);
	    
	    //correlation shifts
	    for (dx=-dMaxx;dx<=dMaxx;dx++)
	    	for (dy=-dMaxy;dy<=dMaxy;dy++)
	    	{
	    		//now calculate correlation value for these shifts
	    		dCC=0;
				//dCC1=0;
				//dCC2=0;
				for(i=0;i<nCorrW;i++)
					for(j=0;j<nCorrW;j++)
					{
						m=i+dx;
						n=j+dy;
						if(m>-1 &&m<nCorrW&& n>-1 &&n<nCorrW)
						{
							val1=ip1.getf(m,n);//-mean1;
							val2=ip2.getf(i,j);//-mean2;
							dCC+=val1*val2;
							//dCC1+=val1*val1;
							//dCC2+=val2*val2;
						}
					}
						
				//if(dCC1!=0.0 && dCC2!=0.0)
					//crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),(float)(dCC/(Math.sqrt(dCC1)*Math.sqrt(dCC2))));
					crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),(float)(dCC));
				//else
					//crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),0);
	    	}
	    
	    //normalization	    
	    ImageProcessor normIP;
	    float dVal;
	    normIP = calcNormCorrelationCoeff(ip1,ip2);
		for(i=0;i<nCorrW;i++)
			for(j=0;j<nCorrW;j++)
			{
				dVal=crosscorr.getf(i,j)/normIP.getf(i,j);
				crosscorr.setf(i,j,dVal);
			}
		
		//crop to original size
		crosscorr.setRoi(new Rectangle((int)((nCorrW-origW)*0.5),(int)((nCorrW-origH)*0.5),origW,origH));
			    
	    return crosscorr.crop();
	}
	
	
	/** Function calculates normalization coefficient (denominator) 
	 *  for correlation based on image overlay area
	 * **/
	public ImageProcessor calcNormCorrelationCoeff(ImageProcessor ip1, ImageProcessor ip2)
	{
		FloatProcessor fpNorm, fpIntegr1, fpIntegr2;
		//FloatProcessor fpCheck1, fpCheck2;
		int imSize = ip1.getWidth(); //width is assumed to be the same as height
		fpNorm = new FloatProcessor(imSize,imSize);
		fpIntegr1 = new FloatProcessor(imSize+1,imSize+1);
		fpIntegr2 = new FloatProcessor(imSize+1,imSize+1);
		//fpCheck1 = new FloatProcessor(imSize,imSize);
		//fpCheck2 = new FloatProcessor(imSize,imSize);
		float dVal1=0, dVal2=0;
		int i,j,dx,dy;
		int imHalfSize;
		
		
		// let's calculate squared integral running sum images
		for (i=0;i<imSize;i++)
			for (j=0;j<imSize;j++)
			{
				dVal1 = ip1.getf(i, j);
				dVal1 = dVal1*dVal1 + fpIntegr1.getf(i,j+1)+ fpIntegr1.getf(i+1,j)-fpIntegr1.getf(i,j);
				fpIntegr1.setf(i+1, j+1, dVal1);
				
				dVal2 = ip2.getf(i, j);
				dVal2 = dVal2*dVal2 + fpIntegr2.getf(i,j+1)+ fpIntegr2.getf(i+1,j)-fpIntegr2.getf(i,j);
				fpIntegr2.setf(i+1, j+1, dVal2);				
			}
		
		//now let's calculate normalization values for each image
		imHalfSize = (int)(imSize*0.5);
		for(dx=(-imHalfSize);dx<imHalfSize;dx++)
			for(dy=(-imHalfSize);dy<imHalfSize;dy++)
			{
				//for both positive values (bottom right square)
				if(dx>=0 && dy>=0)
				{
					 dVal1=fpIntegr1.getf(imSize,imSize)-fpIntegr1.getf(imSize,dy)-fpIntegr1.getf(dx,imSize)+fpIntegr1.getf(dx,dy);
				     dVal2=fpIntegr2.getf(imSize-dx,imSize-dy);
				}
				//for both negative values (top left square)
				if(dx<=0 && dy<=0)
				{
					dVal1=fpIntegr1.getf(imSize+dx,imSize+dy); 
					dVal2=fpIntegr2.getf(imSize,imSize)-fpIntegr2.getf(imSize,-dy)-fpIntegr2.getf(-dx,imSize)+fpIntegr2.getf(-dx,-dy);
				}
				//for x positive and y negative (bottom left square)
				if(dx>0 && dy<0)
				{       
					dVal1=fpIntegr1.getf(imSize,imSize+dy)-fpIntegr1.getf(dx,imSize+dy);
				    dVal2=fpIntegr2.getf(imSize-dx,imSize)-fpIntegr2.getf(imSize-dx,-dy);
				}
				//for x negative and y positive (top right square)
				if(dx<0 && dy>0)
				{
					dVal1=fpIntegr1.getf(imSize+dx,imSize)-fpIntegr1.getf(imSize+dx,dy);
					dVal2=fpIntegr2.getf(imSize,imSize-dy)-fpIntegr2.getf(-dx,imSize-dy);
				} 
				
				fpNorm.setf(dx+imHalfSize,dy+imHalfSize, (float)((Math.sqrt(dVal1*dVal2))));
				//fpCheck1.setf(dx+imHalfSize,dy+imHalfSize, dVal1);
				//fpCheck2.setf(dx+imHalfSize,dy+imHalfSize, dVal2);
				
			}
		//new ImagePlus("xcorrcheck1", fpCheck1.duplicate()).show();
		//new ImagePlus("xcorrcheck2", fpCheck2.duplicate()).show();
		return fpNorm;
	}
	
	/** Function calculates normalization coefficient (denominator) 
	 *  for correlation based on whole image area
	 * **/
	public float calcSimpleCorrelationCoeff(ImageProcessor ip1, ImageProcessor ip2)
	{
		int imSize = ip1.getWidth(); //width is assumed to be the same as height
		float dVal1=0;
		float dVal2=0;
		float dVal;
		int i,j;
		
		for (i=0;i<imSize;i++)
			for (j=0;j<imSize;j++)
			{
				dVal=ip1.getf(i, j);
				dVal1+=dVal*dVal;
				dVal=ip2.getf(i, j);
				dVal2+=dVal*dVal;				
			}
		return (float)Math.sqrt(dVal1*dVal2);
	}
	

	/** Function calculates shift in X and Y between images using
	 * cross correlation calculated through FFT transform
	 * */
	public double [] calcShiftFFTCorrelationDouble(ImageProcessor ip1, ImageProcessor ip2, int nNorm)
	{
		int [] xyshift;
		//int i;
		double [] xyshiftd;
		
		ImageProcessor crossCorrIm = calcFFTCorrelationImage(ip1, ip2, nNorm);
		//new ImagePlus("zz",crossCorrIm).show();
		xyshift = getmaxpositions(crossCorrIm);
		xyshiftd= new double[2];

		xyshiftd[0]=(double)xyshift[0]-Math.round(crossCorrIm.getWidth()*0.5);
		xyshiftd[1]=(double)xyshift[1]-Math.round(crossCorrIm.getHeight()*0.5);
		
		return xyshiftd;
		
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
    
    public static ImageProcessor padzeros(ImageProcessor ip) 
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
   
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(0);
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
    public static void swapQuadrants(ImageStack stack) {
        FHT fht = new FHT(new FloatProcessor(1, 1));
        fht.setShowProgress(false);
        for (int i=1; i<=stack.getSize(); i++)
            fht.swapQuadrants(stack.getProcessor(i));
    }
    
    /** Complex to Complex Inverse Fourier Transform
    *   @author Joachim Wesner
    */
    public static void c2c2DFFT(float[] rein, float[] imin, int maxN, float[] reout, float[] imout) 
    {
            FHT fht = new FHT(new FloatProcessor(maxN,maxN));
            fht.setShowProgress(false);
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
    public static void cplxFHT(int row, int maxN, float[] re, float[] im, boolean reim, float[] fht) 
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
    
    public static ImageStack doComplexInverseTransform(ImageStack stack, int width, int height ) 
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
    
    public static ImageStack unpad(ImageStack stack, int width, int height) 
    {
      
 
        if (width==0 || height==0 || (width==stack.getWidth()&&height==stack.getHeight()))
            return stack;
        StackProcessor sp = new StackProcessor(stack, null);
        ImageStack stack2 = sp.crop(0, 0, width, height);
        return stack2;
    }
    
    /** Gets x and y coordinates of maximum intensity pixel on image. */	
    int [] getmaxpositions(ImageProcessor ipp)
    {
    	int [] results = new int [2];
    	results[0]=0;
    	results[1]=0;

    	double s;
    	double smax=Double.MIN_VALUE;
    	
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
