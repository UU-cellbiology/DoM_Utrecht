package DOM;


import jaolho.data.lma.LMA;

import java.awt.Color;

import java.util.Stack;


import ij.IJ;

import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;

import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.TypeConverter;


public class SMLAnalysis {
	
	ImageStatistics imgstat;
	GaussianBlur lowpassGauss = new GaussianBlur(); //low pass prefiltering
	//Convolver      colvolveOp = new Convolver(); //convolution filter
	float []		 fConKernel;  				//convolution kernel (Gaussian mexican hat)
	float []		 fLPKernel;  				//low-pass Gaussian kernel 
	
	ResultsTable ptable = Analyzer.getResultsTable(); // table with particle's approximate center coordinates
	
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();
	
	
	//default constructor 
	SMLAnalysis()
	{
		ptable.setPrecision(5);
		
	}
	
	
	
	// Particle finding routine based on spots enhancement with
	// 2D PSF Gaussian approximated convolution/backgrounds subtraction, thresholding
	// and particle filtering
	void detectParticles(ImageProcessor ip, SMLDialog fdg, int nFrame, Overlay SpotsPositions_)
	{
		int nThreshold;
		FloatProcessor dupip = null ; //duplicate of image
		ImageProcessor dushort; //duplicate of image
		ByteProcessor dubyte = null; //tresholded image
		TypeConverter tc; 
		
		dupip = (FloatProcessor) ip.duplicate().convertToFloat();
				
		SMLblur1Direction(dupip, fdg.dPSFsigma*0.5, 0.0002, true, (int)Math.ceil(5*fdg.dPSFsigma*0.5));
		SMLblur1Direction(dupip, fdg.dPSFsigma*0.5, 0.0002, false, 0);
		
		//new ImagePlus("gassconvoluted", dupip.convertToFloat().duplicate()).show();
		//low-pass filtering by gaussian blurring
		//lowpassGauss.blurGaussian(dupip, fdg.dPSFsigma*0.5, fdg.dPSFsigma*0.5, 0.0002);
		
		//convolution with gaussian PSF kernel		
		SMLconvolveFloat(dupip, fConKernel, fdg.nKernelSize, fdg.nKernelSize);
		
		//new ImagePlus("convoluted", dupip.duplicate()).show();
		tc = new TypeConverter(dupip, true);
		dushort =  tc.convertToShort();

		  
		
		 
		//thresholding
		
		//old straightforward thresholding
		//imgstat = ImageStatistics.getStatistics(dupip, 22, null); //6 means MEAN + STD_DEV, look at ij.measure.Measurements
		//nThreshold = (int)(imgstat.mean + 3.0*imgstat.stdDev);

		//new smart thresholding
		nThreshold = getThreshold(dushort);
		
		dushort.threshold(nThreshold);
		//convert to byte
		dubyte  = (ByteProcessor) dushort.convertToByte(false);
		//new ImagePlus("threshold", dubyte.duplicate()).show();
		
		//morphological operations on thresholded image	
		//dubyte.dilate(2, 0);
		//dubyte.erode(2, 0);
		
		dubyte.dilate();		
		//new ImagePlus("dilated", dubyte.duplicate()).show();
		dubyte.erode();
		//new ImagePlus("erosion", dubyte.duplicate()).show();

		//dupip.invert();
		
		labelParticles(dubyte, ip, nFrame, fdg.dPixelSize, fdg.nAreaCut, fdg.dPSFsigma, SpotsPositions_, fdg.bShowParticles, fdg.bIgnoreFP);//, fdg.dSymmetry/100);
		//labelParticles(dubyte, duconv, nFrame, fdg.dPixelSize, fdg.nAreaCut, fdg.dPSFsigma, fdg.dSymmetry/100);
		

		/*ImagePlus imp = new ImagePlus("result_decon", duconv);  
		imp.show();
		imp = new ImagePlus("result_threshold", dubyte);  
		imp.show();*/

		//ptable.show("Results");

		//ip.setPixels(pixels)
		//ImagePlus imp2 = new ImagePlus("result2", dupip);  

		//imp2.show();  		
	}
	
	//function that finds centroid x,y and area
	//of spots after thresholding
	//based on connected components	labeling Java code
	//implemented by Mariusz Jankowski & Jens-Peer Kuska
	//and in March 2012 available by link
    //http://www.izbi.uni-leipzig.de/izbi/publikationen/publi_2004/IMS2004_JankowskiKuska.pdf
	
	void labelParticles(ImageProcessor ipBinary, ImageProcessor ipRaw,  int nFrame, double dPixelSize_, int nAreaCut, double dPSFsigma_, Overlay SpotsPositions__, boolean bShow, boolean bIgnore)//, double dSymmetry_)
	{
		int width = ipBinary.getWidth();
		int height = ipBinary.getHeight();
		int nFitRadius = 3; //radius in number of SD around center point to fit Gaussian
		int dBorder; // radius in pixels around center point to fit Gaussian
		int nIntBoxSize = 3; //half length of box for spot integration (in SD around) 
		
		
		int nArea;
		double nFalsePositive;
		int i,j, nCount;
		double dVal, dInt;
		double dIMax, dIMin;
		double [][] spotInt;
	
		double [] dFitParams;
		double [] dFitErrors;
		double dErrCoeff;
		LMA SMLlma;
		SMLTwoDGaussian GFit = new SMLTwoDGaussian(); 

		double xCentroid, yCentroid, xSD, ySD;
		double dIntAmp, dIntNoise;
		double dNoiseAvrg, dNoiseSD;
		double dSNR;
		boolean bBorder;

		int lab = 1;
		int [] pos ;		
		
		Stack<int[]> sstack = new Stack<int[]>( );
		Stack<int[]> stackPost = new Stack<int[]>( );
		int [][] label = new int[width][height] ;
		
		OvalRoi spotROI;

		dBorder= (int)(dPSFsigma_*nFitRadius);
		spotInt = new double [(2*dBorder+1)*(2*dBorder+1)][3];
		dFitParams = new double [6];
		dFitErrors = new double [6];
		
		
		/*ImageProcessor afterFit; //duplicate of image
		afterFit = ipRaw.duplicate();
		double dGValue;
		double [] dCorrsxy = new double [2];*/
		
		for (int r = 1; r < width-1; r++)
			for (int c = 1; c < height-1; c++) {
				
				if (ipBinary.getPixel(r,c) == 0.0) continue ;
				if (label[r][c] > 0.0) continue ;
				/* encountered unlabeled foreground pixel at position r, c */
				/* it means it is a new spot! */
				/* push the position in the stack and assign label */
				sstack.push(new int [] {r, c}) ;
				stackPost.push(new int [] {r, c}) ;
				label[r][c] = lab ;
				nArea = 0;
				dIMax = -1000;
				xCentroid = 0; yCentroid = 0;
				//xMax = 0; yMax = 0;
				dInt = 0;
				bBorder = false;

				/* start the float fill */
				while ( !sstack.isEmpty()) 
				{
					pos = (int[]) sstack.pop() ;
					i = pos[0]; j = pos[1];
					
					//remove all spots at border
					if(i==0 || j==0 || i==(width-1) || j==(height-1))
						bBorder = true;
					nArea ++;
					dVal = ipRaw.getPixel(i,j);
					if (dVal > dIMax)
						dIMax = dVal;
					dInt += dVal;
					xCentroid += dVal*i;
					yCentroid += dVal*j;
					
					
					
					if (ipBinary.getPixel(i-1,j-1) > 0 && label[i-1][j-1] == 0) {
						sstack.push( new int[] {i-1,j-1} );
						stackPost.push( new int[] {i-1,j-1} );
						label[i-1][j-1] = lab ;
					}
					
					if (ipBinary.getPixel(i-1,j) > 0 && label[i-1][j] == 0) {
						sstack.push( new int[] {i-1,j} );
						stackPost.push( new int[] {i-1,j} );
						label[i-1][j] = lab ;
					}
					
					if (ipBinary.getPixel(i-1,j+1) > 0 && label[i-1][j+1] == 0) {
						sstack.push( new int[] {i-1,j+1} );
						stackPost.push( new int[] {i-1,j+1} );
						label[i-1][j+1] = lab ;
					}
					
					if (ipBinary.getPixel(i,j-1) > 0 && label[i][j-1] == 0) {
						sstack.push( new int[] {i,j-1} );
						stackPost.push( new int[] {i,j-1} );
						label[i][j-1] = lab ;
					}
					
					if (ipBinary.getPixel(i,j+1) > 0 && label[i][j+1] == 0) {
						sstack.push( new int[] {i,j+1} );
						stackPost.push( new int[] {i,j+1} );
						label[i][j+1] = lab ;
					}
					if (ipBinary.getPixel(i+1,j-1) > 0 && label[i+1][j-1] == 0) {
						sstack.push( new int[] {i+1,j-1} );
						stackPost.push( new int[] {i+1,j-1} );
						label[i+1][j-1] = lab ;
					}
				
					if (ipBinary.getPixel(i+1,j)>0 && label[i+1][j] == 0) {
						sstack.push( new int[] {i+1,j} );
						stackPost.push( new int[] {i+1,j} );
						label[i+1][j] = lab ;
					}
					
					if (ipBinary.getPixel(i+1,j+1) > 0 && label[i+1][j+1] == 0) {
						sstack.push( new int[] {i+1,j+1} );
						stackPost.push( new int[] {i+1,j+1} );
						label[i+1][j+1] = lab ;
					}
					
				} /* end while */
				if(!bBorder && nArea > nAreaCut)
				{					
					xCentroid /= dInt;
					yCentroid /= dInt;

						if ( (xCentroid>dBorder) && (yCentroid>dBorder) && (xCentroid<(width-1-dBorder)) && (yCentroid<(height-1-dBorder)) )
						{
							//fitting of spot
							
							//creating array of intensity values for fitting
							dIMin = 1000000;
							nCount = 0;
							for(i = (int) (Math.round(xCentroid)- dBorder); i <= Math.round(xCentroid)+ dBorder; i++)
								for(j = (int) (Math.round(yCentroid)- dBorder); j <= Math.round(yCentroid)+ dBorder; j++)
								{
									dVal = ipRaw.getPixel(i,j);
									if (dVal<dIMin)
										dIMin = dVal;
									spotInt[nCount][0] =dVal;
									spotInt[nCount][1] =(double)i;
									spotInt[nCount][2] =(double)j;
									nCount++;
									
								}
							dFitParams[0]= dIMin;
							dFitParams[1]= dIMax - dIMin;
							dFitParams[2]= xCentroid;
							dFitParams[3]= yCentroid;
							dFitParams[4]= dPSFsigma_;
							dFitParams[5]= dPSFsigma_;
							SMLlma = new LMA(GFit, dFitParams, spotInt);
							SMLlma.fit();
							
							dFitErrors = SMLlma.getStandardErrorsOfParameters();
							// scaling coefficient for parameters errors estimation 
							// (Standard deviation of residuals)
							dErrCoeff = Math.sqrt(SMLlma.chi2/(nCount-6));
							for (i=0;i<6;i++)
								dFitErrors[i] *= dErrCoeff; 
							nFalsePositive = 0.0;
							//iterations didn't converge. suspicious point
							if (SMLlma.iterationCount== 101)
								nFalsePositive = 0.5;
							//spot is too big
							if((SMLlma.parameters[4] > 2.0*dPSFsigma_) || (SMLlma.parameters[4]<0.5*dPSFsigma_)||(SMLlma.parameters[5]<0.5*dPSFsigma_)||(SMLlma.parameters[5] > 2.0*dPSFsigma_))
								nFalsePositive = 1.0;
							//localization precision is bigger than PSF size
							if((dFitErrors[2] > dPSFsigma_) || (dFitErrors[3] > dPSFsigma_))
								nFalsePositive = 1.0;
					
							//calculating integrated spot intensity and estimating SNR
							dIntAmp =0;
							dSNR = 0;
							dIntNoise = 0;
							dNoiseAvrg = 0;
							dNoiseSD = 0;
							if(nFalsePositive<0.5)
							{
								xCentroid = Math.round(SMLlma.parameters[2]);
								yCentroid = Math.round(SMLlma.parameters[3]);
								xSD = Math.round(SMLlma.parameters[4]*nIntBoxSize);
								ySD = Math.round(SMLlma.parameters[5]*nIntBoxSize);
								if( ((xCentroid-xSD-1)>0) && ((yCentroid-ySD-1)>0) && ((xCentroid+xSD+1)<(width-1)) && ((yCentroid+ySD+1)<(height-1)))
								{
									//integrated intensity in 3SD * 3SD region
									for(i=(int) (xCentroid-xSD); i<=(xCentroid+xSD); i++)
										for(j=(int) (yCentroid-ySD); j<=(yCentroid+ySD); j++)
											dIntAmp += ipRaw.getPixel(i,j);
									
									//averaged noise around spot
									j = (int)(yCentroid-ySD-1);
									for(i=(int) (xCentroid-xSD-1); i<=(xCentroid+xSD+1); i++)							
										dIntNoise += ipRaw.getPixel(i,j);
									j = (int)(yCentroid+ySD+1);
									for(i=(int) (xCentroid-xSD-1); i<=(xCentroid+xSD+1); i++)							
										dIntNoise += ipRaw.getPixel(i,j);
									i=(int) (xCentroid-xSD-1);
									for(j=(int) (yCentroid-ySD); j<=(yCentroid+ySD); j++)
										dIntNoise += ipRaw.getPixel(i,j);
									i=(int) (xCentroid+xSD+1);
									for(j=(int) (yCentroid-ySD); j<=(yCentroid+ySD); j++)
										dIntNoise += ipRaw.getPixel(i,j);
									dNoiseAvrg = dIntNoise/(4*xSD+4*ySD+8);
									dIntAmp = dIntAmp -  (dNoiseAvrg*((2*xSD+1)*(2*ySD+1))); 
									
									//SD of noise
									j = (int)(yCentroid-ySD-1);
									for(i=(int) (xCentroid-xSD-1); i<=(xCentroid+xSD+1); i++)							
										dNoiseSD += Math.pow(dNoiseAvrg - ipRaw.getPixel(i,j),2);
									j = (int)(yCentroid+ySD+1);
									for(i=(int) (xCentroid-xSD-1); i<=(xCentroid+xSD+1); i++)							
										dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
									i=(int) (xCentroid-xSD-1);
									for(j=(int) (yCentroid-ySD); j<=(yCentroid+ySD); j++)
										dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
									i=(int) (xCentroid+xSD+1);
									for(j=(int) (yCentroid-ySD); j<=(yCentroid+ySD); j++)
										dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
									dNoiseSD = Math.sqrt(dNoiseSD/(4*xSD+4*ySD+8));
									dSNR =  SMLlma.parameters[1]/dNoiseSD;
								}
							}

							if(SMLlma.parameters[2]>0 && SMLlma.parameters[2]<width && SMLlma.parameters[3]>0 && SMLlma.parameters[3]<height)
							{
								if(!(bIgnore && nFalsePositive>0.3))
								{
									ptable_lock.lock();
									ptable.incrementCounter();
	
									ptable.addValue("Amplitude_fit",SMLlma.parameters[1]);
									
									ptable.addValue("X (px)",SMLlma.parameters[2]);							
									ptable.addValue("Y (px)",SMLlma.parameters[3]);
									ptable.addValue("X (nm)",SMLlma.parameters[2]*dPixelSize_);							
									ptable.addValue("Y (nm)",SMLlma.parameters[3]*dPixelSize_);
									ptable.addValue("Z (nm)",0);
									ptable.addValue("False positive?", nFalsePositive);
									ptable.addValue("X_loc_error, pix", dFitErrors[2]);
									ptable.addValue("Y_loc_error, pix", dFitErrors[3]);
	
									ptable.addValue("BGfit",SMLlma.parameters[0]);							
									ptable.addValue("IntegratedInt",dIntAmp);
									ptable.addValue("SNR", dSNR);
	
									ptable.addValue("chi2_fit",SMLlma.chi2);
									ptable.addValue("Frame Number", nFrame+1);					
									ptable.addValue("Iterations_fit",SMLlma.iterationCount);
									ptable.addValue("SD_X_fit, pix",SMLlma.parameters[4]);
									ptable.addValue("SD_Y_fit, pix",SMLlma.parameters[5]);
									ptable.addValue("Amp_loc_error",dFitErrors[1]);
									
									ptable_lock.unlock();
									if(bShow)
									{
										spotROI = new OvalRoi((int)(xCentroid-2*dPSFsigma_),(int)(yCentroid-2*dPSFsigma_),(int)(4.0*dPSFsigma_),(int)(4.0*dPSFsigma_));
										if(nFalsePositive<0.5)
											spotROI.setStrokeColor(Color.green);
										else
											if(nFalsePositive<1)
												spotROI.setStrokeColor(Color.yellow);
											else
												spotROI.setStrokeColor(Color.red);
												
										spotROI.setPosition(nFrame+1);
										SpotsPositions__.add(spotROI);
									}
								}
						}
					}
				
				}
				lab++ ;
			} // end for cycle
		
		/*ImagePlus ipFit = new ImagePlus("after_fit", afterFit);  
		ipFit.show();*/
		return;// label ;

		
	}
	


	// function calculating convolution kernel of 2D Gaussian shape
	// with background subtraction for spots enhancement
	void initConvKernel(SMLDialog fdg)
	{
		int i,j; //counters
		int nBgPixCount; //number of kernel pixels for background subtraction 
		float GaussSum; //sum of all integrated Gaussian function values inside 'spot circle'
		float fSpot; // spot circle radius
		float fSpotSqr; // spot circle radius squared
		
		//intermediate values to speed up calculations
		float fIntensity;
		float fDivFactor;
		float fDist;
		float fPixVal;
		float [][] fKernel;
		
		float nCenter = (float) ((fdg.nKernelSize - 1.0)*0.5);; // center coordinate of the convolution kernel
		
		//kernel matrix
		fKernel = new float [fdg.nKernelSize][fdg.nKernelSize];
		//kernel string 
		fConKernel = new float [fdg.nKernelSize*fdg.nKernelSize];
		
		//Gaussian spot region
		if (3*fdg.dPSFsigma > nCenter)
			fSpot = nCenter;
		else
			fSpot = (float) (3.0*fdg.dPSFsigma);
		
		fSpotSqr = fSpot*fSpot;
		
		//intermediate values to speed up calculations
		fIntensity = (float) (fdg.dPSFsigma*fdg.dPSFsigma*0.5*Math.PI);
		fDivFactor = (float) (1.0/(Math.sqrt(2)*fdg.dPSFsigma));
		GaussSum = 0;
		nBgPixCount = 0;
		
		//first run, filling array with gaussian function Approximation values (integrated over pixel)
		//and calculating number of pixels which will serve as a background
		for (i=0; i<fdg.nKernelSize; i++)
		{
			for (j=0; j<fdg.nKernelSize; j++)
			{
				fDist = (i-nCenter)*(i-nCenter) + (j-nCenter)*(j-nCenter);
				
				//background pixels
				if (fDist > fSpotSqr)
					nBgPixCount++;
				
				//gaussian addition
				fPixVal  = errorfunction((i-nCenter-0.5)*fDivFactor) - errorfunction((i-nCenter+0.5)*fDivFactor);
				fPixVal *= errorfunction((j-nCenter-0.5)*fDivFactor) - errorfunction((j-nCenter+0.5)*fDivFactor);
				fPixVal *= fIntensity;
				fKernel[i][j] = fPixVal;
				GaussSum += fPixVal;													
			}				
		}
		
		//background subtraction coefficient
		fDivFactor = (float)(1.0/(double)nBgPixCount);
		
		//second run, normalization and background subtraction
		for (i=0; i<fdg.nKernelSize; i++)
		{
			for (j=0; j<fdg.nKernelSize; j++)
			{
				fDist = (i-nCenter)*(i-nCenter) + (j-nCenter)*(j-nCenter);
				//normalization
				fConKernel[i+j*(fdg.nKernelSize)] = fKernel[i][j] / GaussSum;
				if (fDist > fSpotSqr)
				{
					//background subtraction
					fConKernel[i+j*(fdg.nKernelSize)] -=fDivFactor;
				}
				
			}		
		}
		
		return;
	}
	
	//returns value of mean intensity+3*SD based on 
	//fitting of image histogram to gaussian function
	int getThreshold(ImageProcessor thImage)
	{
		ImageStatistics imgstat;
		double  [][] dHistogram;		
		int nUpCount, nDownCount;
		double dUpRef, dDownRef;
		int nCount, nMaxCount;
		int i,j;
		double dMean, dSD;
		LMA fitlma;
		
		//mean, sd, min, max
		imgstat = ImageStatistics.getStatistics(thImage, 38, null);
		
		dUpRef   = (imgstat.mean + 3.0*imgstat.stdDev);
		dDownRef = (imgstat.mean - 3.0*imgstat.stdDev);
		if(dUpRef>imgstat.max)
			dUpRef=imgstat.max;
		if(dDownRef<imgstat.min)
			dDownRef=imgstat.min;
		
		nDownCount =   (int) ((dDownRef - imgstat.min)/imgstat.binSize);
		nUpCount =   (int) ((dUpRef - imgstat.min)/imgstat.binSize);
		if(nUpCount>255)
			nUpCount = 255;

		dHistogram = new double [2][nUpCount-nDownCount +1];
		j=0;
		nMaxCount = 0;
		dMean = imgstat.mean;
		for (i=nDownCount; i<=nUpCount; i++)
		{
			
			nCount=imgstat.histogram[i];
			
			dHistogram[0][j] = dDownRef + j*imgstat.binSize;
			dHistogram[1][j] = (double)nCount;
			if(nMaxCount < nCount)
			{
				nMaxCount= nCount;
				dMean = dDownRef + j*imgstat.binSize;
			}
			j++;
		}
		fitlma = new LMA(new SMLOneDGaussian(), new double[] {(double)nMaxCount, dMean, imgstat.stdDev}, dHistogram);
		fitlma.fit();
		dMean = fitlma.parameters[1];
		dSD = fitlma.parameters[2];
		
		if (fitlma.iterationCount== 101 || dMean<imgstat.min || dMean> imgstat.max ||  dSD < imgstat.min || dSD> imgstat.max)
			//fit somehow failed
			return (int)(imgstat.mean + 3.0*imgstat.stdDev);
		else
			return (int)(dMean + 3.0*dSD);
		
		
	}
	
	//implements error function calculation with precision ~ 10^(-4) 
	//used for calculation of Gaussian convolution kernel
	float errorfunction (double d)
	{
		float t = (float) (1.0/(1.0+Math.abs(d)*0.47047));
		float ans = (float) (1.0 - t*(0.3480242 - 0.0958798*t+0.7478556*t*t)*Math.exp((-1.0)*(d*d))); 
		
		if (d >= 0)
			return ans;
		else
			return -ans;
		
	}
	
	//show results table
	void showTable()
	{
            ptable.show("Results");

	}
	
	boolean SMLconvolveFloat(ImageProcessor dupip, float[] kernel, int kw, int kh)
	{
		
		int width = dupip.getWidth();
		int height = dupip.getHeight();
		
		
		int x1 = 0;
		int y1 = 0;
		int x2 = x1 + width;
		int y2 = y1 + height;
		int uc = kw/2;    
		int vc = kh/2;
		float[] pixels = (float[])dupip.getPixels();
		float[] pixels2 = (float[])dupip.getSnapshotPixels();
		if (pixels2==null)
			pixels2 = (float[])dupip.getPixelsCopy();
       
		double sum;
		int offset, i;
		boolean edgePixel;
		int xedge = width-uc;
		int yedge = height-vc;
		
		for(int y=y1; y<y2; y++) {
			for(int x=x1; x<x2; x++) {
				sum = 0.0;
				i = 0;
				edgePixel = y<vc || y>=yedge || x<uc || x>=xedge;
				for(int v=-vc; v <= vc; v++) {
					offset = x+(y+v)*width;
					for(int u = -uc; u <= uc; u++) {
						if (edgePixel) {
 							if (i>=kernel.length) // work around for JIT compiler bug on Linux
 								IJ.log("kernel index error: "+i);
							sum += SMLgetPixel(x+u, y+v, pixels2, width, height)*kernel[i++];
						} else
							sum += pixels2[offset+u]*kernel[i++];
					}
		    	}
				pixels[x+y*width] = (float)(sum);
			}
    	}
	
   		return true;
	}
	
	float SMLgetPixel(int x, int y, float[] pixels, int width, int height) {
		if (x<=0) x = 0;
		if (x>=width) x = width-1;
		if (y<=0) y = 0;
		if (y>=height) y = height-1;
		return pixels[x+y*width];
	}
	
	 void SMLblur1Direction( final FloatProcessor ippp, final double sigma, final double accuracy,
            final boolean xDirection, final int extraLines) {
        
		final float[] pixels = (float[])ippp.getPixels();
        final int width = ippp.getWidth();
        final int height = ippp.getHeight();
        final int length = xDirection ? width : height;     //number of points per line (line can be a row or column)
        final int pointInc = xDirection ? 1 : width;        //increment of the pixels array index to the next point in a line
        final int lineInc = xDirection ? width : 1;         //increment of the pixels array index to the next line
        final int lineFromA = 0 - extraLines;  //the first line to process
        final int lineFrom;
        if (lineFromA < 0) lineFrom = 0;
        else lineFrom = lineFromA;
        final int lineToA = (xDirection ? height : width) + extraLines; //the last line+1 to process
        final int lineTo;
        if (lineToA > (xDirection ? height:width)) lineTo = (xDirection ? height:width);
        else lineTo = lineToA;
        final int writeFrom = 0;    //first point of a line that needs to be written
        final int writeTo = xDirection ? width : height;
 
        final float[][] gaussKernel = lowpassGauss.makeGaussianKernel(sigma, accuracy, length);
        final int kRadius = gaussKernel[0].length;             //Gaussian kernel radius after upscaling
        final int readFrom = (writeFrom-kRadius < 0) ? 0 : writeFrom-kRadius; //not including broadening by downscale&upscale
        final int readTo = (writeTo+kRadius > length) ? length : writeTo+kRadius;
        final int newLength = length;
       
           
          final float[] cache1 = new float[newLength];  //holds data before convolution (after downscaling, if any)
                  
          int pixel0 = 0;
          for (int line=lineFrom; line<lineTo; line ++, pixel0+=lineInc) 
          {
                                    int p = pixel0 + readFrom*pointInc;
                                    for (int i=readFrom; i<readTo; i++ ,p+=pointInc)
                                        cache1[i] = pixels[p];
                                    SMLconvolveLine(cache1, pixels, gaussKernel, readFrom, readTo, writeFrom, writeTo, pixel0, pointInc);
                               
                                    
         }
            
        return;
    }
	 
	    void SMLconvolveLine( final float[] input, final float[] pixels, final float[][] kernel, final int readFrom,
	            final int readTo, final int writeFrom, final int writeTo, final int point0, final int pointInc) {
	        final int length = input.length;
	        final float first = input[0];                 //out-of-edge pixels are replaced by nearest edge pixels
	        final float last = input[length-1];
	        final float[] kern = kernel[0];               //the kernel itself
	        final float kern0 = kern[0];
	        final float[] kernSum = kernel[1];            //the running sum over the kernel
	        final int kRadius = kern.length;
	        final int firstPart = kRadius < length ? kRadius : length;
	        int p = point0 + writeFrom*pointInc;
	        int i = writeFrom;
	        for (; i<firstPart; i++,p+=pointInc) {  //while the sum would include pixels < 0
	            float result = input[i]*kern0;
	            result += kernSum[i]*first;
	            if (i+kRadius>length) result += kernSum[length-i-1]*last;
	            for (int k=1; k<kRadius; k++) {
	                float v = 0;
	                if (i-k >= 0) v += input[i-k];
	                if (i+k<length) v+= input[i+k];
	                result += kern[k] * v;
	            }
	            pixels[p] = result;
	        }
	        final int iEndInside = length-kRadius<writeTo ? length-kRadius : writeTo;
	        for (;i<iEndInside;i++,p+=pointInc) {   //while only pixels within the line are be addressed (the easy case)
	            float result = input[i]*kern0;
	            for (int k=1; k<kRadius; k++)
	                result += kern[k] * (input[i-k] + input[i+k]);
	            pixels[p] = result;
	        }
	        for (; i<writeTo; i++,p+=pointInc) {    //while the sum would include pixels >= length 
	            float result = input[i]*kern0;
	            if (i<kRadius) result += kernSum[i]*first;
	            if (i+kRadius>=length) result += kernSum[length-i-1]*last;
	            for (int k=1; k<kRadius; k++) {
	                float v = 0;
	                if (i-k >= 0) v += input[i-k];
	                if (i+k<length) v+= input[i+k];
	                result += kern[k] * v;
	            }
	            pixels[p] = result;
	        }
	    }
	
	
}