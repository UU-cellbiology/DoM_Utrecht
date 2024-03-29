package fiji.plugin.DOM;


import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Stack;


import ij.IJ;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import jaolho.data.lma.LMA;

/** Comparator for sorting of pixels of thresholded mask
 * by intensity order 
 * **/
class ArrayComparator implements Comparator<int []> {
   
    public int compare(int [] o1, int [] o2) {

    	if(o1[2]==o2[2])
        	return 0;
        if(o1[2]>o2[2])
        	return -1;
        else
        	return 1;               
        
    }
}
public class SMLAnalysis {
	
	ImageStatistics imgstat;
	/** low pass prefiltering */
	GaussianBlur lowpassGauss = new GaussianBlur(); 
	/** convolution kernel (Gaussian mexican hat) */
	float []		 fConKernel;  				
	/** low-pass Gaussian kernel  */
	float []		 fLPKernel;  				
	/** table with results of detection and fitting */
	ResultsTable ptable;// = ResultsTable.getResultsTable();
	
	public Overlay PreviewOverlay_;
	
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();

    final static int[] DIR_X_OFFSET = new int[] {  0,  1,  1,  1,  0, -1, -1, -1 };
    final static int[] DIR_Y_OFFSET = new int[] { -1, -1,  0,  1,  1,  1,  0, -1 };
	
	/** default class constructor 
	 * initializes Results Table **/
	SMLAnalysis()
	{
		ptable = ResultsTable.getResultsTable();
		ptable.setPrecision(2);
		PreviewOverlay_=null;
		
	}
	
	
	
	/** Particle finding routine based on spots enhancement with
	 * 2D PSF Gaussian approximated convolution/backgrounds subtraction, thresholding
	 *  and particle filtering 
	 * **/
	void detectParticles(ImageProcessor ip, SMLDialog fdg, int nFrame, Overlay SpotsPositions_,Roi RoiActive_)
	{
		float [] nThreshold;

		/** float duplicate of image **/
		FloatProcessor dupip = null ; 

		/** tresholded image **/
		ByteProcessor dubyte = null; 
		/** watershed image of maxima **/
		//ByteProcessor dubyteWS = null; 
		//MaximumFinder maxF = new MaximumFinder();
	
		//TypeConverter tc; 
		
		dupip = (FloatProcessor) ip.duplicate().convertToFloat();
				
		SMLblur1Direction(dupip, fdg.dPSFsigma*0.5, 0.0002, true, (int)Math.ceil(5*fdg.dPSFsigma*0.5));
		SMLblur1Direction(dupip, fdg.dPSFsigma*0.5, 0.0002, false, 0);
		
		//new ImagePlus("gassconvoluted", dupip.convertToFloat().duplicate()).show();
		//low-pass filtering by gaussian blurring
		//lowpassGauss.blurGaussian(dupip, fdg.dPSFsigma*0.5, fdg.dPSFsigma*0.5, 0.0002);
		
		//convolution with gaussian PSF kernel		
		SMLconvolveFloat(dupip, fConKernel, fdg.nKernelSize, fdg.nKernelSize);
		
		//new ImagePlus("convoluted", dupip.duplicate()).show();
		//tc = new TypeConverter(dupip, true);
		//dushort =  tc.convertToShort();
		//new ImagePlus("convoluted", dushort.duplicate()).show();
		  
		
		 
		//thresholding
		
		//old straightforward thresholding
		//imgstat = ImageStatistics.getStatistics(dupip, 22, null); //6 means MEAN + STD_DEV, look at ij.measure.Measurements
		//nThreshold = (int)(imgstat.mean + 3.0*imgstat.stdDev);

		//new smart thresholding
		nThreshold = getThreshold(dupip);
		
		//dushort.threshold(nThreshold);
		//convert to byte
		//dubyte  = (ByteProcessor) dushort.convertToByte(false);
		dubyte  =thresholdFloat(dupip,(float)(nThreshold[0]+fdg.dSNR*nThreshold[1]));
		//dubyteWS= maxF.findMaxima(dupip, 3.0*nThreshold[1], nThreshold[0], MaximumFinder.SEGMENTED, true, false);
		//dubyteWS= maxF.findMaxima(dupip, 3.0*nThreshold[1], ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, false, false);
		//dubyteWS= maxF.findMaxima(dupip, 3.0*nThreshold[1], MaximumFinder.SEGMENTED, false);
		
		//new ImagePlus("convoluted", dupip.duplicate()).show();
		
		//new ImagePlus("watershed", dubyteWS.duplicate()).show();
		//new ImagePlus("threshold", dubyte.duplicate()).show();
		
		//dubyte.copyBits(dubyteWS, 0, 0, Blitter.AND);
		//dubyteWS=quickAND(dubyte,dubyteWS);
		//new ImagePlus("ANDresult", dubyte.duplicate()).show();
		//morphological operations on thresholded image	
	
		
		
		//cleaning up image a bit
		if(fdg.nKernelSize>5)
		{
			dubyte.dilate();		
			//new ImagePlus("dilated", dubyte.duplicate()).show();
			dubyte.erode();
			//new ImagePlus("erosion", dubyte.duplicate()).show();
		}		
			
		//new ImagePlus("threshold_eroded", dubyte.duplicate()).show();

		//dupip.invert();
		
		labelParticles(dubyte, ip, nFrame, fdg, SpotsPositions_, RoiActive_);//, fdg.bIgnoreFP);//, fdg.dSymmetry/100);
		

	}
	
	/** function that finds centroids x,y and area
	 * of spots after thresholding
	 * based on connected components labeling Java code
	 * implemented by Mariusz Jankowski & Jens-Peer Kuska
	 * and in March 2012 available by link
	 * http://www.izbi.uni-leipzig.de/izbi/publikationen/publi_2004/IMS2004_JankowskiKuska.pdf
	**/
	void labelParticles(ImageProcessor ipBinary, ImageProcessor ipRaw,  int nFrame, SMLDialog dlg, Overlay SpotsPositions__, Roi RoiAct)
	{
		int width = ipBinary.getWidth();
		int height = ipBinary.getHeight();

		int dBorder; // radius in pixels around center point to fit Gaussian
 				
		//double dPSFsigma_=dlg.dPSFsigma;
		int nArea;

		int i,j,k,m;
		int vNeighbor;

		double dVal, dInt;
		double dIMax, dIMin;

		double xCentroid, yCentroid;
		boolean bBorder;
		boolean bInRoi;
		boolean bIsLocMax;

		int lab = 1;
		int [] pos ;		
		/** stack for segmenting
		 * thresholded spots by flood fill algorithm **/
		Stack<int[]> sstack = new Stack<int[]>( );
		/** Extra storage of thresholded mask region
		 * in case it is bigger than typical spot size.
		 * Used to find local maxima within large spots
		 * **/
		ArrayList<int[]> stackPost = new ArrayList<int[]>( );
		
		int [][] label = new int[width][height] ;		
		
		OvalRoi spotROI;

		dBorder= (int)(dlg.dPSFsigma*DOMConstants.FITRADIUS);
				
		
		for (int r = 1; r < width-1; r++)
			for (int c = 1; c < height-1; c++) {
				
				if (ipBinary.getPixel(r,c) == 0.0) continue ;
				if (label[r][c] > 0.0) continue ;
				/* encountered unlabeled foreground pixel at position r, c */
				/* it means it is a new spot! */
				/* push the position in the stack and assign label */
				sstack.push(new int [] {r, c}) ;
				label[r][c] = lab ;
				stackPost.clear();
				stackPost.add(new int [] {r, c, ipRaw.getPixel(r,c)}) ;
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
						stackPost.add( new int[] {i-1,j-1,ipRaw.getPixel(i-1,j-1)} );						
						label[i-1][j-1] = lab ;
					}
					
					if (ipBinary.getPixel(i-1,j) > 0 && label[i-1][j] == 0) {
						sstack.push( new int[] {i-1,j} );
						stackPost.add( new int[] {i-1,j,ipRaw.getPixel(i-1,j)} );
						label[i-1][j] = lab ;
					}
					
					if (ipBinary.getPixel(i-1,j+1) > 0 && label[i-1][j+1] == 0) {
						sstack.push( new int[] {i-1,j+1} );
						stackPost.add( new int[] {i-1,j+1,ipRaw.getPixel(i-1,j+1)} );
						label[i-1][j+1] = lab ;
					}
					
					if (ipBinary.getPixel(i,j-1) > 0 && label[i][j-1] == 0) {
						sstack.push( new int[] {i,j-1} );
						stackPost.add( new int[] {i,j-1,ipRaw.getPixel(i,j-1)} );
						label[i][j-1] = lab ;
					}
					
					if (ipBinary.getPixel(i,j+1) > 0 && label[i][j+1] == 0) {
						sstack.push( new int[] {i,j+1} );
						stackPost.add( new int[] {i,j+1,ipRaw.getPixel(i,j+1)} );
						label[i][j+1] = lab ;
					}
					if (ipBinary.getPixel(i+1,j-1) > 0 && label[i+1][j-1] == 0) {
						sstack.push( new int[] {i+1,j-1} );
						stackPost.add( new int[] {i+1,j-1,ipRaw.getPixel(i+1,j-1)} );
						label[i+1][j-1] = lab ;
					}
				
					if (ipBinary.getPixel(i+1,j)>0 && label[i+1][j] == 0) {
						sstack.push( new int[] {i+1,j} );
						stackPost.add( new int[] {i+1,j,ipRaw.getPixel(i+1,j)} );
						label[i+1][j] = lab ;
					}
					
					if (ipBinary.getPixel(i+1,j+1) > 0 && label[i+1][j+1] == 0) {
						sstack.push( new int[] {i+1,j+1} );
						stackPost.add( new int[] {i+1,j+1,ipRaw.getPixel(i+1,j+1)} );
						label[i+1][j+1] = lab ;
					}
					
				} /* end while */
				
				// nice small thresholded area, everything is fine
				if(!bBorder && nArea > dlg.nAreaCut && nArea < dlg.nAreaMax)
				{					
					xCentroid /= dInt;
					yCentroid /= dInt;

						if ( (xCentroid>dBorder) && (yCentroid>dBorder) && (xCentroid<(width-1-dBorder)) && (yCentroid<(height-1-dBorder)) )
						{
							bInRoi = true;
							if(RoiAct!=null)
							{
								if(!RoiAct.contains((int)xCentroid, (int)yCentroid))
									bInRoi=false;
							}
							if(bInRoi)
							{
								/** normal detection **/
								if(PreviewOverlay_==null)
								{
									//find minimum value around peak position
									dIMin = 1000000;
									//nCount = 0;
									for(i = (int) (Math.round(xCentroid)- dBorder); i <= Math.round(xCentroid)+ dBorder; i++)
										for(j = (int) (Math.round(yCentroid)- dBorder); j <= Math.round(yCentroid)+ dBorder; j++)
										{
											dVal = ipRaw.getPixel(i,j);
											if (dVal<dIMin)
												dIMin = dVal;																	
										}
									
									ptable_lock.lock();
									ptable.incrementCounter();
									ptable.addValue("Frame Number", nFrame+1);
									ptable.addValue("X_centroid_(px)",xCentroid);	
									ptable.addValue("Y_centroid_(px)",yCentroid);
									ptable.addValue("MaxInt",dIMax);
									ptable.addValue("MinInt",dIMin);
											
									ptable_lock.unlock();
									if(dlg.bShowParticles)
									{
										spotROI = new OvalRoi(0.5+xCentroid-2*dlg.dPSFsigma,0.5+yCentroid-2*dlg.dPSFsigma,4.0*dlg.dPSFsigma,4.0*dlg.dPSFsigma);
										spotROI.setStrokeColor(Color.yellow);	
										spotROI.setPosition(nFrame+1);
										SpotsPositions__.add(spotROI);										
									}
								}
								/** preview **/
								else
								{
									spotROI = new OvalRoi(0.5+xCentroid-2*dlg.dPSFsigma,0.5+yCentroid-2*dlg.dPSFsigma,4.0*dlg.dPSFsigma,4.0*dlg.dPSFsigma);
									spotROI.setStrokeColor(Color.yellow);	
									//spotROI.setPosition(nFrame+1);
									PreviewOverlay_.add(spotROI);																			
								}
							}

					}
				
				}
				////probably many particles in the thresholded area	
				//let's look for the local maxima and use them
				if(nArea >= dlg.nAreaMax && !bBorder)
				{
					//sort array by intensity values
					Collections.sort(stackPost, new ArrayComparator());
					
					while ( !stackPost.isEmpty()) 
					{												
						pos = (int[]) stackPost.get(0);
						k = pos[0]; m = pos[1];
						dVal = pos[2];
						xCentroid = dVal*k;
						yCentroid = dVal*m;
						dInt = dVal;						
						bIsLocMax=true;
						int nEqual = 0;
						//"soft" maximum
						for (int d=0; d<8; d++) 
						{   							
							vNeighbor = ipRaw.getPixel(k+DIR_X_OFFSET[d], m+DIR_Y_OFFSET[d]);
							//if (vNeighbor < dVal) 
							if (vNeighbor <= dVal) 
							{
								xCentroid += vNeighbor*(k+DIR_X_OFFSET[d]);
								yCentroid += vNeighbor*(m+DIR_Y_OFFSET[d]);
								dInt+=vNeighbor;
								if(vNeighbor == dVal)
									nEqual++;
	                        }
							else
							{
								bIsLocMax=false;
								stackPost.remove(0);
								break;
							}
						}
						
						if(nEqual>4&&bIsLocMax)
						{
							bIsLocMax=false;
							stackPost.remove(0);							
						}
						if(bIsLocMax)
						{
							xCentroid /= dInt;
							yCentroid /= dInt;

							if ( (xCentroid>dBorder) && (yCentroid>dBorder) && (xCentroid<(width-1-dBorder)) && (yCentroid<(height-1-dBorder)) )
							{
								bInRoi = true;
								if(RoiAct!=null)
								{
									if(!RoiAct.contains((int)xCentroid, (int)yCentroid))
										bInRoi=false;
								}
								if(bInRoi)
								{
									/** normal detection **/
									if(PreviewOverlay_==null)
									{
										//find minimum value around peak position
										dIMin = 1000000;
										//nCount = 0;
										for(i = (int) (Math.round(xCentroid)- dBorder); i <= Math.round(xCentroid)+ dBorder; i++)
											for(j = (int) (Math.round(yCentroid)- dBorder); j <= Math.round(yCentroid)+ dBorder; j++)
											{
												dVal = ipRaw.getPixel(i,j);
												if (dVal<dIMin)
													dIMin = dVal;																	
											}
										
										ptable_lock.lock();
										ptable.incrementCounter();
										ptable.addValue("Frame Number", nFrame+1);
										ptable.addValue("X_centroid_(px)",xCentroid);	
										ptable.addValue("Y_centroid_(px)",yCentroid);
										ptable.addValue("MaxInt",dIMax);
										ptable.addValue("MinInt",dIMin);
												
										ptable_lock.unlock();
										if(dlg.bShowParticles)
										{
											spotROI = new OvalRoi(0.5+xCentroid-2*dlg.dPSFsigma,0.5+yCentroid-2*dlg.dPSFsigma,4.0*dlg.dPSFsigma,4.0*dlg.dPSFsigma);
											spotROI.setStrokeColor(Color.yellow);	
											spotROI.setPosition(nFrame+1);
											SpotsPositions__.add(spotROI);										
										}
									}
									/** preview **/
									else
									{
										spotROI = new OvalRoi(0.5+xCentroid-2*dlg.dPSFsigma,0.5+yCentroid-2*dlg.dPSFsigma,4.0*dlg.dPSFsigma,4.0*dlg.dPSFsigma);
										spotROI.setStrokeColor(Color.yellow);	
										//spotROI.setPosition(nFrame+1);
										PreviewOverlay_.add(spotROI);											
									}
								}
							}
							//remove max and its neighbors from the list
							int nListLength = stackPost.size();
							i=0;
							while(nListLength>0 && i<nListLength)
							{
								pos =stackPost.get(i);
								//if(nMaxPos[0]>=xCentroid-RoiRad && nMaxPos[0]<=xCentroid+RoiRad && nMaxPos[1]>=yCentroid-RoiRad && nMaxPos[1]<=yCentroid+RoiRad)
								if(Math.abs(xCentroid-pos[0])<3.0*dlg.dPSFsigma && Math.abs(yCentroid-pos[1])<3.0*dlg.dPSFsigma)
								{
									stackPost.remove(i);
									nListLength--;
								}
								else
								{
									i++;
								}
							}
						}//end of if(bIsLocMax)
						
					}
					
				}
				lab++ ;
			} // end for cycle
		
		return;// label ;

		
	}
	
	/** function fitting particles in specific frame
	 * 
	 * @param ipRaw original image processor
	 * @param fdg   dialog item containing all parameters
	 * @param particles_ particles to fit table
	 * @param nFrame  current frame
	 * @param SpotsPositions__  Image Overlay to add particle ROIs
	 */
	void fitParticles(ImageProcessor ipRaw, SMLDialog fdg, double[][] particles_, int nFrame, Overlay SpotsPositions__)
	{
		int width = ipRaw.getWidth();
		int height = ipRaw.getHeight();
		int x,y, i,j, nCount, nParticlesCount, nParticlesNumber;		
		double xCentroid, yCentroid, xSD, ySD;
		double dIntAmp, dIntNoise;
		double dIntAverage,dIntDev;
		double dIMax, dIMin, dInt;
		double dNoiseAvrg, dNoiseSD;
		double dSNR;
		int dBorder; // radius in pixels around center point to fit Gaussian
		double nFalsePositive;				
		
		double [][] spotInt;
		double [][] spotXpos;
		double [][] spotYpos;
		double [] dFitParams;
		double dErrCoeff;
		LevenbergMarquardt SMLlma = new LevenbergMarquardt();
		//LMA SMLlma;
		//SMLTwoDGaussian GFit = new SMLTwoDGaussian(); 
		OvalRoi spotROI;
		
		//initialization
		dBorder= (int)(fdg.dPSFsigma*DOMConstants.FITRADIUS);
		spotInt = new double [(2*dBorder+1)][(2*dBorder+1)];
		spotXpos = new double [(2*dBorder+1)][(2*dBorder+1)];
		spotYpos = new double [(2*dBorder+1)][(2*dBorder+1)];
		dFitParams = new double [6];

		nParticlesNumber = particles_[0].length;
				
		for(nParticlesCount=0;nParticlesCount<nParticlesNumber;nParticlesCount++)
		{
					
			//creating array of intensity values for fitting
			nCount = 0;
			dIMin = Double.MAX_VALUE;
			dIMax = -100;
			
			dIntAverage =0;
			
			for(x = 0, i = (int) (Math.round(particles_[0][nParticlesCount])- dBorder); i <= Math.round(particles_[0][nParticlesCount])+ dBorder; ++x, i++)
				for(y = 0, j = (int) (Math.round(particles_[1][nParticlesCount])- dBorder); j <= Math.round(particles_[1][nParticlesCount])+ dBorder; ++y, j++)
				{					
					dInt = ipRaw.getPixel(i,j);
					spotInt[y][x] = dInt;
					spotXpos[y][x] =(double)i;
					spotYpos[y][x] =(double)j;
					if(dInt>dIMax)
						dIMax=dInt;
					if(dInt<dIMin)
						dIMin=dInt;
					nCount++;
					dIntAverage+=dInt;
				}
			//calculate deviation from average intensity 
			dIntAverage =dIntAverage/nCount;
			dIntDev=0;
			for(x = 0, i = (int) (Math.round(particles_[0][nParticlesCount])- dBorder); i <= Math.round(particles_[0][nParticlesCount])+ dBorder; ++x, i++)
				for(y = 0, j = (int) (Math.round(particles_[1][nParticlesCount])- dBorder); j <= Math.round(particles_[1][nParticlesCount])+ dBorder; ++y, j++)
				{
					dIntDev+=Math.pow(ipRaw.getPixel(i,j)-dIntAverage,2);		
				}			
			
		
				//initial values of fitting parameters				
				dFitParams[0]= dIMin; //minimum, background level
				dFitParams[1]= dIMax - dIMin; // intensity amplitude, max-min
				dFitParams[2]= particles_[0][nParticlesCount];//x center
				dFitParams[3]= particles_[1][nParticlesCount];//y center
				dFitParams[4]= fdg.dPSFsigma;
				dFitParams[5]= fdg.dPSFsigma;
				nFalsePositive = 0.0;
				//SMLlma = new LMA(GFit, dFitParams, spotInt);
				//SMLlma.setMaxIterations(fdg.nIterations);
				double[] fitted_parameters = new double[6];
				double[] fit_errors = new double[6];
				double chi2_fit = 0.0;
		
				// fitting
				fitted_parameters = SMLlma.run(spotInt, spotXpos, spotYpos, 2*dBorder+1, 2*dBorder+1, dFitParams, fdg.nIterations, 0.001);
				fit_errors = SMLlma.calculateStandardErrors(spotInt, spotXpos, spotYpos, 2*dBorder+1, 2*dBorder+1, fitted_parameters);
				chi2_fit = SMLlma.calculateChi2(spotInt, spotXpos, spotYpos, 2*dBorder+1, 2*dBorder+1, fitted_parameters);

				dErrCoeff = Math.sqrt(chi2_fit/(nCount-6));
				for(i =0;i<6;i++)
					fit_errors[i] *= dErrCoeff; 

				//checking for false positives
				
				//spot is too big or too small
				xSD = Math.abs(fitted_parameters[4]);
				ySD = Math.abs(fitted_parameters[5]);
				if((xSD > 1.3*fdg.dPSFsigma) || (xSD<0.70*fdg.dPSFsigma)||(ySD<0.70*fdg.dPSFsigma)||(ySD > 1.3*fdg.dPSFsigma))
					nFalsePositive = 1.0;
				
				//localization precision is bigger than PSF size
				if((fit_errors[2] > fdg.dPSFsigma) || (fit_errors[3] > fdg.dPSFsigma))
					nFalsePositive = 1.0;

				//fitting went out of fitting region
				if( Math.sqrt(Math.pow(particles_[0][nParticlesCount] -fitted_parameters[2],2) +  
						      Math.pow(particles_[1][nParticlesCount] -fitted_parameters[3],2))>fdg.dPSFsigma*DOMConstants.FITRADIUS)
					nFalsePositive = 1.0;
				
				
				dIntAmp =0; dIntNoise = 0;			
				dNoiseAvrg = 0;	dNoiseSD = 0;
				dSNR = 0;
				xCentroid = Math.round(fitted_parameters[2]+0.5);// 0.5 is for correction of pixel shift
				yCentroid = Math.round(fitted_parameters[3]+0.5);

				//calculating integrated spot intensity and estimating SNR
				xSD = Math.round(fdg.dPSFsigma*DOMConstants.FITRADIUS);
				ySD = Math.round(fdg.dPSFsigma*DOMConstants.FITRADIUS);

				if( ((xCentroid-xSD-1)>0) && ((yCentroid-ySD-1)>0) && ((xCentroid+xSD+1)<(width-1)) && ((yCentroid+ySD+1)<(height-1)))
				{
					//integrated intensity in 3SD * 3SD region
					for(i =(int) (xCentroid-xSD); i<=(int)(xCentroid+xSD); i++)
						for(j =(int) (yCentroid-ySD); j<=(int)(yCentroid+ySD); j++)
						{
							dIntAmp += ipRaw.getPixel(i,j);
						}
							
						

					//averaged noise around spot
				
					j = (int)(yCentroid-ySD-1);
					for(i =(int) (xCentroid-xSD-1); i<=(int)(xCentroid+xSD+1); i++)	
					{
						dIntNoise += ipRaw.getPixel(i,j);
						
					}
					j = (int)(yCentroid+ySD+1);
					for(i =(int) (xCentroid-xSD-1); i<=(int)(xCentroid+xSD+1); i++)
					{
						dIntNoise += ipRaw.getPixel(i,j);
					
					}
					i=(int) (xCentroid-xSD-1);
					for(j =(int) (yCentroid-ySD); j<=(int)(yCentroid+ySD); j++)
					{
						dIntNoise += ipRaw.getPixel(i,j);
						
					}
					i=(int) (xCentroid+xSD+1);
					for(j =(int) (yCentroid-ySD); j<=(int)(yCentroid+ySD); j++)
					{
						dIntNoise += ipRaw.getPixel(i,j);
						
					}
					dNoiseAvrg = dIntNoise/(4*xSD+4*ySD+8);
					dIntAmp = dIntAmp - (dNoiseAvrg*((2*xSD+1)*(2*ySD+1)));

					//SD of noise
					j = (int)(yCentroid-ySD-1);
					for(i =(int) (xCentroid-xSD-1); i<=(int)(xCentroid+xSD+1); i++)							
						dNoiseSD += Math.pow(dNoiseAvrg - ipRaw.getPixel(i,j),2);
					j = (int)(yCentroid+ySD+1);
					for(i =(int) (xCentroid-xSD-1); i<=(int)(xCentroid+xSD+1); i++)							
						dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
					i=(int) (xCentroid-xSD-1);
					for(j =(int) (yCentroid-ySD); j<=(int)(yCentroid+ySD); j++)
						dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
					i=(int) (xCentroid+xSD+1);
					for(j =(int) (yCentroid-ySD); j<=(int)(yCentroid+ySD); j++)
						dNoiseSD += Math.pow(dNoiseAvrg -ipRaw.getPixel(i,j),2);
					dNoiseSD = Math.sqrt(dNoiseSD/(4*xSD+4*ySD+8));
					if (nFalsePositive==0)						
						dSNR =  fitted_parameters[1]/dNoiseSD;
					else
						dSNR =  (dIMax - dNoiseAvrg)/dNoiseSD;
				}
				//}
				if(fitted_parameters[2]>0 && fitted_parameters[2]<width && fitted_parameters[3]>0 && fitted_parameters[3]<height && fitted_parameters[1]>0)
					{
						if(!(fdg.bIgnoreFP && nFalsePositive>0.3))
						{
							ptable_lock.lock();
							ptable.incrementCounter();
							
							ptable.addValue("X_(px)",fitted_parameters[2]+0.5);							
							ptable.addValue("Y_(px)",fitted_parameters[3]+0.5);
							ptable.addValue("Frame_Number", nFrame+1);
							ptable.addValue("X_(nm)",(fitted_parameters[2]+0.5)*fdg.dPixelSize);							
							ptable.addValue("X_loc_error(nm)", fit_errors[2]*fdg.dPixelSize);							
							ptable.addValue("Y_(nm)",(fitted_parameters[3]+0.5)*fdg.dPixelSize);
							ptable.addValue("Y_loc_error(nm)", fit_errors[3]*fdg.dPixelSize);
							ptable.addValue("Z_(nm)",0);
							ptable.addValue("Z_loc_error(nm)",0);							
							ptable.addValue("Amplitude_fit",fitted_parameters[1]);
							ptable.addValue("Amp_error",fit_errors[1]);
							ptable.addValue("BGfit",fitted_parameters[0]);
							ptable.addValue("BGfit_error",fit_errors[0]);
							ptable.addValue("SD_X_(nm)",Math.abs(fitted_parameters[4])*fdg.dPixelSize);
							ptable.addValue("SD_X_error(nm)",fit_errors[4]*fdg.dPixelSize);
							ptable.addValue("SD_Y_(nm)",Math.abs(fitted_parameters[5])*fdg.dPixelSize);
							ptable.addValue("SD_Y_error(nm)",fit_errors[5]*fdg.dPixelSize);
							ptable.addValue("False_positive", nFalsePositive);							
							ptable.addValue("IntegratedInt",dIntAmp);
							ptable.addValue("SNR", dSNR);	
							//ptable.addValue("chi2_fit",chi2_fit);
							ptable.addValue("R2_fit",(1-(chi2_fit/dIntDev)));												
							//ptable.addValue("Iterations_fit",fdg.nIterations);
							ptable.addValue("Iterations_fit",SMLlma.iteration_count);

							//case of importing from MTrackJ
							//keeping track information
							if (particles_.length==5)
							{
								ptable.addValue("Track_ID",particles_[2][nParticlesCount]);
								ptable.addValue("Particle_ID",particles_[3][nParticlesCount]);
								ptable.addValue("Track_Length",particles_[4][nParticlesCount]);
							}
				
							
							ptable_lock.unlock();
							if(fdg.bShowParticles)
							{
								spotROI = new OvalRoi(xCentroid-2*fdg.dPSFsigma,yCentroid-2*fdg.dPSFsigma,4.0*fdg.dPSFsigma,4.0*fdg.dPSFsigma);
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
				

			
		}//particles iteration
		//for(nParticlesCount=0;nParticlesCount<nParticlesNumber;nParticlesCount++)
		
	}

	/** function calculating convolution kernel of 2D Gaussian shape
	 * with background subtraction for spots enhancement (kind of mexican hat) **/
	void initConvKernel(SMLDialog fdg)
	{
		int i,j; //counters
		/** number of kernel pixels for background subtraction  **/
		int nBgPixCount;
		/** sum of all integrated Gaussian function values inside 'spot circle' **/
		float GaussSum; 
		/** spot circle radius **/
		float fSpot; 
		/**  spot circle radius squared **/
		float fSpotSqr; 
		
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
		for(i =0; i<fdg.nKernelSize; i++)
		{
			for(j =0; j<fdg.nKernelSize; j++)
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
		for(i =0; i<fdg.nKernelSize; i++)
		{
			for(j =0; j<fdg.nKernelSize; j++)
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
	
	/** Returns value of mean intensity [0] and SD [1] based on 
	 * fitting of image histogram to gaussian function.
	 * Estimates optimal bin size of histogram.
	 * @param thImage input imageprocessor
	 * @return
	 */
	float []  getThreshold(ImageProcessor thImage)
	{
		ImageStatistics imgstat;
		float [] results;
		double  [][] dHistogram;
		double  [] dHistCum;
		double  [][] dNoiseFit;
		int nHistSize;
		int nCount, nMaxCount;
		int nDownCount, nUpCount;
		int i,nPeakPos; 
		double dRightWidth, dLeftWidth, dWidth;
		double dMean, dSD;
		double [] dFitErrors;
		double dErrCoeff;
		double dSum=0.0;
		double dSumFit=0.0;
		LMA fitlma;
		int nBinSizeEst = 256;
		
		nBinSizeEst =getBinOptimalNumber(thImage);
		
		
		
		//mean, sd, min, max
		thImage.setHistogramSize(nBinSizeEst);
		imgstat = ImageStatistics.getStatistics(thImage, Measurements.MEAN+Measurements.STD_DEV+Measurements.MIN_MAX, null);
		dMean = imgstat.mean;
		nPeakPos = 0;
		
		nHistSize = imgstat.histogram.length;		
		dHistogram = new double [2][nHistSize];
		nMaxCount = 0;
	
		//determine position and height of maximum count in histogram (mode)
		//and height at maximum
		for(i =0; i<nHistSize; i++)
		{
			
			nCount=imgstat.histogram[i];
			dHistogram[0][i]=imgstat.min + i*imgstat.binSize;			
			dHistogram[1][i] = (double)nCount;
			dSum += (double)nCount;
			if(nMaxCount < nCount)
			{
				nMaxCount= nCount;
				dMean = imgstat.min + i*imgstat.binSize;
				nPeakPos = i;
			}			
		}
		//estimating width of a peak
		//going to the left
		i=nPeakPos;
		while (i>0 && imgstat.histogram[i]>0.5*nMaxCount)
		{
			i--;			
		}
		if(i<0)
			i=0;
		dLeftWidth = i;
		//going to the right
		i=nPeakPos;
		while (i<nHistSize && imgstat.histogram[i]>0.5*nMaxCount)
		{
			i++;			
		}
		if(i==nHistSize)
			i=nHistSize-1;
		dRightWidth = i;
		//FWHM in bins
		dWidth = (dRightWidth-dLeftWidth);
		dSD = dWidth*imgstat.binSize/2.35;
		//fitting range +/- 3*SD
		dLeftWidth = nPeakPos - 3*dWidth/2.35;
		if(dLeftWidth<0)
			dLeftWidth=0;
		dRightWidth = nPeakPos + 3*dWidth/2.35;
		if(dRightWidth>nHistSize)
			dRightWidth=nHistSize;
		nUpCount = (int)dRightWidth;
		nDownCount = (int)dLeftWidth;
		//preparing histogram range for fitting
		dNoiseFit = new double [2][nUpCount-nDownCount+1];
		for(i =nDownCount;i<=nUpCount;i++)
		{
			dNoiseFit[0][i-nDownCount] = dHistogram[0][i];
			dNoiseFit[1][i-nDownCount] = dHistogram[1][i];
			dSumFit +=dHistogram[1][i];
		}
		//fitting range is too small, less than 80%
		if((dSumFit/dSum)<0.8)
		{
			//build cumulative distribution
			dHistCum = new double [nHistSize];
			dSumFit = 0;
			for(i =0; i<nHistSize; i++)
			{
				dSumFit+=dHistogram[1][i];
				dHistCum[i] = dSumFit/dSum;	
			}
			nUpCount = 0;
			nDownCount = 0;
			for(i =0; i<nHistSize; i++)
			{
				//10% quantile
				if(dHistCum[i]<0.11)
					nDownCount= i;
				//90% quantile
				if(dHistCum[i]<0.91)
					nUpCount= i;			
			}
			//preparing histogram range for fitting
			dNoiseFit = new double [2][nUpCount-nDownCount+1];
			for(i =nDownCount;i<=nUpCount;i++)
			{
				dNoiseFit[0][i-nDownCount] = dHistogram[0][i];
				dNoiseFit[1][i-nDownCount] = dHistogram[1][i];			
			}
			
		}

		fitlma = new LMA(new SMLOneDGaussian(), new double[] {(double)nMaxCount, dMean, dSD}, dNoiseFit);
		fitlma.fit();

		dMean = fitlma.parameters[1];
		dSD = fitlma.parameters[2];
		
		dFitErrors = fitlma.getStandardErrorsOfParameters();
		// scaling coefficient for parameters errors estimation 
		// (Standard deviation of residuals)
		dErrCoeff = Math.sqrt(fitlma.chi2/(nUpCount-nDownCount+1-3));
		for(i =0;i<3;i++)
			dFitErrors[i] *= dErrCoeff;
		for(i =0;i<3;i++)
			dFitErrors[i] *= 100/fitlma.parameters[i]; 
		
		//if (dFitErrors[1]> 20 || dMean<imgstat.min || dMean> imgstat.max ||  dSD < imgstat.min || dSD> imgstat.max)
		if (dMean<imgstat.min || dMean> imgstat.max ||  dSD < imgstat.min || dSD> imgstat.max)
			//fit somehow failed
		{	results = new float [] {(float) imgstat.mean,(float) imgstat.stdDev};}
		else
		{	results = new float [] {(float) dMean,(float) dSD};}

		return results;
		
	}
	
	/** implements error function calculation with precision ~ 10^(-4) 
	  * used for calculation of Gaussian convolution kernel 
	  * **/
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
	
	
	/** function returns optimal bin number for the image histogram
	 *  according to the Freedman-Diaconis rule (check wiki) **/
	int getBinOptimalNumber(ImageProcessor ip)
	{
	
		int width, height;
		int pixelCount;
	
		//first let's calculate median
		//this part is taken from ImageJ source code
		//and adapted for the whole image
		
		width=ip.getWidth();
		height=ip.getHeight();
		pixelCount=width*height;
		

		float[] pixels2 = new float[pixelCount];

		System.arraycopy((float[])ip.getPixels(),0,pixels2,0,pixelCount);
	
		Arrays.sort(pixels2);
		//int middle = pixels2.length/2;
		int qi25 = Math.round(pixelCount*0.25f);
		int qi75 = Math.round(pixelCount*0.75f);

		float IQR = pixels2[qi75]-pixels2[qi25];
		double h= 2*IQR*Math.pow((double)pixelCount, -1.0/3.0);
		
		return (int)Math.round((pixels2[pixelCount-1]-pixels2[0])/h);
			
	}

	/** 
	 * function thresholds FloatProcessor and returns byte version
	 */
	ByteProcessor thresholdFloat(FloatProcessor ip, float dThreshold)
	{
		int width=ip.getWidth();
		int height= ip.getHeight();
		int i;
		ByteProcessor bp = new ByteProcessor(width, height);
		float[] flPixels = (float[])ip.getPixels();
//		byte[] btPixels = (byte[])bp.getPixels();
		
		for (int y=0; y<height; y++) 
		{
			i = y*width;
			for (int x=0; x<width; x++) 
			{
				if(flPixels[i]>=dThreshold)
				{
					//btPixels[i]=(byte)255;
					bp.set(x, y, 255);
				}
				else
				{
					//btPixels[i]=(byte)0;
					bp.set(x, y, 0);
				}				
				i++;
			}
		}
				
		return bp;
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
        final int pointInc = xDirection ? 1 : width;        //increment of the pixels array index to the next poin a line
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