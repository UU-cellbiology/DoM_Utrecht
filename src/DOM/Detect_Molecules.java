package DOM;

// Java
import java.util.Arrays;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.IOException;

// ImageJ
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import ij.gui.YesNoCancelDialog; // There is no Yes/No dialog?

// OpenCL
//import static org.jocl.CL.*;
import org.jocl.*;

/**
 *
 */
public class Detect_Molecules implements PlugIn {
	
	ImagePlus imp;
	ImageProcessor ip;
	Overlay SpotsPositions;
	
	SMLThread [] smlthreads;
	SMLFitThread [] smlfitthreads;
	
	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLProgressCount smlcount = new SMLProgressCount(0);
	//SMLProgressCount smlfitcount = new SMLProgressCount(0);
	int nDetPartNum, nCountThread;
	int nNumberofDetectedColumns;
	
	Roi RoiActive;
	
	//launch search of particles 
	public void run(String arg) {
		
		long startTime;
		long detectionTime=0;
		long fullTime;
		// ... the code being measured ...    
		
		
		int nFreeThread = -1;
		int nSlice = 0;
		int nPositionStart,nPositionEnd;
		boolean bContinue = true;
		boolean bListEnd;
		int nStackSize;
		//double[] nFrameN;
		double[][] detparticles;
		double[][] particles_;
		int nUniqFrames;
		
		IJ.register(Detect_Molecules.class);
		
		
		particles_= new double [1][1];
		//checking whether there is any open images
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
		
		//creating new overlay
		SpotsPositions = new Overlay();
		
		if (!dlg.findParticles()) return;
		
		//generating gaussian mexican-hat convolution kernel with size of PSF
		sml.initConvKernel(dlg);
		//putting limiting criteria on spot size
		//dlg.nAreaCut = dlg.nKernelSize + 1;
		dlg.nAreaCut = (int) (dlg.dPSFsigma * dlg.dPSFsigma);
		
		nStackSize = imp.getStackSize();

		
		// using only necessary number of threads	 
		if(dlg.nThreads < nStackSize) 
			smlthreads = new SMLThread[dlg.nThreads];
		else
			smlthreads = new SMLThread[nStackSize];
		
		//getting active Roi
		RoiActive = imp.getRoi();
		
		//let's start measuring time		
		startTime = System.nanoTime();
		
		//let's log stuff
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Image title: \"" + imp.getTitle() + "\"");
		IJ.log("Detection SD of PSF: " + String.format("%.2f",dlg.dPSFsigma) + " pixels");
		IJ.log("Image pixel size: " + String.format("%.2f",dlg.dPixelSize) + " nm");
		IJ.log("Kernel size: (detection) " + String.format("%d",dlg.nKernelSize) + " pixels, (fitting) "+String.format("%d",(2*(int)(3.0*dlg.dPSFsigma))+1)+" pixels");
		nFreeThread = -1;
		nSlice = 0;
		bContinue = true;
		if(dlg.nDetectionType==0 || dlg.nDetectionType==1)
		{
			sml.ptable.reset(); // erase particle table
			IJ.showStatus("Detecting molecules...");		
			while (bContinue)
			{
				//check whether reached the end of stack
				if (nSlice >= nStackSize) bContinue = false;
				else
				{
					imp.setSliceWithoutUpdate(nSlice+1);
					ip = imp.getProcessor().duplicate();
				}
				
				if (bContinue)
				{
					//filling free threads in the beginning
					if (nSlice < smlthreads.length)
						nFreeThread = nSlice;				
					else
					//looking for available free thread
					{
						nFreeThread = -1;
						while (nFreeThread == -1)
						{
							for (int t=0; t < smlthreads.length; t++)
							{
								if (!smlthreads[t].isAlive())
								{
									nFreeThread = t;
									break;
								}
							}
							if (nFreeThread == -1)
							{
								try
								{
									Thread.currentThread();
									Thread.sleep(1);
								}
								catch(Exception e)
								{
										IJ.error(""+e);
								}
							}
						}
					}
					smlthreads[nFreeThread] = new SMLThread();
					smlthreads[nFreeThread].init(ip, sml, dlg, nSlice, SpotsPositions,nStackSize, smlcount, RoiActive);
					smlthreads[nFreeThread].start();
					//IJ.showProgress(nSlice, nStackSize);
				} //end of if (bContinue)
				nSlice++;
			} // end of while (bContinue)
			
			for (int t=0; t<smlthreads.length;t++)
			{
				try
				{
					smlthreads[t].join();				
				}
				catch(Exception e)
				{
					IJ.error(""+e);
				}
			}
			//detection only
			IJ.showStatus("Detection completed.");
			detectionTime = System.nanoTime() - startTime;
			IJ.log("Detection time:" + String.format("%.2f",((double)Math.abs(detectionTime))*0.000000001) + " s");
		}

		if(dlg.nDetectionType==0)
		{
			//in case of single image, update ROIs
			if(nStackSize == 1 && dlg.bShowParticles)	
			{
				int nRois = SpotsPositions.size();
				for(int i = 0;i<nRois;i++)
				{
					SpotsPositions.get(i).setPosition(0);
				}			
			}
			if(dlg.bShowParticles)
			{
				imp.setOverlay(SpotsPositions);
				imp.updateAndRepaintWindow();
				imp.show();
			}			
			

			
			sml.showTable();
		}
		//detection and fitting
		if(dlg.nDetectionType >= 1)
		{
			//creating new overlay
			SpotsPositions = new Overlay();
			//sorting results by frame
			IJ.showStatus("Sorting Results Table: Preparation...");

			//check if there is non-empty table
			nNumberofDetectedColumns = sml.ptable.getLastColumn()+1; 
			if(nNumberofDetectedColumns==0)
			{
				 IJ.error("There is no valid Results table with detections!");
				 return;
			}
			Sort_Results.sorting_external_silent(sml, 0, true);
			//getting results from table
			double[] nFrameN = sml.ptable.getColumnAsDoubles(0);
			nDetPartNum = nFrameN.length;
			//imported from MTrackJ
			if(nNumberofDetectedColumns==9)
			{
				detparticles = new double [5][nDetPartNum];
			}
			else
			{
				detparticles = new double [2][nDetPartNum];
			}
			detparticles[0] = sml.ptable.getColumnAsDoubles(1);
			detparticles[1] = sml.ptable.getColumnAsDoubles(2);
			if(nNumberofDetectedColumns==9)
			{			
				detparticles[2] = sml.ptable.getColumnAsDoubles(6);
				detparticles[3] = sml.ptable.getColumnAsDoubles(7);
				detparticles[4] = sml.ptable.getColumnAsDoubles(8);
			}
			//erase particle table
			sml.ptable.reset();
						
			IJ.showStatus("Fitting intensities...");
			IJ.showProgress(0, nStackSize);			
			
			nUniqFrames = uniqFrames(nFrameN, nDetPartNum);
			
			/*************************/
			/* GPU Fitting from here */
			/*************************/
			
			if(dlg.bUseGPUAcceleration)
			{
				// Fixed parameters for now
				boolean GPU_USE_DOUBLE_PRECISION = false; // NOTE: double precision not supported
				boolean GPU_AUTOMATIC_MODE_ENABLED = false; // NOTE: not yet implemented
				boolean GPU_PROFILING_MODE_ENABLED = false; // for debugging only
				boolean GPU_DEBUG_MODE_ENABLED = true; // for debugging only
				
				// Fixed parameters in ParallelLMA.java
				//int GPU_MAX_ITERATIONS = 100;
				//int GPU_BATCH_SIZE = 1024;
				//int GPU_LWG_SIZE = 128;
				//
				//double GPU_DEFAULT_X_SIGMA = 2.0f;
				//double GPU_DEFAULT_Y_SIGMA = 2.0f;
				//double GPU_DEFAULT_LAMBDA = 0.001f;
				
				// create new GPU base object and load OpenCL kernel file
				GPUBase gpu = null;
				try
				{
					gpu = new GPUBase(GPU_USE_DOUBLE_PRECISION, GPU_AUTOMATIC_MODE_ENABLED, GPU_PROFILING_MODE_ENABLED, GPU_DEBUG_MODE_ENABLED);
				}
				catch(CLException cle)
				{
					YesNoCancelDialog continue_diag = new YesNoCancelDialog(null, "Problem setting up OpenCL context for GPU acceleration", "Could not set up an OpenCL context on your system. Either\nno suitable device was found for GPU acceleration or\nthe device selection dialog was canceled.\n\nWould you like to continue fitting on the CPU instead?");
					if(continue_diag.yesPressed())
					{
						dlg.bUseGPUAcceleration = false;
					}
					else
					{
						return;
					}
				}
				catch(Throwable e) // catch everything throwable, include errors and exceptions
				{
					YesNoCancelDialog continue_diag = new YesNoCancelDialog(null, "GPU acceleration not supported", "OpenCL is not supported on your computer. Make sure\nyou installed the appropriate OpenCL drivers for your device\nand installed the jar-libraries bundled with this plug-in\ninto the right folder.\n\nWould you like to continue fitting on the CPU instead?");
					if(continue_diag.yesPressed())
					{
						dlg.bUseGPUAcceleration = false;
					}
					else
					{
						return;
					}
				}
				
				// continue only if gpu was constructed successfully
				if(gpu != null)
				{
					// open OpenCL kernel file
					BufferedReader opencl_kernel_reader = null;
					try
					{
						InputStream resource_stream = gpu.getClass().getResourceAsStream("/FittingKernelLMA.cl");
						if(resource_stream == null)
						{
							IJ.error("Could not load OpenCL kernel as resource from JAR file");
							return;
						}
						opencl_kernel_reader = new BufferedReader(new InputStreamReader(resource_stream));
					}
					catch(Exception e)
					{
						IJ.error("Could not load OpenCL kernel");
						return;
					}
					
					// load contents of OpenCL kernel file
					String opencl_kernel_program = "";
					String line = null;
					try
					{
						while((line = opencl_kernel_reader.readLine()) != null)
						{
							// add line to program
							opencl_kernel_program += line + '\n';
						}
					}
					catch(IOException e)
					{
						// early exit
						IJ.error("Could not load contents of OpenCL kernel file");
						return;
					}
					
					// close file
					try
					{
						opencl_kernel_reader.close();
					}
					catch(IOException e)
					{
						// ignore, do nothing
					}
					
					// load OpenCL kernel program
					boolean compile_success = gpu.loadProgramFromString(opencl_kernel_program);
					if(!compile_success)
					{
						IJ.error("Could not load program file into OpenCL");
						return;
					}
					
					// fit spots on GPU
					sml.ptable.reset();
					ParallelLMA.run(gpu, dlg.bUseMLE, dlg.nBatchSize, dlg.nGroupSize, dlg.nIterations, imp, detparticles, nFrameN, dlg.dPSFsigma, dlg.dPixelSize, dlg.bIgnoreFP, sml.ptable);
					
					//detparticles[0][n] = x-coordinate of nth particle
					//detparticles[1][n] = y-coordinate of nth particle
					//nFrameN[n] = frame number of nth particle
					//imp = image plus (including stack)
					//imp.setSliceWithoutUpdate(nSlice+1); set current slice
					//ip = imp.getProcessor().duplicate(); copy ip of current slice
					//smlfitthreads[nFreeThread].init(ip, sml, dlg, particles_, nSlice, SpotsPositions, nStackSize, smlcount); cpu fitting procedure
					//this.sml.fitParticles(this.ip, this.dlg, this.particles, this.nFrame, this.SpotsPositions); called inside sml fit thread
					// particles(_) is subset of detected_particles for current frame nFrame
				}
			}
			
			/**************************************************/
			/* GPU fitting till here, CPU fitting starts here */
			/**************************************************/
			
			if(!dlg.bUseGPUAcceleration)
			{
				/// using only necessary number of threads
				if(dlg.nThreads < nUniqFrames)
					smlfitthreads = new SMLFitThread[dlg.nThreads];
				else
					smlfitthreads = new SMLFitThread[nUniqFrames];
				
				SMLProgressCount smlcount = new SMLProgressCount(0);
				nFreeThread = -1;
				nSlice = 0;
				bContinue = true;
				nPositionStart = 0;
				nPositionEnd = 0;
				nCountThread = -1;
				while (bContinue)
				{
					//check whether reached the end of stack
					if (nSlice >= nStackSize) bContinue = false;
					else
					{
						nPositionStart = nPositionEnd;
						//check whether reached the end of particles table
						if(nPositionStart<nDetPartNum)
						{
							//no particles in this frame
							if((int)(nFrameN[nPositionStart])>(nSlice+1))
							{
								smlcount.SMLProgressCountIncreaseValue((int)(nFrameN[nPositionStart]-nSlice-1));
								IJ.showProgress(smlcount.nSliceLeft-1, nStackSize);
								nSlice = (int) (nFrameN[nPositionStart]-1);
							}
							//finding all particles in current frame
							bListEnd = false;
							while (!bListEnd)
							{
								nPositionEnd++;
								if(nPositionEnd == nDetPartNum)
									bListEnd = true;
								else
									if((int)(nFrameN[nPositionEnd])!=nSlice+1)
										bListEnd = true;
							}
							//making array for them
							if(nNumberofDetectedColumns==9)
								particles_= new double [5][nPositionEnd-nPositionStart];
							else
								particles_= new double [2][nPositionEnd-nPositionStart];
														
							particles_[0] = Arrays.copyOfRange(detparticles[0], nPositionStart, nPositionEnd);
							particles_[1] = Arrays.copyOfRange(detparticles[1], nPositionStart, nPositionEnd);
							if(nNumberofDetectedColumns==9)
							{
								particles_[2] = Arrays.copyOfRange(detparticles[2], nPositionStart, nPositionEnd);
								particles_[3] = Arrays.copyOfRange(detparticles[3], nPositionStart, nPositionEnd);
								particles_[4] = Arrays.copyOfRange(detparticles[4], nPositionStart, nPositionEnd);
							}
							//copying image
							imp.setSliceWithoutUpdate(nSlice+1);
							ip = imp.getProcessor().duplicate();
							nCountThread++;
							
						}
						//reached the end of particles table
						else
							bContinue = false;
						
						
					}
					
					if (bContinue)
					{
						//filling free threads in the beginning
						if (nCountThread < smlfitthreads.length)
							nFreeThread = nCountThread;
						else
							//looking for available free thread
						{
							nFreeThread = -1;
							while (nFreeThread == -1)
							{
								for (int t=0; t < smlfitthreads.length; t++)
								{
									if (!smlfitthreads[t].isAlive())
									{
										nFreeThread = t;
										break;
									}
								}
								if (nFreeThread == -1)
								{
									try
									{
										Thread.currentThread();
										Thread.sleep(1);
									}
									catch(Exception e)
									{
										IJ.error(""+e);
									}
								}
							}
						}
						smlfitthreads[nFreeThread] = new SMLFitThread();
						smlfitthreads[nFreeThread].init(ip, sml, dlg, particles_, nSlice, SpotsPositions, nStackSize, smlcount);
						smlfitthreads[nFreeThread].start();
						
					} //end of if (bContinue)
					nSlice++;
				} // end of while (bContinue)
				
				for (int t=0; t<smlfitthreads.length;t++)
				{
					try
					{
						smlfitthreads[t].join();
					}
					catch(Exception e)
					{
						IJ.error(""+e);
					}
				}
			}


			//in case of single image, update ROIs
			if(nStackSize == 1 && dlg.bShowParticles)	
			{
				int nRois = SpotsPositions.size();
				for(int i = 0;i<nRois;i++)
				{
					SpotsPositions.get(i).setPosition(0);
				}			
			}
			if(dlg.bShowParticles)
			{
				imp.setOverlay(SpotsPositions);
				imp.updateAndRepaintWindow();
				imp.show();
			}	
			sml.showTable();
			IJ.showProgress(nStackSize, nStackSize);
			IJ.showStatus("Fitting completed.");
			fullTime = System.nanoTime() - startTime;
			if (dlg.nDetectionType==2)				
				IJ.log("Fitting time: " + String.format("%.2f",((double)Math.abs(fullTime))*0.000000001)+ " s");
			else
			{
				IJ.log("Fitting time: " + String.format("%.2f", ((double)Math.abs(fullTime-detectionTime))*0.000000001)+ " s");
				IJ.log("Total time:" + String.format("%.2f",((double)Math.abs(fullTime))*0.000000001) + " s");
			}

		}

	}

	private int uniqFrames(double[] nFrameN, int nLength) {
		int i, nCount=1;
		int nCurrFrame;
		boolean bContinue = true;
		
		i=0;
		nCurrFrame = (int)nFrameN[i];		
		while (bContinue)
		{
			i++;
			if(i==nLength)
				bContinue = false;
			else
			{
				if((int)nFrameN[i]!=nCurrFrame)
				{
					nCurrFrame=(int)nFrameN[i];
					nCount++;
					
				}
			}
			
		}
		
		return nCount;
	}

}

class SMLThread extends Thread 
{
	private Overlay SpotsPositions;
	private ImageProcessor ip;
	private SMLDialog dlg;
	private int nFrame;
	private SMLAnalysis sml;
	private SMLProgressCount smlcount;
	private int nStackSize;
	private Roi RoiActive;
	
	public void init(ImageProcessor ip, SMLAnalysis sml, SMLDialog dlg, int nFrame, Overlay SpotsPositions, int nStackSize, SMLProgressCount smlcount, Roi RoiActive)
	{
		this.sml    = sml;
		this.ip     = ip;
		this.dlg    = dlg;
		this.nFrame = nFrame;
		this.SpotsPositions = SpotsPositions;
		this.nStackSize = nStackSize;
		this.smlcount = smlcount;
		this.RoiActive = RoiActive;
	}
	
	public void run()
	{
		this.sml.detectParticles(this.ip, this.dlg, this.nFrame, this.SpotsPositions, this.RoiActive);	
		smlcount.SMLProgressCountIncrease();
		IJ.showProgress(smlcount.nSliceLeft-1, nStackSize);
	}
}

class SMLFitThread extends Thread 
{
	private Overlay SpotsPositions;
	private ImageProcessor ip;
	private SMLDialog dlg;
	private double[][] particles;
	private int nFrame;
	private SMLAnalysis sml;
	private SMLProgressCount smlcount;
	private int nStackSize;
	
	public void init(ImageProcessor ip, SMLAnalysis sml, SMLDialog dlg, double[][] particles, int nFrame, Overlay SpotsPositions, int nStackSize, SMLProgressCount smlcount)
	{
		this.sml    = sml;
		this.ip     = ip;
		this.dlg    = dlg;
		this.nFrame = nFrame;
		this.particles = particles;
		this.SpotsPositions = SpotsPositions;
		this.nStackSize = nStackSize;
		this.smlcount = smlcount;		
	}
	
	public void run()
	{
		this.sml.fitParticles(this.ip, this.dlg, this.particles, this.nFrame, this.SpotsPositions);	
		smlcount.SMLProgressCountIncrease();
		IJ.showProgress(smlcount.nSliceLeft-1, nStackSize);
	}
	
}


class SMLProgressCount
{
	 
	public int nSliceLeft;
	public SMLProgressCount (int ini_value)
	{
		nSliceLeft = ini_value;
	}
	public void SMLProgressCountIncrease()
	{
		synchronized (this)
		{ 
			nSliceLeft++;
		}
		
	}
	public void SMLProgressCountIncreaseValue(int inc_value)
	{
		synchronized (this)
		{ 
			nSliceLeft+=inc_value;
		}
		
	}
	
}


