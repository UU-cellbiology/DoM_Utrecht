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
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

// OpenCL
import static org.jocl.CL.*;
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
	
	//launch search of particles 
	public void run(String arg) {
		
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
		

		nFreeThread = -1;
		nSlice = 0;
		bContinue = true;
		if(dlg.nDetectionType==0 || dlg.nDetectionType==1 || dlg.nDetectionType==2)
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
					smlthreads[nFreeThread].init(ip, sml, dlg, nSlice, SpotsPositions,nStackSize, smlcount);
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
			
			//add width and height of original image to table
			sml.ptable.addValue("Original_image_size",imp.getWidth());
			sml.ptable.addValue("Original_image_size",imp.getHeight());
			
			sml.showTable();
		}
		//detection and fitting
		if(dlg.nDetectionType >= 1)
		{
			//creating new overlay
			SpotsPositions = new Overlay();
			//sorting results by frame
			IJ.showStatus("Sorting Results Table: Preparation...");
			Sort_Results.sorting_external_silent(sml, 0, true);
			//getting results from table
			double[] nFrameN = sml.ptable.getColumnAsDoubles(0);
			nDetPartNum = nFrameN.length;
			detparticles = new double [2][nDetPartNum];
			detparticles[0] = sml.ptable.getColumnAsDoubles(1);
			detparticles[1] = sml.ptable.getColumnAsDoubles(2);
			//detparticles[2] = sml.ptable.getColumnAsDoubles(3);
			//detparticles[3] = sml.ptable.getColumnAsDoubles(4);
			//erase particle table
			sml.ptable.reset();
						
			IJ.showStatus("Fitting intensities...");
			IJ.showProgress(0, nStackSize);
			
			//case of BaLM single molecule intensity analysis
			if(dlg.nDetectionType == 2)
			{
				//clearing all molecules that are not at their final bleach step
				detparticles = BaLMfilter(detparticles, nFrameN);
				nDetPartNum = detparticles[0].length;
				nFrameN = new double[nDetPartNum];
				nFrameN = detparticles[2];
			}
			nUniqFrames = uniqFrames(nFrameN, nDetPartNum);
			
			/*************************/
			/* CPU Fitting from here */
			/*************************/
			
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
							particles_= new double [2][nPositionEnd-nPositionStart];
							particles_[0] = Arrays.copyOfRange(detparticles[0], nPositionStart, nPositionEnd);
							particles_[1] = Arrays.copyOfRange(detparticles[1], nPositionStart, nPositionEnd);
							//particles_[2] = Arrays.copyOfRange(detparticles[2], nPositionStart, nPositionEnd);
							//particles_[3] = Arrays.copyOfRange(detparticles[3], nPositionStart, nPositionEnd);
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
			
			/**************************************************/
			/* CPU fitting till here, GPU fitting starts here */
			/**************************************************/
			
			else
			{
				// Fixed parameters for now
				boolean GPU_USE_DOUBLE_PRECISION = false; // NOTE: double precision not supported
				boolean GPU_AUTOMATIC_MODE_ENABLED = false; // NOTE: not yet implemented
				boolean GPU_PROFILING_MODE_ENABLED = false; // for debugging only
				boolean GPU_DEBUG_MODE_ENABLED = false; // for debugging only
				
				int GPU_MAX_ITERATIONS = 100;
				int GPU_BATCH_SIZE = 1024;
				int GPU_LWG_SIZE = 128;
				
				// create new GPU base object and load OpenCL kernel file
				GPUBase gpu = null;
				try
				{
					gpu = new GPUBase(GPU_USE_DOUBLE_PRECISION, GPU_AUTOMATIC_MODE_ENABLED, GPU_PROFILING_MODE_ENABLED, GPU_DEBUG_MODE_ENABLED);
				}
				catch(CLException cle)
				{
					IJ.error("Could not set up OpenCL environment. MAybe user cancelled device selection dialog, or the OpenCL is not supported on this system");
				}
				
				// open OpenCL kernel file
				BufferedReader opencl_kernel_reader = null;
				try
				{
					InputStream resource_stream = this.getClass().getClassLoader().getResourceAsStream("FittingKernelLMA");
					if(resource_stream == null)
					{
						IJ.error("Could not load OpenCL kernel as resource from JAR file");
					}
					opencl_kernel_reader = new BufferedReader(new InputStreamReader(resource_stream));
				}
				catch(Exception e)
				{
					IJ.error("Could not load OpenCL kernel");
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
				}
				
				// TODO: implement GPU fitting procedure
				IJ.error("TODO: provide GPU implementation");
				
				//detparticles[0][n] = x-coordinate of nth particle
				//detparticles[1][n] = y-coordinate of nth particle
				//nFrameN[n] = frame number of nth particle
				//imp = image plus (including stack)
				//imp.setSliceWithoutUpdate(nSlice+1); set current slice
				//ip = imp.getProcessor().duplicate(); copy ip of current slice
				//smlfitthreads[nFreeThread].init(ip, sml, dlg, particles_, nSlice, SpotsPositions, nStackSize, smlcount); cpu fitting procedure
				//this.sml.fitParticles(this.ip, this.dlg, this.particles, this.nFrame, this.SpotsPositions); called inside sml fit thread
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

	private double[][] BaLMfilter(double[][] detparticles_, double[] nFrameN_) 
	{
		int nDetPartNum, nCount, nPosBegin, nPosEnd, nPosOneBegin;
		int nFrameOne, nFrameTwo, nArrayOneL, nArrayTwoL;
		double[][] filteredparticles_;
		double[][] firstframe_;
		double[][] secondframe_;
		int[] marks;
		int i,j;
		boolean bFlag, bContinue;
		
		nPosBegin = 0;
		nPosEnd = 0;
		nCount = 0;
		nFrameOne = (int) nFrameN_[0];		
		nDetPartNum = nFrameN_.length;
		//whether particle is last bleaching step or not
		marks = new int [nDetPartNum];
		
		//filling first frame
		bFlag = true;
		while(bFlag)
		{
			nPosEnd++;
			if(nPosEnd==nDetPartNum )
				bFlag=false;
			else
				if((int)(nFrameN_[nPosEnd])!=nFrameOne)
					bFlag=false;
			
		}
		firstframe_ = new double[2][nPosEnd-nPosBegin];
		firstframe_[0] = Arrays.copyOfRange(detparticles_[0], nPosBegin, nPosEnd); //x coordinates
		firstframe_[1] = Arrays.copyOfRange(detparticles_[1], nPosBegin, nPosEnd); //y coordinates
		bContinue = true;
		while (bContinue)
		{
			if(nPosEnd == nDetPartNum)
				bContinue = false;
			else
			{
				//next frame does not contain any particles
				if((int)nFrameN_[nPosEnd]>nFrameOne+1)
				{
					//mark all particles as filtered
					for(i=nPosBegin; i<nPosEnd; i++)
					{
						marks[i]=1;
						nCount++;
					}
					nPosBegin=nPosEnd;
					nFrameOne = (int) nFrameN_[nPosBegin];	
					bFlag = true;
					while(bFlag)
					{
						nPosEnd++;
						if(nPosEnd==nDetPartNum )
							bFlag=false;
						else
							if((int)(nFrameN_[nPosEnd])!=nFrameOne)
								bFlag=false;
					}
					//fill frame one array
					firstframe_ = new double[2][nPosEnd-nPosBegin];
					firstframe_[0] = Arrays.copyOfRange(detparticles_[0], nPosBegin, nPosEnd); //x coordinates
					firstframe_[1] = Arrays.copyOfRange(detparticles_[1], nPosBegin, nPosEnd); //y coordinates
				}//if((int)nFrameN_[nPosEnd]>nFrameOne+1)
				//next frame does contains particles
				else
				{
					//let's store those particles
					nPosOneBegin = nPosBegin;
					nPosBegin=nPosEnd;
					nFrameTwo=nFrameOne+1;
					bFlag = true;
					while(bFlag)
					{
						nPosEnd++;
						if(nPosEnd==nDetPartNum )
							bFlag=false;
						else
							if((int)(nFrameN_[nPosEnd])!=nFrameTwo)
								bFlag=false;						
					}
					//fill frame two array
					secondframe_ = new double[2][nPosEnd-nPosBegin];
					secondframe_[0] = Arrays.copyOfRange(detparticles_[0], nPosBegin, nPosEnd); //x coordinates
					secondframe_[1] = Arrays.copyOfRange(detparticles_[1], nPosBegin, nPosEnd); //y coordinates
					nArrayOneL = firstframe_[0].length;
					nArrayTwoL = secondframe_[0].length;
					//looking for presence of particle at second frame
					for(i=0;i<nArrayOneL;i++)
					{
						bFlag = true;
						for(j=0;j<nArrayTwoL;j++)
							if(Math.sqrt(Math.pow(firstframe_[0][i]-secondframe_[0][j], 2)+Math.pow(firstframe_[1][i]-secondframe_[1][j], 2))<=dlg.dPSFsigma)
								bFlag = false;
						if (bFlag)
						{
							marks[nPosOneBegin+i]=1;
							nCount++;
						}
					}
					//done. let's switch arrays
					firstframe_ = new double[2][nArrayTwoL];
					firstframe_ = secondframe_;
					nFrameOne = nFrameTwo;
					
				}
			}
			
		}
		filteredparticles_ = new double[3][nCount];
		nCount = 0;
		for(i=0;i<nDetPartNum;i++)
		{
			if(marks[i]>0)
			{
				filteredparticles_[0][nCount] = detparticles_[0][i];
				filteredparticles_[1][nCount] = detparticles_[1][i];
				//filteredparticles_[2][nCount] = detparticles_[2][i];
				//filteredparticles_[3][nCount] = detparticles_[3][i];
				filteredparticles_[3][nCount] = nFrameN_[i];
				nCount++;
			}
		}
		return filteredparticles_;
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
	
	public void init(ImageProcessor ip, SMLAnalysis sml, SMLDialog dlg, int nFrame, Overlay SpotsPositions, int nStackSize, SMLProgressCount smlcount)
	{
		this.sml    = sml;
		this.ip     = ip;
		this.dlg    = dlg;
		this.nFrame = nFrame;
		this.SpotsPositions = SpotsPositions;
		this.nStackSize = nStackSize;
		this.smlcount = smlcount;
	}
	
	public void run()
	{
		this.sml.detectParticles(this.ip, this.dlg, this.nFrame, this.SpotsPositions);	
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


