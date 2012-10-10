package DOM;


import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class Detect_Molecules implements PlugIn {
	
	ImagePlus imp;
	ImageProcessor ip;
	Overlay SpotsPositions;
	
	SMLThread [] smlthreads;
	
	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	SMLProgressCount smlcount = new SMLProgressCount(0);
	
	//launch search of particles 
	public void run(String arg) {
		
		int nFreeThread = -1;
		int nSlice = 0;
		boolean bContinue = true;
		int nStackSize;
		
		IJ.register(Detect_Molecules.class);
		
		sml.ptable.reset(); // erase particle table
		
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
		
		IJ.showStatus("Detecting molecules...");
		nFreeThread = -1;
		nSlice = 0;
		bContinue = true;
		
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
				smlthreads[nFreeThread].init(ip, sml, dlg, nSlice, SpotsPositions,nStackSize,smlcount);
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
		IJ.showStatus("Detection completed.");
		sml.showTable();

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

	
}

