package fiji.plugin.DOM;

import ij.IJ;
import ij.plugin.PlugIn;

public class LoadDriftCorrection implements PlugIn
{
	SMLAnalysis sml = new SMLAnalysis();
	/** x coordinates of particles in nm*/
	double [] x;
	/** y coordinates of particles in nm*/
	double [] y;
	/** frame numbers of particle*/
	double [] f;

	@Override
	public void run(String arg) {
		
		String sCorrection;
		String [] sCorrSplittedRows;
		String [] sCorrSplittedVals;
		String delimsn = "[\n]+";
		
		String delimst = "[\t]+";
		String delimstab = "[\t]+";
		String delimscomma = ",";
		int nTotalFramesTable;
		double [] xdrift;
		double [] ydrift;
		int [] fdrift;
		int i, nParticlesCount,nFrameN;
		double pxscale;
		int nMaxFrameNTable,nMaxFrameNDrift, nMaxFrame;
		boolean bStop;
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		sCorrection = IJ.openAsString("");
		sCorrSplittedRows = sCorrection.split(delimsn);
		if(sCorrSplittedRows[0].contains(","))
			{delimst=delimscomma;}
		else
			{delimst=delimstab;}
		
		sCorrSplittedVals =sCorrSplittedRows[0].split(delimst);

		if(!(sCorrSplittedVals[0].equals("Frame_Number")&&sCorrSplittedVals[1].equals("X_drift_(nm)")))
		{
			IJ.error("Not able to detect a valid 'Drift Correction' results table, try another file.");
			return;
		}
		//loading values
		nTotalFramesTable = sCorrSplittedRows.length-1;
		xdrift = new double [nTotalFramesTable];
		ydrift = new double [nTotalFramesTable];
		fdrift = new int [nTotalFramesTable];

		for (i=1;i<(nTotalFramesTable+1);i++)
		{
			sCorrSplittedVals =sCorrSplittedRows[i].split(delimst);
			fdrift[i-1]=(int)Double.parseDouble(sCorrSplittedVals[0]);
			xdrift[i-1]=Double.parseDouble(sCorrSplittedVals[1]);
			ydrift[i-1]=Double.parseDouble(sCorrSplittedVals[2]);
		}
	
		// let's sort result table by frame
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_FrameN, true);
		sml.showTable();
		//get frame number
		f   = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		//get coordinates
		x   = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);		
		y   = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
		
		nParticlesCount=f.length;
		nMaxFrameNDrift=fdrift[fdrift.length-1];
		nMaxFrameNTable = (int)f[nParticlesCount-1];
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		if(nMaxFrameNDrift>nMaxFrameNTable)
		{
			nMaxFrame=nMaxFrameNTable;
			
			IJ.log(" WARNING! Drift correction table contains more frames than Results Table! Truncated correction applied!");
		}
		if(nMaxFrameNDrift<nMaxFrameNTable)
		{
			nMaxFrame=nMaxFrameNDrift;
			IJ.log(" WARNING! Results Table contains more frames than Drift correction table! Truncated correction applied!");
		}		
		else
		{ nMaxFrame=nMaxFrameNDrift;  }
		
		//let's restore scale
		pxscale=  sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0);
		//lock table
		sml.ptable_lock.lock();
		i=0;
		bStop=false;
		nFrameN=(int)f[i]-1;
		while(!bStop)
		{	
			x[i]+=xdrift[nFrameN];
			y[i]+=ydrift[nFrameN];
			sml.ptable.setValue(DOMConstants.Col_Xnm, i, x[i]);
			sml.ptable.setValue(DOMConstants.Col_Ynm, i, y[i]);
			sml.ptable.setValue(DOMConstants.Col_X, i, x[i]*pxscale);
			sml.ptable.setValue(DOMConstants.Col_Y, i, y[i]*pxscale);
			
			i++;
			if(i==nParticlesCount)
				bStop=true;
			else
			{
				nFrameN=(int)f[i]-1;
				
				if(nFrameN+1>nMaxFrame)
					bStop=true;
			}
		}
		
		sml.ptable_lock.unlock();
		IJ.log("Drift correction loaded and applied.");
		
		sml.showTable();
	}

}
