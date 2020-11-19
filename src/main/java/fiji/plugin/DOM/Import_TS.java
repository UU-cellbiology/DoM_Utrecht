package fiji.plugin.DOM;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;

public class Import_TS implements PlugIn
{
	SMLAnalysis sml_import = new SMLAnalysis();
	double dPxSize;
	public ResultsTable ptable = Analyzer.getResultsTable();
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();
	
	@Override
	public void run(String arg) 
	{
		// TODO Auto-generated method stub
		
		String filename;	
		String separator;
		String myLine;
		String[] values;
		String[] headings;
		int ncols;
		int i;
		long filesize;
		int nDoMcolsN;
		long nStringCount=0;
		long progressbytes=0;
		int iFrame=1000,iXnm=1000,iYnm=1000,iZnm=1000,iSigma=1000, iSigmaX=1000,iSigmaY=1000,iPhotons=1000,iBG=1000,iBGSD=1000,iErrXY=1000, iErrXY_Z=1000,iErrZ=1000, iChi=1000;
		
		IJ.register(Import_TS.class);
		

		/** whether z values are calculated and present in the table **/
		boolean zValsPresent=false;
		
		
		// TODO Auto-generated method stub
		filename=IJ.getFilePath("Open ThunderSTORM CSV results table");
		if (filename==null)
	       	return;
		
		GenericDialog dgImportTS = new GenericDialog("Import of ThunderSTORM results");
		dgImportTS.addNumericField("Camera's pixel size", Prefs.get("SiMoLoc.TSImportdPixelSize", 65), 2, 4, "nm");	
		dgImportTS.showDialog();
		
		if (dgImportTS.wasCanceled())
            return;
		dPxSize = dgImportTS.getNextNumber();
		Prefs.set("SiMoLoc.TSImportdPixelSize", dPxSize);
		
		
		
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Loading ThunderSTORM results table:"+filename);
		
		File file = new File(filename);
		filesize=file.length();
		Scanner scanner=null;
		try {
			scanner = new Scanner(file);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		//determine separator	
		myLine=scanner.nextLine();	
		values = myLine.split(",");			
		int i_comma=values.length;
		values = myLine.split("\t");			
		int i_tab=values.length;
		
		if(i_comma>i_tab)
		{
			separator=",";
			ncols=i_comma;
		}
		else
		{
			separator="\t";
			ncols=i_tab;
		}
		
		//read header column names
		headings=myLine.split(separator);
		progressbytes+=myLine.getBytes().length;
		//find corresponding column titles (in case it is in random order)
		//y[nm] can be confused with uncertainty [nm], therefore reverse order
		for (i=(headings.length-1);i>=0;i--)
		{

			if(headings[i].contains("frame"))
				iFrame=i;
			if(headings[i].contains("x [nm]"))
				iXnm=i;
			if(headings[i].contains("y [nm]"))
				iYnm=i;
			if(headings[i].contains("z [nm]"))
				iZnm=i;
			if(headings[i].contains("sigma [nm]"))
				iSigma=i;
			if(headings[i].contains("sigma1 [nm]"))
				iSigmaX=i;
			if(headings[i].contains("sigma2 [nm]"))
				iSigmaY=i;
			if(headings[i].contains("intensity [photon]"))
				iPhotons=i;
			if(headings[i].contains("offset [photon]"))
				iBG=i;
			if(headings[i].contains("bkgstd [photon]"))
				iBGSD=i;
			if(headings[i].contains("chi2"))
				iChi=i;
			if(headings[i].contains("uncertainty_xy [nm]"))
				iErrXY_Z=i;
			if(headings[i].contains("uncertainty [nm]"))
				iErrXY=i;
			if(headings[i].contains("uncertainty_z [nm]"))
				iErrZ=i;
			
		}				
		
		if(iZnm<100)
		{
			zValsPresent=true;
			iErrXY=iErrXY_Z;
		}
		
		String [] headingsDoM={"X_(px)","Y_(px)","Frame_Number","X_(nm)","X_loc_error(nm)","Y_(nm)", "Y_loc_error(nm)", 
				"Z_(nm)","Z_loc_error(nm)","Amplitude_fit", "Amp_error", "BGfit", "BGfit_error","SD_X_(nm)","SD_X_error(nm)",
					"SD_Y_(nm)","SD_Y_error(nm)","False_positive","IntegratedInt","SNR", "R2_fit","Iterations_fit"};
		
		nDoMcolsN=headingsDoM.length;
		
		ptable.setPrecision(5);
		IJ.showStatus("Reading ThunderSTORM Results table file...");
		ptable.reset(); // erase Results table
		ptable_lock.lock();
		double [] DoMvalues = new double [nDoMcolsN];

		while ( scanner.hasNext())
		{    

			myLine=scanner.nextLine();
			//read a line from TS file
		    values = myLine.split(separator);
		    if(values.length==0)
		    	break;
		    
		    //IJ.log(myLine);
		    
		    //assign corresponding values
		    DoMvalues[DOMConstants.Col_FrameN]=Double.parseDouble(values[iFrame]);
		    DoMvalues[DOMConstants.Col_Xnm]=Double.parseDouble(values[iXnm]);
		    DoMvalues[DOMConstants.Col_Ynm]=Double.parseDouble(values[iYnm]);
		    
		    DoMvalues[DOMConstants.Col_X]=Double.parseDouble(values[iXnm])/dPxSize;
		    DoMvalues[DOMConstants.Col_Y]=Double.parseDouble(values[iYnm])/dPxSize;
		    
		    //put straight values. in TS export I will ask where data comes from
		    DoMvalues[DOMConstants.Col_loc_errX]=Double.parseDouble(values[iErrXY]);
		    DoMvalues[DOMConstants.Col_loc_errY]=Double.parseDouble(values[iErrXY]);
		    
		    if(zValsPresent)
		    {
			    DoMvalues[DOMConstants.Col_Znm]=Double.parseDouble(values[iZnm]);
			    DoMvalues[DOMConstants.Col_loc_errZ]=Double.parseDouble(values[iErrZ]);
			    DoMvalues[DOMConstants.Col_SD_X]=Double.parseDouble(values[iSigmaX]);
			    DoMvalues[DOMConstants.Col_SD_Y]=Double.parseDouble(values[iSigmaY]);
		    }
		    else
		    {
			    DoMvalues[DOMConstants.Col_SD_X]=Double.parseDouble(values[iSigma]);
			    DoMvalues[DOMConstants.Col_SD_Y]=Double.parseDouble(values[iSigma]);		    	
		    }
		    //chi square is not always present somehow
		    if(iChi<1000)
		    {
		    	DoMvalues[DOMConstants.Col_chi]=Double.parseDouble(values[iChi]);
		    }

		    DoMvalues[DOMConstants.Col_IntegrInt]=Double.parseDouble(values[iPhotons]);
		    DoMvalues[DOMConstants.Col_BGfit]=Double.parseDouble(values[iBG]);
		    
		    //save background STD as fitting error. Ask at TS export if it was the case
		    DoMvalues[DOMConstants.Col_BGfit_error]=Double.parseDouble(values[iBGSD]);
		   
		    
		    
		    //put to table
		    ptable.incrementCounter();
		    for(i=0;i<nDoMcolsN;i++)
		    {
		    	ptable.addValue(headingsDoM[i], DoMvalues[i]);
		    }
		    nStringCount++;

		  
		    progressbytes+=myLine.getBytes().length;
		    IJ.showProgress((double)progressbytes/(double)filesize);
		    
		    
		}
		ptable_lock.unlock();
		IJ.showProgress(0.75);
		IJ.showStatus("Updating Results table window...Please wait...");
		ptable.show("Results");		
		IJ.showStatus("Loading ThunderSTORM Results table...done.");
		IJ.showProgress(1.1);
		IJ.log(Long.toString(nStringCount)+" rows loaded.");	
	}
}
