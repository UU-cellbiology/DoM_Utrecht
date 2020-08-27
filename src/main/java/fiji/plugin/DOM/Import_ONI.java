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

public class Import_ONI implements PlugIn {

	SMLAnalysis sml_import = new SMLAnalysis();
	int nChannel;
	public ResultsTable ptable = Analyzer.getResultsTable();
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();
	
	@Override
	public void run(String arg) {
		
		
		String filename;	
		String separator;
		String myLine;
		String[] values;
		String[] headings;
		int nDoMcolsN;
		int nFileChannel;
		//String[] headingsDoM;
		int ncols;
		int i;
		long filesize;
		long nStringCount=0;
		long progressbytes=0;
		int iFrame=1,iChannel=0,iXnm=2,iYnm=3,iZnm=4,iXpx=7,iYpx=8,iXprnm=5,iYprnm=6, iPhotons=12,iBG=13;
		
		
		IJ.register(Import_ONI.class);
		// TODO Auto-generated method stub
		filename=IJ.getFilePath("Open ONI results table");
		if (filename==null)
	       	return;
		
		GenericDialog dgImportONI = new GenericDialog("Import ONI format");
		String [] chChannel = new String [] {
				"0","1"};
		dgImportONI.addChoice("Load channel:", chChannel, Prefs.get("SiMoLoc.ImportONI", "0"));
		dgImportONI.showDialog();
		if (dgImportONI.wasCanceled())
            return;
		
		nChannel = dgImportONI.getNextChoiceIndex();
		Prefs.set("SiMoLoc.ImportONI", chChannel[nChannel]);

		
		
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Loading ONI results table:"+filename);
		IJ.log("Loading channel number:"+Integer.toString(nChannel));
		
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
		for (i=0;i<headings.length;i++)
		{
			if(headings[i].compareTo("Channel")==0)
				iChannel=i;
			if(headings[i].compareTo("Frame")==0)
				iFrame=i;
			if(headings[i].compareTo("X (nm)")==0)
				iXnm=i;
			if(headings[i].compareTo("Y (nm)")==0)
				iYnm=i;
			if(headings[i].compareTo("Z (nm)")==0)
				iZnm=i;
			if(headings[i].compareTo("X (pix)")==0)
				iXpx=i;
			if(headings[i].compareTo("Y (pix)")==0)
				iYpx=i;
			if(headings[i].compareTo("X precision (nm)")==0)
				iXprnm=i;
			if(headings[i].compareTo("Y precision (nm)")==0)
				iYprnm=i;
			if(headings[i].compareTo("Photons")==0)
				iPhotons=i;
			if(headings[i].compareTo("Background")==0)
				iBG=i;
			
		}
		String [] headingsDoM={"X_(px)","Y_(px)","Frame_Number","X_(nm)","X_loc_error(nm)","Y_(nm)", "Y_loc_error(nm)", 
				"Z_(nm)","Z_loc_error(nm)","Amplitude_fit", "Amp_error", "BGfit", "BGfit_error","SD_X_(nm)","SD_X_error(nm)",
					"SD_Y_(nm)","SD_Y_error(nm)","False_positive","IntegratedInt","SNR", "R2_fit","Iterations_fit"};
		
		nDoMcolsN=headingsDoM.length;
		
		ptable.setPrecision(5);
		IJ.showStatus("Reading ONI Results table file...");
		ptable.reset(); // erase Results table
		ptable_lock.lock();
		double [] DoMvalues = new double [nDoMcolsN];
		//read the rest of the file and load it to results table
		while ( scanner.hasNext())
		{    
			
			myLine=scanner.nextLine();
			//read a line from ONI
		    values = myLine.split(separator);
		    
		    nFileChannel=Integer.parseInt(values[iChannel]);
		    if (nFileChannel==nChannel)
		    {		    
			    //assign corresponding values
			    DoMvalues[DOMConstants.Col_FrameN]=Double.parseDouble(values[iFrame])+1;
			    DoMvalues[DOMConstants.Col_Xnm]=Double.parseDouble(values[iXnm]);
			    DoMvalues[DOMConstants.Col_Ynm]=Double.parseDouble(values[iYnm]);
			    DoMvalues[DOMConstants.Col_Znm]=Double.parseDouble(values[iZnm]);
			    DoMvalues[DOMConstants.Col_loc_errX]=Double.parseDouble(values[iXprnm]);
			    DoMvalues[DOMConstants.Col_loc_errY]=Double.parseDouble(values[iYprnm]);
			    DoMvalues[DOMConstants.Col_IntegrInt]=Double.parseDouble(values[iPhotons]);
			    DoMvalues[DOMConstants.Col_BGfit]=Double.parseDouble(values[iBG]);
			    DoMvalues[DOMConstants.Col_X]=Double.parseDouble(values[iXpx]);
			    DoMvalues[DOMConstants.Col_Y]=Double.parseDouble(values[iYpx]);
			    //channel is stored as Col_BGfit_error
			    DoMvalues[DOMConstants.Col_BGfit_error]=(double)nChannel;
			    //put to table
			    ptable.incrementCounter();
			    for(i=0;i<nDoMcolsN;i++)
			    {
			    	ptable.addValue(headingsDoM[i], DoMvalues[i]);
			    }
			    nStringCount++;
		    }
		  
		    progressbytes+=myLine.getBytes().length;
		    IJ.showProgress((double)progressbytes/(double)filesize);
		    
		    
		}
		ptable_lock.unlock();
		IJ.showProgress(0.75);
		IJ.showStatus("Updating Results table window...Please wait...");
		ptable.show("Results");		
		IJ.showStatus("Loading ONI Results table...done.");
		IJ.showProgress(1.1);
		IJ.log(Long.toString(nStringCount)+" rows loaded.");
	
	}

}
