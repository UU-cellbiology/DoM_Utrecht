package DOM;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;

public class LoadLargeRT implements PlugIn {


	protected ImagePlus image;
	public ResultsTable ptable = Analyzer.getResultsTable();
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();



	@Override
	public void run(String arg0) 
	{
		
		String filename;		
		String separator;
		String myLine;
		String[] values;
		String[] headings;
		int ncols;
		int i;
		long filesize;
		long nStringCount=0;
		long progressbytes=0;
		
		
		//get filename
		
		filename=IJ.getFilePath("Open large Result table");
		
		if (filename==null)
	       	return;
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Loading results table:"+filename);
		
		//open file for parsing
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
		

		//reduce the size of the table by putting less precision
		ptable.setPrecision(2);
		IJ.showStatus("Reading large Results table file...");
		ptable.reset(); // erase Results table
		ptable_lock.lock();
		//read the rest of the file and load it to results table
		while ( scanner.hasNext())
		{    
			
			myLine=scanner.nextLine();
	
		    values = myLine.split(separator);
		    ptable.incrementCounter();									
			
		    for(i=1;i<ncols;i++)
		    {
		    	ptable.addValue(headings[i], Double.parseDouble(values[i]));
		    }
		    //IJ.log(values[2]);
		    nStringCount++;
		    progressbytes+=myLine.getBytes().length;
		    IJ.showProgress((double)progressbytes/(double)filesize);
		    
		}
		ptable_lock.unlock();
		IJ.showProgress(0.75);
		IJ.showStatus("Updating Results table window...Please wait...");
		ptable.show("Results");		
		IJ.showStatus("Loading large Results table...done.");
		IJ.showProgress(1.1);
		IJ.log(Long.toString(nStringCount)+" rows loaded.");
	}

	
}
