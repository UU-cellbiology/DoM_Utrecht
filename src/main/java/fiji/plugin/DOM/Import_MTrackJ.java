package fiji.plugin.DOM;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

public class Import_MTrackJ implements PlugIn {

	double [] trackid;
	double [] patid;
	//double [] tracklength;
	//long nPatNumber;

	SMLAnalysis smlImport = new SMLAnalysis();
	@Override
	public void run(String arg) {
				
		String filename;
		long nPatNumber;
		IJ.register(Import_MTrackJ.class);
		
		//if the function is called from macro and filename already provided		
		boolean bCheck = !arg.isEmpty();	
		if(bCheck)
			filename = arg;
		//otherwise let's ask user for the filename
		else
		{
			OpenDialog dgLoadMTrackJ = new OpenDialog("Load MtrackJ mdf file","", "*.mdf");
			
	        String path = dgLoadMTrackJ.getDirectory();
	        if (path==null)
	        	return;
	
	        filename = path+dgLoadMTrackJ.getFileName();
	        IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
			IJ.log("Loading MTrackJ file: "+filename);
	        smlImport.ptable.reset();
	        nPatNumber = 0;
	        int nTrackN = 0;
	        int nParticleN = 0;
			try {
				
		        BufferedReader br = new BufferedReader(new FileReader(filename));
		        String line;
				while ((line = br.readLine()) != null) 
					{
					   // process the line.
					  String[] line_array = line.split(" ");
					  if(line_array[0].equals("Track"))
					  {
						  nTrackN++;
					      nParticleN = 0;
					  }
					  if(line_array[0].equals("Point"))
					  {						  
						  nParticleN++;
						  nPatNumber++;
						  smlImport.ptable.incrementCounter();
						  smlImport.ptable.addValue("Abs_frame", (int)Float.parseFloat(line_array[5]));
						  smlImport.ptable.addValue("X_(px)", Float.parseFloat(line_array[2]));
						  smlImport.ptable.addValue("Y_(px)", Float.parseFloat(line_array[3]));
						  smlImport.ptable.addValue("Channel", (int)Float.parseFloat(line_array[6]));
						  smlImport.ptable.addValue("Slice", (int)Float.parseFloat(line_array[5]));
						  smlImport.ptable.addValue("Frame", (int)Float.parseFloat(line_array[5]));
						  smlImport.ptable.addValue("Track_ID", nTrackN);
						  smlImport.ptable.addValue("Particle_ID", nParticleN);
					  }
					}
	
		        br.close();
			}
			//catching errors in file opening
			catch (FileNotFoundException e) {
				IJ.error(""+e);
			}	        
			catch (IOException e) {
				IJ.error(""+e);
			}	        
			
			
			//calculate_Tracks_Lengths();
			add_Tracks_Lengths(smlImport,calculate_Tracks_Lengths(smlImport));	
			smlImport.showTable();
			IJ.log("Done loading "+Long.toString(nPatNumber)+" particles.");
		}

	}


	//function calculating the number of detected particles per track
	public static double [] calculate_Tracks_Lengths(final SMLAnalysis sml_import)
	{
		int nCount,i;
		int nCurrTrack, nMaxVal, nIniPosition;
		long nPatNumber;
		
		double [] trackid   = sml_import.ptable.getColumnAsDoubles(6);		
		double [] patid     = sml_import.ptable.getColumnAsDoubles(7);
		nPatNumber = patid.length;
		
		double [] tracklength = new double [(int) nPatNumber];
		

		nMaxVal = 0;
		nCurrTrack = (int) trackid[0];
		nIniPosition = 0;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			if(trackid[nCount]==nCurrTrack)
			{
				if(patid[nCount]>nMaxVal)
					nMaxVal =(int)patid[nCount]; 
			}
			else
			{
				//storing value of length at array
				for(i=nIniPosition;i<nCount;i++)
				{
					tracklength[i] = nMaxVal;
				}
				nMaxVal = (int) patid[nCount];
				nIniPosition = nCount;
				nCurrTrack = (int) trackid[nCount]; 
			}
			//end of the list
			if(nCount == nPatNumber - 1)
			{
				for(i=nIniPosition;i<nPatNumber;i++)
				{
					tracklength[i] = nMaxVal;
				}
			}
			
		}
		return tracklength;
	
	}
	//function adding the number of detected particles per track to results table
	public static void add_Tracks_Lengths(final SMLAnalysis sml_import, final double [] tracklength)
	{
		int nCount;
		//adding to final table
		//init new column
		sml_import.ptable.setValue("Track_Length", 0, tracklength[0]);
		
		//add all stuff
		for(nCount = 1; nCount<tracklength.length; nCount++)
			sml_import.ptable.setValue(8, nCount, tracklength[nCount]);
		
	}


}
