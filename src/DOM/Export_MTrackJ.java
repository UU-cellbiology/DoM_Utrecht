package DOM;



import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import ij.IJ;
import ij.ImagePlus;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

public class Export_MTrackJ implements PlugIn {
	
	SMLAnalysis sml = new SMLAnalysis();
	SMLLinker smlLink;
	ImagePlus imp;	
	
	@Override
	public void run(String arg) {
		
		DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		DecimalFormat df3 = new DecimalFormat ("#.###", symbols);
		DecimalFormat df1 = new DecimalFormat ("#.0", symbols);
		String filename;
		String sFileWrite;
		int i, nCurrTrack, nAbsTracksNumber;


		IJ.register(Export_MTrackJ.class);
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.columnExists(20))
		{
			IJ.error("Not able to detect a valid 'Particles Table' with linked tracks information, please make one.");
			return;
		}
		//if the function is called from macro and filename already provided		
		boolean bCheck = !arg.isEmpty();	
		if(bCheck)
			filename = arg;
		//otherwise let's ask user for the filename
		else
		{
			//let's compose suggested filename
			imp = IJ.getImage();
			if(imp!=null)
			{
				filename = imp.getTitle() + "_tracks";
			}
			else
			{
				filename = "tracks";
			}
			
			SaveDialog sd = new SaveDialog("Save tracks in MTrackJ format", filename, ".mdf");
	        String path = sd.getDirectory();
	        if (path==null)
	        	return;
	
	        filename = path+sd.getFileName();
		}
		
		//sorting tracks to be sure in the Results table organization
		
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_ParticleID, true);
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_TrackID, true);
		//reading Results table to linker object 
		smlLink = new SMLLinker(sml);
		
		//exporting tracks		
        try {
			File file = new File(filename);
			FileWriter writer = new FileWriter(file);
			writer.write("MTrackJ 1.2.0 Data File\n");
			writer.write("Assembly 1\n");
			writer.write("Cluster 1\n");
			i=0;				
			
			//nAbsTracksNumber = smlLink.trackid[0];
			nCurrTrack = (int)smlLink.trackid[0];
			sFileWrite = "Track "+nCurrTrack + "\n";
			writer.write(sFileWrite);
			while(i<smlLink.nPatNumber)
			{
				
				nAbsTracksNumber = (int)smlLink.trackid[i];
				
				//new track
				if(nAbsTracksNumber != nCurrTrack)
				{
					nCurrTrack = nAbsTracksNumber;
					sFileWrite = "Track "+nCurrTrack + "\n";
					writer.write(sFileWrite);
				}
								
				sFileWrite = "Point "+ (int)smlLink.patid[i] + " " + df3.format(smlLink.x[i])+" "+df3.format(smlLink.y[i])+" 1.0 " + df1.format(smlLink.f[i])+ " 1.0\n";
				writer.write(sFileWrite);				
				i++;
			}
			
			writer.write("End of MTrackJ Data File\n");
			writer.close();
		} catch (IOException e) {	
			IJ.error(e.getMessage());
			//e.printStackTrace();
		}
        IJ.showStatus("Exporting tracks... Done.");
	}

}
