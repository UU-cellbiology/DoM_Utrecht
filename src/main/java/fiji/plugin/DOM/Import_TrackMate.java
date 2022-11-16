package fiji.plugin.DOM;



import java.io.FileInputStream;
import java.io.InputStream;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class Import_TrackMate implements PlugIn {
	
	
	SMLAnalysis sml_import = new SMLAnalysis();

	@Override
	public void run(String arg) {
		
		String sFilename;
		long nPatNumber=0;
		
		IJ.register(Import_TrackMate.class);
		sFilename=IJ.getFilePath("Open TrackMate XML file");
		if (sFilename==null)
	       	return;
		if(sFilename.endsWith(".xml"))
    	{
			GenericDialog dgImportTM = new GenericDialog("Import of Trackmate results");
			dgImportTM.addNumericField("Pixel size of original image", Prefs.get("SiMoLoc.TMImportdPixelSize", 1.0), 2, 4, "units");	
			dgImportTM.showDialog();
			
			if (dgImportTM.wasCanceled())
	            return;
			double dPxSize = dgImportTM.getNextNumber();
			Prefs.set("SiMoLoc.TMImportdPixelSize", dPxSize);
			
			IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
			IJ.log("Loading TrackMate file:"+sFilename);
			IJ.log("Original image pixel size:"+Double.toString(dPxSize));
			try
            {
				int nTrack = 0;
            	int nFrame;
				XMLInputFactory inputFactory = XMLInputFactory.newInstance();
            	InputStream in = new FileInputStream(sFilename);
            	XMLStreamReader streamReader = inputFactory.createXMLStreamReader(in);
            	int nPoint = 0;
            	while (streamReader.hasNext()) 
            	{
                    if (streamReader.isStartElement()) 
                    {
                        switch (streamReader.getLocalName()) 
                        {
                        case "Tracks": 
                        {
                        	IJ.log("Loading total of: " + streamReader.getAttributeValue(0) + " tracks... "); 
                        	sml_import.ptable.reset();

                            break;
                        }
                        case "particle": 
                        {
                        	nPoint =1;
                        	nTrack = nTrack + 1;
                            break;
                        }
                        case "detection": 
                        {
                        	nPatNumber++;
                        	nFrame = Integer.parseInt(streamReader.getAttributeValue(0))+1;
  						    sml_import.ptable.incrementCounter();
  						    sml_import.ptable.addValue("Abs_frame", nFrame);
  						    sml_import.ptable.addValue("X_(px)", Float.parseFloat(streamReader.getAttributeValue(1))/dPxSize);
  						    sml_import.ptable.addValue("Y_(px)", Float.parseFloat(streamReader.getAttributeValue(2))/dPxSize);
  						    sml_import.ptable.addValue("Channel", 1);
  						    sml_import.ptable.addValue("Slice", nFrame);
  						    sml_import.ptable.addValue("Frame", nFrame);
  						    sml_import.ptable.addValue("Track_ID", nTrack);
  						    sml_import.ptable.addValue("Particle_ID", nPoint);
                        	nPoint++;
                            break;
                        }
                        }
                    } 
                    streamReader.next();           
            	}
            }
            catch (Exception em)
            {
            	 IJ.log("Error reading "+sFilename);
            	 IJ.log(em.getMessage());
            }
			
			Import_MTrackJ.add_Tracks_Lengths(sml_import,Import_MTrackJ.calculate_Tracks_Lengths(sml_import));		
			sml_import.showTable();
			IJ.log("Done loading "+Long.toString(nPatNumber)+" particles.");
    	}
		else
		{
			IJ.error("Not a TrackMate XML file");
		}
	}

}
