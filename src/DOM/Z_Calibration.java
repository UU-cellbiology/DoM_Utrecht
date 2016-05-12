package DOM;
//in detect molecules on line 329
//add width and height of original image to the results table
//sml.ptable.addValue("Original_image_size",imp.getWidth());
//sml.ptable.addValue("Original_image_size",imp.getHeight());


import java.awt.AWTEvent;
import java.awt.Button;
import java.awt.Color;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.AWTEventListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
//import ij.gui.ImageCanvas;
//import ij.gui.NonBlockingGenericDialog;
//import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;
import ij.plugin.frame.Fitter;

public class Z_Calibration implements PlugIn{

	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	int nCols = 9, maxTracklength, index;
	
	double [] falsepos, xloc, yloc, sdx, sdy, trackid, particleid, tl,  tracklength, sortedTracklength, sdDif;
	double[][] table = new double[nCols][], longestTrack, longestTrackT, temp;
	
	double maxLocError = 0.5;
	
	ImagePlus plotIp;
	
	public void run(String arg) 
	{
		IJ.register(Z_Calibration.class);
	
		if (!dlg.zCalibration()) return;
		//String [] zcUseOptions = new String [] {"Image stack","Particle table","Polynomial coefficients"};
	//1.6 compartibility issue!!!  EUGENE COMMENT
		/*switch(dlg.sZcUse){
			case "Image stack": zCal_imageStack(); break;
			case "Particle table": zCal_particleTable(); break;
			case "Polynomial coefficients": zCal_polCoef(); break;
		}*/
	}
		
	public void zCal_imageStack(){
		//msg("zCal_imageStack: still have to implement this");
		
		NonBlockingGenericDialog dlg = new NonBlockingGenericDialog("title");
		Panel p = new Panel();
		p.add(new Button("po"));
		
		dlg.add(p);
		dlg.pack();
		dlg.showDialog();
		
		
	}
	
	public void zCal_particleTable(){
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.columnExists(13))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		//check that the particles are linked
		if ( !sml.ptable.columnExists(18)){
			IJ.error("Particle table is probably not linked.");
			return;
		}
		
		falsepos = sml.ptable.getColumnAsDoubles(6);
		xloc = sml.ptable.getColumnAsDoubles(7);
		yloc = sml.ptable.getColumnAsDoubles(8);
		sdx = sml.ptable.getColumnAsDoubles(15);
		sdy = sml.ptable.getColumnAsDoubles(16);
		trackid = sml.ptable.getColumnAsDoubles(18);
		particleid = sml.ptable.getColumnAsDoubles(19);
		tracklength = sml.ptable.getColumnAsDoubles(20);
		
		sortedTracklength = tracklength.clone();
		Arrays.sort(sortedTracklength);
		maxTracklength = (int)sortedTracklength[sortedTracklength.length-1];
		
		table[0] = falsepos;
		table[1] = xloc;
		table[2] = yloc;
		table[3] = sdx;
		table[4] = sdy;
		table[5] = trackid;
		table[6] = particleid;
		table[7] = tracklength;
		table[8] = new double[table[0].length];
		//insert relative z-values (based on particleid and distance between z planes)
		for (int i=0; i<table[0].length; i++) {
			  table[8][i] = table[6][i] * dlg.zCalDistBetweenPlanes;
		}
		table = transpose(table);

		//sort table base on tracklength, then trackid, then particleid
		Arrays.sort(table, new Comparator<double[]>() {
		    public int compare(double[] a, double[] b) {
		    	if(Double.compare(a[7], b[7]) != 0){
		    		return Double.compare(a[7], b[7]);
		    	}
		    	else if(Double.compare(a[5], b[5]) != 0){
		    		return Double.compare(a[5], b[5]);
		    	}
		    	else{
		    		return Double.compare(a[6], b[6]);
		    	}
		    }
		});
		
		longestTrack = new double[maxTracklength][];//=> rows are rows of original particle table
		System.arraycopy(table,(table.length - maxTracklength),longestTrack,0,maxTracklength);
		longestTrackT = transpose(longestTrack); //=> rows are columns of original particle table
		
		Plot p = new Plot("Width and height of PSF","z-position (nm)","size of PSF (px)");
		//plot sdx
		p.addPoints(longestTrackT[8], longestTrackT[3], Plot.CONNECTED_CIRCLES);
		//plot sdy
		p.setColor(Color.red);
		p.addPoints(longestTrackT[8], longestTrackT[4], Plot.CONNECTED_CIRCLES);
		p.setColor(Color.black);
		p.setLegend("width	height", Plot.AUTO_POSITION);
		ImagePlus ip = p.getImagePlus();
		ip.show();
		
		int approxZero = dlgSetZero();
		ip.close();
		if(approxZero < 0)
			return;
		
		for(index = 0; index < maxTracklength; index++){
			if(longestTrackT[8][index] > approxZero)
				break;
		}
		
		int zRange = 400; //nm
		int nFrames = (int)(zRange/dlg.zCalDistBetweenPlanes)/2;
		temp = new double[nFrames*2+1][];
		
		System.arraycopy(longestTrack,index - nFrames,temp,0,2*nFrames+1);
		
		//find z-offset
		double dif, minDif = 1000, zOffset = 0;
		for(index = 0; index<temp.length; index++){
			dif = Math.abs(temp[index][3]-temp[index][4]);
			if(dif<minDif){
				zOffset = temp[index][8];
				minDif = dif;
			}
		}
		
		//subtract z-offset from all z-values
		for(index = 0; index < maxTracklength; index++){
			longestTrackT[8][index] = longestTrackT[8][index] - zOffset;  
		}
		
		longestTrack = transpose(longestTrackT);
		
		sdDif = new double[maxTracklength];
		for(int i = 0; i<maxTracklength;i++){
			sdDif[i] = longestTrackT[3][i] - longestTrackT[4][i];
		}
		longestTrack = transpose(longestTrackT);
		
		p = new Plot("Difference in width and height","z-position (nm)","width-height (px)");
		
		p.addPoints(longestTrackT[8], sdDif, Plot.CONNECTED_CIRCLES);
		
		ip = p.getImagePlus();
		ip.show();
		
		//fitPars = fitRangeMin, fitRangeMax, polDegree
		int[] fitPars = dlgFitRange();
		ip.close();
		if(fitPars[0] == 0 & fitPars[1] == 0)
			return;
		
		int indexMin = findIndex(longestTrackT[8],fitPars[0]);
		int indexMax = findIndex(longestTrackT[8],fitPars[1]);
		
		
		CurveFitter cf = polyFit(Arrays.copyOfRange(sdDif,indexMin,indexMax),Arrays.copyOfRange(longestTrackT[8],indexMin,indexMax),dlg.fitPolynomialDegree);
		//get plot window (for closing later)
		ImagePlus impFitPlot = IJ.getImage();
		
		double[] fitCoef = cf.getParams();
		//for(int i=0;i<fitCoef.length;i++){
		//	msg(fitCoef[i]);
		//}

		if(dlgStoreCal()){
			Prefs.set("SiMoLOc.ZC_polDegree", dlg.fitPolynomialDegree);
			switch(dlg.fitPolynomialDegree){
				case 1: Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", 0);
						Prefs.set("SiMoLOc.ZC_polCoef3", 0);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;
				case 2: Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", fitCoef[2]);
						Prefs.set("SiMoLOc.ZC_polCoef3", 0);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;
				case 3: Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", fitCoef[2]);
						Prefs.set("SiMoLOc.ZC_polCoef3", fitCoef[3]);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;
			}
			msg("Calibration stored!");
		}
		
		impFitPlot.close();
		//msg(fitRange[0]+" "+fitRange[1]);
	}
	
	Plot plot;
	
	public void zCal_polCoef(){
		plot = new Plot("title","x-label","y-label");
		plot.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, new double[]{1,2,3,4,5,6,7,8,9,10}, 2);
		
		final ImagePlus imp = new ImagePlus();//plot.getImagePlus();
		plot.setImagePlus(imp);
		plot.draw();
		
		final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("title");
		nb.addCheckbox("clabel", true);
		//nb.addImage(imp);
		ImageCanvas ic = new ImageCanvas(imp);
		final Panel jp = new Panel();
		jp.add(ic);
		nb.add(jp);
		
		Button b = new Button("test");
		b.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent arg0) {
				Random rand = new Random();
				double[] y = new double[10];
				 for ( int j = 0 ; j < y.length ; j++ ) { 
		              y[j] = rand.nextInt(101);
		           }
				plot = new Plot("title","x-label","y-label");
				plot.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, y, 2);				
				plot.setImagePlus(imp);
				plot.draw();
				for(int i=0;i<jp.getComponentCount();i++){
					if(jp.getComponent(i) instanceof ImageCanvas){
						jp.remove(i);
					}
				}
				
				ImageCanvas ic = new ImageCanvas(imp);
				
				jp.add(ic);
				
				nb.pack();
				
				///WHERE IS THIS METHOD??? EUGENE
				//nb.revalidate();
				
				
			}
			
		});
		nb.add(b);
		
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length == 0){
			IJ.log("no listeneres yet");
			Toolkit.getDefaultToolkit().addAWTEventListener(new AWTEventListener(){

				@Override
				public void eventDispatched(AWTEvent arg0) {
					// TODO Auto-generated method stub
					MouseEvent me = (MouseEvent) arg0;
					if(arg0.getSource() instanceof ImageCanvas){
						if(me.getID() == MouseEvent.MOUSE_CLICKED){
							
							IJ.log("x = "+String.valueOf(plot.descaleX(me.getX()))+" y = "+String.valueOf(plot.descaleY(me.getY())));
						}
					}
					
				}
			
			}
			,AWTEvent.MOUSE_EVENT_MASK);

		}
		else{
			IJ.log("listerner already in place");
		}
		
		
		nb.showDialog();
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length != 0){
			AWTEventListener l = Toolkit.getDefaultToolkit().getAWTEventListeners()[0];
			Toolkit.getDefaultToolkit().removeAWTEventListener(l);
			IJ.log("Listener removed");
		}
	}
	
	
	class nbdlgZCal extends NonBlockingGenericDialog{
		
		String title = "Z-calibration NonBlocking Dialog";
		
		
		public nbdlgZCal(String title) {
			super(title);
		}

		NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Z-calibration NonBlocking Dialog");
		
		public void toggleFitOptions(boolean on){
			if(on){
				
			}
		}
		
		
	}
	
	public CurveFitter polyFit(double[]x, double[]y, int degree){
		
		CurveFitter cf = new CurveFitter(x, y);
		int polDegree;
		
		switch(degree){
			case 1: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
			case 2: polDegree = ij.measure.Calibration.POLY2;break;
			case 3: polDegree = ij.measure.Calibration.POLY3;break;
			default: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
		}
		cf.doFit(polDegree); 
		Fitter.plot(cf);
		
		IJ.log(cf.getResultString());
		return cf;
	}
	
	
	public int findIndex(double[] array, int value){
		for(int i=0; i<array.length; i++){
			if((int)array[i] == value)
				return i;
		}
		return -1;
	}
	
	
	public int dlgSetZero() {
		GenericDialog dlg = new GenericDialog("Type approximate zero position");
		dlg.addNumericField("Zero at approximate: ", 0, 0);
		
		
		dlg.setResizable(false);
		dlg.showDialog();
		if (dlg.wasCanceled())
            return -1;
		
		return (int)dlg.getNextNumber();
	}
	
	public int[] dlgFitRange() {
		GenericDialog dlg = new GenericDialog("Fit Range");
		dlg.addNumericField("Fit range min (nm): ", -400, 0);
		dlg.addNumericField("Fit range max (nm): ", 400, 0);
		
		dlg.setResizable(false);
		dlg.showDialog();
		if (dlg.wasCanceled())
            return new int[]{0,0};
		
		return new int[]{(int)dlg.getNextNumber(),(int)dlg.getNextNumber()};
	}
	
	public boolean dlgStoreCal() {
		GenericDialog dlg = new GenericDialog("Store calibration?");
		dlg.addMessage("Store the new calibration?");

		dlg.setResizable(false);
		dlg.showDialog();
		return dlg.wasOKed();
	}
	
	double[][] transpose (double[][] array) {
	  if (array == null || array.length == 0)//empty or unset array, nothing do to here
	    return array;

	  int width = array.length;
	  int height = array[0].length;

	  double[][] array_new = new double[height][width];

	  for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
	      array_new[y][x] = array[x][y];
	    }
	  }
	  return array_new;
	}

	void msg(String message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(message);
		gd.showDialog();
	}
	void msg(double message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(String.valueOf(message));
		gd.showDialog();
	}
	void msg(int message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(String.valueOf(message));
		gd.showDialog();
	}

}
