package DOM;

import java.util.Arrays;
import java.util.Comparator;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

class tableComparator implements Comparator<double []> {
    private int columnToSortOn;
    private boolean ascending;
   
    //contructor to set the column to sort on.
    tableComparator(int columnToSortOn, boolean ascending) {
      this.columnToSortOn = columnToSortOn;
      this.ascending = ascending;
    }

    // Implement the abstract method which tells
    // how to order the two elements in the array.
    public int compare(double [] o1, double [] o2) {
        double[] row1 = (double[])o1;
        double[] row2 = (double[])o2;
        int res;
        if(row1[columnToSortOn]==row2[columnToSortOn])
        	return 0;
        if(row1[columnToSortOn]>row2[columnToSortOn])
        	res = 1;
        else
        	res = -1;
        if (ascending)
        	return res;
        else 
        	return (-1)*res;
               
        
    }
}

public class Sort_Results implements PlugIn {

	SMLAnalysis sml = new SMLAnalysis();
	
	@Override
	public void run(String arg) {
		
		int nColN=0;
		boolean bOrder;
		//asking user for sorting criteria
		GenericDialog dgSortParticles = new GenericDialog("Sort Particles");
		String [] SortColumn = new String [] {
				"Amplitude_fit","X_(px)","Y_(px)","False_positive","X_loc_error_(px)","Y_loc_error_(px)","IntegratedInt","SNR","Frame Number"};
		int Colindex;
		String [] SortOrder = new String [] {
				"Ascending","Descending"};
		int Sortindex;
		
		dgSortParticles.addChoice("Sort by column:", SortColumn, Prefs.get("SiMoLoc.SortColumn", "SNR"));
		dgSortParticles.addChoice("Sorting order:", SortOrder, Prefs.get("SiMoLoc.SortOrder", "Ascending"));
		
		dgSortParticles.showDialog();
		if (dgSortParticles.wasCanceled())
            return;
		Colindex = dgSortParticles.getNextChoiceIndex();
		Sortindex = dgSortParticles.getNextChoiceIndex();
		Prefs.set("SiMoLoc.SortColumn", SortColumn[Colindex]);
		Prefs.set("SiMoLoc.SortOrder", SortOrder[Sortindex]);

		switch (Colindex)
		{
			case 0:nColN=0; break;
			case 1:nColN=1; break;
			case 2:nColN=2; break;
			case 3:nColN=6; break;
			case 4:nColN=7; break;
			case 5:nColN=8; break;
			case 6:nColN=10; break;
			case 7:nColN=11; break;
			case 8:nColN=13; break;		
		}
		if (Sortindex == 0)
			bOrder = true;
		else 
			bOrder = false;
		sorting(nColN,bOrder);
		
	}
	
	public void sorting(int nSortingColumn, boolean ascending) 
	{
		int colNumber;
		int len;
		int i,j;
		
		colNumber = sml.ptable.getLastColumn()+1;

		double [] s = sml.ptable.getColumnAsDoubles(0);
		len = s.length;
		double [][] data = new double[len][colNumber];
		
		IJ.showStatus("Sorting Results Table: Preparation...");
		for (i=0; i<colNumber;i++)
		{
			s = sml.ptable.getColumnAsDoubles(i);
			for(j=0; j<len;j++)
			data[j][i]= s[j];
		}
		IJ.showStatus("Sorting Results Table: Sorting...");
		Arrays.sort(data, new tableComparator(nSortingColumn, ascending));
		
		IJ.showStatus("Sorting Results Table: Updating Table..."); 		
		for (i=0;i<colNumber;i++)
			for(j=0;j<len;j++)
				sml.ptable.setValue(i, j, data[j][i]);		
		sml.showTable();
		
	}
	
	//sorting function for external calls
	public static void sorting_external_silent(SMLAnalysis smlext, int nSortingColumn, boolean ascending) 
	{
		int colNumber;
		int len;
		int i,j;
		
		colNumber = smlext.ptable.getLastColumn()+1;

		double [] s = smlext.ptable.getColumnAsDoubles(0);
		len = s.length;
		double [][] data = new double[len][colNumber];
		
		//IJ.showStatus("Sorting Results Table: Preparation...");
		for (i=0; i<colNumber;i++)
		{
			s = smlext.ptable.getColumnAsDoubles(i);
			for(j=0; j<len;j++)
			data[j][i]= s[j];
		}
		//IJ.showStatus("Sorting Results Table: Sorting...");
		Arrays.sort(data, new tableComparator(nSortingColumn, ascending));
		
		//IJ.showStatus("Sorting Results Table: Updating Table..."); 		
		for (i=0;i<colNumber;i++)
			for(j=0;j<len;j++)
				smlext.ptable.setValue(i, j, data[j][i]);		
		//smlext.showTable();
		
	}
	
}
