package DOM;

import java.awt.Color;
import java.util.ArrayList;

import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;

public class SMLLinker {

	ImagePlus implnk;
	SMLDialog settings;
	SMLAnalysis smllink;
	double [] x;
	double [] y;	
	double [] f;
	double [] fp;
	double [] trackid;
	double [] patid;
	double [] nLength;
	long nPatNumber;
	int nUniqFrames;
	Overlay ovTracks;
	double dFPThreshold;
	
	ArrayList<double[]> TrackCoords;
	
	
	SMLLinker( SMLAnalysis sml_, SMLDialog dlg_, Overlay ovTracks_, ImagePlus imp_)
	{
		settings = dlg_;
		smllink = sml_;
		ovTracks = ovTracks_;
		implnk = imp_;
			
		//frame number
		x   = smllink.ptable.getColumnAsDoubles(1);		
		y   = smllink.ptable.getColumnAsDoubles(2);
		//coordinates
		f   = smllink.ptable.getColumnAsDoubles(13);
		//false positive mark
		fp  = smllink.ptable.getColumnAsDoubles(6);
		
		//total particles number
		nPatNumber = f.length;
		//number of unique frames
		nUniqFrames = uniqFrames();
		trackid = new double [(int) nPatNumber];
		patid = new double [(int) nPatNumber];
		//threshold for particles linking (include false positives or not)
		switch (dlg_.nLinkFP)
		{ //only true positives
			case 0:
				dFPThreshold = 0.3;
				break;
		 //true and half positives
			case 1:
			dFPThreshold = 0.7;
				break;
		//all particles
			case 2:
				dFPThreshold = 1.1;
				break;
			default:
				dFPThreshold = 0.3;
				break;
				
		}
		
	}
	
	void Link_Closest()
	{
		
		double [][] x_un;
		double [][] y_un;
		double [][] fp_un;
		
		int [][] framesstat;
		int nCurrFrame, nCount, nFrameCount, nCurrentParticle, nCompareFrame, nCompareParticle,nFrameCompareCount;
		double dDistance;
		int i;
		int nCurrentTrack;
		int nCurrentParticleinTrack;
		boolean bContinue;
		int nLinksNumber;
		int nPenalty, nCompareAbsCount;
		double x_compare, y_compare;
		
		TrackCoords = new ArrayList<double[]>();
	
		
		//frame number[0] and number of particles in each frame [1] and cumulative number of particles in list till this frame [2] 
		framesstat = new int [nUniqFrames][3];
		//calculate number of particles in each frame and
		//store it for each frame in framesstat array
		nCurrFrame = (int) f[0];
		nCount = 0;
		nFrameCount = 0;
		framesstat[nFrameCount][0]=nCurrFrame;
		for (i=0;i<nPatNumber;i++)
		{
			if((int)f[i]==nCurrFrame)
				nCount++;
			else
			{	
				//number of particles at this frame
				framesstat[nFrameCount][1]=nCount;
				nCurrFrame = (int)f[i];
				//cumulative number of particles
				if(nFrameCount==0)
					framesstat[nFrameCount][2]=nCount;
				else
					framesstat[nFrameCount][2]=nCount+framesstat[nFrameCount-1][2];
				nFrameCount++;
				framesstat[nFrameCount][0]=nCurrFrame;
				nCount = 1;	
			}
			
		}
		//last element
		framesstat[nFrameCount][1]=nCount;
		framesstat[nFrameCount][0]=nCurrFrame;
		
		//resorting data to new arrays for convenience
		
		//1) allocating arrays with uneven rows 
		//  for particles from each frame separately
		x_un   = new double [nUniqFrames][];
		y_un   = new double [nUniqFrames][];
		fp_un  = new double [nUniqFrames][];
				
		for(i=0;i<nUniqFrames;i++)
		{
			x_un[i]   = new double[framesstat[i][1]];
			y_un[i]   = new double[framesstat[i][1]];
			fp_un[i]  = new double[framesstat[i][1]];
		}
		//2) putting all data to those arrays
		nFrameCount = 0;
		nCurrentParticle = -1;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			nCurrentParticle++;
			if(nCurrentParticle == framesstat[nFrameCount][1])
			{
				nCurrentParticle=0;
				nFrameCount++;		
			}			
			x_un[nFrameCount][nCurrentParticle]   = x[nCount];
			y_un[nFrameCount][nCurrentParticle]   = y[nCount];
			fp_un[nFrameCount][nCurrentParticle]  = fp[nCount];					
		}
		
		
		
		
		
		//all right!
		//now let's link particles
		nFrameCount = 0;
		nCurrentParticle = -1;
		nCurrentTrack = 0;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{
			nCurrentParticle++;
			if(nCurrentParticle == framesstat[nFrameCount][1])
			{
				nCurrentParticle=0;
				nFrameCount++;
			}
					
			//checking whether particle was already linked
			if(fp_un[nFrameCount][nCurrentParticle]<dFPThreshold && nFrameCount<(nUniqFrames-1) && trackid[nCount]<1) 
			{ 	
				nCurrentParticleinTrack = 1;
				nCurrentTrack++;
				trackid[nCount]=nCurrentTrack;
				patid[nCount]=nCurrentParticleinTrack;
				nCurrFrame = framesstat[nFrameCount][0];
				TrackCoords.clear();
				x_compare = x[nCount];
				y_compare = y[nCount];
				TrackCoords.add(new double[] {x_compare,y_compare});
				
				bContinue = true;
				nFrameCompareCount = nFrameCount+1;
				nCompareFrame = framesstat[nFrameCompareCount][0];
				nCompareParticle = -1;
				nLinksNumber = 0;
				nPenalty = 0;
				while(bContinue) 
				{
					nCompareParticle++;
					if(nCompareParticle == framesstat[nFrameCompareCount][1]) //end of particles list in current frame
					{
						if (nFrameCompareCount == nUniqFrames - 1)	//reached the end of total particle list, stop						
							{bContinue = false;}
						else
						{//jumping to next frame
						 // but let's check whether we found something while scanning previous one?
							//no, we didn't let's increase penalty
							if(nCompareFrame-nCurrFrame>nLinksNumber+nPenalty)
							{
								nPenalty++;
								if(nPenalty>settings.nLinkFrameGap)
								{
									bContinue = false;
								}
								else
								{									
									nFrameCompareCount++;								
									nCompareFrame = framesstat[nFrameCompareCount][0];
									nCompareParticle = 0;										
								}
							}
							//yes, we did, let's jump then
							else
							{
								nLinksNumber = nPenalty+nLinksNumber;
								nPenalty = 0;
								nFrameCompareCount++;									
								nCompareFrame = framesstat[nFrameCompareCount][0];
								nCompareParticle = 0;
							}													
						}
							
					}
					//absolute number of particle in the Particle Table
					nCompareAbsCount = framesstat[nFrameCompareCount-1][2]+nCompareParticle;
					//ok, let's see if all checks are passed
					if(bContinue && fp_un[nFrameCompareCount][nCompareParticle]<dFPThreshold && trackid[nCompareAbsCount]<1)
					{//seems so. let's measure distances between molecules then
						
						dDistance = Math.sqrt(Math.pow(x_compare-x_un[nFrameCompareCount][nCompareParticle], 2) +Math.pow(y_compare-y_un[nFrameCompareCount][nCompareParticle], 2)); 
						//they are closer then one pixel
						if(dDistance<settings.dLinkDistance)
						{
							//found new link!
							//let's store values							
							nLinksNumber++;
							//mark it as already used 
							trackid[nCompareAbsCount] = nCurrentTrack;
							nCurrentParticleinTrack++;
							patid[nCompareAbsCount] = nCurrentParticleinTrack;
							if(settings.nLinkTrace != 0) //change reference position
							{
								x_compare = x_un[nFrameCompareCount][nCompareParticle];
								y_compare = y_un[nFrameCompareCount][nCompareParticle];
							}
														
							TrackCoords.add(new double[] {x_un[nFrameCompareCount][nCompareParticle],y_un[nFrameCompareCount][nCompareParticle]});							

							//exit from current frame (a bit lame one)
							nCompareParticle = framesstat[nFrameCompareCount][1] -1;
						}
					}
				}//end of "while" cycle looking for next linkers
				
			}//end of "checking whether particle was already linked" condition
		}//end of particles check		
		//last frame
		//pick up leftovers
		for(nCurrentParticle=0;nCurrentParticle<framesstat[nUniqFrames-1][1];nCurrentParticle++)
		{
			nCompareAbsCount = framesstat[nUniqFrames-2][2]+nCurrentParticle;
			if (trackid[nCompareAbsCount]<1 && fp_un[nUniqFrames-1][nCurrentParticle]<dFPThreshold)
			{
				nCurrentTrack++;
				trackid[nCompareAbsCount]=nCurrentTrack;
				patid[nCompareAbsCount]=1;
			}
		}
		
		smllink.ptable.setValue("Track_ID", 0, trackid[0]);
		smllink.ptable.setValue("Particle_ID", 0, patid[0]);
		//adding linking information
		for(nCount = 1; nCount<nPatNumber; nCount++)
		{
			smllink.ptable.setValue(18, nCount, trackid[nCount]);
			smllink.ptable.setValue(19, nCount, patid[nCount]);
		}	
		
		
		
	}
	

	//function calculating the number of detected particles per track
	void calculate_Tracks_Lengths()
	{
		int nCount,i;
		int nCurrTrack, nMaxVal, nIniPosition;
		
		nLength = new double [(int) nPatNumber];
		
		trackid   = smllink.ptable.getColumnAsDoubles(18);		
		patid     = smllink.ptable.getColumnAsDoubles(19);
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
					nLength[i] = nMaxVal;
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
					nLength[i] = nMaxVal;
				}
			}
			
		}
		
	
	}
	//function adding the number of detected particles per track to results table
	void add_Tracks_Lengths()
	{
		int nCount;
		//adding to final table
		smllink.ptable.setValue("Track_Det_N", 0, nLength[0]);
		for(nCount = 1; nCount<nPatNumber; nCount++)
			smllink.ptable.setValue(20, nCount, nLength[nCount]);
		
	}
	
	void mark_Last_Detections(int nLastParts)
	{
		double [] nMark = new double [(int) nPatNumber];

		int nCount;
		int nCurrTrack, nCurrLength;
		boolean bMark;
						
		//lengths of tracks
		nLength = smllink.ptable.getColumnAsDoubles(20);
		
		nCurrTrack = (int) trackid[0];
		nCurrLength = (int) nLength[0];
		if(nCurrLength>=nLastParts)
		{
			bMark = true;
			if(nCurrLength==nLastParts)
				nMark[0] = 1;
		}
		else 
			bMark = false;
		for(nCount = 0; nCount<nPatNumber; nCount++)
		{		
			if(trackid[nCount]==nCurrTrack)
			{
				if(bMark)
					if (patid[nCount]+nLastParts>nCurrLength)
						nMark[nCount] = 1;			
			}
			else
			{
				nCurrTrack = (int) trackid[nCount]; 
				nCurrLength = (int) nLength[nCount];
				if(nCurrLength>=nLastParts)
				{
					bMark = true;
					if(nCurrLength==nLastParts)
						nMark[nCount] = 1;			
				}
				else 
					bMark = false;
			}
			
		}
		//adding to final table
		
		smllink.ptable.setValue("Last_detections", 0, nMark[0]);

		for(nCount = 1; nCount<nPatNumber; nCount++)
			smllink.ptable.setValue(21, nCount, nMark[nCount]);				

	}
	
	
	//function counting number of unique frames	
	int uniqFrames()
	{
		int i, nCount=1, nLength = f.length;
		int nCurrFrame;
		boolean bContinue = true;
		
		i=0;
		nCurrFrame = (int)f[i];		
		while (bContinue)
		{
			i++;
			if(i==nLength)
				bContinue = false;
			else
			{
				if((int)f[i]!=nCurrFrame)
				{
					nCurrFrame=(int)f[i];
					nCount++;
					
				}
			}
			
		}
		
		return nCount;
		
	}
	
	void addTrack(Color trackColor, int nFrame)
	{
		int i;
		float [][] xy_roi;
		double dx, dy;
		Roi polyLine;
		Roi startPoint;
		if(TrackCoords.size()>1)
		{
			xy_roi = new float[2][TrackCoords.size()];
			for (i=0;i<TrackCoords.size();i++)
			{
				xy_roi[0][i]= (float)(TrackCoords.get(i)[0]);
				xy_roi[1][i]= (float)(TrackCoords.get(i)[1]);
			}
			polyLine = new PolygonRoi(xy_roi[0], xy_roi[1], xy_roi[0].length, Roi.POLYLINE);		
			polyLine.setStrokeColor(trackColor);
			polyLine.setStrokeWidth(0.5);
			polyLine.setPosition(nFrame);
			ovTracks.add(polyLine);
		}
		
		//starting point of the track
		dx=TrackCoords.get(0)[0];
		dy=TrackCoords.get(0)[1];
		dx-=0.5;
		dy-=0.5;
		startPoint = new OvalRoi(dx,dy, 1.0,1.0);
		startPoint.setStrokeColor(trackColor);
		startPoint.setStrokeWidth(0.5);
		startPoint.setPosition(nFrame);				
		ovTracks.add(startPoint);
		
	}
	Color SwitchingColorPalette(int nNumber)	
	{
				
		switch(nNumber%8){
		case 0:
			return Color.blue;
		case 1:
			return Color.cyan;
		case 2:
			return Color.green;
		case 3:
			return Color.magenta;
		case 4:
			return Color.orange;
		case 5:
			return Color.pink;
		case 6:
			return Color.red;
		case 7:
			return Color.yellow;
		default:
			return Color.white;		
		
		}				
	}
	
	void addTracksToOverlay()
	{
		int nCount=0;
		int i, nPrevFrame, nCurrFrame;
		int nTrackCount=0;
		int nTrackLength, nTrackNumber;
		
		//table was sorted by track number
		//so let's update values
		//frame number
		x   = smllink.ptable.getColumnAsDoubles(1);		
		y   = smllink.ptable.getColumnAsDoubles(2);
		//coordinates
		f   = smllink.ptable.getColumnAsDoubles(13);
		trackid = smllink.ptable.getColumnAsDoubles(18);
		
		//go through all table
		while (nCount<nPatNumber)
		{
			if(nLength[nCount]<2)
				//do not show not connected particles
				nCount++;
			else
			{
				//first point of track
				TrackCoords.clear();
				TrackCoords.add(new double[] {x[nCount], y[nCount]});
				nPrevFrame = (int) f[nCount];
				nTrackLength = (int) nLength[nCount];
				nTrackNumber = (int)trackid[nCount];
				addTrack(SwitchingColorPalette(nTrackNumber), nPrevFrame);
				//remaining points
				for (nTrackCount=nCount+1; nTrackCount<nCount+nTrackLength;nTrackCount++)
				{
					nCurrFrame = (int) f[nTrackCount];
					
					//in case there is a gap in the track
					if(nCurrFrame>nPrevFrame+1)
						for(i=nPrevFrame+1;i<nCurrFrame; i++)
							addTrack(SwitchingColorPalette(nTrackNumber), i);
					//current point of track
					TrackCoords.add(new double[] {x[nTrackCount],y[nTrackCount]});
					addTrack(SwitchingColorPalette(nTrackNumber), nCurrFrame);
					nPrevFrame = nCurrFrame;
				}
				//jump to next track
				nCount+=nTrackLength;
			}
				
			
		}
		
		
	}
}


