/*
 * Title:        EdgeCloudSim - Idle/Active Load Generator implementation
 * 
 * Description: 
 * IdleActiveLoadGenerator implements basic load generator model where the
 * mobile devices generate task in active period and waits in idle period.
 * Task interarrival time (load generation period), Idle and active periods
 * are defined in the configuration file.
 * 
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 * Copyright (c) 2017, Bogazici University, Istanbul, Turkey
 */

package edu.boun.edgecloudsim.task_generator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import edu.boun.edgecloudsim.core.SimSettings;
import edu.boun.edgecloudsim.core.SimSettings.APP_TYPES;
import edu.boun.edgecloudsim.utils.EdgeTask;
import edu.boun.edgecloudsim.utils.SimLogger;
import edu.boun.edgecloudsim.utils.SimUtils;

public class IdleActiveLoadGenerator extends LoadGeneratorModel{

	public IdleActiveLoadGenerator(int _numberOfMobileDevices, double _simulationTime, String _simScenario) {
		super(_numberOfMobileDevices, _simulationTime, _simScenario);
	}

	@Override
	public void initializeModel() {
		
		double[] arrList = new double[26000];
		
		File file = new File("/home/c00303945/ResearchWork/Fall_2018/EdgeCloudSim/AnnSimulation/Latest/EdgeCloudSim_A/denseArrival1.dat");
        
	    BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	       
	    try {
	           int index = 0;
	           String line = null;
	           while ((line=br.readLine())!=null) {
	              
	               String[] lineArray = line.split(",");
	               arrList[index]= Double.parseDouble(lineArray[2]);
	               index++;
	           							} 
	       }catch (FileNotFoundException e) {
	           e.printStackTrace();
	       } catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	    
	    
		taskList = new ArrayList<EdgeTask>();
		
		//exponential number generator for file input size, file output size and task length
		ExponentialDistribution[][] expRngList = new ExponentialDistribution[SimSettings.APP_TYPES.values().length][3];
		
		//NormalDistribution 
		NormalDistribution[][] nrmlRngList = new NormalDistribution[SimSettings.APP_TYPES.values().length][3]; 
		
		//create random number generator for each place
		for(int i=0; i<SimSettings.APP_TYPES.values().length; i++) {
			if(SimSettings.getInstance().getTaskLookUpTable()[i][0] ==0)
				continue;
			
			expRngList[i][0] = new ExponentialDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][5]);
			expRngList[i][1] = new ExponentialDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][6]);
			expRngList[i][2] = new ExponentialDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][7]);
		}
		
		
		//create random number generator for each place
	    for(int i=0; i<SimSettings.APP_TYPES.values().length; i++) {
					if(SimSettings.getInstance().getTaskLookUpTable()[i][0] ==0)
						continue;
					
					nrmlRngList[i][0] = new NormalDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][5],10);
					nrmlRngList[i][1] = new NormalDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][6],10);
					nrmlRngList[i][2] = new NormalDistribution(SimSettings.getInstance().getTaskLookUpTable()[i][7],100);
				}
		
		
		
		//Each mobile device utilizes an app type (task type)
		for(int i=0; i<numberOfMobileDevices; i++) {
			APP_TYPES randomTaskType = null;
			double taskTypeSelector = SimUtils.getRandomDoubleNumber(0,100);
			double taskTypePercentage = 0;
			for (SimSettings.APP_TYPES taskType : SimSettings.APP_TYPES.values()) {
				taskTypePercentage += SimSettings.getInstance().getTaskLookUpTable()[taskType.ordinal()][0];
				if(taskTypeSelector <= taskTypePercentage){
					randomTaskType = taskType;
					break;
				}
			}
			if(randomTaskType == null){
				SimLogger.printLine("Impossible is occured! no random task type!");
				continue;
			}
			
			double poissonMean = SimSettings.getInstance().getTaskLookUpTable()[randomTaskType.ordinal()][2];
			double activePeriod = SimSettings.getInstance().getTaskLookUpTable()[randomTaskType.ordinal()][3];
			double idlePeriod = SimSettings.getInstance().getTaskLookUpTable()[randomTaskType.ordinal()][4];
			double activePeriodStartTime = SimUtils.getRandomDoubleNumber(10, 10+activePeriod);  //start from 10th seconds
			double virtualTime = activePeriodStartTime; //start from 10th seconds

			ExponentialDistribution rng = new ExponentialDistribution(poissonMean);
			int index = 0;
			while(virtualTime < simulationTime) {
				double interval = rng.sample();

				if(interval <= 0){
					SimLogger.printLine("Impossible is occured! interval is " + interval + " for device " + i + " time " + virtualTime);
					continue;
				}
				//SimLogger.printLine(virtualTime + " -> " + interval + " for device " + i + " time ");
				virtualTime += interval;
				
				if(virtualTime > activePeriodStartTime + activePeriod){
					activePeriodStartTime = activePeriodStartTime + activePeriod + idlePeriod;
					virtualTime = activePeriodStartTime;
					continue;
				}
				
				taskList.add(new EdgeTask(i,randomTaskType, virtualTime+arrList[index], nrmlRngList));
				index++;
			}
			
		}
	}

}
