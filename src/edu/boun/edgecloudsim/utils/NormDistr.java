package edu.boun.edgecloudsim.utils;

import edu.boun.edgecloudsim.core.SimSettings;

public class NormDistr {
	
	protected double mean;
	protected double stdev;
	protected int hostID;
	protected SimSettings.APP_TYPES taskType;
	
	@SuppressWarnings("unused")
	private NormDistr() {
		
	};
	
	public NormDistr(double _mean, double _stdev) {
		
		mean = _mean;
		stdev = _stdev;
	}

	public double getMean() {
		return mean;
	}

	public void setMean(double mean) {
		this.mean = mean;
	}

	public double getStdev() {
		return stdev;
	}

	public void setStdev(double stdev) {
		this.stdev = stdev;
	}

	public int getHostID() {
		return hostID;
	}

	public void setHostID(int hostID) {
		this.hostID = hostID;
	}

	public SimSettings.APP_TYPES getTaskType() {
		return taskType;
	}

	public void setTaskType(SimSettings.APP_TYPES taskType) {
		this.taskType = taskType;
	}
	
	
	

}
