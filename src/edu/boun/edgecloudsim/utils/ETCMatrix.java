package edu.boun.edgecloudsim.utils;

import java.util.HashMap;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ETCMatrix {
	/**
	 *  The matrix is holding a mean and a standard deviation 
	 *  for each type of task on every base station
	 */
	protected NormDistr[][] etcMatrix = null;
	
	/**
	 * The number of Virtual Machines in the simulation and
	 * the number of Task Types
	 */
	protected int dataCenterNum = 0;
	protected int taskTypeNum = 0;
	
	protected HashMap<String, NormDistr> distributions;

	private String[] arrays;
	
	/**
	 * A private constructor to ensure that only 
	 * an correct initialized matrix could be created
	 */
	
	@SuppressWarnings("unused")
	private ETCMatrix() {
		
	};
	
	/**
	 * A parameterized constructor
	 * @param vmTotalNum takes the total number of VMs in the simulation 
	 * @param cloudletTotalnum takes the total number of Cloudlet Types in the simulation
	 */
	
	public ETCMatrix(int _dataCenterNum, int _taskTypeNum, HashMap<String, NormDistr> _distributions) {
		
		this.dataCenterNum = _dataCenterNum;
		this.taskTypeNum = _taskTypeNum;
		this.etcMatrix = new NormDistr[dataCenterNum][taskTypeNum];
		this.distributions = _distributions;
		
		for(String key: distributions.keySet()) {
			
			arrays = key.split("\\.");
			
			int row = Integer.parseInt(arrays[0]);
			int column = Integer.parseInt(arrays[1]);
			NormDistr newDistribution = distributions.get(key);
			this.etcMatrix[row][column] = newDistribution;
			
			
			}
		}
	
	public int getDataCnum() {
		return dataCenterNum;
	}
	

	
	/**
	 * Returns the 
	 * @param dataCenterID
	 * @param cloudletType
	 * @return
	 */
	
	public NormDistr getDistribution(int dataCenterID, int taskType) {
			
			if (taskType > taskTypeNum || dataCenterID > dataCenterNum) {
				throw new ArrayIndexOutOfBoundsException("The Virtual Machine or the Task Type does not exist in this ETC");
			}
			
			return etcMatrix[dataCenterID][taskType];
		}

		/**
		 * Inputs the LognormalDistributions from the HashMap into
		 * the matrix (2d-array)
		 */
	

	public double getMu(int dataCenter, int taskType) {
		
		if(etcMatrix[dataCenter][taskType] == null) {
			return 0;
		}
		
		NormDistr distr = etcMatrix[dataCenter][taskType];
	
		return distr.mean;
	}
	
	public double getSigma(int dataCenter, int taskType) {
		
		if(etcMatrix[dataCenter][taskType] == null) {
			return 0;
		}
		
		NormDistr distr = etcMatrix[dataCenter][taskType];
	
		return distr.stdev;
	}
	
	public double getWorseCaseTime(int dataCenter, int taskType) {
		
		if(etcMatrix[dataCenter][taskType] == null) {
			return 0;
		}
		
		NormDistr distr = etcMatrix[dataCenter][taskType];
		return distr.mean+distr.stdev;
		
	}
	
	public double getProbability(int dataCenter, int taskType, double deadLine) {
		
		if(etcMatrix[dataCenter][taskType] == null) {
			return 0.0;
		}
		NormDistr distr = etcMatrix[dataCenter][taskType];
		
		if(distr.mean == 0 || distr.stdev == 0) {
			return 1.0;
		}
		
		NormalDistribution newDistr = new NormalDistribution(distr.mean, distr.stdev);
		
		return newDistr.cumulativeProbability(deadLine);
		
		
	}
	
	
	public void printMatrix() {
		
		for(int i= 0; i < dataCenterNum; i++) {
			for(int j = 0; j < taskTypeNum; j++) {
				NormDistr distr = etcMatrix[i][j];
				if(distr == null) {
					break;
				}
				else {
					System.out.println("Task type " + j + " on DataCenter " + i + " has mean of " + distr.getMean() + " and stdv of "+ distr.getStdev());
				}
			}
				
			
		}
		
	}


}
